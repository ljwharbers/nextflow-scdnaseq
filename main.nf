// Get genome attributes
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_index = WorkflowMain.getGenomeAttribute(params, 'fasta_index')
params.bwaindex = WorkflowMain.getGenomeAttribute(params, 'bwa')

// Get penalty to use
params.penalty = params.multipcf ? params.multipcf_penalty : params.segmentation_alpha

// Print pipeline info
log.info """\

		S I N G L E - C E L L   D N A - S E Q 
		=====================================
		Input directory	    : ${params.indir}
		Output directory    : ${params.outdir}
		Reference genome    : ${params.genome ? params.genome : params.fasta.replaceAll(".+\\/", "")}
		"""
		.stripIndent()

// Align fastq files with bwa-mem2
process ALIGN {
    label "process_highcpu"
    tag "bwa mem on ${sample}"
    
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'
    
    input:
        tuple val(sample), path(reads)
        path(index)

    output:
        tuple val(sample), path("${sample}.bam"), emit: bam
        
    script:
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa mem \\
            -t $task.cpus \\
            \$INDEX \\
            $reads \\
            | samtools sort -@ $task.cpus -o ${sample}.bam
        """
}

// Deduplicate bamfiles using MarkDuplicates
process MARKDUPLICATES {
    tag "MarkDuplicatesSpark on ${sample}"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0"
	
    publishDir "${params.outdir}/bamfiles", mode:'copy'

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${sample}.dedup.bam"),     emit: bam
        tuple val(sample), path("${sample}.dedup.bai"),     emit: bai
        path "${sample}.dedup.metrics", emit: metrics

    script:
        def avail_mem = 3072
        if (!task.memory) {
            log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }

        """
        gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
            MarkDuplicates \\
            --INPUT ${bam} \\
            --OUTPUT ${sample}.dedup.bam \\
            --METRICS_FILE ${sample}.dedup.metrics \\
            --CREATE_INDEX true \\
            --TMP_DIR . 
        """
}


// Get list of chromosomes to use
process GETCHROMS {
    tag "Get list of chromosomes"
    label 'process_single'
    
    container "https://depot.galaxyproject.org/singularity/grep:2.14--1"
    
    input:
        path fasta_index
        val sex
        
    output:
        stdout
        
    script:
        def grep_expression
        grep_expression = sex == "male" ? '^(chr)?[0-9]?[0-9]?[XY]?\\t' : '^(chr)?[0-9]?[0-9]?[X]?\\t'
        """
        cat ${fasta_index} | grep -P '${grep_expression}' | cut -f1 | tr '\n' ' '
        """
}

// CNA calling using ASCAT.sc
process ASCATSC {
    tag "Copy number calling using ASCAT.sc using binsize: ${binsize} and penalty: ${penalty}"
    label 'process_high'
    
    container "library://ljwharbers/scdnaseq/ascatsc_renv:0.0.1"
    
    publishDir "${params.outdir}/ASCATsc/${binsize}/${penalty}", mode:'copy'
    
    input:
        path bams
        path bam_bai
        val genome
        val sex
        val chroms
        tuple val (binsize), val(penalty) 
        val multipcf
        val min_ploidy
        val max_ploidy
        val min_purity
        val max_purity
        
    output:
        path "ASCAT.sc_results.rds"
        path "ASCAT.sc_profiles.rds"
        path "plots/genomewideheatmap.png"
        path "plots/profiles/*"
        
    script:
        def multipcf_arg
        multipcf_arg = params.multipcf ? "--multipcf" : ""
        
        def chromosomes
        chromosomes = chroms.collect { it[0] }.join(",")
        
        """
        ASCATsc_run.R --bams ${bams} --outdir \$PWD --genome ${genome} --binsize ${binsize} \\
        --sex ${sex} --chroms ${chroms} --penalty ${penalty} --threads ${task.cpus}  \\
        --min_ploidy ${min_ploidy} --max_ploidy ${max_ploidy} --min_purity ${min_purity} \\
        --max_purity ${max_purity} ${multipcf_arg} 
        """
}

process FLAGSTAT {
    tag "samtools flagstat on ${sample}"
    label 'process_single'

    container "https://depot.galaxyproject.org/singularity/samtools:1.8--4"

    input:
        tuple val(sample), path(bam)

    output:
        path "${sample}.flagstat" 
        
    script:
        """
        samtools flagstat ${bam} > ${sample}.flagstat
        """
}

// Run Fastqc on initial fastq files
process  FASTQC {
	label "process_low"
	tag "FASTQC on ${sample}"

	container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
	
	input:
		tuple val(sample), path(reads)
	
	output:
		path "fastqc_${sample}_logs"
	
	script:
		"""
		mkdir fastqc_${sample}_logs
		fastqc -o fastqc_${sample}_logs -f fastq -q ${reads}
		"""
}

// Get report from MultiQC
process MULTIQC {
	label "process_medium"
	tag "MultiQC on all samples"
	
	container "https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0"
	
	publishDir "${params.outdir}/preprocessing", mode:'copy'

	input:
		path '*'

	output:
		path 'multiqc_report.html'

	script:
		"""
		multiqc .
		"""
}

// Define workflow
workflow {
    if (params.paired) {
        Channel
            .fromFilePairs(params.indir + "/*R{1,2}*{fastq,fastq.gz,fq,fq.gz}")
            .set { reads }
    } else {
        Channel
            .fromPath(params.indir + "/*{fastq,fastq.gz,fq,fq.gz}")
            .map { file -> tuple( file.baseName.replaceAll("_R.+", ""), file) }
            .set { reads }
    }

    // PREPROCESS
    aligned = ALIGN(reads, params.bwaindex) // BWA-MEM2
    deduped = MARKDUPLICATES(aligned) // GATK MARKDUPLICATES
    
    // CNA Calling
    chroms = GETCHROMS(params.fasta_index, params.sex)
    ascat_bams = deduped.bam.map { tuple( it[1] ) }
    ascat_bai = deduped.bai.map { tuple ( it[1] )} 
    
    // Make combinations of binsizes and penalties
    binsize = channel.from(params.binsize)
    penalty = channel.from(params.penalty)
    binsize.combine(penalty).set { binsize_penalty }
    
    ASCATSC(ascat_bams.collect(), ascat_bai.collect(),         
            params.genome, params.sex,            
            chroms, binsize_penalty, params.multipcf,
            params.min_ploidy, params.max_ploidy,
            params.min_purity, params.max_purity)

    // QC
    fastqc_ch = FASTQC(reads) // FASTQC
    flagstat_ch = FLAGSTAT(deduped.bam) // FLAGSTAT
    MULTIQC(flagstat_ch.mix(fastqc_ch, deduped.metrics).collect()) // MultiQC
    
}