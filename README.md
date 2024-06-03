## Nextflow pipeline for processing of GPSeq data

### Requirements
To be able to run this pipeline you need nextflow (version 23.04 or higher) and singularity (tested on version 3.8.6) installed. This does not work on Mac.

The easiest way to install these tools is with conda package manager.  
For example using the following lines (assuming you have conda installed):  
`conda create -n nextflow -c conda-forge -c bioconda nextflow=23.10.0 singularity=3.8*`

If you don't have conda installed yet you can install and initialize it in the following way:  
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

You can then activate the environment by running:  
`conda activate nextflow`

### Testing the pipeline
First clone the pipeline (this will create a folder in your current working directory)  
```
git clone https://github.com/ljwharbers/nextflow-scdnaseq
cd nextflow-scdnaseq
```

You can test the pipeline by typing:  
`nextflow run main.nf -profile test`  

This will run the pipeline with default parameters and a test samplesheet and dataset that is included in the repository. If this runs without any issues you can carry on and run it in your own dataset.

### Running the pipeline with your own dataset
To run the pipeline with your own dataset, there are a few steps to take.

1. Specify the directory where all the fastq files are located in the command line with the following tag: `--indir {path/to/fastq/dir}`
2. Adjust `nextflow.config`.  There are some default parameters used and specified in the configuration file and, depending on your most common usecase, it is advisable to change some of these defaults.
   * If you will mostly run GPSeq on human, you can write the path to your own local reference file and bowtie2 index in the config file under `fasta` and `bwt2index`. However, I recommend using the iGenomes reference files (described further down).
   * Check `max_memory` and `max_cpus`. It is important that these do not go above your system values.
   * Go over other parameters defined within the `params { }` section in the config file and change whatever you feel fit.  

Following this you can either change the default parameters in the `nextflow.config` file or supply the parameters related to your own dataset in the command you type. I suggest that you change parameters that won't change much between runs in the `nextflow.config` file, while you specify parameters such as `indir` and `outdir` through the command line.

An example command to run this on your own data could be:  
`nextflow run main.nf --indir path/to/fastq/dir --outdir path/to/results --fasta path/to/reference.fa --fasta_index path/to/reference.fa.fai --bwaindex /path/to/bwa/index/folder`

If you do not specify a fasta file and bowtie2 index, you can specify the reference genome you want to use and it will download it from an AWS s3 bucket. For example in the following way:  
`nextflow run main.nf --indir path/to/fastq/dir --outdir path/to/results --genome GRCh38`

Downloading the fasta file and index might be slow so you can also download the files that you would need through using this tool: https://ewels.github.io/AWS-iGenomes/ Note: you need `aws` tool for this. Once you've downloaded the reference and index files you need you can change the `igenomes_base` parameter in `nextflow.config` and it will take the fasta/index files from there instead of downloading it through nextflow.

Finally, if you want to resume canceled or failed runs, you can add the tag `-resume`. Usually, I always use this tag irregardless of if I run something for the first time or if I am resuming a run.