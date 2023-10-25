## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

####################################################################
packages = c("data.table", "argparser", "pbapply", "ASCAT.sc")
invisible(
    sapply(
        packages, function(x) {
            suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
        })
    )

source("plotProfile.R")
source("plotHeatmap.R")
####################################################################

## Parse arguments
parser = arg_parser("Run ASCAT.sc on a directory of bamfiles")
parser = add_arguments(parser, "--bamdir", help = "Path to counts", nargs = 1)
parser = add_arguments(parser, "--outdir", help = "Path to output directory", nargs = 1, default = "./")
parser = add_arguments(parser, "--genome", help = "Build used, currently only supports hg19/38 or GRCh37/38", nargs = 1)
parser = add_arguments(parser, "--binsize", help = "binsize", nargs = 1)
parser = add_arguments(parser, "--chroms", help = "list of chromosomes to use", default = c(1:22, "X"), nargs = Inf)
parser = add_arguments(parser, "--sex",  help = "Sex of samples, either 'male' or 'female'")
parser = add_argument(parser, "--multipcf", help = "whether to use multipcf or not", nargs = 1)
parser = add_argument(parser, "--penalty", help = "Segmentation alpha or multipcf penalty to use", nargs = 1)
parser = add_arguments(parser, "--threads", help = "Number of threads to use", nargs = 1, type = "integer")
argv = parse_args(parser)

####################################################################
## Set params
bamdir  = argv$bamdir
stopifnot(argv$genome %in% c("hg19", "GRCh37", "hg38", "GRCh38"))
genome = ifelse(argv$genome %in% c("hg19", "GRCh37"), "hg19", "hg38")
binsizes = argv$binsizes
chroms = argv$chroms
sex = argv$sex
prefix = ifelse(grepl("chr", chroms[1]), "chr", "")
multipcf = argv$multipcf
penalty = argv$penalty
threads = argv$threads

####################################################################
bams = dir(bamdir, pattern = "bam$", full = TRUE, rec = FALSE)
workdir = argv$outdir
setwd(workdir)
####################################################################

####################################################################
res = run_sc_sequencing(tumour_bams = bams,
                        build = genome,
                        allchr = chroms, ## Needs a "chr" instead of "" if reference genome has 'chr' prefix
                        sex = rep(sex, length(bams)),
                        chrstring_bam = prefix,  ## Needs a "chr" instead of "" if reference genome has 'chr' prefix
                        purs = 1, ## purity grid values
                        ploidies = seq(1.6, 5, 0.01), ## average ploidy grid values
                        maxtumourpsi = 5, ## maximum tumour ploidy
                        binsize = binsize, ## bin size - reduce if enough sequencing reads (look at dpb (depth per bin) value in plots can go down to 100bp or even lower)
                        projectname = "",
                        MC.CORES = threads, ##number of cores available
                        multipcf = multipcf, ##use multipcf for multi-track segmentation if multi-sample sequencing
                        segmentation_alpha = penalty,
                        predict_refit = FALSE)
####################################################################
saveRDS(res, "ASCAT.sc_results.rds") ## save results for reuse later on
####################################################################

####################################################################
# Get bins
bins = pblapply(names(res$lSe), function(chr) {
  data.table(chr = chr, 
             start = res$lSe[[chr]]$start,
             end = res$lSe[[chr]]$end)
}, cl = 38)
bins = rbindlist(bins)

# Get profiles
profiles = pblapply(res$allProfiles, function(cell) {
  dt = as.data.table(cell)
  rep(dt$total_copy_number, as.numeric(dt$num.mark)) |> as.numeric()
}, cl = 38)
profiles = do.call(cbind, profiles) |> data.table()

# Get raw segments
raw = pblapply(colnames(profiles), function(cell) {
  dt = rbindlist(res$allTracks.processed[[cell]]$lCTS)
  # return(dt$smoothed * res$allSolutions[[dt$file[1]]]$ploidy)
  return(res$allSolutions[[dt$file[1]]]$ploidy * 2 ^ (dt$smoothed))
  
})
raw = do.call(cbind, raw) |> data.table()
setnames(raw, colnames(profiles))

# Save profiles
saveRDS(list(bins, profiles, raw), "ASCAT.sc_profiles.rds")

# Plot profiles
pblapply(colnames(profiles), function(cell) {
  plotProfile(profiles[[cell]], raw[[cell]], bins) |> save_and_plot(paste0("./plots/", cell), height = 6, width = 14, output = "png")
}, cl = 38)

# Plot Heatmap
plotHeatmap(profiles, bins, linesize = 1.75) |> save_and_plot(paste0("./plots/genomewideheatmap"), height = 4, width = 20, output = "png")
