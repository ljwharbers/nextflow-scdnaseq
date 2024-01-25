#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Script for running ASCAT.sc through commandline

packages = c("data.table", "argparser", "pbapply", "ASCAT.sc",
             "ggplot2", "cowplot", "dplyr", "ggdendro", "patchwork",
             "RColorBrewer")

invisible(
    sapply(
        packages, function(x) {
            suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
        })
    )

## Parse arguments
parser = arg_parser("Run ASCAT.sc on a directory of bamfiles")
parser = add_argument(parser, "--bams", help = "Path to bam files", nargs = Inf)
parser = add_argument(parser, "--outdir", help = "Path to output directory", nargs = 1, default = "./")
parser = add_argument(parser, "--genome",
                      help = "Build used, currently only supports hg19/38 (or GRCh37/38)", nargs = 1)
parser = add_argument(parser, "--binsize", help = "binsize (can be multiple, space separated)",
                      nargs = 1, type = "numeric")
parser = add_argument(parser, "--chroms", help = "list of chromosomes to use", default = c(1:22, "X"), nargs = "+")
parser = add_argument(parser, "--sex",  help = "Sex of samples, either 'male' or 'female'")
parser = add_argument(parser, "--min_ploidy",
                      help = "Min Ploidy value to use in grid search (default = 1.5)",
                      nargs = 1, default = 1.5, type = "numeric")
parser = add_argument(parser, "--max_ploidy",
                      help = "Max Ploidy value to use in grid search (default = 5)",
                      nargs = 1, default = 5, type = "numeric")
parser = add_argument(parser, "--min_purity",
                      help = "Min Purity value to use in grid search (default = 1)",
                      nargs = 1, default = 1, type = "numeric")
parser = add_argument(parser, "--max_purity",
                      help = "Max Purity value to use in grid search (default = 1)",
                      nargs = 1, default = 1, type = "numeric")
parser = add_argument(parser, "--psi", help = "Maximum tumour Psi (default = 5)",
                      default = 5, type = "numeric")
parser = add_argument(parser, "--multipcf", help = "Use multipcf", flag = TRUE)
parser = add_argument(parser, "--penalty",
                      help = "Segmentation alpha or multipcf penalty to use (multiple allowed, space separated)",
                      nargs = 1, type = "numeric")
parser = add_argument(parser, "--threads", help = "Number of threads to use", nargs = 1, type = "integer")
argv = parse_args(parser)


## Set params
bams  = argv$bams
workdir = argv$outdir
stopifnot(argv$genome %in% c("hg19", "GRCh37", "hg38", "GRCh38"))
genome = ifelse(argv$genome %in% c("hg19", "GRCh37"), "hg19", "hg38")
binsize = argv$binsize
chroms = argv$chroms
sex = argv$sex
prefix = ifelse(grepl("chr", chroms[1]), "chr", "")
ploidies = seq(argv$min_ploidy, argv$max_ploidy, 0.01)
purs = seq(argv$min_purity, argv$max_purity, 0.01)
psi = argv$psi
multipcf = ifelse(argv$multipcf, "TRUE", "FALSE")
penalty = argv$penalty
threads = argv$threads

# Run ASCAT.sc
res = run_sc_sequencing(tumour_bams = bams,
                    build = genome,
                    allchr = chroms, ## Needs a "chr" instead of "" if reference genome has 'chr' prefix
                    sex = rep(sex, length(bams)),
                    chrstring_bam = prefix,  ## Needs a "chr" instead of "" if reference genome has 'chr' prefix
                    purs = purs, ## purity grid values
                    ploidies = ploidies, ## average ploidy grid values
                    maxtumourpsi = psi, ## maximum tumour ploidy
                    binsize = binsize, ## bin size - reduce if enough sequencing reads (look at dpb (depth per bin) value in plots can go down to 100bp or even lower)
                    projectname = "",
                    MC.CORES = threads, ##number of cores available
                    multipcf = multipcf, ##use multipcf for multi-track segmentation if multi-sample sequencing
                    segmentation_alpha = penalty,
                    predict_refit = FALSE)


####################################################################
saveRDS(res, paste0(workdir, "/ASCAT.sc_results.rds")) ## save full results
####################################################################



# plotProfile = function(segments, raw, bins, sc=TRUE, cn=NULL){
plotProfile = function(segments, raw, bins, sc=TRUE, linesize = 1){
  # Set theme
  theme_set(theme_cowplot())

  # Process bins
  bins = bins
  segments = segments
  raw = raw
  bins[, bin := seq_along(chr)]
  bins[, end_cum := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]
  
  # Make chr_bounds
  chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
  chr_bounds = chr_bounds %>% 
    mutate(mid = round(min + (max-min) / 2,0),
           end_bp=cumsum(as.numeric(chrlen_bp)), 
           start_bp = end_bp - chrlen_bp, 
           mid_bp = round((chrlen_bp / 2) + start_bp, 0))
  
  
  
  if(sc) {
    # dt = cbind(bins, segments, raw * cn)
    dt = cbind(bins, segments, raw)
    setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "cn", "raw"))
    
    dt[, cn := ifelse(cn < 11, cn, 11)]
    dt[, col := ifelse(cn < 11, as.character(cn), "10+")]
    dt[, col := factor(col, levels = c(as.character(0:10), "10+"))]
    
    colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
               "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
    names(colors) = c(as.character(0:10), "10+")
    
    # save plot
    plot = ggplot(dt, aes(x = bin)) +
      geom_point(aes(y = raw, color = col), size = 0.7) +
      geom_point(aes(y = cn), size = linesize) +
      scale_color_manual(values = colors, drop = F) +
      scale_y_continuous(labels=scales::comma_format(accuracy = 1), breaks = scales::pretty_breaks(6), limits = c(0, 12)) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(y = "Copy Number", x = "") +
      geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
      geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", inherit.aes = F) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  else {
  dt = cbind(bins, segments, raw)
  setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "segment", "raw"))
  
  # Assign colors
  dt[segment > log2(2.5/2), col := "Amplification"]
  dt[segment < log2(1.5/2), col := "Deletion"]
  dt[is.na(col), col := "Neutral"]

  colors = c(brewer.pal(3, "Set1")[1:2], "#c1c1c1")
  names(colors) = c("Amplification", "Deletion", "Neutral")
  
  # save plot
  plot = ggplot(dt, aes(x = bin)) +
    geom_point(aes(y = raw, color = col), size = 0.7) +
    geom_point(aes(y = segment), size = linesize) +
    scale_color_manual(values = colors, drop = F) +
    scale_y_continuous(limits = c(-4, 4), labels=scales::comma_format(accuracy = 1), breaks = scales::pretty_breaks(6)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(y = "Log2 ratio", x = "") +
    geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
    geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", inherit.aes = F) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  }
  return(plot)
}

# Function to plot genomewide copy number profile heatmaps
plotHeatmap = function(profiles, bins, order, dendrogram = TRUE, linesize = .5, rownames = FALSE, annotation = NULL) {
  # Check that order is not specified while dendrogram is also requested
  if(!missing(order)) dendrogram = FALSE
  # Load libraries
  sapply(packages, require, character.only = T)
  
  # Make sure they are data.tables
  setDT(bins)
  setDT(profiles)
  
  # Set cowplot theme
  theme_set(theme_cowplot())
  # Get cumulative locations
  bins[, bin := seq_along(chr)]
  bins[, end_cum := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]
  
  # Make chr_bounds
  chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
  chr_bounds = chr_bounds %>% 
    mutate(mid = round(min + (max-min) / 2,0),
           end_bp=cumsum(as.numeric(chrlen_bp)), 
           start_bp = end_bp - chrlen_bp, 
           mid_bp = round((chrlen_bp / 2) + start_bp, 0))
  
  #Colors
  colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
             "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
  
  names(colors) = c(as.character(0:10), "10+")

  dt = data.table(cbind(bins, profiles))
  
  # Set Theme depending on rownames
  if(rownames) {
    custom_theme = theme(axis.ticks.y = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title = element_blank())
  } else {
    custom_theme = theme(axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title = element_blank())
  }
  

  if(dendrogram) {
    # Clustering
    hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
    dhc = as.dendrogram(hc)
    
    # Rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    
    # Plot Dendrogram
    dendro = ggplot(ggdendro::segment(ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0.004, 0.004)) +
      theme_dendro() 
    
    # Prepare for heatmap
    dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
    dt_melt[, value := factor(value)]
    dt_melt[as.numeric(value) > 10, value := "10+"]
    dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
    
    # Set sample order
    dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
    
    # Plot heatmap
    heatmap = ggplot(dt_melt) +
      geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), linewidth = linesize) + 
      coord_flip() +
      scale_color_manual(values = colors, drop = F) +
      labs(color = "Copy Number") + 
      scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
      geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, linewidth = .8) +
      custom_theme
    
    
    if(!is.null(annotation)) {
      setnames(annotation, c("sample", "variable", "value"))
      
      # Plot annotation
      annotation[, sample := factor(sample, levels = ddata$labels$label)]
      annot_plt = ggplot(annotation, aes(x = 1, y = sample, fill = value)) +
        geom_bar(stat = "identity", width = 1) +
        scale_fill_npg() +
        theme_void() +
        theme(legend.position = "right")
      
      combined = dendro + annot_plt + heatmap + plot_layout(ncol = 3, widths = c(0.2, 0.05, 2), guides = "collect")
      
    } else {
      combined = cowplot::plot_grid(dendro, heatmap, ncol = 2, rel_widths = c(0.1, 1))
    }
    
    return(combined)
  }
  
  # Prepare for heatmap
  dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
  dt_melt[, value := factor(value)]
  dt_melt[as.numeric(value) > 10, value := "10+"]
  dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
  
  if(!missing(order)) {
    # Set sample order
    dt_melt[, variable := factor(variable, levels = order)]
  } else {
    # Clustering
    hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
    dhc = as.dendrogram(hc)
    
    dt_melt[, variable := factor(variable, levels = labels(dhc))]
  }
  
  # Plot heatmap
  heatmap = ggplot(dt_melt) +
    geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), linewidth = linesize) +    
    coord_flip() +
    scale_color_manual(values = colors, drop = F) +
    labs(color = "Copy Number") + 
    scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
    geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, linewidth = .8) +
    custom_theme
  
  
  return(heatmap)
}


# Function to save a plot in different file formats and show it within the notebook
save_and_plot <- function(x, bname, width, height, dpi=300, 
                          use.cowplot=FALSE, ncol=1, nrow=1, base_height=4, base_width=NULL, 
                          base_aspect_ratio=1, plot=FALSE, onefile = TRUE, pointsize = 8, output=c("cairo_ps", "cairo_pdf", "postscript", "png")){
  
  # Create dir if not exists
  if(!dir.exists(dirname(bname))) {
    dir.create(dirname(bname), recursive = T)
  }
  
  if( !use.cowplot ){
    if("png" %in% output){
      png(filename=file.path(paste0(bname, ".png")), units="in", 
          height=height, width=width, res=dpi)
      print(x)
    }
    if("cairo_ps" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      cairo_ps(filename=file.path(paste0(bname, "_cairo_ps.eps")), 
             onefile = TRUE, height=height, width=width, family="Helvetica", 
             pointsize=pointsize, antialias="none")
    print(x)
    }
    if("cairo_pdf" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      cairo_pdf(filename=file.path(paste0(bname, "_cairo_pdf.pdf")), 
              onefile = TRUE, height=height, width=width, family="Helvetica", 
              pointsize=pointsize, antialias="none")
    print(x)
    }
    if("postscript" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      postscript(file=file.path(paste0(bname, "_postscript.eps")), 
               onefile = TRUE, paper="special", height=height, width=width, 
               family="Helvetica", pointsize=pointsize, horizontal=FALSE)
    print(x)
    }
    # while(length(dev.list())>0) invisible(dev.off())
    # svg(file=file.path(paste0(bname, ".svg")), 
    #            onefile = TRUE, height=height, width=width, 
    #            family="Helvetica", pointsize=8)
    # print(x)
    # while(length(dev.list())>0) invisible(dev.off())
  }else{
    if("png" %in% output){
      save_plot(x, filename=file.path(paste0(bname, ".png")), ncol = 
                  ncol, nrow = nrow, base_height = base_height, base_width = base_width, 
                base_aspect_ratio = base_aspect_ratio, dpi=dpi)
    }
    if("cairo_ps" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      save_plot(x, filename=file.path(paste0(bname, "_cairo_ps.eps")), 
                ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                  base_width, base_aspect_ratio = base_aspect_ratio, device=cairo_ps)
    }
    if("cairo_pdf" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      save_plot(x, filename=file.path(paste0(bname, "_cairo_pdf.pdf")), 
              ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device=cairo_pdf)
    }
    if("postscript" %in% output){
      while(length(dev.list())>0) invisible(dev.off())
      save_plot(x, filename=file.path(paste0(bname, "_postscript.eps")), 
              ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device="ps")
    }
    while(length(dev.list())>0) invisible(dev.off())
    
    # save_plot(x, filename=file.path(paste0(bname, ".svg")), 
    #           ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
    #             base_width, base_aspect_ratio = base_aspect_ratio, device="svg")
    # while(length(dev.list())>0) invisible(dev.off())
  }
  if( plot )
    print(x)
  dev.off()
}

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
saveRDS(list(bins, profiles, raw), paste0(workdir, "/ASCAT.sc_profiles.rds"))

# Plot profiles
invisible(pblapply(colnames(profiles), function(cell) {
  plotProfile(profiles[[cell]], raw[[cell]], bins) |> save_and_plot(paste0(workdir, "/plots/profiles/", cell), height = 6, width = 14, output = "png")
}, cl = 38))

# Plot Heatmap
invisible(plotHeatmap(profiles, bins, linesize = 1.75) |> save_and_plot(paste0(workdir, "/plots/genomewideheatmap"), height = 4, width = 20, output = "png"))
