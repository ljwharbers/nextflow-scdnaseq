require(ggplot2)
require(cowplot)
require(dplyr)
require(RColorBrewer)

theme_set(theme_cowplot())

# plotProfile = function(segments, raw, bins, sc=TRUE, cn=NULL){
plotProfile = function(segments, raw, bins, sc=TRUE, linesize = 1){
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
