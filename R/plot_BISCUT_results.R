### Author: Juliann Shih, jshih@broadinstitute.org
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: July 24, 2023

### License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
### Dependencies: tested using R 4.1 and Python 3.9
#'
#' Visualize BISCUT peaks
#' 
#' Summary plots are generated automatically with a BISCUT run. If you like, you can use this function
#' to generate custom plots for arbitrary sets of peaks.
#'
#' @import ggplot2
#' @param plot_data A data.table (or data.frame) with peak plotting information. For format, see the
#'   text files in the \code{summary_plots} subdirectory of a BISCUT output directory.
#' @return A ggplot, or NULL if the plot_data is empty.
#' @export
plot_BISCUT_results <- function(plot_data, title, xlab, abslocs = get_chromosome_coordinates()) {
  if(! is.data.frame(plot_data)) {
    stop('Expected plot_data to be data.frame')
  }
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  df = plot_data
  gg = ggplot() + theme_classic() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                       limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                       breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                       labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
    geom_rect(data=df, mapping=aes(xmin=ymin, xmax=ymax, ymin=-xmax, ymax=-xmin), color=df$color, fill=df$color) +
    geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
    geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
    geom_vline(xintercept= 0) +
    scale_x_continuous(name=xlab,
                       limits = c(-(max(abs(df$ymin),abs(df$ymax))+1), max(abs(df$ymin),abs(df$ymax))+1))
  return(gg)
}