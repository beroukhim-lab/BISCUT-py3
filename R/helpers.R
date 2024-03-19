#' Filter out largest/smallest tel/cent events
#' @param df breakpoints data.frame
#' @param telcent_thres Threshold for filtering out events that are very small or vary large. Events
#'   that cover a fraction of the arm that is less than telcent_thres or greater than
#'   (1-telcent_thres) are removed.
filter_big_small <- function(df, telcent_thres = 1e-3) {
  df <- df[df$percent>=telcent_thres,]
  df <- df[df$percent<=(1-telcent_thres),]
  df <- df[order(df$percent),]
  return(df)
} 

#' Read in breakpoint files from a directory
#' 
#' @param breakpoint_file_dir directory containing breakpoint files
#' @param telcent_thres Threshold for filtering out events that are very small or vary large. Events
#'   that cover a fraction of the arm that is less than telcent_thres or greater than
#'   (1-telcent_thres) are removed.
load_breakpoint_files = function(breakpoint_file_dir, telcent_thres = 1e-3) {
  if(! is.character(breakpoint_file_dir) || length(breakpoint_file_dir) != 1) {
    stop('breakpoint_file_dir should be a 1-length character (specifically, a directory path).')
  }
  if(! dir.exists(breakpoint_file_dir)) {
    stop('breakpoint_file_dir ', breakpoint_file_dir, 'could not be found.')
  }
  breakpoint_files <- list.files(path = breakpoint_file_dir, pattern = '_(amp|del)_(cent|tel).txt$',
                                 full.names = T)
  concat_background <- function(bfiles){
    d <- data.frame()
    for (f in bfiles)
    {
      df <- read.csv(f, sep="\t", header =T)
      ampdel <- ifelse(grepl('_amp', f, fixed=TRUE), 'amp', 'del')
      telcent <- ifelse(grepl('_tel', f, fixed=TRUE), 'tel', 'cent')
      df['amp_del'] <- ampdel
      df['tel_cent'] <- telcent
      d <- rbind(d,df)
    }
    output <- d
  }
  breakpoints <- concat_background(breakpoint_files)
  tel <- breakpoints[breakpoints$tel_cent=='tel', c(1:5)]
  tel <- filter_big_small(tel, telcent_thres = telcent_thres)
  cent <- breakpoints[breakpoints$tel_cent=='cent', c(1:5)]
  cent <- filter_big_small(cent, telcent_thres = telcent_thres)
  return(list(tel = tel, cent = cent))
}

#' Get chromosome arms
#' 
#' Returns a character vector of all chromosome arms. Acrocentromeric chromosomes are represented by
#' just chromosome names (e.g., "14"), while the arms of other chromosomes get suffixes (e.g., "2p",
#' "2q").
#' 
#' @export
get_chromosome_arms = function() {
  return(all_arms)
}

#' Get a table of chromosomal p and q arm coordinates
#' 
#' @param genome Which genome build to use. Currently, only "hg19" is available.
#' @export
get_chromosome_coordinates = function(genome = 'hg19') {
  if(! identical(genome, 'hg19')) {
    stop('Currently, only hg19 coordinates are available.')
  }
  return(validate_chr_coordinates(default_abslocs))
}

#' Reads and validates chromosome coordinates
#' 
#' Takes either a file path or a data.frame, validates the data.frame, and returns a version with
#' just the necessary columns.
#' @param chromosome_coordinates File path or data.frame.
validate_chr_coordinates = function(chromosome_coordinates) {
  if(is.character(chromosome_coordinates)) {
    if(length(chromosome_coordinates) != 1) {
      stop('chromosome coordinates should be a data.frame of chromosome coordinates or a file path.')
    }
    if(! file.exists(chromosome_coordinates)) {
      stop('Specified chromosome_coordinates file does not exist.')
    }
    chromosome_coordinates = read.table(chromosome_coordinates, sep='\t', header=T)
  } else if(! is(chromosome_coordinates, 'data.frame')) {
    stop('chromosome coordinates should be a data.frame of chromosome coordinates or a file path.')
  }
  
  required_cols = c('chromosome_info', 'size', 'p_start', 'p_end', 'q_start', 'q_end', 'centromere')
  missing_cols = setdiff(required_cols, colnames(chromosome_coordinates))
  if(length(missing_cols) > 0) {
    stop('chromosome_coordinates is missing required columns: ', paste(missing_cols, collapse = ', '), '.')
  }
  if(! all(unique(sapply(chromosome_coordinates[, required_cols], is.numeric)))) {
    stop('The following columns in chromosome_coordinates must all be numeric:\n',
         paste(required_cols, collapse = ", "), '.')
  }
  if(length(colnames(chromosome_coordinates)[colnames(chromosome_coordinates) %in% required_cols]) != length(required_cols)) {
    stop('Some column names in chromosome_coordinates are repeated.')
  }
  chromosome_coordinates = chromosome_coordinates[, required_cols]
  
  # Produce calculated columns
  chromosome_coordinates$offset = c(0, cumsum(as.numeric(chromosome_coordinates$size[1:(nrow(chromosome_coordinates) - 1)])))
  chromosome_coordinates$middle_of_arm = ceiling(chromosome_coordinates$offset + chromosome_coordinates$size/2)
  chromosome_coordinates$q_end1 = cumsum(as.numeric(chromosome_coordinates$size))
  chromosome_coordinates$q_start1 = chromosome_coordinates$centromere + chromosome_coordinates$offset
  return(as.data.frame(chromosome_coordinates))
}


#' Get a table of gene locations
#' 
#' @param genome Which genome build to use. Currently, only "hg19" is available.
#' @export
get_gene_locations = function(genome = 'hg19') {
  if(! identical(genome, 'hg19')) {
    stop('Currently, only hg19 coordinates are available.')
  }
  return(default_genelocs)
}

#' Reads and validates gene locations
#' 
#' @param gene_locations File path or data.frame.
#' @return validated gene locations data.frame
validate_gene_locations = function(gene_locations) {
  if(is.character(gene_locations)) {
    if(length(gene_locations) != 1) {
      stop('gene_locations should be a data.frame of gene_locations or a file path.')
    }
    if(! file.exists(gene_locations)) {
      stop('Specified gene_locations file does not exist.')
    }
    gene_locations = read.table(gene_locations, sep='\t', header=T)
  } else if(! is(gene_locations, 'data.frame')) {
    stop('gene_locationsshould be a data.frame of gene_locations or a file path.')
  }
  required_cols = c("Chr", "Start", "End", "Gene")
  missing_cols = setdiff(required_cols, colnames(gene_locations))
  if(length(missing_cols) > 0) {
    stop('gene_locations is missing required columns: ', paste(missing_cols, collapse = ', '), '.')
  }
  if(length(colnames(gene_locations)[colnames(gene_locations) %in% required_cols]) != length(required_cols)) {
    stop('Some column names in gene_locations are repeated.')
  }
  if(any(duplicated(gene_locations$Gene))) {
    stop('Some genes appear more than once in the gene_locations Gene column.')
  }
  return(as.data.frame(gene_locations))
}





