#' Wrapper for filter_BISCUT_knowngenes (Python)
#' 
#' @param results_dir Output directory from do_biscut().
#' @param genes Vector of gene names
#' @export
extract_gene_results = function(results_dir, genes) {
  if(! rlang::is_scalar_character(results_dir)) {
    stop('results_dir should be 1-length character (the file path of a BISCUT output directory)')
  }
  if(! dir.exists(results_dir)) {
    stop('results_dir directory not found.')
  }
  results_dir = normalizePath(results_dir)
  if(! is.character(genes)) {
    stop('genes should be type character.')
  }
  reticulate::source_python(system.file('python/combine_BISCUT_results.py', package = 'BISCUT'))
  filter_BISCUT_knowngenes(results_dir, genes)
}

