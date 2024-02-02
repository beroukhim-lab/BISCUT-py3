#' Prepare breakpoint files for BISCUT
#' 
#' @param segment_file Copy number segmentation file path. See README for formatting details.
#' @param output_dir Where to create an output directory and save breakpoint files.
#' @param chromosome_coordinates The coordinates of p and q arms of the human chromosomes. By
#'   default, hg19-based coordinates are used. Please note that these coordinates are not
#'   biologically well-defined regions but rather what appeared to be the boundaries of regions that
#'   can be considered telomeric and (particularly) centromeric from SNP-array based TCGA copy
#'   number profiles. You might want to use your custom set of coordinates based on your copy number
#'   profiles (ex. IGV view of start/end of event near telomere/centromere). From our experience,
#'   choosing the telomeric/centromeric boundaries wide enough so that the events are not missed due
#'   to imperfections of assaying technologies is more important than genome reference differences
#'   (ex., hg19 vs. hg38).
#' @param arms Which chromosome arms to process. Defaults to all arms (see ?get_chromosome_arms()). You 
#' may choose to run a subset of arms for testing purposes.
#' @param amplitude_threshold Copy number changes less than this threshold are considered noise and
#'   ignored. By default is 0.2 for impurity-corrected input (ISAR correction). We suggest using a
#'   smaller threshold (ex., 0.1) if the input is not corrected for impurity.
#' @param cores How many cores to use for processing chromosome arms in parallel. As always with the
#'  \code{parallel} library, use of multiple cores is not supported on Windows systems.
#' @export
#' @return Returns NULL if breakpoint files are successfully written. Otherwise, attempts to return
#'   error information.
make_breakpoint_files = function(segment_file = NULL, output_dir = NULL, 
                                 chromosome_coordinates = get_chromosome_coordinates(), 
                                 arms = get_chromosome_arms(),
                                 amplitude_threshold = .2,
                                 cores = 1) {
  
  
  if(! is.character(segment_file) || length(segment_file) != 1) {
    stop('segment_file should be 1-length character (specifying the copy number segementation file)')
  }
  if(! file.exists(segment_file)) {
    stop('Specified segment_file does not exist.')
  }
  segment_file = normalizePath(segment_file)
  
  
  # Validate choice of output directory
  if(! is.character(output_dir) || length(output_dir) != 1 || nchar(output_dir) == 0) {
    stop('output_dir should be 1-length character (specifically, a path where breakpoint files should be saved).')
  }
  if(! dir.exists(dirname(output_dir))) {
    stop('Parent directory for input output_dir (', dirname(output_dir), ' does not exist.')
  }
  if(dir.exists(output_dir)) {
    stop('Specified output_dir directory ', output_dir, ' already exists.')
  }
  
  # Validate amplitude_threshold
  if(! rlang::is_double(amplitude_threshold, finite = T)|| length(amplitude_threshold) != 1 || amplitude_threshold <= 0 || amplitude_threshold >= 1) {
    stop('amplitude_threshold must be a scalar numeric on the interval (0, 1).')
  }
  
  chromosome_coordinates = validate_chr_coordinates(chromosome_coordinates)
  
  if(! rlang::is_scalar_integerish(cores) || ! is.finite(cores) || cores < 1) {
    stop('cores should be positive integer')
  }
  max_cores = parallel::detectCores()
  if(cores > max_cores) {
    message('Running with ', max_cores, " detected cores.")
    cores = max_cores
  }
  
  arms = unique(arms)
  bad_arms = setdiff(arms, get_chromosome_arms())
  if(length(bad_arms) > 0) {
    stop('Unrecognized arms: ', paste0(bad_arms, collapse = ', '), '.')
  }
  
  if(! identical(arms, get_chromosome_arms())) {
    message('Running with a subset of chromosome arms. Adding "some_arms" to output directory as a reminder.')
    output_dir = gsub('/*$', '_some_arms', output_dir)
  }
  
  if(! dir.create(output_dir)) {
    stop('Could not create output_dir.')
  }
  
  # To avoid confusion, arranging to delete all output if this function exits without finishing.
  finished_processing = FALSE
  on.exit({
    if(! finished_processing) {
      unlink(output_dir, recursive = T)
      message('Creating breakpoint files did not complete. Deleted incomplete output directory.')
    }
  })
  output_dir = normalizePath(output_dir)
  output_prefix = basename(output_dir)
  
  reticulate::source_python(system.file('python/BISCUT_preprocessing.py', package = 'BISCUT'))
  
  # Get all combinations of amp/del and arms.
  run_args = expand.grid(aneu = c('amp', 'del'), arm = arms)
  
  # Advise regarding cores option
  if(cores == 1 && max_cores > 1) {
    message("Running with 1 core. Restart with more cores (argument \"cores = ...\") if you want to reduce runtime.")
  }
  
  # Run the preprocess_arm Python function
  message("Processing chromosome arms...")
  is_success = pbapply::pblapply(1:nrow(run_args),
                                  function(x) {
                                    preprocess_arm(arm = run_args$arm[x],
                                                   aneu = run_args$aneu[x],
                                                   segment_file = segment_file,
                                                   output_dir = output_dir,
                                                   output_prefix = output_prefix,
                                                   chr_info = chromosome_coordinates,
                                                   threshold = amplitude_threshold)
                                  }, cl = cores)
  which_failed = is_success[! sapply(is_success, rlang::is_true)]
  if(length(which_failed) > 0) {
    warning('For some reason, creation of some or all breakpoint files failed.\nDeleted incomplete output directory.')
    return(which_failed)
  }
  finished_processing = TRUE
  message('Breakpoint files created.')
}


