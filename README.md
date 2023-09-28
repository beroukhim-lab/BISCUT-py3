[![DOI](https://zenodo.org/badge/558464969.svg)](https://zenodo.org/badge/latestdoi/558464969)

# BISCUT-py3
This is a (slightly) more user-friendly version of `BISCUT`,
**B**reakpoint **I**dentification of **S**ignificant **C**ancer **U**ndiscovered **T**argets. The method detects genomic loci under selective pressures by analyzing partial-SCNAs. The paper where we describe the method and the results of applying it to ~11,000 TCGA tumors:
 - [Shih, J., Sarmashghi, S., Zhakula-Kostadinova, N. et al. Cancer aneuploidies are shaped primarily by effects on tumour fitness. Nature (2023). https://doi.org/10.1038/s41586-023-06266-3][biscut-paper]

This version has been tested using `R 4.1.2` and `Python 3.9.7`. 

## How to run BISCUT
### Preparing telomeric/centromeric copy number segments
The first step is to run `BISCUT_preprocessing.py` on the input segmentation file. 
The columns of the input should be (in order)
```
Sample	Chromosome	Start	End	Num_Probes	Segment_Mean
```
 - This is what a SNP-array based (relative) copy-number segmentation file looks like. `Segment_Mean` is assumed
   to be in log2ratio scale. If you have your data in absolute allelic copy number format, to obtain `Segment_Mean`
   for each segment, simply add the allelic copy numbers and then divide that by the ploidy of the sample and take log2.
   This way the values will be centered aournd zero; positive values indicate copy number > ploidy (amps) and negative values
   indicate copy number < ploidy (dels). 
 - If the number of probes per segment does not apply to your input data (ex.,WGS-based copy numbers),
   set them all to `10` (anything larger than the current threshold which is `4`)
 
There are a few variables that are hard-coded in the beginning of this file.

 - **amplitude_threshold** - Copy number changes less than this threshold are considered noise and ignored. By default is `0.2` for impurity-corrected input (ISAR correction). We suggest using a smaller threshold (ex., `0.1`) if the input is not corrected for impurity.
 - **chromosome_coordinates** - The coordinates of p and q arms of the human chromosomes. By default, hg19-based coordinates provided under `\docs` directory are used. Please note these coordinates are not some clean-cut biologically well-defined regions but rather what appeared to be the boundaries of regions that can be considered telomeric and (particularly) centromeric from SNP-array based TCGA copy number profiles. You might want to use your custom set of coordinates based on your copy number profiles (ex. IGV view of start/end of event near telomere/centromere). From our experience, choosing the telomeric/centromeric boundaries wide enough so that the events are not missed due to imperfections of assaying technologies is more important than genome reference differences (ex., hg19 vs. hg38).
 - **tumor_type** - It's used for naming input and outputs. Together with the next variable they set the name of the input seg file.
 - **seg_file_suffix** - Any part of the input seg file name after the tumor type.
 - **n_proc** - Number of processors to use for multi-processing over chromosome arms.
 - **date_suffix** - It is used throughout the pipeline to uniquely name outputs.

### Finding peaks and summarizing the results
After running the pre-processing step and preparing telomeric/centromeric segments, it's time to run `BISCUT_peak_finding.R` which runs functions from multiple scripts including `Python` functions and for that you need `reticulate` package. Similar to the first step, there are a number of variables in the beginning that need to be set.

 - **date_suffix** - Should be identical to one used above for the pre-processing.
 - **tumor_type** - Should be identical to one used above for the pre-processing.
 - **ci** - Confidence interval threshold around peaks. Higher threshold results in wider intervals around peaks. We use `0.95` as default.
 - **qval_thres** - The threshold on q-values computed for the significance of 
 identified peaks after correcting for multiple hypothesis testing. We use `0.05` as default.
 - **telcent_thres** - Threshold to filter out very short or whole-arm events. It's a complimentary way besides the `chromosome coordiantes` to correct for noise in calling the copy number breakpoints. We use `0.001`, but you can increase it if you see peaks called too close to telomere/centromere that you think are not real.
  - **USE_BACKGROUND** - This is a switch where if set `TRUE`, pre-computed background length distributions (specified by `tel_background_file` and `cent_background_file`) are used for peak calling. By default is `FALSE` which means a new background is computed from the current data. Please refer to the section on per-lineage analysis below for more info on this.
  - **tel_background_file** - The pre-computed background distribution of telomeric events. If `USE_BACKGROUND=TRUE`, BISCUT uses the background from all TCGA tumors (provided under `\docs`) unless you set this variable to a different background. See section on per-lineage analysis for more info.
  - **cent_background_file** - Same as `tel_background_file` except it's for the centromeric events. 
  - **n_cores** - Number of cpu cores to use for multi-processing peak finding over chromosome arms.
  - **genelocs_file** - The gene locations (hg19-based) to generate list of genes per peaks in the output
  - **abslocs_file** - The coordinates of p and q arms of human chromosomes (hg19). Same files as `chromosome coordiantes` used in pre-processing.
  

### How lineage-specific analysis should be done
Originally BISCUT code worked with a list of tumor types and iterated over them. However, it was confusing how to select/calculate the right background distribution. In this version, we have made the whole pipeline (pre-processing and peak-finding) to run on a single lineage or pan-cancer (specified by `tumor_type`).

 - If you would like to do a **pan-cancer analysis**, you just need to provide all CN segments from all samples in a single seg file as input to the pre-processing step, and then run the peak finding using either the TCGA background (`USE_BACKGROUND=TRUE`), or a new background generated from your own data (`USE_BACKGROUND=FALSE`).
  - If you want to to do a **per-lineage analysis**, pick a lineage one at a time, subset your seg file to only those samples, run the pipeline on that seg file using either a pan-cancer or lineage-specific background. For the pan-cancer background, set `USE_BACKGROUND=TRUE` and set `tel_background_file` and `cent_background_file` accordingly based on whether you want to use TCGA background (provided under `/docs`) or the background you generated in the pan-cancer analysis. If you have enough samples and you want to use a lineage-specific background, set `USE_BACKGROUND=FALSE`. In any case, move to the next lineage and repeat the process. This can be streamlined by for example a bash script where you loop over tumor types and provide `tumor_type` as input to `BISCUT_preprocessing.py` and `BISCUT_peak_finding.R`.

### Tutorial
Here we provide the data and steps to generate the results for TCGA pancancer analysis described in the BISCUT manuscript. You can download the segmented copy number data that has been corrected for sample impurity from this [link][pancan-data]. Next, run `BISCUT_preprocessing.py` on that data to generate `breakpoint_files_*date*` directory which contains the breakpoints for telomeric and centromeric copy number alterations across all chromosome arms. Once you have successfully generated the breakpoint files, use them as input to run `BISCUT_peak_finding.R` to identify loci under selection. This will generate a directory `results_*date*` which containes all the results. There are different parameters used during the computation that are descibed above but if you use the default hard-coded parameters, you will be able to reproduce the manuscript results (with slight differences due to the random seed). In particular, `results_*date*/all_BISCUT_results.txt` formatted to provide all genes in loci under selection where each row is a gene, and `results_*date*/all_cols/PANCAN_BISCUT_results_cols_*ci*.txt` formatted so that each column is a locus and genes are listed underneath. You can also look under `results_*date*/PANCAN` for plots of length distribution for each arm, alteration direction, and event type (telomeric or centromeric) compared with the background to review the indentified loci.


[pancan-data]: https://drive.google.com/drive/folders/1-VZ_A0uodEs04Jg-Gphkl3HU5VdEyICW?usp=share_link
[biscut-paper]: https://www.nature.com/articles/s41586-023-06266-3
