# Run interactively from BISCUT project directory
default_genelocs_file = system.file('extdata/hg19_genes_refseq_telcentfiltered_020518.txt', package = 'BISCUT')
default_genelocs = read.csv(default_genelocs_file,sep='\t',header=T)

default_abslocs_file = system.file('extdata/SNP6_hg19_chromosome_locs_200605.txt', package = 'BISCUT')
default_abslocs = read.table(default_abslocs_file, sep='\t', header=T)

# Prepare default background data. We will not filter out the (0, .001) and (.999, 1) events
# because the user can control that filtering via telcent_thres.
tel_background_file = system.file('extdata/PANCAN_tels_200605.txt.gz', package = 'BISCUT')
cent_background_file = system.file('extdata/PANCAN_cents_200605.txt.gz', package = 'BISCUT')
tel <- read.csv(tel_background_file, sep='\t')
cent <- read.csv(cent_background_file,sep='\t')
default_background_data = list(tel = tel, cent = cent)

acro_chr_arms = c(13, 14, 15, 21, 22)
other_chr = setdiff(1:22, acro_chr_arms)
other_chr_arms = c(paste0(other_chr, 'p'), paste0(other_chr, 'q'))
all_arms = stringr::str_sort(c(acro_chr_arms, other_chr_arms), numeric = TRUE)


usethis::use_data(default_genelocs, default_abslocs, default_background_data,
                  all_arms, internal = TRUE, overwrite = TRUE)


