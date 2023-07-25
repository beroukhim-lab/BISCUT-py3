### Author: Juliann Shih, jshih@broadinstitute.org
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: July 24, 2023

### License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
### Dependencies: tested using R 4.1 and Python 3.9
### See README for guide on how to run this package

library(pastecs)
library(ismev)
library(extRemes)
library(ggplot2)
library(fitdistrplus)
library(GenomicRanges)
library(gridExtra)


plot_BISCUT_results <- function(thedate, ci, abslocs) {
  for(ttpath in list.dirs(paste('results',thedate,toString(ci),sep='_'),recursive = F)) {
    tt = strsplit(ttpath,'/')[[1]][2]
    print(tt)
    
    #del
    if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del.txt',sep=''))$size>0 ) {
      df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del.txt',sep=''),sep='\t',header=F,as.is=T)
      
      del <- ggplot() + theme_classic() + ggtitle('Deletions') + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                           limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                           breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                           labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
        geom_rect(data=df_del, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_del$V5, fill=df_del$V5) +
        geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
        geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
        geom_vline(xintercept= 0) +
        scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                           limits = c(-(max(abs(df_del$V3),abs(df_del$V4))+1), max(abs(df_del$V3),abs(df_del$V4))+1)) 
      
      ggsave(plot=del, paste(ttpath,'/summary/',tt,'_BISCUT_results_del.pdf',sep=''),height=11,width=8)
    } else {
      del <- ggplot() +theme_classic()
    }
    
    #amp
    if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp.txt',sep=''))$size>0 ) {
      df_amp <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp.txt',sep=''),sep='\t',header=F,as.is=T)
      
      amp <- ggplot() + theme_classic() + ggtitle('Amplifications') + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                           limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                           breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                           labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
        geom_rect(data=df_amp, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_amp$V5, fill=df_amp$V5) +
        geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
        geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
        geom_vline(xintercept= 0) +
        scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                           limits = c(-(max(abs(df_amp$V3),abs(df_amp$V4))+1), max(abs(df_amp$V3),abs(df_amp$V4))+1)) 
      
      ggsave(plot=amp, paste(ttpath,'/summary/',tt,'_BISCUT_results_amp.pdf',sep=''),height=11,width=8)
    } else {
      amp <- ggplot() +theme_classic()
    }
    
    #neg
    if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg.txt',sep=''))$size>0 ) {
      df_neg <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg.txt',sep=''),sep='\t',header=F,as.is=T)
      
      neg <- ggplot() + theme_classic() + ggtitle('Negative Selection') + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                           limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                           breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                           labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
        geom_rect(data=df_neg, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_neg$V5, fill=df_neg$V5) +
        geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
        geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
        geom_vline(xintercept= 0) +
        scale_x_continuous(name="Significance\nessential-like <----------------------------------------------> toxic-like",
                           limits = c(-(max(abs(df_neg$V3),abs(df_neg$V4))+1), max(abs(df_neg$V3),abs(df_neg$V4))+1)) 
      
      ggsave(plot=neg, paste(ttpath,'/summary/',tt,'_BISCUT_results_neg.pdf',sep=''),height=11,width=8)
    } else {
      neg <- ggplot() +theme_classic()
    }
    
    #pos
    if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos.txt',sep=''))$size>0 ) {
      df_pos <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos.txt',sep=''),sep='\t',header=F,as.is=T)
      
      pos <- ggplot() + theme_classic() + ggtitle('Positive Selection') + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                           limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                           breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                           labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
        geom_rect(data=df_pos, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_pos$V5, fill=df_pos$V5) +
        geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
        geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
        geom_vline(xintercept= 0) +
        scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                           limits = c(-(max(abs(df_pos$V3),abs(df_pos$V4))+1), max(abs(df_pos$V3),abs(df_pos$V4))+1)) 
      
      ggsave(plot=pos, paste(ttpath,'/summary/',tt,'_BISCUT_results_pos.pdf',sep=''),height=11,width=8)
    } else {
      pos <- ggplot() +theme_classic()
    }
    
    #pos
    if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_negpos.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_negpos.txt',sep=''))$size>0 ) {
      df_negpos <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_negpos.txt',sep=''),sep='\t',header=F,as.is=T)
      
      negpos <- ggplot() + theme_classic() + ggtitle('Positive+Negative Selection') + theme(plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                           limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                           breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                           labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
        geom_rect(data=df_negpos, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_negpos$V5, fill=df_negpos$V5) +
        geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
        geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
        geom_vline(xintercept= 0) +
        scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                           limits = c(-(max(abs(df_negpos$V3),abs(df_negpos$V4))+1), max(abs(df_negpos$V3),abs(df_negpos$V4))+1)) 
      
      ggsave(plot=negpos, paste(ttpath,'/summary/',tt,'_BISCUT_results_negpos.pdf',sep=''),height=11,width=8)
    } else {
      negpos <- ggplot() +theme_classic()
    }
    
    for(telcent in c('tel','cent','telcent')) {
      #dels
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_del_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        del <- ggplot() + theme_classic() + ggtitle('Deletions') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_rect(data=df_del, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_del$V5, fill=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                             limits = c(-(max(abs(df_del$V3),abs(df_del$V4))+1), max(abs(df_del$V3),abs(df_del$V4))+1)) 
        
        ggsave(plot=del, paste(ttpath,'/summary/',tt,'_BISCUT_results_del_',telcent,'.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_del_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_del_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_del_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        delss <- ggplot() + theme_classic() + ggtitle('Deletions') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000),
                            limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                            breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                            labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_segment(data=df_del, mapping=aes(x=V2, y=-V1, xend=V4, yend=-V3), color=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                             limits = c(-(max(abs(df_del$V2),abs(df_del$V4))+1), max(abs(df_del$V2),abs(df_del$V4))+1)) 
        
        ggsave(plot=delss, paste(ttpath,'/summary/',tt,'_BISCUT_results_del_',telcent,'_jagged.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      #dels
    
      #amps
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp_',telcent,'.txt',sep=''))$size>0) {
        df_amp <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_amp_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        amp <- ggplot() + theme_classic() + ggtitle('Amplifications') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_rect(data=df_amp, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_amp$V5, fill=df_amp$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                             limits = c(-(max(abs(df_amp$V3),abs(df_amp$V4))+1), max(abs(df_amp$V3),abs(df_amp$V4))+1))
        ggsave(plot=amp, paste(ttpath,'/summary/',tt,'_BISCUT_results_amp_',telcent,'.pdf',sep=''),height=11,width=8)
      } else {
        amp <- ggplot() + theme_classic()
      }
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_amp_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_amp_',telcent,'.txt',sep=''))$size>0 ) {
        df_amp <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_amp_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        ampss <- ggplot() + theme_classic() + ggtitle('Amplifications') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000),
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_segment(data=df_amp, mapping=aes(x=V2, y=-V1, xend=V4, yend=-V3), color=df_amp$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\nnegative selection <----------------------------------------------> positive selection",
                             limits = c(-(max(abs(df_amp$V2),abs(df_amp$V4))+1), max(abs(df_amp$V2),abs(df_amp$V4))+1)) 
        
        ggsave(plot=ampss, paste(ttpath,'/summary/',tt,'_BISCUT_results_amp_',telcent,'_jagged.pdf',sep=''),height=11,width=8)
      } else {
        amp <- ggplot() +theme_classic()
      }
      #amps
      
      #negs
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_neg_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        del <- ggplot() + theme_classic() + ggtitle('Negative Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_rect(data=df_del, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_del$V5, fill=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V3),abs(df_del$V4))+1), max(abs(df_del$V3),abs(df_del$V4))+1)) 
        
        ggsave(plot=del, paste(ttpath,'/summary/',tt,'_BISCUT_results_neg_',telcent,'.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_neg_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_neg_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_neg_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        delss <- ggplot() + theme_classic() + ggtitle('Positive Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000),
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_segment(data=df_del, mapping=aes(x=V2, y=-V1, xend=V4, yend=-V3), color=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V2),abs(df_del$V4))+1), max(abs(df_del$V2),abs(df_del$V4))+1)) 
        
        ggsave(plot=delss, paste(ttpath,'/summary/',tt,'_BISCUT_results_neg_',telcent,'_jagged.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      #negs
      
      #pos
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_pos_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        del <- ggplot() + theme_classic() + ggtitle('Positive Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_rect(data=df_del, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_del$V5, fill=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V3),abs(df_del$V4))+1), max(abs(df_del$V3),abs(df_del$V4))+1)) 
        
        ggsave(plot=del, paste(ttpath,'/summary/',tt,'_BISCUT_results_pos_',telcent,'.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_pos_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_pos_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_pos_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        delss <- ggplot() + theme_classic() + ggtitle('Positive Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000),
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_segment(data=df_del, mapping=aes(x=V2, y=-V1, xend=V4, yend=-V3), color=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V2),abs(df_del$V4))+1), max(abs(df_del$V2),abs(df_del$V4))+1)) 
        
        ggsave(plot=delss, paste(ttpath,'/summary/',tt,'_BISCUT_results_pos_',telcent,'_jagged.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      #pos
      
      #posneg
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_posneg_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_posneg_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_plotting_posneg_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        del <- ggplot() + theme_classic() + ggtitle('Positive + Negative Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000), 
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_rect(data=df_del, mapping=aes(xmin=V3, xmax=V4, ymin=-V2, ymax=-V1), color=df_del$V5, fill=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V3),abs(df_del$V4))+1), max(abs(df_del$V3),abs(df_del$V4))+1)) 
        
        ggsave(plot=del, paste(ttpath,'/summary/',tt,'_BISCUT_results_posneg_',telcent,'.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      if(file.exists(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_posneg_',telcent,'.txt',sep='')) && file.info(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_posneg_',telcent,'.txt',sep=''))$size>0 ) {
        df_del <- read.csv(paste(ttpath,'/summary/',tt,'_BISCUT_results_for_jagged_plotting_posneg_',telcent,'.txt',sep=''),sep='\t',header=F,as.is=T)
        
        delss <- ggplot() + theme_classic() + ggtitle('Positive + Negative Selection') + theme(plot.title = element_text(hjust = 0.5)) +
          scale_y_continuous(name="Chromosomes", expand = c(0.01,10000),
                             limits = c(-max(abslocs$q_end1[1:length(abslocs$q_end1)-1]),0),
                             breaks=c(-(abslocs$middle_of_arm[1:length(abslocs$middle_of_arm)-1])), 
                             labels = c(abslocs$chromosome_info[1:length(abslocs$chromosome_info)-1])) +
          geom_segment(data=df_del, mapping=aes(x=V2, y=-V1, xend=V4, yend=-V3), color=df_del$V5) +
          geom_hline(yintercept=c(0,-abslocs$q_end1[1:length(abslocs$q_end1)-1]),linetype='solid',color='#808080') +
          geom_hline(yintercept=c(0,-abslocs$q_start1[1:length(abslocs$q_start1)-1]),linetype='dotted',color='#808080') +
          geom_vline(xintercept= 0) +
          scale_x_continuous(name="Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like",
                             limits = c(-(max(abs(df_del$V2),abs(df_del$V4))+1), max(abs(df_del$V2),abs(df_del$V4))+1)) 
        
        ggsave(plot=delss, paste(ttpath,'/summary/',tt,'_BISCUT_results_posneg_',telcent,'_jagged.pdf',sep=''),height=11,width=8)
      } else {
        del <- ggplot() +theme_classic()
      }
      #posneg
    }
  }
}