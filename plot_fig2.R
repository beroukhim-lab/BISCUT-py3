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
library(gridExtra)
library(stringr)
library(gtable)
library(cowplot)

plot_fig2 <- function(thedate, ci, abslocs) {
  dir.create(file.path(paste('results',thedate,toString(ci),sep='_'),'arms/plots'),showWarnings=T)
  dir.create(file.path(paste('results',thedate,toString(ci),sep='_'),'genes/plots'),showWarnings=T)
  armfiles <- list.files(file.path(paste('results',thedate,toString(ci),sep='_'),'arms','files'),pattern='*fig2.txt', full.names=T, recursive= T)
  genefiles <- list.files(file.path(paste('results',thedate,toString(ci),sep='_'),'genes','files'),pattern='*fig2.txt', full.names=T)
  for(ttpath in c(armfiles,genefiles)) {
    if(file.info(ttpath)$size<=1) {next('empty df')} else {
      results <- read.csv(ttpath,sep='\t',header=T,as.is=T,na.strings='NA')
    }
    if(nrow(results)== 0) {next('empty df')}
    print(ttpath)
    print(results)
    results$thing <- paste(results$type, results$telcent, results$code, sep=' ')
    print(results$thing)
    # this if statement is for arms, not genes.
    # if(strsplit(ttpath,'_')[[1]][4]=='amp' | strsplit(ttpath,'_')[[1]][4]=='del') {
       
    if(ttpath %in% armfiles) {
      arm = strsplit(basename(ttpath),'_')[[1]][1]
      dir.create(file.path(paste('results',thedate,toString(ci),sep='_'),'arms/plots/', arm),showWarnings=T)
      if (substr(arm,nchar(arm),nchar(arm))=='p') {
        armno = as.numeric(substr(arm,1,str_length(arm)-1))
        st = abslocs[armno,'p_start']
        en = abslocs[armno,'p_end']
        p1 <- ggplot() + theme_classic() + ggtitle('') + theme(plot.title = element_text(hjust = 0.5), axis.title.y=element_blank(), 
                                                               axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
                                                               
                                                               plot.margin=margin(5,10,5,50)) +
          geom_rect(data=results, mapping=aes(xmin=Peak.Start, xmax=Peak.End, ymin=ymin, ymax=ymax), fill=results$colors) +
          geom_text(aes(label=thing), data=results, hjust=1, y=results$ymin+0.45, x=0)+
          xlim(st,en) + 
          xlab(' ') + scale_y_continuous(expand = c(0,0))+
          coord_cartesian(clip = 'off')
      } else {
        if(substr(arm,nchar(arm),nchar(arm))=='q') {armno = as.numeric(substr(arm,1,str_length(arm)-1))} else {armno = as.numeric(arm)}
        st = abslocs[armno, 'q_start']
        en = abslocs[armno,'q_end']
        p1 <- ggplot() + theme_classic() + ggtitle('') + theme(plot.title = element_text(hjust = 0.5), axis.title.y=element_blank(), 
                                                               axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
                                                               
                                                               plot.margin=margin(5,10,5,50)) +
          geom_rect(data=results, mapping=aes(xmin=Peak.Start, xmax=Peak.End, ymin=ymin, ymax=ymax), fill=results$colors) +
          geom_text(aes(label=thing), data=results, hjust=1, y=results$ymin+0.45, x=st-7e6)+
          xlim(st,en) + 
          xlab(' ') + scale_y_continuous(expand = c(0,0))+
          coord_cartesian(clip = 'off')
      }
      #print(arm)

     } else{
       p1 <- ggplot() + theme_classic() + ggtitle('') + theme(plot.title = element_text(hjust = 0.5), axis.title.y=element_blank(), 
                                                              axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
                                                              
                                                              plot.margin=margin(5,10,5,50)) +
         geom_rect(data=results, mapping=aes(xmin=Peak.Start, xmax=Peak.End, ymin=ymin, ymax=ymax), fill=results$colors) +
         geom_text(aes(label=thing), data=results, hjust=1, y=results$ymin+0.45, x=min(results$Peak.Start)- (max(results$Peak.End)-min(results$Peak.Start))*0.02)+
         #
         xlab(' ') + scale_y_continuous(expand = c(0,0))+
         coord_cartesian(clip = 'off')
     }

    p2 <- ggplot() + theme_classic() + ggtitle('Significance Score') + theme(plot.title = element_text(hjust = 0),axis.title.y=element_blank(), 
                                                                                                        axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
                                                                                                        axis.line.y=element_blank(), 
                                                                                                        plot.margin=margin(5,20,5,10)) +
      geom_rect(data=results, mapping=aes(xmin=0, xmax=combined_sig, ymin=ymin, ymax=ymax), fill=results$colors) + 
      xlab('-log10 BY FDR * KS stat') +
      #coord_cartesian(clip = 'off') +
      scale_y_continuous(expand = c(0,0))
      #geom_text(aes(label=type),data=results, x=max(results$combined_sig, na.rm=T) + 0.1*max(results$combined_sig, na.rm=T), y=results$ymin+0.45, hjust=0) +
    
    resize_widths <- function(g, widths = rep(1, length(idpanels))){
      idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
      g$widths <- grid:::unit.list(g$widths)
      g$widths[idpanels] <- unit.c(do.call(unit, list(widths, 'null')))
      g
    }
    
    g1 <- ggplotGrob(p1)
    g2 <- ggplotGrob(p2)
    g0 <- plot_grid(p1,p2,rel_widths=c(5,1.5),align='h')
    #g0 <- resize_widths(g0, c(3,1))
    #print(g0)
    #p0 <- arrangeGrob(p1, p2, nrow = 1, widths=c(5,1.5))
    #grid.arrange(p0)
    #if(strsplit(ttpath,'_')[[1]][4]=='amp' | strsplit(ttpath,'_')[[1]][4]=='del') {
    if(ttpath %in% armfiles) {
      ggsave(plot=g0, file.path(paste('results',thedate,toString(ci),sep='_'),'arms','plots',arm,str_replace(basename(ttpath),'txt','pdf')),height=4,width=11)
      
    } else {
      ggsave(plot=g0, file.path(paste('results',thedate,toString(ci),sep='_'),'genes','plots',str_replace(basename(ttpath),'txt','pdf')),height=4,width=11)
    }
  }
}

