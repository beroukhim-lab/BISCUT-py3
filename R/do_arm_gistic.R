#' Apply BISCUT method to one chromosome arm amp/del tel/cent
#' 
#' @keywords internal
do_arm_gistic <- function(arm, direc, telcent, mode, ci,qval_thres, telcent_thres, breakpoint_file_dir, n,
                          emp_bg, results_dir, abslocs, genelocs, seed = NULL) {
  if(! is.null(seed)) {
    set.seed(seed)
  }
  chromosome = as.integer(sub('[pq]$', '', arm))
  pq = substr(arm, nchar(arm), nchar(arm))
  pq = ifelse(pq %in% c('p', 'q'), pq, 'q') # acrocentromeric chromosomes lack p/q in arm name; always q
  
  # Using abslocs to define fractionize/unfractionize
  fractionize <- function(x, chromosome, pq, telcent) {
    p_arm_length <- (abslocs[chromosome,'p_end']-abslocs[chromosome,'p_start']+1)
    q_arm_length <- (abslocs[chromosome,'q_end']-abslocs[chromosome,'q_start']+1)
    if(pq == 'p' && telcent=='tel') {
      return ((x-abslocs[chromosome,'p_start']+1)/p_arm_length) 
    } else if (pq=='q' && telcent=='tel') {
      return (1-((x-abslocs[chromosome,'q_start']+1)/q_arm_length))
    } else if (pq=='q' &&telcent=='cent') {
      return ((x-abslocs[chromosome,'q_start']+1)/q_arm_length) 
    } else if (pq=='p' && telcent=='cent') {
      return (1-((x-abslocs[chromosome,'p_start']+1)/p_arm_length) )
    }
  }
  
  
  unfractionize <- function(y, chromosome, pq, telcent) {
    p_arm_length <- (abslocs[chromosome,'p_end']-abslocs[chromosome,'p_start']+1)
    q_arm_length <- (abslocs[chromosome,'q_end']-abslocs[chromosome,'q_start']+1)
    if(pq == 'p' && telcent=='tel') {
      return (round(y*p_arm_length-1+abslocs[chromosome,'p_start']))
    } else if (pq=='q' && telcent=='tel') {
      return(round(q_arm_length-1+abslocs[chromosome,'q_start']-(y*q_arm_length)))
      #return (1-((x-abslocs[chromosome,'q_start']+1)/q_arm_length))
    } else if (pq=='q' &&telcent=='cent') {
      return (round(y*q_arm_length-1+abslocs[chromosome,'q_start']))
    } else if (pq=='p' && telcent=='cent') {
      return (round(p_arm_length-1+abslocs[chromosome,'p_start']-(y*p_arm_length)))
    }
  }
  
  tablify <- function(df, chromosome, pq, telcent, remove_ones = F) {
    if((pq=='p'&&telcent=='tel') || (pq=='q' && telcent=='cent')) {
      tabs <- as.data.frame(table(unfractionize(df$percent,chromosome,pq,telcent)))
      rownames(tabs) <- tabs$Var1
      if(remove_ones) {tabs <- tabs[tabs$Freq>1,]}
    } else {
      tabs <- as.data.frame(table(unfractionize(df$percent,chromosome,pq,telcent)))
      rownames(tabs) <- tabs$Var1
      if(remove_ones) {tabs <- tabs[tabs$Freq>1,]}
      tabs <- tabs[order(-as.integer(tabs$Var1)),]
    }
    return(tabs)
  }
  
  curr_file_pattern = paste0('_', arm, '_', direc, '_', telcent, '.txt$')
  probefilename = list.files(path = breakpoint_file_dir, 
                             pattern = curr_file_pattern,
                             full.names = TRUE)
  if(length(probefilename) > 1) {
    stop('Unexpectedly found more than one file ending in ', curr_file_pattern, ' in breakpoint_file_dir.')
  } else if (length(probefilename) == 0) {
    stop('Could not find breakpoint file ending in ', curr_file_pattern, ' in breakpoint_file_dir.')
  }
  probes <- read.table(probefilename,sep='\t',header=T,stringsAsFactors = F,fill=T)
  probes <- filter_big_small(probes)
  tabprobes <- tablify(probes,chromosome,pq,telcent,F)
  probelist <- as.numeric(as.character(tabprobes$Var1))
  
  emp = emp_bg$emp
  betafit = emp_bg$fit
  alpha = summary(betafit)$estimate[1]
  beta = summary(betafit)$estimate[2]
  
  # Read in file again
  df <- read.table(probefilename,sep='\t',header=T,stringsAsFactors = F,fill=T)
  df$percent = as.numeric(df$percent)
  
  df1 <-df
  df1 <- filter_big_small(df1)
  
  ## if there's a coverage desert, spread the breakpoints until the next PROBE.
  tabs <- tablify(df1,chromosome,pq,telcent,T)
  
  newc <- df1$percent
  intnewc <- unfractionize(newc,chromosome,pq,telcent)
  
  for (k in rownames(tabs)) {
    z = as.integer(match(k,intnewc)) #index of first instance
    num_affected <- tabs[k,'Freq']
    if(num_affected+z > length(intnewc)) {
      newseq = seq(from=fractionize(intnewc[z],chromosome,pq,telcent), to=1-telcent_thres, length.out=num_affected+1)
    }
    else {
      if((pq=='p'&&telcent=='tel') || (pq=='q' && telcent=='cent')) {
        nextprobe <- probelist[match(intnewc[z],probelist)+1]
      } else{
        nextprobe <- probelist[match(intnewc[z],probelist)+1]
      }
      intnewseq = seq(from=intnewc[z], to=nextprobe,length.out=num_affected+1)
      newseq = fractionize(intnewseq, chromosome,pq,telcent)
    }
    # Previously, newseq[1:length(newseq-1)] gave warning due to off-by-one on indexing. However, actual behavior was the same.
    newc[z:(z+num_affected-1)] = newseq[1:(length(newseq)-1)] 
  }
  df1$percent <- newc
  
  x = qbeta(seq(0,1,length=length(rownames(df1))),alpha,beta)
  
  decide_selection <- function(dis) {
    if(abs(max(dis)) > abs(min(dis))) return('positive')
    else return('negative')
  }
  
  find_peaks = function(leftright, df1, lims, x, prior_peaks, prefix, iteration) {
    # Search for peak, put in peak list
    # If found:
    #   results_left = Search left (list)
    #   results_right = Search right (list)
    #   append results to peak list
    # Return peak list
    # 
    
    # Will stop early when insufficient samples
    if(length(rownames(df1))<4) {
      return(NULL)
    }
    df1$rownamez <- seq(1,length(rownames(df1)))
    lowlim = lims[1]
    highlim = lims[2]
    abs_dis = truncdist::ptrunc(df1$percent,spec='beta',alpha,beta, a=lowlim,b=highlim)*length(df1$percent) - df1$rownamez
    dis = abs_dis/length(df1$percent)
    df1$dis =dis
    selection <- decide_selection(dis)
    if(selection=='positive') {
      peak_index <- which(dis==max(dis))
    } else { #negative selection
      peak_index <- which(dis==min(dis))
    }
    
    if(is.null(leftright)) {
      parent_peak = NA
      prefix = substr(selection,1,1)
    } else {
      curr_prefix = paste0(substr(leftright,1,1),substr(selection,1,1),sep='')
      parent_peak = prefix
      prefix = paste0(prefix, '-', curr_prefix)
    }
    curemp <- emp[(emp<=highlim)&(emp>=lowlim)]
    
    # Get significance with ks.test; suppress the approximate p-value warning.
    theks <- withCallingHandlers(ks.test(df1$percent, curemp),
                                 warning = function(w) {
                                   if (grepl("p-value will be approximate in the presence of ties", conditionMessage(w))) {
                                     invokeRestart("muffleWarning")
                                   }})
    p_ks <- theks$p
    stat_ks <- theks$statistic[['D']]
    
    iteration_name <- paste(arm,direc,telcent,prefix,toString(iteration),sep='_')
    parent_name <- paste(arm,direc,telcent,parent_peak,toString(iteration - 1),sep='_')

    #find number of peak
    peak<- df1$percent[peak_index]
    
    #TOGGLE FOR LINE
    vector_of_interest <- df1$percent
    
    if(p_ks > qval_thres) {
      return(NULL)
    }
    
    make_right_dist <- function(rightdf, lims) {
      rightnumrows = length(rightdf$percent)
      
      rlims = c(min(rightdf$percent),lims[2])
      rightx = qbeta(seq(pbeta(rlims[1],alpha,beta),pbeta(rlims[2],alpha,beta),length=rightnumrows),alpha,beta)
      return(list(rlims,rightx))
    }
    
    make_left_dist <- function(leftdf, lims) {
      leftnumrows = length(leftdf$percent)
      llims = c(lims[1],max(leftdf$percent))
      leftx = qbeta(seq(pbeta(llims[1],alpha,beta),pbeta(llims[2],alpha,beta),length=leftnumrows),alpha,beta)
      return(list(llims,leftx))
    }
    
    calc_bounds <-function(p) {
      #number of events (telomere-bounded deletions)
      right_numevents <- length(df1[df1$percent>p,]$percent)
      left_numevents <- length(df1[df1$percent<p,]$percent)
      
      bound_helper <- function(num_events, low, high, min_or_max) {
        if(num_events == 0) {
          # Matches previous behavior of returning -Inf when there are no events
          return(quantile(-Inf,c((1-ci)/2, 1-((1-ci)/2))))
        }
        unif_boot <- runif(min=pbeta(low,alpha,beta), max = pbeta(high,alpha,beta), n = num_events * n)
        extremes <- numeric(num_events)
        for (i in seq_len(num_events)) {
          extremes[i] <- min_or_max(unif_boot[((i - 1) * n + 1):(i*n)])
        }
        gevdist = qbeta(extremes, alpha, beta)
        return(quantile(gevdist,c((1-ci)/2, 1-((1-ci)/2))))
      }
      return(list(bound_helper(num_events = right_numevents, low = p, high = highlim, min_or_max = min), 
                  bound_helper(num_events = left_numevents, low = lowlim, high = p, min_or_max = max)))
    }
    
    calc_left_boundary <- function(p,interval=1e-4) {
      left_boundary <- p
      right_high <- calc_bounds(left_boundary)[[1]][[2]]
      #if left_boundary is less than or eq to right_high, then drift left until it's not
      if (p<=right_high) {
        while (T) {
          #print(left_boundary)
          if (calc_bounds(left_boundary)[[1]][[2]] < p | left_boundary<=0) { #if goes out of bounds
            break
          }
          left_boundary <- left_boundary - interval
          
        }
      }
      
      return(left_boundary)
    }
    
    calc_right_boundary <- function(p,interval=1e-4) {
      right_boundary <- p
      left_high <- calc_bounds(right_boundary)[[2]][[1]]
      #if right_boundary is greater than or eq to left_high, then drift right until it's not
      if (p>=left_high) {
        while (T) {
          if (calc_bounds(right_boundary)[[2]][[1]] > p | right_boundary>=1) { #if goes out of bounds
            break
          }
          right_boundary <- right_boundary + interval
          
        }
      }
      
      return(right_boundary)
    }
    
    
    calc_both_boundaries <- function(peak) {
      if( (selection=='positive' && telcent == 'tel') || (selection=='negative' && telcent=='cent')) {
        #print(paste('original peak', peak,sep=' '))
        right_peak <- peak
        left_peak <- peak
        
        # calculate right-most peak
        left_bound_a <- calc_bounds(right_peak)[[2]]
        right_bound_a <- calc_bounds(right_peak)[[1]]
        #while not at ends, and not in between good ok bounds
        counter<- match(right_peak,vector_of_interest)
        
        while (counter!=1 &&
               counter!=length(vector_of_interest) && 
               vector_of_interest[counter-1]<left_bound_a[[2]]) {     #next one to the left is far enough, so move right
          temp_right_peak <- vector_of_interest[counter+1]
          left_bound_a <- calc_bounds(temp_right_peak)[[2]]
          if (right_peak >=left_bound_a[[2]]) {break}
          right_peak<- temp_right_peak
          counter <- counter+1
        }
        
        # calculate left-most peak
        left_bound_b <- calc_bounds(left_peak)[[2]]
        right_bound_b <- calc_bounds(left_peak)[[1]]
        counter <- match(left_peak,vector_of_interest)
        while (counter!=1 && 
               counter!=length(vector_of_interest) &&
               vector_of_interest[counter+1]<=right_bound_b[[2]]) {# next one to the right is close enough, so move on left
          temp_left_peak <- vector_of_interest[counter-1]
          right_bound_b<- calc_bounds(temp_left_peak)[[1]]
          if (left_peak>right_bound_b[[2]]) {break}
          left_peak <- temp_left_peak
          counter <- counter-1
        }
        
        left_boundary <- calc_left_boundary(left_peak)
        right_boundary <- right_peak
        
        
      } else { 
        
        #print(paste('original peak', peak,sep=' '))
        right_peak <- peak
        left_peak <-peak
        
        
        # calculate left-most peak
        left_bound_a <- calc_bounds(left_peak)[[2]]
        right_bound_a <- calc_bounds(left_peak)[[1]]
        counter<- match(left_peak,vector_of_interest)
        
        # NOT to be checking next one to the left is far enough, so move right; now it's
        # checking that the next one to the right is far enough, so move left.
        while (counter!=1 && 
               counter!=length(vector_of_interest) && 
               vector_of_interest[counter+1]>right_bound_a[[1]]) { #next one to the right is far enough, so "true" peak should be left.
          
          temp_left_peak <- vector_of_interest[counter-1]
          right_bound_a <- calc_bounds(temp_left_peak)[[1]]
          if (left_peak <=right_bound_a[[1]]) {break} #oriingal left peak is smaller than new lower right_bound of temp left peak
          left_peak<- temp_left_peak
          counter<-counter-1
        }
        
        # calculate right-most peak
        left_bound_b <- calc_bounds(right_peak)[[2]]
        right_bound_b <- calc_bounds(right_peak)[[1]]
        counter <- match(right_peak,vector_of_interest)
        while (counter!=1 && 
               counter!=length(vector_of_interest) && 
               vector_of_interest[counter-1]>=left_bound_b[[1]]) {# next one to the left is close enough, so "true" peak should be right.
          temp_right_peak <- vector_of_interest[counter+1]
          left_bound_b <- calc_bounds(temp_right_peak)[[2]]
          if (right_peak<left_bound_b[[1]]) {break} #original right peak is smaller than the new lower left bound of temp right peak
          right_peak <- temp_right_peak
          counter <- counter+1
        }
        
        right_boundary <- calc_right_boundary(right_peak)
        left_boundary <- left_peak
      }
      return(c(left_boundary,right_boundary))
    }
    
    boundaries <- calc_both_boundaries(peak)
    
    # putting restrictions here on left_boundary and right_boundary
    if(pq=='p') {
      if (telcent=='tel') { #p arm, simple, telomere, so 0.1-0.4 is exactly as is.
        left_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
      } else { #p arm, simple, centromere, so 0.1-0.4 is actually 0.9-0.6
        left_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
      }
      
    } else {
      if (telcent=='tel') { #q arm, simple, telomere, 0.1-0.4 is actually 0.9-0.6
        left_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
      } else { #q arm, simple, centromere, so 0.1-0.4 is exactly as is.
        left_boundary_coord <- unfractionize(boundaries[1],chromosome,pq,telcent)
        right_boundary_coord <- unfractionize(boundaries[2],chromosome,pq,telcent)
      }
    }
    
    if(right_boundary_coord<left_boundary_coord) {
      low <- right_boundary_coord
      high <- left_boundary_coord
      right_boundary_coord <- high
      left_boundary_coord <- low
    }
    
    if(is.na(right_boundary_coord) || is.na(left_boundary_coord)) {
      return(NULL)
    }
    
    #print(probelist)
    sorted_probelist<-probelist[order(probelist)]
    left_closest_index <- which(abs(sorted_probelist-left_boundary_coord)==min(abs(sorted_probelist-left_boundary_coord)))
    if(left_boundary_coord!= sorted_probelist[left_closest_index][1]) {
      left_boundary_coord <- sorted_probelist[left_closest_index-1][1]
    } else {
      left_boundary_coord <- sorted_probelist[left_closest_index][1]
    }
    
    right_closest_index <- which(abs(sorted_probelist-right_boundary_coord)==min(abs(sorted_probelist-right_boundary_coord)))
    #print(right_closest_index)
    if(right_boundary_coord!= sorted_probelist[right_closest_index][1]) {
      right_boundary_coord <- sorted_probelist[right_closest_index+1][1]
    } else {
      right_boundary_coord <- sorted_probelist[right_closest_index][1]
    }
    
    if((pq=='p'&&telcent=='tel')||(pq=='q'&&telcent=='cent')) {
      left_boundary <- fractionize(left_boundary_coord,chromosome,pq,telcent)
      right_boundary <- fractionize(right_boundary_coord,chromosome,pq,telcent)
      
    } else {
      left_boundary <- fractionize(right_boundary_coord,chromosome,pq,telcent)
      right_boundary <- fractionize(left_boundary_coord,chromosome,pq,telcent)
    }
    

    
    chrgenelocs <- genelocs[genelocs$Chr == chromosome, ]
    
    filteredgenelocs <- chrgenelocs[chrgenelocs$End >= left_boundary_coord,]
    filteredgenelocs <- filteredgenelocs[filteredgenelocs$Start <= right_boundary_coord,]
    if(length(rownames(filteredgenelocs))==0) {
      rownames(chrgenelocs) <- NULL
      leftsearch <- chrgenelocs$End
      leftdistance <- min(abs(leftsearch-left_boundary_coord))
      leftindex <- which(abs(leftsearch-left_boundary_coord)==leftdistance)
      rightsearch <- chrgenelocs$Start
      rightdistance <- min(abs(rightsearch-right_boundary_coord))
      rightindex <- which(abs(rightsearch-right_boundary_coord)==rightdistance)
      if(leftdistance<rightdistance) {
        filteredgenelocs <- chrgenelocs[leftindex,]
      } else if (rightdistance < leftdistance) {
        filteredgenelocs <- chrgenelocs[rightindex,]
      } else {filteredgenelocs <- chrgenelocs[c(rightindex,leftindex),]
      }
      filteredgenelocs$Gene <- paste('[',filteredgenelocs$Gene,']',sep='')
    }
    filteredgenelocs$Start.1 <- fractionize(filteredgenelocs$Start,chromosome, pq,telcent)
    filteredgenelocs$End.1 <- fractionize(filteredgenelocs$End,chromosome, pq,telcent)
    filteredgenelocs$ks_stat <- stat_ks
    filteredgenelocs$ks_p <- p_ks
    filteredgenelocs$log10_ks_p <- -1*log10(p_ks)
    filteredgenelocs$n_events <- length(rownames(df1))
    filteredgenelocs$Peak.Start <- floor(left_boundary_coord)
    filteredgenelocs$Peak.End <- ceiling(right_boundary_coord)
    filteredgenelocs$Peak.Start.1 <- left_boundary
    filteredgenelocs$Peak.End.1 <- right_boundary
    filteredgenelocs$arm <- arm
    filteredgenelocs$direction <- direc
    filteredgenelocs$telcent <- telcent
    filteredgenelocs$negpos <- substr(selection,1,1)
    filteredgenelocs$iter <- iteration
    filteredgenelocs$conf <- ci
    filteredgenelocs <- filteredgenelocs[order(filteredgenelocs$Start),] 
    
    #NEW 102119, REMOVE BIGGER THAN 0.5 PEAKS
    if(is.na(right_boundary) || is.na(left_boundary)) {
      return(NULL)
    }
    if (right_boundary-left_boundary > 0.5) {
      return(NULL)
    }
    
    #NEW 210311, stop if boundaries are out of range
    if ((right_boundary>= 1) || (left_boundary<=0)) {
      return(NULL)
    }
    
    check_overlaps <- function(prior_peaks, new_peak) {
      statusquo = FALSE
      for(i in prior_peaks) {
        statusquo <- max(new_peak[1],i[1]) <= min(new_peak[2],i[2])
        if(statusquo){break}
      }
      return(statusquo)
    }
    
    #080918; got toggled, left and right.
    if(mode=='overlap') {
      rightdf <- df1[vector_of_interest>left_boundary,]
      leftdf <- df1[vector_of_interest<right_boundary,]
      #print(prior_peaks)
      #print(left_boundary, right_boundary)
      if(check_overlaps(prior_peaks,c(left_boundary,right_boundary))) {
        return(NULL)
      }
    } else {
      rightdf <- df1[vector_of_interest>=right_boundary,]
      leftdf <- df1[vector_of_interest<=left_boundary,]
    }
    
    rightx <- make_right_dist(rightdf,lims)[[2]]
    rlims <- make_right_dist(rightdf,lims)[[1]]
    
    leftx <- make_left_dist(leftdf,lims)[[2]]
    llims <- make_left_dist(leftdf,lims)[[1]]
    
    filteredgenelocs$n_right <- length(rightx)
    filteredgenelocs$n_left <- length(leftx)
    filteredgenelocs$peak_id = iteration_name
    filteredgenelocs$parent = parent_name
    filteredgenelocs$code = prefix
    
    prior_peaks[[prefix]] <- c(left_boundary, right_boundary)
    
    # Produce peak plots
    df1$x = x
    peak_pos = df1$percent[peak_index]
    title <- paste(arm, direc,telcent,selection,prefix,toString(ci),'iteration',toString(iteration),sep=' ')
    title = paste0(title,'\nks p-value = ',toString(signif(p_ks, 3)))
    
    p1 = ggplot(df1, aes(x = percent, y = rownamez)) + geom_point(size = .75) + xlab('pSCNA Length') + ylab('Ranked Tumors') + ggtitle(title) +
      geom_point(aes(x = x, y = rownamez), shape = '.', size = 1) + 
      geom_vline(xintercept = peak_pos) + theme_classic()
    
    p2 = ggplot(data= df1, aes(x = percent, y = dis)) + geom_point(size = .75) +
      xlab('pSCNA Length') + ylab('Distance from background') + 
      scale_x_continuous(limits = c(lowlim, highlim)) + 
      geom_vline(xintercept = peak_pos) +
      geom_vline(xintercept = left_boundary, lty = 'dashed') + 
      geom_vline(xintercept = right_boundary, lty = 'dashed') + theme_classic()
    
    peak_plot_file = paste0(results_dir,'/peak_plots/', arm,'_',direc,'_',telcent,'_',prefix,'_',
                            toString(ci),'_iter',toString(iteration),'.pdf')
    peak_plots = cowplot::plot_grid(p1, p2, nrow = 2) 
    cowplot::save_plot(peak_plot_file, peak_plots)
    
    iteration = iteration + 1
    left_results = find_peaks(leftright = 'left', df1 = leftdf, lims = llims, x = leftx, prior_peaks = prior_peaks, prefix = prefix, iteration = iteration)
    right_results = find_peaks(leftright = 'right', df1 = rightdf,lims = rlims, x = rightx, prior_peaks = prior_peaks, prefix = prefix, iteration = iteration)
    
    
    return(c(left_results, list(as.data.table(filteredgenelocs)), right_results))
  }
  return(rbindlist(find_peaks(leftright = NULL, df1 = df1, lims = c(0, 1), x = x, prior_peaks = list(), prefix = NULL, iteration = 1)))
}
