#!/usr/bin/env python

### Author: Juliann Shih, jshih@broadinstitute.org
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: July 24, 2023

### License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
### Dependencies: tested using R 4.1 and Python 3.9
### See README for guide on how to run this package

__author__ = 'jshih'

import pandas as pd
import os
import multiprocessing as mp


amplitude_threshold = 0.2
chromosome_coordinates = 'docs/SNP6_hg19_chromosome_locs_200605.txt'
tumor_type = 'PANCAN'
seg_file_suffix = '_ISAR.seg.txt'
n_proc = 8 
date_suffix = '2022_07_26'


#for given length n that is 2 or greater, give list of tuples between each element
def junctions(n):
    res = []
    for i in range(0,n-1):
        res.append((i,i+1))
    return res

def do_arm(arm, tumortype, aneu, date, threshold=amplitude_threshold):
    info = pd.read_csv(chromosome_coordinates,sep='\t',index_col='chromosome_info').transpose().to_dict()
    if str(arm) in ['13','14','15','21','22']:
        coord = (info[int(arm)]['q_start'], info[int(arm)]['q_end'])
    elif arm.endswith('q'):
        coord = (info[int(arm[:-1])]['q_start'], info[int(arm[:-1])]['q_end'])
    elif arm.endswith('p'):
        coord = (info[int(arm[:-1])]['p_start'], info[int(arm[:-1])]['p_end'])
    else:
        coord = (info[int(arm)]['p_start'], info[int(arm)]['p_end'], info[int(arm)]['q_start'], info[int(arm)]['q_end'])
    coord=[int(i) for i in coord]

    psedftel = get_percent_start_end(aneu, arm, coord, 'docs/' + tumortype + seg_file_suffix, threshold, 'tel', False)
    psedfcent = get_percent_start_end(aneu, arm, coord, 'docs/' + tumortype + seg_file_suffix, threshold, 'cent', False)

    psedftel.to_csv('breakpoint_files_' + date + '/' + tumortype + '/' + tumortype + '_' + arm + '_' + aneu +'_tel'+ '.txt', sep='\t')
    psedfcent.to_csv('breakpoint_files_' + date + '/' + tumortype + '/' + tumortype + '_' + arm + '_' + aneu +'_cent'+ '.txt', sep='\t')


def preprocess_seg(arm, coord,segdfloc,js=4):
    segdf = pd.read_csv(segdfloc, sep='\t', index_col='Sample')
    segdf = segdf[segdf['Num_Probes']>=js]
    if arm.endswith('p'):
        chr = int(arm[:-1])
        segdf = segdf[segdf['Chromosome'] == chr]
        segdf = segdf[segdf['Start'] < coord[1]]  # segment starts in 8pmatch(left_peak,df1$end)
        segdf['End'] = segdf['End'].clip(upper = coord[1])  # sets 'End' to the limit of 8p if the segment goes thru centromere
        segdf['Start'] = segdf['Start'].clip(lower = coord[0])
        segdf = segdf[(segdf['Start'] <= coord[1]) & (segdf['End'] <= coord[1])]
        segdf = segdf[(segdf['Start'] >= coord[0]) & (segdf['End'] >= coord[0])]
        return segdf
    elif arm.endswith('q') or str(arm) in ['13', '14', '15', '21', '22']:  # if acrocentric or ends with q
        chr = int(arm) if str(arm) in ['13', '14', '15', '21', '22'] else int(arm[:-1])
        segdf = segdf[segdf['Chromosome'] == chr]
        segdf = segdf[segdf['End'] > coord[0]]  # segment ends in q arm
        segdf['Start'] = segdf['Start'].clip(lower = coord[0])  # sets 'Start' to the limit of q arm if the segment goes thru centromere
        segdf['End'] = segdf['End'].clip(upper = coord[1])
        segdf = segdf[(segdf['Start'] <= coord[1]) & (segdf['End'] <= coord[1])]
        segdf = segdf[(segdf['Start'] >= coord[0]) & (segdf['End'] >= coord[0])]
        return segdf
    else:
        return pd.concat([preprocess_seg(arm + 'p', [coord[0], coord[1]],segdfloc), preprocess_seg(arm + 'q', [coord[2], coord[3]],segdfloc)])

def get_percent_start_end(aneu, arm, coord, segdfloc,ncutoff,telcent,integer):
    if aneu=='amp': direction=1
    elif aneu =='del': direction= -1
    elif aneu=='non-del': direction= -2
    elif aneu=='non-amp': direction =2
    len_arm = float(coord[1]-coord[0]+1) if len(coord)==2 else float(coord[1]-coord[0]+1)+float(coord[3]-coord[2]+1)
    percent_start_end = {} # will be sample : [%del, loc start, loc end] s.t. loc start to loc end contains no more than nopercent% non-deletion
    segdf = preprocess_seg(arm, coord,segdfloc)
    allthesamples = list(set(segdf.index))
    for sample in allthesamples:
        miniseg = segdf[segdf.index==sample].reset_index()
        segments_in_arm=[] #all segments in the arm

        for i in miniseg.index:
            if aneu=='del':
                if miniseg.loc[i,'Segment_Mean']<=(-1)*ncutoff :
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],direction))
                else:
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],0))

            elif aneu=='amp':
                if miniseg.loc[i,'Segment_Mean']>=ncutoff:
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],direction))
                else:
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],0))

            elif aneu =='non-amp':
                if miniseg.loc[i,'Segment_Mean']<ncutoff :
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],direction))
                else:
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],0))
            elif aneu =='non-del':
                if miniseg.loc[i,'Segment_Mean']>(-1)*ncutoff :
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],direction))
                else:
                    segments_in_arm.append((miniseg.loc[i,'Start'],miniseg.loc[i,'End'],0))


        alt_segments_in_arm= list(filter(lambda x:x[2]==direction,segments_in_arm)) #alt segments
        # no segments, or doesn't start in telomere/centromere
        if len(alt_segments_in_arm)==0:
            percent_start_end[sample] = [0, 0, 0]
        elif arm.endswith('p') and telcent == 'tel' and alt_segments_in_arm[0][0]!=coord[0]: #looking for tel
            percent_start_end[sample] = [0, 0, 0]
        elif arm.endswith('p') and telcent =='cent' and alt_segments_in_arm[-1][1]!=coord[1]: #p, doing centromere, needs to end at centromere
            percent_start_end[sample] = [0, 0, 0]
        elif (not arm.endswith('p')) and telcent == 'tel' and alt_segments_in_arm[-1][1]!=coord[1]:
            percent_start_end[sample] = [0, 0, 0]
        elif (not arm.endswith('p')) and telcent == 'cent' and alt_segments_in_arm[0][0]!=coord[0]:
          #      ((arm.endswith('q') or arm in ['13','14','15','21','22']) and alt_segments_in_arm[-1][1]!=coord[1]): #no seg
            percent_start_end[sample]=[0,0, 0]
        elif len(alt_segments_in_arm)==1: #uses integers for now
            percent_start_end[sample]=[(alt_segments_in_arm[0][1]-alt_segments_in_arm[0][0]+1), alt_segments_in_arm[0][0], alt_segments_in_arm[0][1]]
        else: #2 or more segments
            joined_segments_in_arm = [] #joined all of thresholded level
            c = segments_in_arm[0][2] #starting c (-1, 0, or 1)
            begin = segments_in_arm[0][0]
            end = segments_in_arm[0][1]
            for seg in segments_in_arm[1:-1]:
                if seg[2]==c: #same direction as last segment
                    end = seg[1]
                else: #if there's a change
                    joined_segments_in_arm.append((begin,end,c)) #append the last several segments, joined together
                    c = seg[2]
                    begin = seg[0]
                    end = seg[1]
            if segments_in_arm[-1][2]==c:
                joined_segments_in_arm.append((begin,segments_in_arm[-1][1],c))
            else:
                joined_segments_in_arm.append((begin,end,c))
                joined_segments_in_arm.append(segments_in_arm[-1])
            #joined_alt_segments_in_arm: first to last alt seg, not including middle stuff
            joined_alt_segments_in_arm = [i for i in joined_segments_in_arm if i[2]==direction]
            #adding fourth value in tuple for amount of segment that is already non-aneuploid
            #print sample
            joined_alt_segments =[(i[0],i[1],i[2],0) for i in joined_alt_segments_in_arm]
            #print segments_in_arm
            #print joined_alt_segments
            percent_start_end[sample] = join_segs(joined_alt_segments, arm, telcent)

        #print sample, percent_start_end[sample]

    psedf =  pd.DataFrame(percent_start_end).transpose().rename(columns={0:'percent',1:'start',2:'end'})


    #UNTOGGLE THIS!!!
    if not integer:
        psedf['start']=psedf['start'].apply(lambda x: min(1,max(0.0,(x-coord[0])/len_arm)))
        psedf['end']=psedf['end'].apply(lambda x: max(0,min(1.0,(x-coord[0])/len_arm)))
        psedf['percent']=psedf.index.map(lambda x: (psedf.loc[x,'percent']-(coord[2]-coord[1]-1))/len_arm if psedf.loc[x,'start']<=coord[1] and len(coord)==4 and psedf.loc[x,'end']>=coord[2] else psedf.loc[x,'percent']/len_arm)

        psedf['percent']=psedf['percent'].apply(lambda x: 0.0 if x<0.0 else 1.0 if x>1.0 else x)

        if not ((arm.endswith('p') and telcent =='tel') or ((not arm.endswith('p')) and telcent=='cent')):
            temp = psedf['start']
            newstart = 1-psedf['end']
            newend = 1-temp
            psedf = psedf.drop(['start','end'],axis=1)
            psedf['start'] = newstart
            psedf['end'] = newend


    #UNTOGGLE THIS!!!

    psedf.index.name = 'Sample'
    return psedf

def join_segs(joined_alt_segments, arm, telcent):
    if (arm.endswith('p') and telcent =='tel') or ((not arm.endswith('p')) and telcent=='cent'):
        # DONE FOR P ARMS NOW
        #if I can do this in one script and not recursively that would be great
        segments = joined_alt_segments

        def calc_lens(segments, finished=False):
            if finished and joined_alt_segments[0][1]==joined_alt_segments[0][0]:
                segments.insert(0,joined_alt_segments[0])
            if len(segments)==1:
                altlen=segments[0][1]-segments[0][0]+1
                notlen=0
            else:
                notlen = 0
                for i in range(1, len(segments)):
                    notlen = notlen + (segments[i][0] - segments[i - 1][1] - 1)
                altlen = sum([i[1] - i[0] + 1 for i in segments[1:]])
            return notlen, altlen

        notlen, altlen= calc_lens(segments)

        while True:
            if notlen <= altlen:
                notlen, altlen = calc_lens(segments, True)
                break
            else: #if there's too much white space
                segments = segments[:-1]
                notlen, altlen = calc_lens(segments, False)
        return [segments[-1][1] - segments[0][0] + 1, segments[0][0], segments[-1][1]]

    else:
        segments = joined_alt_segments
        def calc_lens(segments, finished=False):
            if finished and joined_alt_segments[-1][1] == joined_alt_segments[-1][0]:
                segments.append(joined_alt_segments[-1])

            if len(segments) == 1:
                altlen = segments[0][1] - segments[0][0] + 1
                notlen = 0
            else:
                notlen = 0
                for i in range(0, len(segments)-1):
                    notlen = notlen + (segments[i+1][0] - segments[i][1] - 1)
                altlen = sum([i[1] - i[0] + 1 for i in segments[:-1]])
            return notlen, altlen

        notlen, altlen = calc_lens(segments)
        while True:
            if notlen <= altlen:
                notlen, altlen = calc_lens(segments, True)
                break
            else:  # if there's too much white space
                segments = segments[1:]
                notlen, altlen = calc_lens(segments, False)
        return [segments[-1][1] - segments[0][0] + 1, segments[0][0], segments[-1][1]]




def take_care_arms(tt,date):
    pool = mp.Pool(n_proc)
    arms_list = ['13','14','15','21','22'] + [str(i)+'q' for i in range(1,23) if i not in [13,14,15,21,22]] + [str(i)+'p' for i in range(1,23) if i not in [13,14,15,21,22]]
    if not os.path.isdir('breakpoint_files_'+date): os.mkdir('breakpoint_files_'+date)
    if not os.path.isdir('breakpoint_files_'+date+'/'+ tt): os.mkdir('breakpoint_files_'+date+'/'+tt)
    results = [pool.apply_async(do_arm, args=(arm,tt,'amp',date)) for arm in arms_list] + [pool.apply_async(do_arm, args=(arm,tt,'del',date)) for arm in arms_list]
    for result in results:
      result.get(9999999)
    pool.close()
    pool.join()

if __name__=='__main__':
    take_care_arms(tumor_type, date_suffix)

