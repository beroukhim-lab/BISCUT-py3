### Author: Juliann Shih, jshih@broadinstitute.org
### Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
### Date last updated: July 24, 2023

### License: GNU GPL2, Copyright (C) 2023 Dana-Farber Cancer Institute
### Dependencies: tested using R 4.1 and Python 3.9
### See README for guide on how to run this package

import pandas as pd
import numpy as np
import os
from operator import itemgetter
import re
import itertools
import fnmatch

darkred = '#a50f15'
lightred = '#fcae91'
lightblue = '#6baed6'
darkblue = '#08519c'

# returns the overlap
def overlap_helper(a, b):
    return [max(a[0], b[0]), min(a[1], b[1])] if min(a[1], b[1]) - max(a[0], b[0]) + 1 > 0 else False

def overlap_helper_simple(a, b):
    return True if min(a[1], b[1]) - max(a[0], b[0]) + 1 > 0 else False

def calc_overlaps(folder, genelocs):
    df = pd.read_csv(folder+'/all_BISCUT_results.txt', sep='\t')
    armgroups = df.groupby('arm')
    #typegroups = df.groupby('type')
    s_list = []
    # for tt, ttdf in typegroups:
    #     armgroups = ttdf.groupby('arm')
    for arm, armdf in armgroups:  # specific arm within specific tumor type
        # if arm == '1p':
        armdf = armdf.drop_duplicates(['Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'code'])
        #print armdf
        for dir, telcent, negpos in itertools.product(['amp', 'del'], ['tel', 'cent'], ['n', 'p']):
            # within is everything that shares same path
            within = armdf[(armdf.direction == dir) & (armdf.telcent == telcent) & (armdf.negpos == negpos)]
            without = armdf[(armdf.direction != dir) | (armdf.telcent != telcent) | (armdf.negpos != negpos)]
            within = within.filter(
                ['Chr', 'arm', 'Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'iter', 'code',
                 'ksby', 'combined_sig'])
            without = without.filter(
                ['Chr', 'arm', 'Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'iter', 'code',
                 'ksby', 'combined_sig'])
            #print dir, telcent, negpos
            #print len(within)
            #print len(without)
            if not within.empty:
                for i in within.index:
                    #print type1
                    newwithout = without
                    #print 'within',within
                    #print 'difftypewithin',difftypewithin
                    #print 'without',newwithout
                    for j in newwithout.index:
                        overlap = overlap_helper([within.loc[i, 'Peak.Start'], within.loc[i, 'Peak.End']],
                                                 [newwithout.loc[j, 'Peak.Start'], newwithout.loc[j, 'Peak.End']])
                        if overlap != False:
                            # print within.loc[i], without.loc[j], overlap
                            irow = within.loc[i]
                            jrow = newwithout.loc[j]
                            #print overlap
                            chrgenelocs = genelocs[genelocs['Chr'] == irow.Chr]
                            peakgenelocs = chrgenelocs[chrgenelocs.Start <= overlap[1]]
                            peakgenelocs = peakgenelocs[peakgenelocs.End >= overlap[0]]
                            peakgenelocs = peakgenelocs.sort_values('Start') #steph edit
                            genes = peakgenelocs['Gene'].tolist()
                            consistent = True if (irow.direction == jrow.direction and irow.negpos == jrow.negpos) or (irow.direction != jrow.direction and irow.negpos != jrow.negpos) else False
                            s = pd.Series(
                                {'Chr': irow.Chr, 'arm': arm, 'Overlap.Start': overlap[0],
                                 'Overlap.End': overlap[1], 'start1': irow['Peak.Start'], 'end1': irow['Peak.End'], 'direction1': irow.direction, 'telcent1': irow.telcent,
                                 'negpos1': irow.negpos, 'iter1': irow.iter, 'code1': irow.code, 'ksby1': irow.ksby,
                                 'combined_sig1': irow.combined_sig, 'direction2': jrow.direction,
                                 'telcent2': jrow.telcent, 'start2': jrow['Peak.Start'], 'end2': jrow['Peak.End'],
                                 'negpos2': jrow.negpos, 'iter2': jrow.iter, 'code2': jrow.code, 'ksby2': jrow.ksby,
                                 'combined_sig2': jrow.combined_sig,
                                 'combined_sig_sum': irow.combined_sig + jrow.combined_sig, 'genes': genes, 'consistent': consistent})
                            #print s
                            s_list.append(s)
    if len(s_list) > 0:
        dfdf0 = pd.DataFrame(s_list)
        dfdf = dfdf0.drop_duplicates(['Overlap.Start', 'Overlap.End', 'arm', 'combined_sig_sum'])
        dfdf = dfdf.sort_values('combined_sig_sum',ascending=False) #steph edit
        cols = ['Chr','arm','Overlap.Start','Overlap.End'] + [i+'1' for i in ['start','end','direction','telcent','negpos','iter','code','ksby','combined_sig']] + \
               [i + '2' for i in ['start','end','direction', 'telcent', 'negpos', 'iter', 'code', 'ksby', 'combined_sig']] + ['combined_sig_sum','genes','consistent']
        dfdf = dfdf[cols]
        dfdf.to_csv(folder + '/BISCUT_overlapping_peaks.txt', sep='\t', index=False)
        

def overlap_significance(folder, overlapsdf, abslocs_file, num_perms = 1000):
    info = pd.read_csv(abslocs_file, sep='\t',index_col='chromosome_info').transpose().to_dict()
    def coords(arm):
        if arm in ['13', '14', '15', '21', '22']:
            coord = (info[int(arm)]['q_start'], info[int(arm)]['q_end'])
        elif arm.endswith('q'):
            coord = (info[int(arm[:-1])]['q_start'], info[int(arm[:-1])]['q_end'])
        elif arm.endswith('p'):
            coord = (info[int(arm[:-1])]['p_start'], info[int(arm[:-1])]['p_end'])
        else:
            coord = (
            info[int(arm)]['p_start'], info[int(arm)]['p_end'], info[int(arm)]['q_start'], info[int(arm)]['q_end'])
        return coord

    def make_permutations(results):
        try:
            newdf = pd.DataFrame(index=results.index)
            for id in results.index:
                # print id
                peak_len = results.loc[[id], 'peak_length'].iloc[0]
                type = id[0]
                dir = id[1]
                telcent = id[2]
                code = id[3]
                len_code = len(code.split('-'))
                results['r_dep'] = results['code'].map(lambda x: all(['r' in i for i in x.split('-')[len_code:]]))
                results['l_dep'] = results['code'].map(lambda x: all(['l' in i for i in x.split('-')[len_code:]]))
                #print results
                tempresults = results[(results.code.str.startswith(code)) & (results.code != code)]
                tempresults = tempresults[tempresults['direction'] == dir]
                tempresults = tempresults[tempresults['telcent'] == telcent]
                r_dep = sum(tempresults[tempresults.r_dep].peak_length)
                l_dep = sum(tempresults[tempresults.l_dep].peak_length)
                #print 'r_dep', r_dep
                #print 'l_dep', l_dep
                if results.loc[[id], 'iter'].iloc[0] == 1:  # first iter
                    if (arm.endswith('p') and telcent == 'tel') or (not arm.endswith('p') and telcent == 'cent'):  # normal
                        rand_start = np.random.randint(arm_start + l_dep, arm_end - r_dep - peak_len)
                        newdf.set_value(id, 'start', rand_start)
                        newdf.set_value(id, 'end', rand_start + peak_len)
                        newdf.set_value(id, 'r_interval_start', rand_start + peak_len + 1)
                        newdf.set_value(id, 'r_interval_end', arm_end)
                        newdf.set_value(id, 'l_interval_start', arm_start)
                        newdf.set_value(id, 'l_interval_end', rand_start - 1)
                    else:  # q
                        rand_start = np.random.randint(arm_start + r_dep, arm_end - l_dep - peak_len)
                        newdf.set_value(id, 'start', rand_start)
                        newdf.set_value(id, 'end', rand_start + peak_len)
                        newdf.set_value(id, 'l_interval_start', rand_start + peak_len + 1)
                        newdf.set_value(id, 'l_interval_end', arm_end)
                        newdf.set_value(id, 'r_interval_start', arm_start)
                        newdf.set_value(id, 'r_interval_end', rand_start - 1)
                # WORK ON THIS
                else:  # more than first iter
                    prev_id = (type, dir, telcent, '-'.join(code.split('-')[:-1]))
                    new_direction = code[-2]
                    lstart = newdf.loc[[prev_id], 'l_interval_start'].iloc[0]
                    lend = newdf.loc[[prev_id], 'l_interval_end'].iloc[0]
                    rstart = newdf.loc[[prev_id], 'r_interval_start'].iloc[0]
                    rend = newdf.loc[[prev_id], 'r_interval_end'].iloc[0]
                    if (arm.endswith('p') and telcent == 'tel') or (
                                not arm.endswith('p') and telcent == 'cent'):  # normal direction
                        if new_direction == 'l':
                            # pass
                            rand_end = np.random.randint(lstart + l_dep + peak_len, lend)
                            newdf.set_value(id, 'end', rand_end)
                            newdf.set_value(id, 'start', rand_end - peak_len)
                            newdf.set_value(id, 'r_interval_start', rand_end + 1)
                            newdf.set_value(id, 'r_interval_end', rend)
                            newdf.set_value(id, 'l_interval_end', rand_end - peak_len - 1)
                            newdf.set_value(id, 'l_interval_start', lstart)
                        else:  # if new_direction is r
                            rand_start = np.random.randint(rstart, rend - peak_len - r_dep)
                            newdf.set_value(id, 'start', rand_start)
                            newdf.set_value(id, 'end', rand_start + peak_len)
                            newdf.set_value(id, 'r_interval_start', rand_start + peak_len + 1)
                            newdf.set_value(id, 'r_interval_end', rend)
                            newdf.set_value(id, 'l_interval_start', lstart)
                            newdf.set_value(id, 'l_interval_end', rand_start - 1)
                    else:  # opposite direction
                        if new_direction == 'r':
                            rand_end = np.random.randint(rstart + r_dep + peak_len, rend)
                            newdf.set_value(id, 'end', rand_end)
                            newdf.set_value(id, 'start', rand_end - peak_len)
                            newdf.set_value(id, 'r_interval_start', rstart)
                            newdf.set_value(id, 'r_interval_end', rand_end - peak_len - 1)
                            newdf.set_value(id, 'l_interval_start', rand_end + 1)
                            newdf.set_value(id, 'l_interval_end', lend)
                        else:  # if new direction is l
                            rand_start = np.random.randint(lstart, lend - peak_len - l_dep)
                            newdf.set_value(id, 'start', rand_start)
                            newdf.set_value(id, 'end', rand_start + peak_len)
                            newdf.set_value(id, 'l_interval_start', rand_start + peak_len + 1)
                            newdf.set_value(id, 'l_interval_end', lend)
                            newdf.set_value(id, 'r_interval_start', rstart)
                            newdf.set_value(id, 'r_interval_end', rand_start - 1)
            return newdf  # return a df with a permuted set of peaks
        except:
            pass


    #groups = overlapsdf.groupby('arm')
    all_results = pd.read_csv(folder+'/all_BISCUT_results.txt', sep='\t')


    # permute each type + arm separately and concat them for the comparisons
    fullcombolist = []
    for i in overlapsdf.index:
        fullcombolist.append((overlapsdf.loc[i,'arm'],overlapsdf.loc[i,'type1']))
        fullcombolist.append((overlapsdf.loc[i,'arm'],overlapsdf.loc[i,'type2']))

    fullcomboset = set(fullcombolist) #len is 418; this is each individual combo that shows up in the overlaps

    fullcombodic = {}
    for arm, type in fullcomboset:
        results = all_results[all_results.arm == arm]
        results = results[results.type == type]
        results = results.drop_duplicates(['Peak.Start', 'Peak.End', 'code'])
        results['peak_id'] = zip(results.type, results.direction, results.telcent, results.code)
        results = results.set_index('peak_id')
        results['peak_length'] = results['Peak.End'] - results['Peak.Start']
        results = results.sort_values(['type','iter'])
        #print results

        arm_start, arm_end = coords(arm)
        iteration_dfs = []
        while len(iteration_dfs) < num_perms:
            n = len(iteration_dfs) + 1
            newdf = make_permutations(results)
            if newdf is not None:
                newdf['perm'] = n
                iteration_dfs.append(newdf)

        permdf = pd.concat(iteration_dfs)
        #print permdf
        fullcombodic[(arm,type)] = permdf
    #print fullcombodic
    for i in overlapsdf.index:
        arm = overlapsdf.loc[i, 'arm']
        type1 = overlapsdf.loc[i, 'type1']
        type2 = overlapsdf.loc[i, 'type2']
        peak1 = (type1, overlapsdf.loc[i, 'direction1'], overlapsdf.loc[i, 'telcent1'], overlapsdf.loc[i, 'code1'])
        peak2 = (type2, overlapsdf.loc[i, 'direction2'], overlapsdf.loc[i, 'telcent2'], overlapsdf.loc[i, 'code2'])
        # if type1==type2: #within same tumor type
        #     comparison_df = fullcombodic[(arm,type1)]
        # else:
        #     comparison_df = pd.concat([fullcombodic[(arm,type1)], fullcombodic[(arm,type2)]])


        df1 = fullcombodic[(arm,type1)].loc[[peak1]].set_index('perm')
        df1['coords1'] = zip(df1.start, df1.end)
        df2 = fullcombodic[(arm,type2)].loc[[peak2]].set_index('perm')
        df2['coords2'] = zip(df2.start, df2.end)

        df1['coords2'] = df2['coords2']
        zz = df1.index.map(lambda x: overlap_helper_simple(df1.loc[x, 'coords1'], df1.loc[x, 'coords2']))
        #print len(zz[zz])
        overlapsdf.set_value(i, 'perm_sig', len(zz[zz]) / float(num_perms))
        #print overlapsdf

    # finaldf = pd.concat(recombine)
    # finaldf = finaldf.sort('combined_sig_sum', ascending=False)
    # finaldf.to_csv('something_011420.txt', sep='\t', index=False)
    overlapsdf.to_csv(results+'/BISCUT_overlaps_with_permuted_sig_210319.txt', sep='\t', index=False)
    overlapsdf[overlapsdf['consistent'] == True].to_csv(folder + '/BISCUT_overlaps_with_permuted_sig_consistent_only_210319.txt', sep='\t',
                                            index=False)
    return overlapsdf

def process_for_ggplot_jagged(results, abslocs_file):
    locs = pd.read_csv(abslocs_file, sep='\t',index_col='chromosome_info')

    #print tt
    df = pd.read_csv(results + '/all_BISCUT_results.txt', sep='\t')
    df['peak_id'] = zip(df.arm, df.direction, df.telcent,
                        df.code)
    sizes = df.groupby('peak_id').size()
    sizes = sizes[sizes <= 50]
    df = df[df['peak_id'].isin(sizes.index)]

    df = df.drop_duplicates(subset=['log10_ksby', 'Peak.Start', 'Peak.End', 'arm', 'direction','telcent', 'negpos'])

    tempfullname = list(set(df['type'].str.cat([df.arm.astype(str),df.s,df.telcent,df.code,df.conf.astype(str),'iter'+df.iter.astype(int).astype(str),'|'+df.log10_ksby.astype(str)],sep='_')))
    #print tempfullname
    fullname = {}
    for fn in tempfullname:
        fullname[fn.split('_|')[0]]=float(fn.split('_|')[1])
    #print fullname

    # for each peak, make a thing [xmin, xmax, ymin, ymax] [0,10,0,20]
    df['pq'] = df['Cytoband'].str[:1]
    # df['pq'] = 'void'

    ampfilltel = []
    delfilltel = []
    ampfillcent = []
    delfillcent = []
    for fa in fullname:
        arm = fa.split('_')[1]
        direc = fa.split('_')[2]
        telcent = fa.split('_')[3]
        prefix = fa.split('_')[4]
        #print fa

        pp = pd.read_csv(results+'/iterations/'+fa+'plotpeaks.txt',sep='\t')
        if arm.endswith('q') or arm in ['13','14','15','21','22']:
            pp['locx'] = sorted(pp['locx'],reverse=True)
        pp = pp.drop_duplicates(subset=['locx','distancey'])
        pp['distancey'] = pp['distancey']*fullname[fa]
        if arm in ['13','14','15','21','22']:
            pp['locx'] = pp['locx'] + locs.loc[int(arm),'offset']
        else:
            pp['locx'] = pp['locx'] + locs.loc[int(arm[:-1]), 'offset']
        #if fa == 'PANCAN_5q_del_cent_p_0.95_iter1': print pp
        pp = pp.drop('fraclocx', axis=1)
        pp = pp.reset_index(drop=True)
        pp = pp.join(pp.loc[1:].reset_index(drop=True), rsuffix='r')[:-1]
        #if fa == 'PANCAN_5q_del_cent_p_0.95_iter1': print pp
        if direc=='amp':
            if prefix[-1]=='p':
                pp['color'] = darkred
            else:
                pp['color'] = lightblue
            if telcent =='tel':
                ampfilltel.append(pp)
            else:
                ampfillcent.append(pp)
        else:
            if prefix[-1] =='p':
                pp['color'] = darkblue
            else:
                pp['color'] = lightred
            if telcent =='tel':
                delfilltel.append(pp)
            else:
                delfillcent.append(pp)
    #print 'printing' + tt
    #if genefilter==10000:
    try:
        pd.concat(ampfilltel).to_csv(results + '/summary/BISCUT_results_for_jagged_plotting_amp_tel.txt',
                                sep='\t', index=False, header=False)
    except:
        pass
    try:
        pd.concat(ampfillcent).to_csv(results + '/summary/BISCUT_results_for_jagged_plotting_amp_cent.txt',
                                sep='\t', index=False, header=False)
    except:
        pass
    try:
        pd.concat(delfilltel).to_csv(results + '/summary/BISCUT_results_for_jagged_plotting_del_tel.txt',
                                sep='\t', index=False, header=False)
    except:
        pass
    try:
        pd.concat(delfillcent).to_csv(results + '/summary/BISCUT_results_for_jagged_plotting_del_cent.txt',
                                sep='\t', index=False, header=False)
    except:
        pass

        centdfs = [i for i in os.listdir(os.path.join(results, 'summary')) if fnmatch.fnmatch(i, '*_for_jagged_plotting_*_cent.txt')]
        teldfs = [i for i in os.listdir(os.path.join(results, 'summary')) if fnmatch.fnmatch(i, '*_for_jagged_plotting_*_tel.txt')]
        togetherdfs = centdfs+teldfs
        #print centdfs
        #print teldfs

        alltogether = []
        if len(centdfs) !=0:
            try:
                centdf = pd.concat([pd.read_csv(os.path.join(results, 'summary',i),sep='\t',header=None) for i in centdfs if os.stat(os.path.join(results, 'summary',i)).st_size!=0 ] )
                alltogether.append(('cent',centdf))
            except: pass
        if len(teldfs) !=0:
            try:
                teldf = pd.concat([pd.read_csv(os.path.join(results, 'summary',i),sep='\t',header=None) for i in teldfs if os.stat(os.path.join(results, 'summary',i)).st_size!=0 ] )
                alltogether.append(('tel',teldf))
            except: pass
        if len(togetherdfs) != 0:
            try:
                togetherdf = pd.concat([pd.read_csv(os.path.join(results, 'summary',i),sep='\t',header=None) for i in togetherdfs if os.stat(os.path.join(results,'summary',i)).st_size!=0 ] )
                alltogether.append(('telcent',togetherdf))
            except:pass
        #print alltogether
        if len(alltogether)>0:
            for tc, df in alltogether:
                tempdf = df
                tempdf['abssmall'] = [min(i) for i in zip(tempdf[1].abs(),tempdf[3].abs())]
                tempdf['absbig'] = [max(i) for i in zip(tempdf[1].abs(),tempdf[3].abs())]

                tempdf = tempdf.reset_index(drop=True)
                #print tempdf
                for i in tempdf.index:
                    #print i
                    #print df.loc[i,4]
                    if tempdf.loc[i,4] == lightred or tempdf.loc[i,4]==darkred: #on the right side, onco-gene liike
                        tempdf.set_value(i,1,tempdf.loc[i,'abssmall'])
                        tempdf.set_value(i,3,tempdf.loc[i,'absbig'])
                    else:
                        tempdf.set_value(i,1,-tempdf.loc[i,'absbig'])
                        tempdf.set_value(i,3,-tempdf.loc[i,'abssmall'])
                tempdf = tempdf.filter(items=[0,1,2,3,4])
                #print tempdf
                posdf = tempdf[tempdf[4].isin([darkblue,darkred])]
                negdf = tempdf[tempdf[4].isin([lightblue,lightred])]
                #print posdf
                #print negdf

                # posdf.to_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results_for_jagged_plotting_pos_'+tc+'.txt',sep='\t',index=False,header=False)
                # negdf.to_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results_for_jagged_plotting_neg_'+tc+'.txt',sep='\t',index=False,header=False)
                # pd.concat([posdf,negdf]).to_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results_for_jagged_plotting_posneg_'+tc+'.txt',sep='\t',index=False,header=False)

                #if genefilter == 10000:
                posdf.to_csv(
                    results + '/summary/BISCUT_results_for_jagged_plotting_pos_' + tc + '.txt',
                    sep='\t', index=False, header=False)
                negdf.to_csv(
                    results + '/summary/BISCUT_results_for_jagged_plotting_neg_' + tc + '.txt',
                    sep='\t', index=False, header=False)
                pd.concat([posdf, negdf]).to_csv(
                    results + '/summary/BISCUT_results_for_jagged_plotting_posneg_' + tc + '.txt',
                    sep='\t', index=False, header=False)

# Jeff M: Function currently not called in workflow, so not editing. At the least, this function
# would need to be updated to be compatible changes made to the output directory structure.
def make_all_peaks(results, qval_thres):
    cancer = pd.read_csv('cancer_genes_030718.txt', sep='\t', index_col='Unnamed: 0')
    for tt in os.listdir(results):
        try:
            BISCUT = pd.read_csv(results + '/' + tt + '/' + tt + '_BISCUT_results.txt', sep='\t')
            BISCUT = BISCUT[BISCUT['ksby'] <= qval_thres]
            BISCUT['log10_ksby'] = BISCUT['log10_ksby'].replace({np.inf: 16})
            groups = BISCUT.groupby(['arm', 'direction', 'code'])

            allgenes = []
            allcancergenes = []
            dicfordf = {}
            for combo, df in groups:
                direction = df['direction'].tolist()[0]
                negpos = df['negpos'].tolist()[0]
                stopgo = 'STOP' if (direction == 'del' and negpos == 'p') or (direction == 'amp' and negpos == 'n') else 'GO'
                if direction == 'amp' and negpos == 'n':
                    color = lightblue
                elif direction == 'del' and negpos == 'p':
                    color = darkblue
                elif direction == 'del' and negpos == 'n':
                    color = lightred
                elif direction == 'amp' and negpos == 'p':
                    color = darkred
                cyto = df.Cytoband.value_counts().index[0]
                # cyto = 'void'
                genes = df['Gene'].tolist()
                peakloc = 'chr' + str(df['Chr'].tolist()[0]) + ':' + str(df['Peak.Start'].tolist()[0]) + '-' + str(
                    df['Peak.End'].tolist()[0])
                if genes[0].startswith('['):
                    genes = [genes[0][1:-1]]
                minicancer = cancer[cancer.index.isin(genes)]
                minicancergenes = minicancer.index.tolist()
                miniscores = [minicancer.loc[i, 'Total_Score'] for i in minicancergenes]
                dicfordf['_'.join(combo)] = {'n_total_genes': len(genes), 'n_driver_genes': len(minicancer.index.tolist()),
                                             'driver_genes': ', '.join(minicancergenes),
                                             'driver_score': sum(minicancer['Total_Score'] / len(genes)),
                                             'log10_ksby': df['log10_ksby'].tolist()[0], 'ks_stat': df['ks_stat'].tolist()[0],
                                             'combined_sig': df['log10_ksby'].tolist()[0] * df['ks_stat'].tolist()[0],
                                             'all_genes': ', '.join(genes),
                                             'driver_gene_points': ', '.join(
                                                 [str(minicancer.loc[i, 'Total_Score']) for i in minicancergenes]),
                                             'STOP_or_GO': stopgo, 'direction': direction, 'negpos': negpos, 'color': color,
                                             'max_driver_gene_points': max(miniscores) if len(miniscores) > 0 else 0,
                                             'cytoband': cyto, 'peak_location': peakloc}
                allgenes = allgenes + genes
                allcancergenes = allcancergenes + minicancergenes
            #print 'in peaks, duplicates', len(allgenes)
            #print 'in peaks, no duplicates', len(list(set(allgenes)))
            #print 'driver genes in peaks, duplicates', len(allcancergenes)
            #print 'driver genes in peaks, no duplicates', len(list(set(allcancergenes)))
            # print dicfordf
            # print pd.DataFrame(dicfordf).transpose()
            newdf = pd.DataFrame(dicfordf).transpose()
            newdf =newdf[['peak_location','cytoband','STOP_or_GO','direction','negpos','color','combined_sig','log10_ksby','ks_stat',
                          'all_genes','n_total_genes','driver_genes','n_driver_genes','driver_gene_points','max_driver_gene_points',
                          'driver_score']]
            newdf = newdf.sort_values('combined_sig',ascending=False)
            newdf.to_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results_all_peaks.txt', sep='\t', index_label='combo')

        except:
            pass

def filter_BISCUT_knowngenes(folder, genes):
    df = pd.read_csv(folder+'/all_BISCUT_results.txt',sep='\t')
    if not os.path.exists(os.path.join(folder,'genes')): os.mkdir(os.path.join(folder,'genes'))
    if not os.path.exists(os.path.join(folder,'genes','files')): os.mkdir(os.path.join(folder,'genes','files'))
    
    # If just 1 gene is fed in from R, will get interpreted as string; convert to list
    if isinstance(genes, str):
      genes = [genes]
    for g in genes:
        try:
            minidf = df[df['Gene']==g]
            minidf = minidf.sort_values(by=['combined_sig','n_events'],ascending=[False,False])
            minidf.to_csv(folder+'/genes/files/'+g+'_BISCUT_results.txt',sep='\t',index=False)
            minidf = minidf.reset_index(drop=True)

            def colors(x):
                if x['direction']=='amp' and x['negpos']=='p':
                    return darkred
                elif x['direction']=='del' and x['negpos']=='n':
                    return lightred
                elif x['direction'] == 'del' and x['negpos'] == 'p':
                    return darkblue
                else:
                    return lightblue
            minidf['colors'] = minidf.apply(lambda x: colors(x),axis=1)
            start = minidf.loc[0,'Start']
            end = minidf.loc[0,'End']
            truncdf = minidf.filter(items=['Peak.Start','Peak.End','n_events','type','direction','telcent','iter','code','combined_sig','colors'])
            truncdf = truncdf.append(
                pd.DataFrame([[start, end, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, '#000000']], columns=truncdf.columns))
            truncdf = truncdf.reset_index(drop=True)
            truncdf['ymax'] = (-truncdf.index)-0.1
            truncdf['ymin'] = (-truncdf.index)-1
            truncdf.to_csv(folder+'/genes/files/'+g+'_BISCUT_fig2.txt',sep='\t',index=False)
        except:
            pass

def filter_BISCUT_arms(folder, arms):
    # If just 1 arm is fed in from R, will get interpreted as string; convert to list
    if isinstance(arms, str):
      arms = [arms]
    df = pd.read_csv(folder+'/all_BISCUT_results.txt',sep='\t')
    if not os.path.exists(os.path.join(folder,'arms')): os.mkdir(os.path.join(folder,'arms'))
    if not os.path.exists(os.path.join(folder,'arms','files')): os.mkdir(os.path.join(folder,'arms','files'))
    for g in arms:
        #try:
        if not os.path.exists(os.path.join(folder, 'arms', 'files',g)): os.mkdir(os.path.join(folder, 'arms', 'files',g))
        arm_dic = {i:{} for i in [('amp','n'),('del','n'),('amp','p'),('del','p')]}
        for c in [('amp','n'),('del','n'),('amp','p'),('del','p')]:
            minidf = df.drop_duplicates(subset=['arm','direction','telcent','negpos','code'])
            minidf = minidf[minidf['arm']==g]
            minidf = minidf[minidf['direction']==c[0]]
            minidf = minidf[minidf['negpos']==c[1]]
            minidf = minidf.sort_values(by=['combined_sig','n_events'],ascending=[False,False])
           #includes multiple iterations
            minidf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_multiter_BISCUT_results.txt',sep='\t',index=False)
            arm_dic[c]['_multiter_BISCUT_results'] = minidf

            minidf = minidf.reset_index(drop=True)
            if minidf.empty:
                continue

            def colors(x):
                if x['direction']=='amp' and x['negpos']=='p':
                    return darkred
                elif x['direction']=='del' and x['negpos']=='n':
                    return lightred
                elif x['direction'] == 'del' and x['negpos'] == 'p':
                    return darkblue
                else:
                    return lightblue
            minidf['colors'] = minidf.apply(lambda x: colors(x),axis=1)
            start = minidf.loc[0,'Start']
            end = minidf.loc[0,'End']
            truncdf = minidf.filter(items=['Peak.Start','Peak.End','n_events','direction','telcent','iter','code','combined_sig','colors'])
            truncdf = truncdf.reset_index(drop=True)
            truncdf['ymax'] = (-truncdf.index)-0.1
            truncdf['ymin'] = (-truncdf.index)-1
            truncdf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_multiter_BISCUT_fig2.txt',sep='\t',index=False)
            arm_dic[c]['_multiter_BISCUT_fig2'] = truncdf
            truncdf = truncdf[~truncdf.code.str[:-2].str.contains(c[1])]
            truncdf = truncdf.reset_index(drop=True)
            truncdf['ymax'] = (-truncdf.index)-0.1
            truncdf['ymin'] = (-truncdf.index)-1
            truncdf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_firstiter_BISCUT_fig2.txt',sep='\t',index=False)
            arm_dic[c]['_firstiter_BISCUT_fig2'] = truncdf


            truncdf = truncdf[truncdf['code']==c[1]]
            truncdf = truncdf.reset_index(drop=True)
            truncdf['ymax'] = (-truncdf.index)-0.1
            truncdf['ymin'] = (-truncdf.index)-1
            truncdf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_BISCUT_fig2.txt',sep='\t',index=False)
            arm_dic[c]['_BISCUT_fig2'] = truncdf

           #filters out multiple iterations
            minidf = minidf[~minidf.code.str[:-2].str.contains(c[1])]
            minidf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_firstiter_BISCUT_results.txt',sep='\t',index=False)
            arm_dic[c]['_firstiter_BISCUT_results'] = minidf


            minidf = minidf[minidf['code']==c[1]]
            minidf.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+c[0]+'_'+c[1]+'_BISCUT_results.txt',sep='\t',index=False)
            arm_dic[c]['_BISCUT_results'] = minidf
        #print g
        #print arm_dic
        for d in ['_multiter_BISCUT_results','_multiter_BISCUT_fig2','_firstiter_BISCUT_fig2','_BISCUT_fig2','_firstiter_BISCUT_results','_BISCUT_results']:
            for c in [('amp','n'),('del','n'),('amp','p'),('del','p')]:
                if d not in arm_dic[c]:
                    arm_dic[c][d] = pd.DataFrame()
            onco = pd.concat([arm_dic[('amp','p')][d], arm_dic[('del','n')][d]])
            ts = pd.concat([arm_dic[('amp','n')][d], arm_dic[('del','p')][d]])
           # print onco
            #print ts
            pre = 'onco'
            for mini in [onco,ts]:
                if not mini.empty:
                    mini = mini.sort_values(by=['combined_sig','n_events'],ascending=[False,False])
                    #ts = ts.sort(columns=['combined_sig','n_events'],ascending=[False,False])
                    mini = mini.reset_index(drop=True)
                    mini['ymax'] = (-mini.index)-0.1
                    mini['ymin'] = (-mini.index)-1
                mini.to_csv(folder+'/arms/files/'+g+'/'+g+'_'+pre+d+'.txt',sep='\t',index=False)
                pre = 'ts'


def all_processing(results_dir, arms, genelocs):
    calc_overlaps(results_dir, genelocs) 
    #process_for_ggplot_jagged(folder, qval_thres, abslocs_file)
    #filter_BISCUT_knowngenes(folder, genes)
    filter_BISCUT_arms(results_dir,arms)
