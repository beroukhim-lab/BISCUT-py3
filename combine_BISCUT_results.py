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

def make_table_results(folder, kspvals, qval_thres):
    ci = folder.split('_')[-1]
    kspvalsdf = pd.read_csv(kspvals,sep='\t',index_col='combo')
    alldfs = []
    #for tt in [i for i in os.listdir(folder) if not i.startswith('.') and not i.endswith('.txt') and not i.endswith('.py') and not i.endswith('.pdf') and not i=='stats']:
    for tt in next(os.walk(folder))[1]:
   # for tt in ['PANCAN']:
        if not os.path.exists(os.path.join(folder,tt,'summary')): os.mkdir(os.path.join(folder,tt,'summary'))
        li = []
        r = re.compile('iter\d.txt')
        #x =  [i for i in os.listdir(folder+'/'+tt) if i.endswith('.txt') and not i.endswith('BISCUT_results.txt')
        #      and not i.endswith('plotpeaks.txt')]
        x = list(filter(r.search,os.listdir(folder+'/'+tt))) #steph edit
        if len(x)>0:
            for f in x: #all files
                #print f
                forindex = '_'.join(f.split('_')[:5]+[f.split('_')[-1][4]])
                df = pd.read_csv(folder+'/'+tt+'/'+f,sep='\t')
                df= df.rename(columns={'ks':'ksp','log10_ks':'log10_ksp'})
                df['ksby'] = kspvalsdf.loc[forindex,'by']
                df['log10_ksby'] = -np.log10(df['ksby'])
                df = df.replace(to_replace={'log10_ksby':{np.inf:16}})
                df['combined_sig'] = df['ks_stat'] * df['log10_ksby']
                df['code'] = forindex.split('_')[4]
                #this is for pruning
                code_parts = forindex.split('_')[4].split('-')
                iter = int(forindex.split('_')[5])
                if iter>1:
                    all_previous= []
                    all_previous_files =[]
                    for j in range(1,iter): #for iter==5, go from 1 to 4
                        previous_code = '-'.join(code_parts[:j])
                        previous_iter = str(j)
                        previous_index = '_'.join(forindex.split('_')[:4])+'_'+previous_code+'_'+previous_iter
                        previous_file = '_'.join(forindex.split('_')[:4])+'_'+previous_code+'_'+ci+'_iter'+previous_iter+'.txt'
                        all_previous.append(previous_index)
                        all_previous_files.append(previous_file)
                        #print forindex, previous_index
                    #print all_previous
                    if any([kspvalsdf.loc[previous_index,'by'] > qval_thres for previous_index in all_previous]) \
                            or any([d.empty for d in [pd.read_csv(folder+'/'+tt+'/'+e) for e in all_previous_files]]):
                        #print forindex
                        df['ksby'] = 1
                        kspvalsdf.at[forindex,'by'] = 1 #steph edit
                #END pruning
                df = df[df['ksby']<=qval_thres] #empty df if not significant.

                # #Pruning part 2
                # if df.empty and iter ==1: #either because not significant or because peak was too big
                #     print forindex
                #     all_downstream_codes = [c for c in kspvals.index if c.startswith('_'.join(forindex.split('_')[:-1])) and c!=forindex]
                #     kspvalsdf.set_value(all_downstream_codes,'by',1)
                # #end pruning part 2

                li.append(df)
                df.to_csv(folder+'/'+tt+'/'+f,sep='\t',index=False)


            bigone = pd.concat(li)
            bigone.to_csv(folder+'/'+tt+'/summary/'+tt+'_BISCUT_results.txt',sep='\t',index=False)
            rnk = bigone.filter(['Gene','combined_sig'])
            aprnk = bigone[(bigone['direction']=='amp')&(bigone['negpos']=='p')].filter(['Gene','combined_sig'])
            dprnk = bigone[(bigone['direction']=='del')&(bigone['negpos']=='p')].filter(['Gene','combined_sig'])
            anrnk = bigone[(bigone['direction']=='amp')&(bigone['negpos']=='n')].filter(['Gene','combined_sig'])
            dnrnk = bigone[(bigone['direction']=='del')&(bigone['negpos']=='n')].filter(['Gene','combined_sig'])
            tsrnk = pd.concat([dprnk,anrnk])
            oncrnk = pd.concat([dnrnk,aprnk])
            names = [folder+'/'+tt+'/summary/'+tt+'_BISCUT_results'+x+'.rnk' for x in ['','_ts-like','_onc-like','_amp-p','_del-p','_amp-n','_del-n']]
            for rn,n in zip([rnk,tsrnk,oncrnk,aprnk,dprnk,anrnk,dnrnk],names):
                rnz = rn.sort_values(by = 'combined_sig',ascending=False) #steph edit
                rnz = rnz.drop_duplicates('Gene')
                rnz.to_csv(n,sep='\t',index=False)
                 #rnz.to_csv('
            alldfs.append(bigone)
    combined = pd.concat(alldfs)
    combined.to_csv(folder+'/all_BISCUT_results.txt',sep='\t',index=False)
    kspvalsdf.to_csv(kspvals, sep='\t')

def calc_overlaps(folder, genelocs_file):
    genelocs = pd.read_csv(genelocs_file,sep='\t')
    df = pd.read_csv(folder+'/all_BISCUT_results.txt', sep='\t')
    armgroups = df.groupby('arm')
    #typegroups = df.groupby('type')
    s_list = []
    # for tt, ttdf in typegroups:
    #     armgroups = ttdf.groupby('arm')
    for arm, armdf in armgroups:  # specific arm within specific tumor type
        # if arm == '1p':
        armdf = armdf.drop_duplicates(['type','Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'code'])
        #print armdf
        for dir, telcent, negpos in itertools.product(['amp', 'del'], ['tel', 'cent'], ['n', 'p']):
            # within is everything that shares same path
            within = armdf[(armdf.direction == dir) & (armdf.telcent == telcent) & (armdf.negpos == negpos)]
            without = armdf[(armdf.direction != dir) | (armdf.telcent != telcent) | (armdf.negpos != negpos)]
            within = within.filter(
                ['Chr', 'arm', 'type', 'Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'iter', 'code',
                 'ksby', 'combined_sig'])
            without = without.filter(
                ['Chr', 'arm', 'type', 'Peak.Start', 'Peak.End', 'direction', 'telcent', 'negpos', 'iter', 'code',
                 'ksby', 'combined_sig'])
            #print dir, telcent, negpos
            #print len(within)
            #print len(without)
            if not within.empty:
                for i in within.index:
                    type1 = within.loc[i,'type']
                    #print type1
                    difftypewithin = within[within['type']!=type1]
                    newwithout = pd.concat([without,difftypewithin], ignore_index=True)
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
                                 'Overlap.End': overlap[1], 'type1': irow.type, 'start1': irow['Peak.Start'], 'end1': irow['Peak.End'], 'direction1': irow.direction, 'telcent1': irow.telcent,
                                 'negpos1': irow.negpos, 'iter1': irow.iter, 'code1': irow.code, 'ksby1': irow.ksby,
                                 'combined_sig1': irow.combined_sig, 'type2': jrow.type, 'direction2': jrow.direction,
                                 'telcent2': jrow.telcent, 'start2': jrow['Peak.Start'], 'end2': jrow['Peak.End'],
                                 'negpos2': jrow.negpos, 'iter2': jrow.iter, 'code2': jrow.code, 'ksby2': jrow.ksby,
                                 'combined_sig2': jrow.combined_sig,
                                 'combined_sig_sum': irow.combined_sig + jrow.combined_sig, 'genes': genes, 'consistent': consistent})
                            #print s
                            s_list.append(s)
    dfdf0 = pd.DataFrame(s_list)
    dfdf = dfdf0.drop_duplicates(['Overlap.Start', 'Overlap.End', 'arm', 'combined_sig_sum'])
    dfdf = dfdf.sort_values('combined_sig_sum',ascending=False) #steph edit
    cols = ['Chr','arm','Overlap.Start','Overlap.End'] + [i+'1' for i in ['type','start','end','direction','telcent','negpos','iter','code','ksby','combined_sig']] + \
           [i + '2' for i in ['type','start','end','direction', 'telcent', 'negpos', 'iter', 'code', 'ksby', 'combined_sig']] + ['combined_sig_sum','genes','consistent']
    dfdf = dfdf[cols]
    #print dfdf
    dfdf.to_csv(folder+'/BISCUT_overlaps_011320.txt', sep='\t', index=False)
    dfdf[dfdf['consistent']==True].to_csv(folder+'/BISCUT_overlaps_consistent_only_210319.txt',sep='\t',index=False)
    return dfdf

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

def make_column_results(folder, qval_thres):
    if not os.path.exists(os.path.join(folder,'all_cols')): os.mkdir(os.path.join(folder,'all_cols'))
    list_of_cols = []
    #for tt in [i for i in os.listdir(folder) if not i.startswith('.') and not i.endswith('.txt') and not i.endswith('.py') and not i.endswith('pdf')]:
    for tt in next(os.walk(folder))[1]:
        if not os.path.exists(os.path.join(folder,tt,'summary')): os.mkdir(os.path.join(folder,tt,'summary'))

        r = re.compile('iter\d.txt')
        #x =  [i for i in os.listdir(folder+'/'+tt) if i.endswith('.txt') and not i.endswith('BISCUT_results.txt') and not i.endswith('0.9.txt') and not i.endswith('plotpeaks.txt')]
        x =  list(filter(r.search,os.listdir(folder+'/'+tt)))
        if len(x)>0:
            li = []
            for f in x:
                #print f
                df = pd.read_csv(folder+'/'+tt+'/'+f,sep='\t')
                df = df.replace(to_replace={'log10_ksby': {np.inf: 16}})
                if df.empty:
                    continue
                if df.loc[0,'ksby'] >qval_thres:
                    continue
                # cyto =df.Cytoband.value_counts().index[0]
                # peakband = str(df.loc[0,'Chr']) + cyto
                peakband = 'void'
                peakloc = 'chr'+str(df.loc[0,'Chr'])+':'+str(df.loc[0,'Peak.Start'])+'-'+str(df.loc[0,'Peak.End'])
                negpos = df.loc[0,'negpos']
                direction = df.loc[0,'direction']
                telcent = df.loc[0,'telcent']
                n_events = df.loc[0,'n_events']
                genes =df['Gene'].tolist()
                log10ksby = df.loc[0,'log10_ksby']
                sig = df.loc[0,'ks_stat']
                combinedsig = df.loc[0,'log10_ksby'] * df.loc[0,'ks_stat']
                code = df.loc[0,'code']
                if (direction=='del' and negpos == 'p') or (direction=='amp' and negpos == 'n'):
                    supposed = 'TS-like'
                else:
                    supposed='onco-like'

                thelist = [peakband, peakloc, combinedsig, log10ksby, sig, n_events, direction, telcent,negpos, code,supposed] +genes
                li.append(thelist)
            li = sorted(li, key=itemgetter(2),reverse=True)
            #print li
            if len(li)==0:
                continue
            maxlen = max([len(x) for x in li])
            newdf = pd.DataFrame(index=range(0,maxlen))
            for k in range(0,len(li)):
                newdf[k] = li[k] + ([np.nan]*(maxlen-len(li[k])))
            #print newdf
            #print maxlen
            leftheaders = ['cytoband','peak_location','combined_sig','log10_ksby','ks_stat','n_events','direction','telomeric or centromeric','selection','code','TS or onco-like','genes']
            newdf.insert(0,'stuff',leftheaders + ([np.nan]*(maxlen-len(leftheaders))))
            newdf.to_csv(folder+'/'+tt+'/summary/'+tt+'_BISCUT_results_cols_'+folder.split('_')[-1]+'.txt',sep='\t',index=False,header=False)
            newdf.to_csv(folder+'/all_cols/'+tt+'_BISCUT_results_cols_'+folder.split('_')[-1]+'.txt',sep='\t',index=False,header=False)
            list_of_cols.append((tt,newdf))
    return list_of_cols

def process_for_ggplot(results, qval_thres, abslocs_file):

    locs = pd.read_csv(abslocs_file,sep='\t',index_col='chromosome_info')
    for tt in [i for i in os.listdir(results) if os.path.isdir(os.path.join(results,i)) and i not in ('all_cols','genes','arms')]:
        try:
            #print tt
            df = pd.read_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results.txt',sep='\t')

            df['peak_id'] = zip(df.arm, df.direction, df.telcent,
                                df.code)

            if tt=='PANCAN':
                sizes = df.groupby('peak_id').size()
                sizes = sizes[sizes <= 50]
                df = df[df['peak_id'].isin(sizes.index)]

            df = df.replace(to_replace={'log10_ksby':{np.inf:16}})

            df = df.drop_duplicates(subset=['log10_ksby','Peak.Start','Peak.End','arm','direction','telcent','negpos'])
            df = df[df['ksby']<=qval_thres]

            df['log10_ksby']=df.apply(lambda x: -1*x['log10_ksby'] if x['direction']=='del' else x['log10_ksby'],axis=1)
            df['combined_sig'] = df['log10_ksby'] * df['ks_stat']

            # for each peak, make a thing [xmin, xmax, ymin, ymax] [0,10,0,20]
            # df['pq'] = df['Cytoband'].str[:1]
            df['pq'] = 'void'
            #amps
            #for tc in ['tel','cent']:
            try:
                df_amp = df[df['direction']=='amp']
                #df_amp = df_amp[df_amp['telcent']==tc]
                coords = []
                for i in df_amp.index:
                    xmin = df_amp.loc[i,'Peak.Start'] + locs.loc[int(df_amp.loc[i,'Chr']),'offset']
                    xmax = df_amp.loc[i,'Peak.End'] + locs.loc[int(df_amp.loc[i,'Chr']),'offset']
                    ymin = min([0,df_amp.loc[i,'combined_sig']])
                    ymax = max([0,df_amp.loc[i,'combined_sig']])
                    if min([ymin,ymax]) < 0: color =lightblue
                    if max([ymin, ymax]) > 0: color = darkred
                    #print [xmin,xmax,ymin,ymax,color]
                    coords.append([xmin,xmax,ymin,ymax,color])
                #if genefilter == 10000:
                pd.DataFrame(coords).to_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results_for_plotting_amp.txt',sep='\t',index=False,header=False)
            except:
                pass

            try:
                df_del = df[df['direction']=='del']
                coords = []
                for i in df_del.index:
                    xmin = df_del.loc[i,'Peak.Start']+ locs.loc[int(df_del.loc[i,'Chr']),'offset']
                    xmax = df_del.loc[i,'Peak.End']+ locs.loc[int(df_del.loc[i,'Chr']),'offset']
                    ymin = min([0,df_del.loc[i,'combined_sig']])
                    ymax = max([0,df_del.loc[i,'combined_sig']])
                    if min([ymin,ymax]) < 0: color =lightred
                    if max([ymin,ymax]) > 0: color =darkblue
                    #print [xmin,xmax,ymin,ymax,color]
                    coords.append([xmin,xmax,ymin,ymax,color])
                #if genefilter == 10000:
                pd.DataFrame(coords).to_csv(results+'/'+tt+'/summary/'+tt+'_BISCUT_results_for_plotting_del.txt',sep='\t',index=False,header=False)
            except:
                pass

            try:
                df_pos = df[df['negpos'] == 'p']
                coords = []
                for i in df_pos.index:
                    xmin = df_pos.loc[i, 'Peak.Start'] + locs.loc[int(df_pos.loc[i, 'Chr']), 'offset']
                    xmax = df_pos.loc[i, 'Peak.End'] + locs.loc[int(df_pos.loc[i, 'Chr']), 'offset']
                    ymin = min([0, df_pos.loc[i, 'combined_sig']])
                    ymax = max([0, df_pos.loc[i, 'combined_sig']])
                    if min([ymin, ymax]) < 0: color = darkblue
                    if max([ymin, ymax]) > 0: color = darkred
                    # print [xmin,xmax,ymin,ymax,color]
                    coords.append([xmin, xmax, ymin, ymax, color])
                # if genefilter == 10000:
                pd.DataFrame(coords).to_csv(
                    results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_plotting_pos.txt',
                    sep='\t', index=False, header=False)
            except:
                pass

            try:
                df_neg = df[df['negpos'] == 'n']
                # df_del = df_del[df_del['telcent']==tc]
                coords = []
                for i in df_neg.index:
                    xmin = df_neg.loc[i, 'Peak.Start'] + locs.loc[int(df_neg.loc[i, 'Chr']), 'offset']
                    xmax = df_neg.loc[i, 'Peak.End'] + locs.loc[int(df_neg.loc[i, 'Chr']), 'offset']
                    ymin = min([0, df_neg.loc[i, 'combined_sig']])
                    ymax = max([0, df_neg.loc[i, 'combined_sig']])
                    if min([ymin, ymax]) < 0: color = lightred
                    if max([ymin, ymax]) > 0: color = lightblue
                    # print [xmin,xmax,ymin,ymax,color]
                    coords.append([xmin, xmax, ymin, ymax, color])
                # if genefilter == 10000:
                pd.DataFrame(coords).to_csv(
                    results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_plotting_neg.txt',
                    sep='\t', index=False, header=False)
            except:
                pass

        except:
            pass


def process_for_ggplot_jagged(results, qval_thres, abslocs_file):
    locs = pd.read_csv(abslocs_file, sep='\t',index_col='chromosome_info')
    #for tt in ['PANCAN']:
    for tt in [i for i in next(os.walk(results))[1] if i not in ['stats','all_cols','arms','genes']]:
        #print tt
        try:
            df = pd.read_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results.txt', sep='\t')

            df['peak_id'] = zip(df.arm, df.direction, df.telcent,
                                df.code)
            sizes = df.groupby('peak_id').size()
            sizes = sizes[sizes <= 50]
            df = df[df['peak_id'].isin(sizes.index)]

            df = df.replace(to_replace={'log10_ksby': {np.inf: 16}})

            df = df.drop_duplicates(subset=['log10_ksby', 'Peak.Start', 'Peak.End', 'arm', 'direction','telcent', 'negpos'])
            df = df[df['ksby'] <= qval_thres]

            tempfullname = list(set(df['type'].str.cat([df.arm.astype(str),df.direction,df.telcent,df.code,df.conf.astype(str),'iter'+df.iter.astype(int).astype(str),'|'+df.log10_ksby.astype(str)],sep='_')))
            #print tempfullname
            fullname = {}
            for fn in tempfullname:
                fullname[fn.split('_|')[0]]=float(fn.split('_|')[1])
            #print fullname

            # for each peak, make a thing [xmin, xmax, ymin, ymax] [0,10,0,20]
            # df['pq'] = df['Cytoband'].str[:1]
            df['pq'] = 'void'

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

                pp = pd.read_csv(results+'/'+tt+'/'+fa+'plotpeaks.txt',sep='\t')
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
                pd.concat(ampfilltel).to_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_amp_tel.txt',
                                        sep='\t', index=False, header=False)
            except:
                pass
            try:
                pd.concat(ampfillcent).to_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_amp_cent.txt',
                                        sep='\t', index=False, header=False)
            except:
                pass
            try:
                pd.concat(delfilltel).to_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_del_tel.txt',
                                        sep='\t', index=False, header=False)
            except:
                pass
            try:
                pd.concat(delfillcent).to_csv(results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_del_cent.txt',
                                        sep='\t', index=False, header=False)
            except:
                pass
        except:
            pass

        centdfs = [i for i in os.listdir(os.path.join(results,tt,'summary')) if fnmatch.fnmatch(i, '*_for_jagged_plotting_*_cent.txt')]
        teldfs = [i for i in os.listdir(os.path.join(results,tt,'summary')) if fnmatch.fnmatch(i, '*_for_jagged_plotting_*_tel.txt')]
        togetherdfs = centdfs+teldfs
        #print centdfs
        #print teldfs

        alltogether = []
        if len(centdfs) !=0:
            try:
                centdf = pd.concat([pd.read_csv(os.path.join(results,tt,'summary',i),sep='\t',header=None) for i in centdfs if os.stat(os.path.join(results,tt,'summary',i)).st_size!=0 ] )
                alltogether.append(('cent',centdf))
            except: pass
        if len(teldfs) !=0:
            try:
                teldf = pd.concat([pd.read_csv(os.path.join(results,tt,'summary',i),sep='\t',header=None) for i in teldfs if os.stat(os.path.join(results,tt,'summary',i)).st_size!=0 ] )
                alltogether.append(('tel',teldf))
            except: pass
        if len(togetherdfs) != 0:
            try:
                togetherdf = pd.concat([pd.read_csv(os.path.join(results,tt,'summary',i),sep='\t',header=None) for i in togetherdfs if os.stat(os.path.join(results,tt,'summary',i)).st_size!=0 ] )
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
                    results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_pos_' + tc + '.txt',
                    sep='\t', index=False, header=False)
                negdf.to_csv(
                    results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_neg_' + tc + '.txt',
                    sep='\t', index=False, header=False)
                pd.concat([posdf, negdf]).to_csv(
                    results + '/' + tt + '/summary/' + tt + '_BISCUT_results_for_jagged_plotting_posneg_' + tc + '.txt',
                    sep='\t', index=False, header=False)
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
                # cyto = df.Cytoband.value_counts().index[0]
                cyto = 'void'
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
    df = pd.read_csv(folder+'/all_BISCUT_results.txt',sep='\t')
    if not os.path.exists(os.path.join(folder,'arms')): os.mkdir(os.path.join(folder,'arms'))
    if not os.path.exists(os.path.join(folder,'arms','files')): os.mkdir(os.path.join(folder,'arms','files'))
    for g in arms:
        #try:
        if not os.path.exists(os.path.join(folder, 'arms', 'files',g)): os.mkdir(os.path.join(folder, 'arms', 'files',g))
        arm_dic = {i:{} for i in [('amp','n'),('del','n'),('amp','p'),('del','p')]}
        for c in [('amp','n'),('del','n'),('amp','p'),('del','p')]:
            minidf = df.drop_duplicates(subset=['arm','type','direction','telcent','negpos','code'])
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
            truncdf = minidf.filter(items=['Peak.Start','Peak.End','n_events','type','direction','telcent','iter','code','combined_sig','colors'])
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

def extract_cols(folder, lop, arms):
    if not os.path.exists(os.path.join(folder,'arms','cols')): os.mkdir(os.path.join(folder,'arms','cols'))
    for a in arms:
        #print a
        list_of_minicols=[]
        for tt, col in lop:
            #print a, tt
            try:
                zz = col.reset_index(drop=True)
                tcol = zz.transpose()
                minitcol = tcol[tcol[0].str.startswith(a)]
                #if minitcol.empty:
     #              break
                #print minitcol
                minitcol.insert(0,'tt',tt)
   #             minitcol = minitcol.reset_
                #print minitcol
                minicol = minitcol.transpose()
                minicol=minicol.reset_index(drop=True)
                list_of_minicols.append(minicol)
            except: pass
 #       li = sorted(list_of_minicols, key=itemgetter(3),reverse=True)
        df = pd.concat(list_of_minicols,axis=1)
        #print df
        df = df.transpose().sort_values(3,ascending=False).transpose()
 #       df = df.transpose()
        leftheaders = ['tumor_type','cytoband','peak_location','combined_sig','log10_ksby','ks_stat','n_events','direction','telcent','selection','code','TS or onco-like','genes']
        df.insert(0,'stuff',leftheaders + ([np.nan]*(df.shape[0]-len(leftheaders))))
        df.to_csv(os.path.join(folder,'arms','cols',a+'_BISCUT_results_cols_'+folder.split('_')[-1]+'.txt'),sep='\t',index=False,header=False)


def all_processing(date, ci, qval_thres, genes, arms, genelocs_file, abslocs_file):
    folder='results_'+date+'_'+str(ci)
    kspval=folder+'/KS_pvalues_'+date+'_'+str(ci)+'.txt'

    make_table_results(folder,kspval, qval_thres)
    calc_overlaps(folder, genelocs_file)
    lop = make_column_results(folder, qval_thres)
    process_for_ggplot(folder, qval_thres, abslocs_file)
    process_for_ggplot_jagged(folder, qval_thres, abslocs_file)
    #make_all_peaks(folder, qval_thres)
    filter_BISCUT_knowngenes(folder, genes)
    filter_BISCUT_arms(folder,arms)
    extract_cols(folder,lop,arms)
