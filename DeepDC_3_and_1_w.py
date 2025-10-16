import os
import subprocess
import argparse
import itertools
import pandas as pd
import numpy as np
import pyBigWig
from pyfaidx import Fasta
from tqdm import tqdm

import utils.DeepSpCas9_main as DeepSpCas9
import utils.DeepCRISPR_main as DeepCRISPR
import utils.DistillatedRuleSet2 as Ruleset2
import utils.DistillatedCRISPRedict as CRISPRedict
import utils.HyenaDNA_main as Hyena
import xgboost as xgb


trans_base_dict={'A':'T','C':'G','G':'C','T':'A','N':'N'}
chr_list = [f'chr{i}' for i in range(1,23)]+['chrX','chrY']
cell_line_withepi = ['K562','A549']
multi_loci = False

def reverse_DNA(seq):
    ### 用于反转DNA
    return ''.join([trans_base_dict[c] for c in seq])[::-1]

def gc_content(seq_list):
    ### 计算GC含量
    content = []
    for seq in seq_list:
        if len(seq)==0:
            content.append(0.0)
        else:
            content.append( (seq.upper().count('C')+seq.upper().count('G')) / len(seq) )
    return content

def find_contextual_single(data,fa=None,col_name = 'sequence',expand_flank=30,upstream=4,downstream=6):
    
    # 自动寻找上下文并且把补全上下文之后的序列整理成有30bp的格式
    ## 因为参考基因组并不会随着细胞系的改变而改变，所以这里可以这么做

    p1,p2=[],[];e1,e2=[],[]
    adj_s1,adj_s2=[],[]
    merge_table=data
    
    for _,row in merge_table.iterrows():
    
        s1 = row[col_name].upper().replace(' ','') # 因为输出是ChopChop的,非常稳定，就是20 bp protospacer + 3 bp PAM

        try:
            c,s,e = row['chr'],int(row['start']),int(row['end'])
            if not c in chr_list:
                p1.append('-1');e1.append('N'*total_len);adj_s1.append(-1)
                continue
        except:
            p1.append('-1');e1.append('N'*upstream+s1+'N'*downstream);adj_s1.append(-1)
            continue

        sequence=fa[c][s-expand_flank:e+expand_flank].seq.upper()

        l1=sequence.find(s1)

        if l1<0: # 没有抓取到位置信息的情况,自动认为是负向的
            l1=sequence.find(reverse_DNA(s1))
            p1.append('0');e1.append( reverse_DNA( fa[c][s-expand_flank+l1-downstream:s-expand_flank+l1+23+upstream].seq.upper() ) )
        else:
            p1.append('1');e1.append( fa[c][s-expand_flank+l1-upstream:s-expand_flank+l1+23+downstream].seq.upper() )
        adj_s1.append(s-expand_flank+l1)
    
    return p1,e1,adj_s1

def gain_epivalue_TPM(region,bw=None):
    
    chromo=region[0]
    start=min(region[1],region[2])
    end=max(region[1],region[2])
    
    if start==end:
        end+=1
        
    delta=end-start
    
    try:
        values = bw.stats(chromo, start, end, type='sum')[0]
        tpm = (values/delta)/(bw.header()['sumData']/bw.header()['nBasesCovered'])*1000
    except Exception as e:
        tpm = 0
    return tpm

def gain_series_TPM(data,epi='DNase',cellline='A549',genome='hg38',lo=None):
    
    epi_value_series=[]
    
    if genome=='hg38': 
        bw=pyBigWig.open(f"epigenome_ref/{cellline}/{epi}.bigWig")
    elif genome=='hg19':
        bw=pyBigWig.open(f"epigenome_ref/{cellline}/hg19/{epi}.bigWig")
    
    for i in range(len(data)):
        c,s,e = data.at[i,'chr'],data.at[i,'start'],data.at[i,'end']
        if s==e:
            e+=1
        else:
            s,e = min(s,e),max(s,e)
        if not lo is None:
            c,s,e = convert_loci([c,s,e],lo=lo)
        
        try:
            epi_value_series.append(gain_epivalue_TPM([c,s,e],bw=bw))
        except:
            epi_value_series.append(0.0)
    
    bw.close()
    
    return epi_value_series

def clean_path(directory, extension):

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        
        if os.path.isfile(file_path):
            if not filename.endswith(extension):
                os.remove(file_path)
                
parser = argparse.ArgumentParser(description="Run ChopChop")
parser.add_argument(
        "--region",required=True,help="The target region in the format 'chr:start-end', e.g., 'chr11:5292000-5293000'."
    )
parser.add_argument(
        "--genome",required=True,help="The genome uesd. Choose in hg19,hg38"
    )
parser.add_argument(
        "--cellline",required=False,default="K562",help="The cell line uesd. Choose in K562,A549"
    )
parser.add_argument(
        "--oc",required=False,default="./result/",help="The output directory. Default is ./result/"
    )
                
assert parser.parse_args().genome in ['hg19','hg38'], 'Genome illegal'
assert parser.parse_args().cellline in ['K562','A549'], 'Cell line illegal'

cmd = f"-J -BED -G {parser.parse_args().genome} -filterGCmin 30 -filterGCmax 70 -t WHOLE -n N -R 4 " + \
       "-3 PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60 " + \
       "-A 290 -DF 300 -a 20 -T 1 -g 20 -scoringMethod DOENCH_2014 -f NN -v 3 -M NGG -BB AGGCTAGTCCGT"

try:
    region = parser.parse_args().region
    #print(region)
    if ',' in region:
        multi_loci = True # 开启多位点标签
        region = region.split(',')[:20] # 最多支持20个位点连续预测
        #print(region)
        os.makedirs(parser.parse_args().oc, exist_ok=True)
        for i,r in enumerate(region):
            chromo,start,end = r.split(':')[0],r.split(':')[1].split('-')[0],r.split('-')[1]
            if not chromo in [f'chr{i}' for i in range(1,23)]+['chrX','chrY']: continue
            if start is None: continue
            if end is None: continue
            if start==end: continue
            cmd_user   = " ".join( [f"--targets {chromo}:{start}-{end} -o {parser.parse_args().oc}",cmd] )
            outputfile = parser.parse_args().oc + f'/_result_.{i}.tsv'
            with open(outputfile, 'w') as f:
                subprocess.run(
                    ["python", "./chopchop/chopchop_py3.py"] + cmd_user.split(),
                    stdout=f,
                    text=True
                )
            clean_path(parser.parse_args().oc,'.tsv')
        merge_df = pd.DataFrame()
        for j in range(0,i+1):
            sub_df = pd.read_csv(parser.parse_args().oc + f'/_result_.{j}.tsv',sep='\t',skiprows=1)
            sub_df['element'] = region[j]
            merge_df = pd.concat( [ merge_df, sub_df] )
        merge_df.to_csv(parser.parse_args().oc + f'/_result_.tsv',sep='\t',index=False)
        table = pd.read_csv(f'{parser.parse_args().oc}/_result_.tsv',sep='\t')
        table = table[(table['MM0']==0) & (table['MM1']<=1)].reset_index(drop=True)
        
        table[['chr','start']] = table['Genomic location'].str.split(':',expand=True).values
        table['start'] = table['start'].astype(int);table['end'] = table['start']+23
        with Fasta(f"./genome_ref/{parser.parse_args().genome}.fa") as fa:
            _,table['seq_adj'],_ = find_contextual_single(table,col_name='Target sequence',fa=fa,upstream=4,downstream=3)

        table['DeepSpCas9_score'] = DeepSpCas9.predict_sequence(table['seq_adj'])
        table['DeepCRISPR_score'] = DeepCRISPR.predict_sequence([seq[4:27] for seq in table['seq_adj']])
        table['Ruleset2_score']   = Ruleset2.predict(table['seq_adj'],model_path='DistillatedRuleSet2.pth')
        table['CRISPRedit_score'] = CRISPRedict.predict(table['seq_adj'],model_path='DistillatedCRISPRedict.pth')
        table = table[['chr','start','end','Strand','seq_adj','DeepSpCas9_score','DeepCRISPR_score','Ruleset2_score','CRISPRedit_score','element']]
        
        merge_pair_df = pd.DataFrame()
        for element,group in table.groupby(['element']):
            permutations = list(itertools.permutations(group.index, 2))
            pair_df = pd.DataFrame({
                'chr':           [group.loc[i, 'chr'] for i, j in permutations],
                'Cutsite_1':     [group.loc[i, 'start'] for i, j in permutations],
                'Cutsite_2':     [group.loc[j, 'start'] for i, j in permutations],
                'seq1':          [group.loc[i, 'seq_adj'] for i, j in permutations],
                'seq2':          [group.loc[j, 'seq_adj'] for i, j in permutations],
                'DeepSpCas9_s1': [group.loc[i, 'DeepSpCas9_score'] for i, j in permutations],
                'DeepCRISPR_s1': [group.loc[i, 'DeepCRISPR_score'] for i, j in permutations],
                'Ruleset2_s1':   [group.loc[i, 'Ruleset2_score'] for i, j in permutations],
                'CRISPRedit_s1': [group.loc[i, 'CRISPRedit_score'] for i, j in permutations],
                'DeepSpCas9_s2': [group.loc[j, 'DeepSpCas9_score'] for i, j in permutations],
                'DeepCRISPR_s2': [group.loc[j, 'DeepCRISPR_score'] for i, j in permutations],
                'Ruleset2_s2':   [group.loc[j, 'Ruleset2_score'] for i, j in permutations],
                'CRISPRedit_s2': [group.loc[j, 'CRISPRedit_score'] for i, j in permutations]
            })
            pair_df['element'] = element[0]
            merge_pair_df = pd.concat([merge_pair_df,pair_df])
        pair_df = merge_pair_df
        pair_df['pair_length'] = pair_df['Cutsite_2'].astype(int) - pair_df['Cutsite_1'].astype(int)
        pair_df_filtered = pair_df[(pair_df['pair_length'] >= 50) & (pair_df['pair_length'] <= 200)].reset_index(drop=True)
        pair_df_filtered.columns=['chr','start','end','seq1','seq2',
                                  'DeepSpCas9_s1','DeepCRISPR_s1','Ruleset2_s1','CRISPRedit_s1',
                                  'DeepSpCas9_s2','DeepCRISPR_s2','Ruleset2_s2','CRISPRedit_s2','element','length']
        assert len(pair_df_filtered)>2, 'No enough sgRNA pair is found'

        with Fasta(f"./genome_ref/{parser.parse_args().genome}.fa") as fa:
            pair_df_filtered['Median_sequence']=pair_df_filtered.apply( lambda x: fa[x['chr']][ x['start']:x['end'] ].seq.upper() , axis=1 )

        pair_df_filtered['GC']=gc_content(pair_df_filtered['Median_sequence'])

        if not pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() == 0:
            pair_df_filtered['DeepSpCas9_s1'] = ( pair_df_filtered['DeepSpCas9_s1']-pair_df_filtered['DeepSpCas9_s1'].min() ) / ( pair_df_filtered['DeepSpCas9_s1'].max() - pair_df_filtered['DeepSpCas9_s1'].min() )*10
        if not pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() == 0:
            pair_df_filtered['DeepSpCas9_s2'] = ( pair_df_filtered['DeepSpCas9_s2']-pair_df_filtered['DeepSpCas9_s2'].min() ) / ( pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() )*10

        if not pair_df_filtered['DeepCRISPR_s1'].max() - pair_df_filtered['DeepCRISPR_s1'].min() == 0:
            pair_df_filtered['DeepCRISPR_s1'] = ( pair_df_filtered['DeepCRISPR_s1']-pair_df_filtered['DeepCRISPR_s1'].min() ) / ( pair_df_filtered['DeepCRISPR_s1'].max() - pair_df_filtered['DeepCRISPR_s1'].min() )*10
        if not pair_df_filtered['DeepCRISPR_s2'].max() - pair_df_filtered['DeepCRISPR_s2'].min() == 0:
            pair_df_filtered['DeepCRISPR_s2'] = ( pair_df_filtered['DeepCRISPR_s2']-pair_df_filtered['DeepCRISPR_s2'].min() ) / ( pair_df_filtered['DeepCRISPR_s2'].max() - pair_df_filtered['DeepCRISPR_s2'].min() )*10

        if not pair_df_filtered['CRISPRedit_s1'].max() - pair_df_filtered['CRISPRedit_s1'].min() == 0:
            pair_df_filtered['CRISPRedit_s1'] = ( pair_df_filtered['CRISPRedit_s1']-pair_df_filtered['CRISPRedit_s1'].min() ) / ( pair_df_filtered['CRISPRedit_s1'].max() - pair_df_filtered['CRISPRedit_s1'].min() )*10
        if not pair_df_filtered['CRISPRedit_s2'].max() - pair_df_filtered['CRISPRedit_s2'].min() == 0:
            pair_df_filtered['CRISPRedit_s2'] = ( pair_df_filtered['CRISPRedit_s2']-pair_df_filtered['CRISPRedit_s2'].min() ) / ( pair_df_filtered['CRISPRedit_s2'].max() - pair_df_filtered['CRISPRedit_s2'].min() )*10

        if not pair_df_filtered['Ruleset2_s1'].max() - pair_df_filtered['Ruleset2_s1'].min() == 0:
            pair_df_filtered['Ruleset2_s1'] = ( pair_df_filtered['Ruleset2_s1']-pair_df_filtered['Ruleset2_s1'].min() ) / ( pair_df_filtered['Ruleset2_s1'].max() - pair_df_filtered['Ruleset2_s1'].min() )*10
        if not pair_df_filtered['Ruleset2_s2'].max() - pair_df_filtered['Ruleset2_s2'].min() == 0:
            pair_df_filtered['Ruleset2_s2'] = ( pair_df_filtered['Ruleset2_s2']-pair_df_filtered['Ruleset2_s2'].min() ) / ( pair_df_filtered['Ruleset2_s2'].max() - pair_df_filtered['Ruleset2_s2'].min() )*10

        if parser.parse_args().cellline in cell_line_withepi: # 如果当前细胞系有表观数据
            try:
                pair_df_filtered['DNase']  =gain_series_TPM( pair_df_filtered, epi='DNase'  , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['ATAC']   =gain_series_TPM( pair_df_filtered, epi='ATAC'   , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['H3K27ac']=gain_series_TPM( pair_df_filtered, epi='H3K27ac', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['H3K4me3']=gain_series_TPM( pair_df_filtered, epi='H3K4me3', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None) 
            except:
                assert 1==0, 'Error in epigenetic feature capture'
        else: # 没有表观数据的情况下就全部置0，只用单端分数进行打分
            pair_df_filtered['DNase']   = 0
            pair_df_filtered['ATAC']    = 0
            pair_df_filtered['H3K27ac'] = 0
            pair_df_filtered['H3K4me3'] = 0    

        xgb_model = xgb.XGBRegressor()
        xgb_model.load_model('./xgboost_model.model')

        ipdata=pair_df_filtered
        ipdata['DeepSpCas9_harmony'] = 1/(1/(ipdata['DeepSpCas9_s1']+1)+1/(ipdata['DeepSpCas9_s2']+1))
        ipdata['DeepCRISPR_harmony'] = 1/(1/(ipdata['DeepCRISPR_s1']+1)+1/(ipdata['DeepCRISPR_s2']+1))
        ipdata['CRISPRedict_harmony']= 1/(1/(ipdata['CRISPRedit_s1']+1)+1/(ipdata['CRISPRedit_s2']+1))
        ipdata['Ruleset2_harmony']   = 1/(1/(ipdata['Ruleset2_s1']+1)+1/(ipdata['Ruleset2_s2']+1))
        ipdata = ipdata.sort_values(by=['DeepSpCas9_harmony','DeepCRISPR_harmony','CRISPRedict_harmony','Ruleset2_harmony','DNase','ATAC','H3K27ac','H3K4me3'],ascending=False)
        #ipdata = ipdata.iloc[:50]
        ipdata['hyena_score']=np.array([Hyena.hyena_inference(s[:16384]) for s in list(ipdata['Median_sequence'].fillna('N'))])
        ipdata = ipdata.fillna(0.0)

        features =  ['DNase','ATAC','H3K27ac','H3K4me3','GC','length',
                     'DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score']

        result = xgb_model.predict(ipdata[features].values)
        ipdata['DeepDC_score'] = result.flatten()
        result = [];profile_result = []
        for element,group in ipdata.groupby('element'):
            group['element'] = element[0]
            group['DeepDC_score'] = (group['DeepDC_score']-group['DeepDC_score'].min())/(group['DeepDC_score'].max()-group['DeepDC_score'].min())
            group = group.sort_values(by=['DeepDC_score','DeepSpCas9_harmony','Ruleset2_harmony'],ascending=False)
            group[['DNase','ATAC','H3K27ac','H3K4me3','length']] = group[['DNase','ATAC','H3K27ac','H3K4me3','length']].astype(int)
            group[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']] = \
                group[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']].astype(float).round(2)
            profile_result.append(group.iloc[:50])
            result.append(group)
        
        result = pd.concat(result)
        profile_result = pd.concat(profile_result)
        
        result[['chr','start','end','seq1','seq2']+features+['DeepDC_score','element']].to_csv(f'{parser.parse_args().oc}/result.score.tsv',sep='\t',index=False)
        profile_result[['chr','start','end','seq1','seq2']+features+['DeepDC_score']].to_csv(f'{parser.parse_args().oc}/result.profile.tsv',sep='\t',index=False)
    else:
        multi_loci = False # 关闭多位点标签,正常走完全程
        chromo,start,end = region.split(':')[0],region.split(':')[1].split('-')[0],region.split('-')[1]
        #print(chromo,start,end)
        assert chromo in [f'chr{i}' for i in range(1,23)]+['chrX','chrY'], 'Chromosome is illegal'
        assert not start is None, 'Region should contain start'
        assert not end is None, 'Region should contain end'
        assert start!=end, 'Region start could not equal to end'
        os.makedirs(parser.parse_args().oc, exist_ok=True)
        cmd_user   = " ".join( [f"--targets {parser.parse_args().region} -o {parser.parse_args().oc}",cmd] )
        outputfile = parser.parse_args().oc + f'/_result_.tsv'   
        with open(outputfile, 'w') as f:
            subprocess.run(
                ["python", "./chopchop/chopchop_py3.py"] + cmd_user.split(),
                stdout=f,
                text=True
            )
        
        clean_path(parser.parse_args().oc,'.tsv')

        table = pd.read_csv(f'{parser.parse_args().oc}/_result_.tsv',sep='\t',skiprows=1)
        table = table[(table['MM0']==0) & (table['MM1']<=1)].reset_index(drop=True)

        table[['chr','start']] = table['Genomic location'].str.split(':',expand=True).values
        table['start'] = table['start'].astype(int);table['end'] = table['start']+23
        with Fasta(f"./genome_ref/{parser.parse_args().genome}.fa") as fa:
            _,table['seq_adj'],_ = find_contextual_single(table,col_name='Target sequence',fa=fa,upstream=4,downstream=3)

        table['DeepSpCas9_score'] = DeepSpCas9.predict_sequence(table['seq_adj'])
        table['DeepCRISPR_score'] = DeepCRISPR.predict_sequence([seq[4:27] for seq in table['seq_adj']])
        table['Ruleset2_score']   = Ruleset2.predict(table['seq_adj'],model_path='models/DistillatedRuleSet2.pth')
        table['CRISPRedit_score'] = CRISPRedict.predict(table['seq_adj'],model_path='models/DistillatedCRISPRedict.pth')
        table = table[['chr','start','end','Strand','seq_adj','DeepSpCas9_score','DeepCRISPR_score','Ruleset2_score','CRISPRedit_score']]

        permutations = list(itertools.permutations(table.index, 2))
        pair_df = pd.DataFrame({
            'chr':           [table.loc[i, 'chr'] for i, j in permutations],
            'Cutsite_1':     [table.loc[i, 'start'] for i, j in permutations],
            'Cutsite_2':     [table.loc[j, 'start'] for i, j in permutations],
            'seq1':          [table.loc[i, 'seq_adj'] for i, j in permutations],
            'seq2':          [table.loc[j, 'seq_adj'] for i, j in permutations],
            'DeepSpCas9_s1': [table.loc[i, 'DeepSpCas9_score'] for i, j in permutations],
            'DeepCRISPR_s1': [table.loc[i, 'DeepCRISPR_score'] for i, j in permutations],
            'Ruleset2_s1':   [table.loc[i, 'Ruleset2_score'] for i, j in permutations],
            'CRISPRedit_s1': [table.loc[i, 'CRISPRedit_score'] for i, j in permutations],
            'DeepSpCas9_s2': [table.loc[j, 'DeepSpCas9_score'] for i, j in permutations],
            'DeepCRISPR_s2': [table.loc[j, 'DeepCRISPR_score'] for i, j in permutations],
            'Ruleset2_s2':   [table.loc[j, 'Ruleset2_score'] for i, j in permutations],
            'CRISPRedit_s2': [table.loc[j, 'CRISPRedit_score'] for i, j in permutations]
        })

        pair_df['pair_length'] = pair_df['Cutsite_2'].astype(int) - pair_df['Cutsite_1'].astype(int)
        pair_df_filtered = pair_df[(pair_df['pair_length'] >= 50) & (pair_df['pair_length'] <= 200)].reset_index(drop=True)
        pair_df_filtered.columns=['chr','start','end','seq1','seq2',
                                  'DeepSpCas9_s1','DeepCRISPR_s1','Ruleset2_s1','CRISPRedit_s1',
                                  'DeepSpCas9_s2','DeepCRISPR_s2','Ruleset2_s2','CRISPRedit_s2','length']
        assert len(pair_df_filtered)>2, 'No enough sgRNA pair is found'

        with Fasta(f"./genome_ref/{parser.parse_args().genome}.fa") as fa:
            pair_df_filtered['Median_sequence']=pair_df_filtered.apply( lambda x: fa[x['chr']][ x['start']:x['end'] ].seq.upper() , axis=1 )

        pair_df_filtered['GC']=gc_content(pair_df_filtered['Median_sequence'])

        if not pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() == 0:
            pair_df_filtered['DeepSpCas9_s1'] = ( pair_df_filtered['DeepSpCas9_s1']-pair_df_filtered['DeepSpCas9_s1'].min() ) / ( pair_df_filtered['DeepSpCas9_s1'].max() - pair_df_filtered['DeepSpCas9_s1'].min() )*10
        if not pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() == 0:
            pair_df_filtered['DeepSpCas9_s2'] = ( pair_df_filtered['DeepSpCas9_s2']-pair_df_filtered['DeepSpCas9_s2'].min() ) / ( pair_df_filtered['DeepSpCas9_s2'].max() - pair_df_filtered['DeepSpCas9_s2'].min() )*10

        if not pair_df_filtered['DeepCRISPR_s1'].max() - pair_df_filtered['DeepCRISPR_s1'].min() == 0:
            pair_df_filtered['DeepCRISPR_s1'] = ( pair_df_filtered['DeepCRISPR_s1']-pair_df_filtered['DeepCRISPR_s1'].min() ) / ( pair_df_filtered['DeepCRISPR_s1'].max() - pair_df_filtered['DeepCRISPR_s1'].min() )*10
        if not pair_df_filtered['DeepCRISPR_s2'].max() - pair_df_filtered['DeepCRISPR_s2'].min() == 0:
            pair_df_filtered['DeepCRISPR_s2'] = ( pair_df_filtered['DeepCRISPR_s2']-pair_df_filtered['DeepCRISPR_s2'].min() ) / ( pair_df_filtered['DeepCRISPR_s2'].max() - pair_df_filtered['DeepCRISPR_s2'].min() )*10

        if not pair_df_filtered['CRISPRedit_s1'].max() - pair_df_filtered['CRISPRedit_s1'].min() == 0:
            pair_df_filtered['CRISPRedit_s1'] = ( pair_df_filtered['CRISPRedit_s1']-pair_df_filtered['CRISPRedit_s1'].min() ) / ( pair_df_filtered['CRISPRedit_s1'].max() - pair_df_filtered['CRISPRedit_s1'].min() )*10
        if not pair_df_filtered['CRISPRedit_s2'].max() - pair_df_filtered['CRISPRedit_s2'].min() == 0:
            pair_df_filtered['CRISPRedit_s2'] = ( pair_df_filtered['CRISPRedit_s2']-pair_df_filtered['CRISPRedit_s2'].min() ) / ( pair_df_filtered['CRISPRedit_s2'].max() - pair_df_filtered['CRISPRedit_s2'].min() )*10

        if not pair_df_filtered['Ruleset2_s1'].max() - pair_df_filtered['Ruleset2_s1'].min() == 0:
            pair_df_filtered['Ruleset2_s1'] = ( pair_df_filtered['Ruleset2_s1']-pair_df_filtered['Ruleset2_s1'].min() ) / ( pair_df_filtered['Ruleset2_s1'].max() - pair_df_filtered['Ruleset2_s1'].min() )*10
        if not pair_df_filtered['Ruleset2_s2'].max() - pair_df_filtered['Ruleset2_s2'].min() == 0:
            pair_df_filtered['Ruleset2_s2'] = ( pair_df_filtered['Ruleset2_s2']-pair_df_filtered['Ruleset2_s2'].min() ) / ( pair_df_filtered['Ruleset2_s2'].max() - pair_df_filtered['Ruleset2_s2'].min() )*10

        if parser.parse_args().cellline in cell_line_withepi: # 如果当前细胞系有表观数据
            try:
                pair_df_filtered['DNase']  =gain_series_TPM( pair_df_filtered, epi='DNase'  , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['ATAC']   =gain_series_TPM( pair_df_filtered, epi='ATAC'   , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['H3K27ac']=gain_series_TPM( pair_df_filtered, epi='H3K27ac', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
                pair_df_filtered['H3K4me3']=gain_series_TPM( pair_df_filtered, epi='H3K4me3', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None) 
            except:
                assert 1==0, 'Error in epigenetic feature capture'
        else: # 没有表观数据的情况下就全部置0，只用单端分数进行打分
            pair_df_filtered['DNase']   = 0
            pair_df_filtered['ATAC']    = 0
            pair_df_filtered['H3K27ac'] = 0
            pair_df_filtered['H3K4me3'] = 0    

        xgb_model = xgb.XGBRegressor()
        xgb_model.load_model('models/xgboost_model.model')

        ipdata=pair_df_filtered
        ipdata['DeepSpCas9_harmony'] = 1/(1/(ipdata['DeepSpCas9_s1']+1)+1/(ipdata['DeepSpCas9_s2']+1))
        ipdata['DeepCRISPR_harmony'] = 1/(1/(ipdata['DeepCRISPR_s1']+1)+1/(ipdata['DeepCRISPR_s2']+1))
        ipdata['CRISPRedict_harmony']= 1/(1/(ipdata['CRISPRedit_s1']+1)+1/(ipdata['CRISPRedit_s2']+1))
        ipdata['Ruleset2_harmony']   = 1/(1/(ipdata['Ruleset2_s1']+1)+1/(ipdata['Ruleset2_s2']+1))
        ipdata = ipdata.sort_values(by=['DeepSpCas9_harmony','DeepCRISPR_harmony','CRISPRedict_harmony','Ruleset2_harmony','DNase','ATAC','H3K27ac','H3K4me3'],ascending=False)
        #ipdata = ipdata.iloc[:50]
        ipdata['hyena_score']=np.array([Hyena.hyena_inference(s[:16384]) for s in list(ipdata['Median_sequence'].fillna('N'))])
        ipdata = ipdata.fillna(0.0)

        features =  ['DNase','ATAC','H3K27ac','H3K4me3','GC','length',
                     'DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score']

        result = xgb_model.predict(ipdata[features].values)
        ipdata['DeepDC_score'] = result.flatten()
        result = ipdata.sort_values(by=['DeepDC_score','DeepSpCas9_harmony','Ruleset2_harmony'],ascending=False)
        result['DeepDC_score'] = (result['DeepDC_score']-result['DeepDC_score'].min())/(result['DeepDC_score'].max()-result['DeepDC_score'].min())
        result[['DNase','ATAC','H3K27ac','H3K4me3','length']] = result[['DNase','ATAC','H3K27ac','H3K4me3','length']].astype(int)

        result[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']] = \
          result[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']].astype(float).round(2)

        result[['chr','start','end','seq1','seq2']+features+['DeepDC_score']].to_csv(f'{parser.parse_args().oc}/result.score.tsv',sep='\t',index=False)
        result[['chr','start','end','seq1','seq2']+features+['DeepDC_score']].iloc[:50].to_csv(f'{parser.parse_args().oc}/result.profile.tsv',sep='\t',index=False)
except:
    assert 1==0, 'Region is illegal'