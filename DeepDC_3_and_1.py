import os
import subprocess
import time,base64
import argparse
import itertools
import pandas as pd
import numpy as np
import secrets,string
import pyBigWig
from pyfaidx import Fasta
from tqdm import tqdm


t0 = time.time()
trans_base_dict={'A':'T','C':'G','G':'C','T':'A','N':'N'}
chr_list = [f'chr{i}' for i in range(1,23)]+['chrX']

def reverse_DNA(seq):
    return ''.join([trans_base_dict[c] for c in seq])[::-1]

def gc_content(seq_list):
    content = []
    for seq in seq_list:
        if len(seq)==0:
            content.append(0.0)
        else:
            content.append( (seq.upper().count('C')+seq.upper().count('G')) / len(seq) )
    return content

def find_contextual_single(data,fa=None,col_name = 'sequence',expand_flank=30,upstream=4,downstream=6):

    p1,p2=[],[];e1,e2=[],[]
    adj_s1,adj_s2=[],[]
    merge_table=data
    
    for _,row in merge_table.iterrows():
    
        s1 = row[col_name].upper().replace(' ','')

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

        if l1<0:
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
        bw=pyBigWig.open(f"./epigenome_ref/{cellline}/{epi}.bigWig")
    elif genome=='hg19':
        bw=pyBigWig.open(f"./epigenome_ref/{cellline}/hg19/{epi}.bigWig")
    
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

def generate_random_string(length=8):
    characters = string.ascii_letters + string.digits
    random_string = ''.join(secrets.choice(characters) for _ in range(length))
    return random_string

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
parser.add_argument(
        "--cuda",required=False,default=None,help="Whether use GPU device. Default is not. If you have GPU device, please input as an index like cuda:0"
    )
                
assert parser.parse_args().genome in ['hg19','hg38'], 'Genome illegal'
assert parser.parse_args().cellline in ['K562','A549'], 'Cell line illegal'

cmd = f"-J -BED -G {parser.parse_args().genome} -filterGCmin 30 -filterGCmax 70 -t WHOLE -n N -R 4 " + \
       "-3 PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60 " + \
       "-A 290 -DF 300 -a 20 -T 1 -g 20 -scoringMethod DOENCH_2014 -f NN -v 3 -M NGG -BB AGGCTAGTCCGT"

try:
    region = parser.parse_args().region
    #print(region)
    chromo,start,end = region.split(':')[0],region.split(':')[1].split('-')[0],region.split('-')[1]
    #print(chromo,start,end)
    assert chromo in [f'chr{i}' for i in range(1,23)]+['chrX','chrY'], 'Chromosome is illegal'
    assert not start is None, 'Region should contain start'
    assert not end is None, 'Region should contain end'
    assert start!=end, 'Region start could not equal to end'
except:
    assert 1==0, 'Region is illegal'
    
objectname = base64.urlsafe_b64encode(int(time.time() * 1000).to_bytes(8, byteorder='big')).decode().rstrip("=")+'_'+generate_random_string()
os.mkdir(f'./result/{objectname}')
cmd_user   = " ".join( [f"--targets {parser.parse_args().region} -o {parser.parse_args().oc}",cmd] )
outputfile = parser.parse_args().oc + f'/{objectname}/_result_.tsv'

with open(outputfile, 'w') as f:
    subprocess.run(
        ["python", "./chopchop/chopchop_py3.py"] + cmd_user.split(),
        stdout=f,
        text=True
    )

clean_path(parser.parse_args().oc,'.tsv')

table = pd.read_csv(f'{parser.parse_args().oc}/{objectname}/_result_.tsv',sep='\t')
table = table[(table['MM0']==0) & (table['MM1']<=1)].reset_index(drop=True)

table[['chr','start']] = table['Genomic location'].str.split(':',expand=True).values
table['start'] = table['start'].astype(int);table['end'] = table['start']+23

with Fasta(f"/cluster2/huanglab/liquan/genome_ref/{parser.parse_args().genome}.fa") as fa:
    _,table['seq_adj'],_ = find_contextual_single(table,col_name='Target sequence',fa=fa,upstream=4,downstream=3)
    
if parser.parse_args().cuda is None:
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
elif 'cuda:' in parser.parse_args().cuda:
    os.environ["CUDA_VISIBLE_DEVICES"] = 'cuda:' in parser.parse_args().cuda[5:]
else:
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    
import sys
import DeepSpCas9_main as DeepSpCas9
import DeepCRISPR_main as DeepCRISPR
import DistillatedRuleSet2 as Ruleset2
import DistillatedCRISPRedict as CRISPRedict
import HyenaDNA_main as Hyena
import xgboost as xgb
    
table['DeepSpCas9_score'] = DeepSpCas9.predict_sequence(table['seq_adj'])
table['DeepCRISPR_score'] = DeepCRISPR.predict_sequence([seq[4:27] for seq in table['seq_adj']])
table['Ruleset2_score']   = Ruleset2.predict(table['seq_adj'],model_path='./DistillatedRuleSet2.pth')
table['CRISPRedit_score'] = CRISPRedict.predict(table['seq_adj'],model_path='./DistillatedCRISPRedict.pth')
table = table[['chr','start','end','Strand','seq_adj','DeepSpCas9_score','DeepCRISPR_score','Ruleset2_score','CRISPRedit_score']]

import itertools

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

#pair_df_filtered.to_csv(f'{parser.parse_args().oc}result.pair.tsv',sep='\t',index=False)

if parser.parse_args().genome=='hg19':
    with Fasta("./genome_ref/hg19.fa") as hg19:
        pair_df_filtered['Median_sequence']=pair_df_filtered.apply( lambda x: hg19[x['chr']][ x['start']:x['end'] ].seq.upper() , axis=1 )
elif parser.parse_args().genome=='hg38':
    with Fasta("./genome_ref/hg38.fa") as hg38:
        pair_df_filtered['Median_sequence']=pair_df_filtered.apply( lambda x: hg38[x['chr']][ x['start']:x['end'] ].seq.upper() , axis=1 )
        
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
    
try:
    pair_df_filtered['DNase']  =gain_series_TPM( pair_df_filtered, epi='DNase'  , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
    pair_df_filtered['ATAC']   =gain_series_TPM( pair_df_filtered, epi='ATAC'   , cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
    pair_df_filtered['H3K27ac']=gain_series_TPM( pair_df_filtered, epi='H3K27ac', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)
    pair_df_filtered['H3K4me3']=gain_series_TPM( pair_df_filtered, epi='H3K4me3', cellline=parser.parse_args().cellline, genome=parser.parse_args().genome, lo=None)      
except:
    assert 1==0, 'Error in epigenetic feature capture'

xgb_model = xgb.XGBRegressor()
xgb_model.load_model('./xgboost_model.model')

ipdata=pair_df_filtered
ipdata['DeepSpCas9_harmony'] = 1/(1/(ipdata['DeepSpCas9_s1']+1)+1/(ipdata['DeepSpCas9_s2']+1))
ipdata['DeepCRISPR_harmony'] = 1/(1/(ipdata['DeepCRISPR_s1']+1)+1/(ipdata['DeepCRISPR_s2']+1))
ipdata['CRISPRedict_harmony']= 1/(1/(ipdata['CRISPRedit_s1']+1)+1/(ipdata['CRISPRedit_s2']+1))
ipdata['Ruleset2_harmony']   = 1/(1/(ipdata['Ruleset2_s1']+1)+1/(ipdata['Ruleset2_s2']+1))
ipdata = ipdata.sort_values(by=['DeepSpCas9_harmony','DeepCRISPR_harmony','CRISPRedict_harmony','Ruleset2_harmony','DNase','ATAC','H3K27ac','H3K4me3'],ascending=False)
ipdata['hyena_score']=np.array([Hyena.hyena_inference(s[:16384]) for s in list(ipdata['Median_sequence'].fillna('N'))])
ipdata = ipdata.fillna(0.0)

features =  ['DNase','ATAC','H3K27ac','H3K4me3','GC','length',
             'DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score']

result = xgb_model.predict(ipdata[features].values)
ipdata['DeepDC_score'] = result.flatten()
result = ipdata
result['DeepDC_score'] = (result['DeepDC_score']-result['DeepDC_score'].min())/(result['DeepDC_score'].max()-result['DeepDC_score'].min())
result[['DNase','ATAC','H3K27ac','H3K4me3','length']] = result[['DNase','ATAC','H3K27ac','H3K4me3','length']].astype(int)

result[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']] = \
  result[['GC','DeepSpCas9_harmony','DeepCRISPR_harmony','Ruleset2_harmony','CRISPRedict_harmony','hyena_score','DeepDC_score']].astype(float).round(2)

result[['chr','start','end','seq1','seq2']+features+['DeepDC_score']].to_csv(f'{parser.parse_args().oc}/{objectname}/result.score.tsv',sep='\t',index=False)
t1 = time.time()
print(f'Using time:{t1-t0:.2f} s')
