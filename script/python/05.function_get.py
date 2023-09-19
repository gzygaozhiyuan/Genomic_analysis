import numpy as np
import pandas as pd
import sys, os
import argparse,time

args = argparse.ArgumentParser('For get each gene function')

args.add_argument('-gp', '--gene_phe', help='Gene to Phenotype file', required=True, type=str)
args.add_argument('-o', '--output', help='Output file', default='output.txt', type=str)

args = args.parse_args()

gene_phe = args.gene_phe
output = args.output
# bg_gene = args.bcgene

# 获取感兴趣基因集
gene_list = []
f1 = open(gene_phe,'r')
for l in f1:
    gene_list.append(l.strip().split('\t')[0])

file_path = '/public3/home/scg9832/script/01.use/02.func'

# GO
bc_gene_go = pd.read_csv(f'{file_path}/Osa_GO.txt', sep='\t', header=0)
bc_gene_go['Gene set name'] = bc_gene_go['Gene set name'].apply(lambda x:x.replace('_',' ').title())
bc_gene_go['GO'] = bc_gene_go['Gene set name'] + '|' + bc_gene_go['Category']
bc_gene_go = bc_gene_go.loc[:,['Gene','GO']]
gene_go_need = bc_gene_go[bc_gene_go['Gene'].isin(gene_list)]

# Gfam
bc_gene_gfam = pd.read_csv(f'{file_path}/Osa_GFam.txt', sep='\t', header=0)
bc_gene_gfam['Gene set name'] = bc_gene_gfam['Gene set name'].apply(lambda x:x.replace('_',' ').title())
bc_gene_gfam['GFam'] = bc_gene_gfam['Gene set name'] + '|' + bc_gene_gfam['Category']
bc_gene_gfam = bc_gene_gfam.loc[:,['Gene','GFam']]
gene_gfam_need = bc_gene_gfam[bc_gene_gfam['Gene'].isin(gene_list)]

# KEGG
bc_gene_kegg = pd.read_csv(f'{file_path}/Osa_ALL.KEGG.txt', sep='\t', header=0)
bc_gene_kegg['Gene set name'] = bc_gene_kegg['Gene set name'].apply(lambda x:x.replace('KEGG_',''))
bc_gene_kegg['Gene set name'] = bc_gene_kegg['Gene set name'].apply(lambda x:x.replace('_',' ').title())
bc_gene_kegg['KEGG'] = bc_gene_kegg['Gene set name']
bc_gene_kegg = bc_gene_kegg.loc[:,['Gene','KEGG']]
gene_kegg_need = bc_gene_kegg[bc_gene_kegg['Gene'].isin(gene_list)]

# Mapman
bc_gene_map = pd.read_csv(f'{file_path}/Osa_ALL.Map.txt', sep='\t',header=0)
bc_gene_map['Mapman'] = bc_gene_map['Gene set name'].apply(lambda x:x.replace('_',' ').title())
bc_gene_map = bc_gene_map.loc[:,['Gene','Mapman']]
gene_map_need = bc_gene_map[bc_gene_map['Gene'].isin(gene_list)]

# Cyc
bc_gene_cyc = pd.read_csv(f'{file_path}/Osa_ALL.Cyc.txt', sep='\t', header=0)
bc_gene_cyc['Cyc'] = bc_gene_cyc['Gene set name'].apply(lambda x:x.replace('_',' ').title())
bc_gene_cyc = bc_gene_cyc.loc[:,['Gene','Cyc']]
gene_cyc_need = bc_gene_cyc[bc_gene_cyc['Gene'].isin(gene_list)]

# rice_fn 
rice_fn = pd.read_csv(f'{file_path}/rice_fn.txt',sep='\t',header=None,names=['chr','st','sp','Gene','ms'])

# All
all_func = pd.merge(gene_go_need,gene_gfam_need,on='Gene',how='outer')
all_func = pd.merge(all_func,gene_kegg_need,on='Gene',how='outer')
all_func = pd.merge(all_func,gene_map_need,on='Gene',how='outer')
all_func = pd.merge(all_func,gene_cyc_need,on='Gene',how='outer')
all_func = pd.merge(all_func,rice_fn,on='Gene',how='outer')

all_func.to_csv('./func_TMP.txt',header=True,index=False,na_rep='NA',sep='\t')

# 重新排列一下，每个基因一行信息
f1 = open('./func_TMP.txt','r')
dd = {}
for l in f1:
    if not l.startswith('Gene'):
        ls = l.strip().split('\t')
        gene = ls[0]
        go = ls[1]
        gfam = ls[2]
        kegg = ls[3]
        map = ls[4]
        cyc = ls[5]
        ms = ls[-1]
        if gene not in dd:
            dd[gene] = {'go':[],'gfam':[],'kegg':[],'map':[],'cyc':[],'ms':[]}
        if go not in dd[gene]['go']:
            dd[gene]['go'].append(go)
        if gfam not in dd[gene]['gfam']:
            dd[gene]['gfam'].append(gfam)
        if kegg not in dd[gene]['kegg']:
            dd[gene]['kegg'].append(kegg)
        if map not in dd[gene]['map']:
            dd[gene]['map'].append(map)
        if cyc not in dd[gene]['cyc']:
            dd[gene]['cyc'].append(cyc)
        if ms not in dd[gene]['ms']:
            dd[gene]['ms'].append(ms)
    else:
        pass

with open(output,'w')as fd:
    fd.write(f"GENE\tGO\tGfam\tKEGG\tMapman\tCyc\tMs\n")
f1 = open(gene_phe,'r')
for l in f1:
    gene = l.strip().split('\t')[0]
    if gene in dd:
        go_ms = ';'.join(dd[gene]['go'])
        gfam_ms = ';'.join(dd[gene]['gfam'])
        kegg_ms = ';'.join(dd[gene]['kegg'])
        map_ms = ';'.join(dd[gene]['map'])
        cyc_ms = ';'.join(dd[gene]['cyc'])
        ms_ms = ';'.join(dd[gene]['ms'])
        with open(output,'a')as fd:
            fd.write(f"{gene}\t{go_ms}\t{gfam_ms}\t{kegg_ms}\t{map_ms}\t{cyc_ms}\t{ms_ms}\n")
    else:
        with open(output,'a')as fd:
            fd.write(f"{gene}\tNA\tNA\tNA\tNA\tNA\tNA\n")

os.system('rm -rf ./func_TMP.txt')