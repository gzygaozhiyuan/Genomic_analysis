import os,argparse
import numpy as np
import pandas as pd

args = argparse.ArgumentParser('Get gene region vcf')

args.add_argument('-g', '--gene_ls', help='Gene list file', required=True, type=str)
args.add_argument('-i', '--vcf', help='VCF file', required=True, type=str)
args.add_argument('-s', '--sample_ls', help='Sample list file', required=True, type=str)
args.add_argument('-t', '--sub_type', help='Sub type', required=True, type=str)
args.add_argument('-d', '--dis', help='Gene upstream and downstream regions [kb]', default=10, type=int)
args.add_argument('-p', '--prefix', help='Output prefix file', default='output', type=str)

args = args.parse_args()

gene_ls_file = args.gene_ls
vcf_file = args.vcf
distence = args.dis
out_prefix = args.prefix
sample_ls_file = args.sample_ls
sub_type = args.sub_type

# 定义一个函数，根据type对pos进行缩放
def scale_pos(row,df,gene_st,gene_sp,distence):
    if row['Type'] == 'Gene':
        return (row['POS']-gene_st)/(gene_sp-gene_st)
    elif row['Type'] == 'Up':
        return ((row['POS']-gene_st+distence*1000)/(distence*1000)-1)*distence
    elif row['Type'] == 'Down':
        return (row['POS']-gene_sp)/(distence*1000)*distence+1
    else:
        return row['POS']

# 将sub_type转换为小写
sub_lower = sub_type.lower()

# ['GJ','XI','cA','cB'],['Japonica','Indica','Aus','Basmati']
sub_d = {'gj':'Japonica','xi':'Indica','ca':'Aus','cb':'Basmati'}

gene_d = {}
rice_fn = open('/public3/home/scg9832/script/01.use/02.func/rice_fn.txt','r')
for l in rice_fn:
    ls = l.strip().split('\t')
    gene_d[ls[-2]] = [int(ls[1]),int(ls[2]),ls[0][3:]]

# 计算整个基因组的pi
os.system(f"vcftools --vcf {vcf_file} --keep {sample_ls_file} --site-pi --out all_{sub_type}")

# 读取结果文件
result = pd.read_csv(f"all_{sub_type}.sites.pi",sep='\t',header=0)

# 整理每个基因的绘图文件
f1 = open(gene_ls_file,'r')
for l in f1:
    gene = l.strip().split('\t')[0]
    gene_st = gene_d[gene][0]
    gene_sp = gene_d[gene][1]
    gene_region_st = gene_st - int(distence)*1000
    gene_region_sp = gene_sp + int(distence)*1000
    gene_chr = int(gene_d[gene][-1])
    # 取出gene_chr的行
    result_need = result[result['CHROM']==gene_chr]
    # 取出特定区域的行
    result_need = result_need[(result_need['POS']>gene_region_st) & (result_need['POS']<gene_region_sp)]
    # 添加一列Sub
    result_need['Sub'] = sub_d[sub_lower]
    # 添加一列Type
    result_need['Type'] = result_need['POS'].apply(lambda x: 'Up' if x < gene_st else 'Down' if x > gene_sp else 'Gene')
    # 对每类Type的POS列进行缩放
    result_need['x'] = result_need.apply(lambda row: scale_pos(row, result_need,gene_st,gene_sp,distence), axis=1)
    result_need['gene_st'] = gene_st
    result_need['gene_sp'] = gene_sp
    # 保存文件
    result_need.to_csv(f"{gene}_{sub_type}.txt",sep='\t',header=True,index=False)

    