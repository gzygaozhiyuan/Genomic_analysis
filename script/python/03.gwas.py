# !/usr/bin/python3
# -*- coding: UTF-8 -*-
# by zhiyuan

import argparse,os

parser = argparse.ArgumentParser(description='One click GWAS by zhiyuan only in xu_lab')

parser.add_argument('-p','--phenotype',type=str,required=True,help='Phenotype file')
parser.add_argument('-s','--sample',type=str,required=True,help='Sample list file')
parser.add_argument('-m','--maf',type=float,required=True,help='maf')
parser.add_argument('-g','--geno',type=float,required=True,help='geno')
parser.add_argument('-sw','--software',type=str,required=True,choices=['emmax','gemma'],help='Software type for gwas,like:emmax or gemma')
parser.add_argument('-pi','--phenotypeid',type=str,required=True,help='Phenotype id')
parser.add_argument('-o','--output',type=str,default='outfile.txt',help='Output file name')
args = parser.parse_args()

phenotype = args.phenotype
software = args.software
sample_ls = args.sample
out = args.output
maf = args.maf
geno = args.geno
phe_id = args.phenotypeid

# ld pruning file
ld1 = f'plink --bfile /public3/home/scg9832/common_data/01.genomic/4.8M/base_filtered_v0.7 --keep {sample_ls} --geno {geno} --maf {maf} --indep-pairwise 50 10 0.2 --out LD'
ld2 = f'plink --bfile /public3/home/scg9832/common_data/01.genomic/4.8M/base_filtered_v0.7 --keep {sample_ls} --geno {geno} --maf {maf} --extract LD.prune.in --recode12 --output-missing-genotype 0 --transpose --out LD.prune'
ld3 = f'plink --bfile /public3/home/scg9832/common_data/01.genomic/4.8M/base_filtered_v0.7 --keep {sample_ls} --geno {geno} --maf {maf}  --extract LD.prune.in --make-bed --out LD.prune'
ld_all = f"# LD pruning\n{ld1}\n{ld2}\n{ld3}\n\n"

# association file
gene1 = f'plink --bfile /public3/home/scg9832/common_data/01.genomic/4.8M/base_filtered_v0.7 --keep {sample_ls} --geno {geno} --maf {maf} --make-bed --out gene'
gene2 = f'plink --bfile gene --recode12 --output-missing-genotype 0 --transpose --out genet'
if software == 'emmax':
    gene_all = f"# Quality control\n{gene1}\n{gene2}\n\n"
elif software == 'gemma':
    gene_all = f"# Quality control\n{gene1}\n\n"

# sorted sample
phe1 = f'python /public3/home/scg9832/script/python/sample_sort_for_gwas.py {phenotype} LD.prune.fam {phe_id}_sorted.txt'
phe2 = "awk '{print $3}' " + f"{phe_id}_sorted.txt > {phe_id}_p.txt"
if software == 'emmax':
    phe_all = f"# Sorted sample\n{phe1}\n\n"
elif software == 'gemma':
    phe_all = f"# Sorted sample\n{phe1}\n{phe2}\n\n"

# pca
p1 = 'plink --bfile LD.prune --pca 5'
p2_e = "awk '{print $1,$2,1,$3,$4,$5,$6,$7}' plink.eigenvec > pca.emmax.cov"
p2_g = "awk '{print $3,$4,$5,$6,$7}' plink.eigenvec > pca.gemma.cov"
if software == 'emmax':
    pca_all = f"# PCA\n{p1}\n{p2_e}\n\n"
elif software == 'gemma':
    pca_all = f"# PCA\n{p1}\n{p2_g}\n\n"

# kinship
k_e = "emmax-kin -v -s -d 10 LD.prune 1>log.kin"
k_g = f"gemma-0.98.5-linux -bfile LD.prune -gk 2 -p {phe_id}_p.txt"
if software == 'emmax':
    kin_all = f"# kinship by emmax\n{k_e}\n\n"
elif software == 'gemma':
    kin_all = f"# kinship by gemma\n{k_g}\n\n"

# GWAS
g_e_1 = f"emmax -v -d 10 -p {phe_id}_sorted.txt -t genet -k LD.prune.hIBS.kinf -c pca.emmax.cov -o {out}"
g_e_2 = f"perl /public3/home/scg9832/script/perl/emmax2plotData.pl {out}.ps > {out}"
g_g_1 = f"gemma-0.98.5-linux -bfile gene -k output/result.sXX.txt -lmm 1 -p {phe_id}_p.txt -c pca.gemma.cov -o {out}"
if software == 'emmax':
    gwas_all = f"# GWAS by emmax\n{g_e_1}\n{g_e_2}\n"
elif software == 'gemma':
    gwas_all = f"# GWAS by gemma\n{g_g_1}\n"

# tt
tt = f"# !/usr/bin/bash\n# by zhiyuan\n# GWAS by {software.title()}\n\n"

with open(f"GWAS_{software}.sh",'w')as fd:
    fd.write(tt)
    fd.write(ld_all)
    fd.write(gene_all)
    fd.write(phe_all)
    fd.write(pca_all)
    fd.write(kin_all)
    fd.write(gwas_all)
