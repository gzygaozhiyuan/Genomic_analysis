import sys, os
import argparse,time

args = argparse.ArgumentParser('For Population haplotype analysis')

args.add_argument('-g', '--gene_ls', help='Gene list file', required=True, type=str)
args.add_argument('-s', '--sample_ls', help='Sample list file prefix', required=True, type=str)
args.add_argument('-gp', '--gene_phe', help='Gene to Phenotype file', required=True, type=str)
args.add_argument('-p', '--phe', help='Phenotype file', required=True, type=str)
args.add_argument('-n', '--num', help='The least num per haplotype', required=True, type=int)

args = args.parse_args()

gene_ls_file = args.gene_ls
sample_ls_file = args.sample_ls
gene_phe_file = args.gene_phe
phe_file = args.phe
num = args.num

st1 = time.time()

a = f'==============================================================================================\n'

print(a)
print('To running haplotype analysis')
print('Author: Gao Zhiyuan\nTime: 2023-08-31\n')
print(f'python 04.hap.py\n\nGene list file: {gene_ls_file}\nSample list file: {sample_ls_file}\nGene phe file: {gene_phe_file}\nPhe file: {phe_file}\nNum: {num}\n\nStart running...\n\n')
print(a)

# 对promoter和cds区分别进行单倍型分析
time_start = time.time()
print(f'Perform haplotype analysis on the promoter region...')
os.system(f'perl /san/hp127T-5/Students.Dir/gao.zhiyuan/scripts/perl/haplotype_MSU.pl --input {gene_ls_file} --region promoter --maf 0.05 --sample {sample_ls_file}')
os.system(f'perl /san/hp127T-5/Students.Dir/gao.zhiyuan/scripts/perl/haplotype_MSU.pl --input {gene_ls_file} --region cds --maf 0.05 --sample {sample_ls_file}')
time_end = time.time()
print(f'Convertion finished! Time used: {time_end - time_start}s')
print(a)

# 获取CDS和promoter文件的路径
file_ls = os.listdir('./')
for each_file in file_ls:
    if each_file.endswith('promoter1000'):
        promoter_path = each_file
    elif each_file.endswith('CDS'):
        cds_path = each_file

os.system('mkdir 01.cds 02.pro 03.cds 04.pro 05.cds 06.pro 07.hap_fin 08.hap_phe_sub')
hap_script_path = '/san/hp127T-5/Students.Dir/gao.zhiyuan/scripts/python/01.quntixiangguan/01.hap'

# 单倍型分析
## 脚本1：将原单倍型结果拆分开，并只保留双等位位点
sc1_cds = f'python {hap_script_path}/01.new_hap_each_SNP.py {phe_file} {gene_phe_file} {cds_path} 01.cds cds'
sc1_pro = f'python {hap_script_path}/01.new_hap_each_SNP.py {phe_file} {gene_phe_file} {promoter_path} 02.pro pro'
## 脚本2：对每个位点做测验t-test
sc2_cds = f'python {hap_script_path}/02.new_hap_each_SNP.py ./01.cds ./03.cds cds'
sc2_pro = f'python {hap_script_path}/02.new_hap_each_SNP.py ./02.pro ./04.pro pro'
## 脚本3：多个表型合并
sc3_cds = f'python {hap_script_path}/03.new_hap_muphe.py {gene_ls_file} ./03.cds ./05.cds'
sc3_pro = f'python {hap_script_path}/03.new_hap_muphe.py {gene_ls_file} ./04.pro ./06.pro'
## 脚本4：合并promoter区与CDS区
sc4 = f'python {hap_script_path}/04.new_hap_hebin.py {gene_ls_file} ./05.cds ./06.pro ./07.hap_fin'
## 脚本5：重新组合单倍型结合上表型
sc5 = f'python {hap_script_path}/05.new_hap_newhap.py ./07.hap_fin ./08.hap_phe_sub {phe_file}'
## 脚本6：进行方差分析
sc6 = f'python {hap_script_path}/06.new_hap_anova.py {gene_phe_file} {str(num)} ./08.hap_phe_sub result_anova.txt'

print('Perform haplotype analysis...')
time_start = time.time()
print('Script 1 running...')
os.system(sc1_cds)
os.system(sc1_pro)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print('Script 2 running...')
time_start = time.time()
os.system(sc2_cds)
os.system(sc2_pro)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print('Script 3 running...')
time_start = time.time()
os.system(sc3_cds)
os.system(sc3_pro)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print('Script 4 running...')
time_start = time.time()
os.system(sc4)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print('Script 5 running...')
time_start = time.time()
os.system(sc5)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print('Script 6 running...')
time_start = time.time()
os.system(sc6)
time_end = time.time()
print(f'Finished! Time used: {time_end - time_start}s\n')
print(a)
sp1 = time.time()
print(f'Finished! All Time used: {st1 - sp1}s\n')