import sys, os
import argparse,time

args = argparse.ArgumentParser('For Population LD decay calculation')

args.add_argument('-i', '--input', help='Input file[can be vcf or plink bfile]', required=True, type=str)
args.add_argument('-p', '--prefix', help='Output file prefix', required=True, type=str)
args.add_argument('-t', '--type', help='Input file type', required=True, choices=['vcf', 'plink-bfile'], type=str)
args.add_argument('-s', '--sub', help='Subpopulation file', required=False, type=str)

args = args.parse_args()

inputfile = args.input
prefix = args.prefix
type = args.type
subfile = args.sub
a = f'==============================================================================================\n'
print(a)
print('To calculate LD decay')
print('Author: Gao Zhiyuan\nTime: 2023-08-29\n')
print(f'python 01.ld.py\n\nInput file: {inputfile}\nOutput file prefix: {prefix}\nInput file type: {type}\nSubpopulation file: {subfile}\n\nStart running...\n\n')
print(a)
if type == 'plink-bfile':
    print('Input file is plink bfile, converting to vcf...')
    time_start = time.time()
    inputfile_z = inputfile.split('/')[-1]
    # 将plink bfile转换为vcf
    os.system(f'plink --bfile {inputfile} --recode vcf-iid --out {inputfile_z}')
    inputfile = inputfile_z + '.vcf'
    time_end = time.time()
    print(f'Convertion finished! Time used: {time_end - time_start}s')
    print(a)

print('Split into multiple subpopulations...')
# 将样本按亚群体分开
time_start = time.time()
sub_ls = []
f1 = open(subfile, 'r')
for l in f1:
    if not l.startswith('ID'):
        ls = l.strip().split('\t')
        sample = ls[0]
        sub = ls[-1]
        with open(f'./{prefix}_{sub}.txt', 'a') as f2:
            f2.write(f'{sample}\n')
        if sub not in sub_ls:
            sub_ls.append(sub)
f1.close()
time_end = time.time()
print(f'Split finished! Time used: {time_end - time_start}s')
print(a)

# 计算整个群体的LD decay
print('Calculating LD decay for whole population...')
time_start = time.time()
os.system(f'PopLDdecay -InVCF {inputfile} -OutStat {prefix}_whole.stat.gz -MaxDist 1000')
time_end = time.time()
print(f'Calculation finished! Time used: {time_end - time_start}s')
print(a)

# 计算亚群体的LD decay
print('Calculating LD decay for each subpopulation...')
for sub in sub_ls:
    print(f'Calculating LD decay for {sub}...')
    time_start = time.time()
    os.system(f'PopLDdecay -InVCF {inputfile} -OutStat {prefix}_{sub}.stat.gz -MaxDist 1000 -SubPop {prefix}_{sub}.txt')
    time_end = time.time()
    print(f'Calculation finished! Time used: {time_end - time_start}s\n')
print(a)

# 制备绘图所需文件
print('Preparing files for plotting...')
time_start = time.time()
os.system(f'Plot_OnePop.pl -inFile {prefix}_whole.stat.gz -output wholeFig --keepR')
for sub in sub_ls:
    os.system(f'Plot_OnePop.pl -inFile {prefix}_{sub}.stat.gz -output {sub}Fig --keepR')
time_end = time.time()
print(f'Preparation finished! Time used: {time_end - time_start}s')
print(a)

