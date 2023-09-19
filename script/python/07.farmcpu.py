import sys, os
import argparse,time

args = argparse.ArgumentParser('Run farmCPU for GWAS')

args.add_argument('-i', '--input', help='Input file[can be vcf or plink bfile]', required=True, type=str)
args.add_argument('-p', '--phe', help='Phenotype file', required=True, type=str)
args.add_argument('-o', '--output', help='Output Path', default="./", type=str)
args.add_argument('-t', '--type', help='Input file type', required=True, choices=['vcf', 'plink-bfile'], type=str)

args = args.parse_args()

inputfile = args.input
phefile = args.phe
output = args.output
type = args.type
a = f'==============================================================================================\n'
print(a)
print('To run farmCPU for GWAS')
print('Author: Gao Zhiyuan\nTime: 2023-08-29\n')
print(f'python 07.farmcpu.py\n\nInput file: {inputfile}\nPhenotype file: {phefile}\nOutput Path: {output}\nInput file type: {type}\n\nStart running...\n\n')
print(a)

# 判断输入文件类型
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

# 制备info文件和geno文件
print("Generating info and geno files...")
time_start = time.time()
os.system(f'python /public3/home/scg9832/script/01.use/03.farmcpu/vcf2num.py {inputfile} tpm')
time_end = time.time()
print(f'Generate finished! Time used: {time_end - time_start}s')
print(a)

# 进行PCA分析
