import sys, os
import argparse,time

args = argparse.ArgumentParser('For Population tree calculation')

args.add_argument('-i', '--input', help='Input file[can be vcf or plink bfile]', required=True, type=str)
args.add_argument('-p', '--prefix', help='Output file prefix', required=True, type=str)
args.add_argument('-t', '--type', help='Input file type', required=True, choices=['vcf', 'plink-bfile'], type=str)

args = args.parse_args()

inputfile = args.input
prefix = args.prefix
type = args.type
a = f'==============================================================================================\n'
print(a)
print('To calculate Tree use iqtree2')
print('Author: Gao Zhiyuan\nTime: 2023-08-29\n')
print(f'python 02.tree.py\n\nInput file: {inputfile}\nOutput file prefix: {prefix}\nInput file type: {type}\n\nStart running...\n\n')
print(a)

# 将plink bfile转换为vcf
if type == 'plink-bfile':
    print('Input file is plink bfile, converting to vcf...')
    time_start = time.time()
    inputfile_z = inputfile.split('/')[-1]
    os.system(f'plink --bfile {inputfile} --recode vcf-iid --out {inputfile_z}')
    inputfile = inputfile_z + '.vcf'
    time_end = time.time()
    print(f'Convertion finished! Time used: {time_end - time_start}s')
    print(a)

# 根据vcf文件得到phy文件
print('Converting vcf to phylip...')
time_start = time.time()
os.system(f'python /public3/home/scg9832/script/python/vcf2phylip.py -i {inputfile}')
time_end = time.time()
print(f'Convertion finished! Time used: {time_end - time_start}s')
print(a)

# 根据phy文件计算进化树
print('Calculating tree...')
time_start = time.time()
os.system(f'iqtree2 -s {inputfile[:-3]}min4.phy -st DNA -T 20 -mem 80G -m GTR -redo -B 1000 -bnni --prefix {prefix}_tree')
time_end = time.time()
print(f'Calculation finished! Time used: {time_end - time_start}s')
print(a)
