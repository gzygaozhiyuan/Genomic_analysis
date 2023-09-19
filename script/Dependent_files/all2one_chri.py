import sys

# 功能：将vcf文件整合为一条染色体的vcf文件

vcf_file = sys.argv[1]
out_file = sys.argv[2]

chr_length = {}

f1 = open(vcf_file, 'r')
for l in f1:
    if l.startswith('##'):
        with open(out_file, 'a') as f2:
            f2.write(l)
        if l.startswith('##contig'):
            ls = l.strip().split(',')
            # 获取染色体号
            chr = ls[0].split('=')[-1]
            print(chr)
            # 获取染色体长度
            length = int(ls[1].split('=')[1][:-1])
            print(length)
            chr_length[chr] = length
    elif l.startswith('#CHROM'):
        with open(out_file, 'a') as f2:
            f2.write(l)
    else:
        ls = l.strip().split('\t')
        chr = ls[0]
        pos = int(ls[1])
        if chr == '1':
            with open(out_file, 'a') as f2:
                f2.write(l)
        elif chr == '2':
            pos = pos + chr_length['1']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '3':
            pos = pos + chr_length['1'] + chr_length['2']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '4':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '5':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '6':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '7':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '8':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6'] + chr_length['7']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '9':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6'] + chr_length['7'] + chr_length['8']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '10':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6'] + chr_length['7'] + chr_length['8'] + chr_length['9']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '11':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6'] + chr_length['7'] + chr_length['8'] + chr_length['9'] + chr_length['10']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
        elif chr == '12':
            pos = pos + chr_length['1'] + chr_length['2'] + chr_length['3'] + chr_length['4'] + chr_length['5'] + chr_length['6'] + chr_length['7'] + chr_length['8'] + chr_length['9'] + chr_length['10'] + chr_length['11']
            ls[0] = '1'
            ls[1] = str(pos)
            with open(out_file, 'a') as f2:
                f2.write('\t'.join(ls) + '\n')
f1.close()