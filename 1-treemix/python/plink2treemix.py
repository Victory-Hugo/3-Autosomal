#!/usr/bin/python

import sys, os, gzip

if len(sys.argv) < 3:
    print("plink2treemix.py [gzipped input file] [gzipped output file]")
    print("ERROR: improper command line")
    exit(1)
# Python 3: use text mode and encoding
infile = gzip.open(sys.argv[1], "rt", encoding="utf-8")
outfile = gzip.open(sys.argv[2], "wt", encoding="utf-8")

pop2rs = dict()
rss = list()
rss2 = set()

# 跳过表头
line = infile.readline()
line = infile.readline()
line_count = 0  # 添加行计数器用于生成唯一SNP标识符

while line:
    line = line.strip().split()
    if len(line) < 8:  # 确保行有足够的列
        line = infile.readline()
        continue
        
    chr_id = line[0]
    rs = line[1]
    pop = line[2]
    mc = line[6]
    total = line[7]
    
    # 创建唯一SNP标识符: chr_rs_linenum
    unique_rs = f"{chr_id}_{rs}_{line_count}"
    line_count += 1
    
    if unique_rs not in rss2:
        rss.append(unique_rs)
    rss2.add(unique_rs)
    
    if pop not in pop2rs:
        pop2rs[pop] = dict()
    if unique_rs not in pop2rs[pop]:
        # 确保mc和total是字符串
        pop2rs[pop][unique_rs] = " ".join([str(mc), str(total)])
    line = infile.readline()

pops = list(pop2rs.keys())
for pop in pops:
    print(pop, end=' ', file=outfile)
print("", file=outfile)

for rs in rss:
    row = []
    for pop in pops:
        # 检查该群体是否有这个SNP的数据
        if rs in pop2rs[pop]:
            tmp = pop2rs[pop][rs].split()
            c1 = int(tmp[0])
            c2 = int(tmp[1])
            c3 = c2-c1
            row.append(f"{c1},{c3}")
        else:
            # 如果该群体没有此SNP的数据，输出0,0
            row.append("0,0")
    # 检查该行是否全为0,0
    if all(x == "0,0" for x in row):
        continue
    print(" ".join(row), file=outfile)