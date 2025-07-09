#!/bin/bash

# 接收输入参数
file=$1  # 输入的二进制PLINK文件前缀

# 设置PCA分析的主成分数
pcaNum=10

# 使用PLINK进行PCA分析
echo "正在对 $file 进行PCA分析..."
plink --bfile $file --pca $pcaNum

# 准备群组信息文件
echo "正在准备群组信息文件..."
# 提取FAM文件的前两列（FID和IID）
cut -d' ' -f1-2 $file.fam > first.two.csv
# 提取样本组名
cut -d' ' -f1 $file.fam > group.name.csv
# 合并信息生成群组文件
paste first.two.csv group.name.csv >  $file.csv
# 生成不重复的群组列表
awk '!seen[$0]++' group.name.csv > group.csv
