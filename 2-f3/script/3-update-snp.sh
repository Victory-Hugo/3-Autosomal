#!/usr/bin/env bash
set -euo pipefail

# 原始文件路径
orig="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.snp"
# 临时输出文件
upd="${orig%.snp}.update.snp"
# 备份文件
bak="${orig}.bak"
# 染色体编号
chrom="1"

# 1) 生成 update.snp：  
#    - 第一列设为空  
#    - 第二列设为 “<染色体>_<物理位置>”（取自原第4列）  
#    - 第三列插入 染色体编号  
#    - 剩余列原样后移  
awk -v CHR="$chrom" 'BEGIN { OFS="\t" }
{
  # 第1列：空字符串
  printf "%s", ""
  # 第2列：CHR_物理位置（原第4列）
  printf "%s%s", OFS, CHR "_" $4
  # 第3列：CHR
  printf "%s%s", OFS, CHR
  # 剩余列：原第3列及其之后
  for(i=3; i<=NF; i++){
    printf "%s%s", OFS, $i
  }
  # 行尾换行
  printf "\n"
}' "$orig" > "$upd"

# 2) 备份原文件
mv -- "$orig" "$bak"

# 3) 用 update.snp 覆盖原文件名
mv -- "$upd" "$orig"

echo "已完成："
echo "  原文件已备份为：$bak"
echo "  新文件替换为：  $orig"
