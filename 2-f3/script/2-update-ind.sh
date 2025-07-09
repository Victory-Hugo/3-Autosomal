#!/bin/bash
/home/luolintao/miniconda3/envs/pyg/bin/python3 \
    "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/2-f3/src/1-update_ind.py" \
  --ind /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.ind \
  --groups /mnt/f/OneDrive/文档（科研）/脚本/Download/6-format-trans/5-geno-ind-snp/conf/ID_pop.txt \
  --output /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data_update.ind

mv /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.ind \
   /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.ind.backup

# 将新的文件替代旧文件，旧文件备份
mv /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data_update.ind \
   /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.ind
