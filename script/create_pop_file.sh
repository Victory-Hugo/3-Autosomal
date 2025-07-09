#!/bin/bash

# 辅助脚本：帮助创建群体文件
echo "=== 创建群体文件辅助脚本 ==="

FAM_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data.fam"
OUTPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/conf/pop.cov"

echo "Step 1: 检查FAM文件中的样本..."
if [ ! -f "$FAM_FILE" ]; then
    echo "错误: FAM文件不存在: $FAM_FILE"
    exit 1
fi

echo "样本总数: $(wc -l < $FAM_FILE)"
echo "前10个样本:"
head -10 "$FAM_FILE"

echo ""
echo "Step 2: 创建群体文件模板..."
echo "# 请根据样本名称规律手动编辑此文件，为每个样本分配正确的群体" > "$OUTPUT_FILE"
echo "# 格式: FID IID CLST" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# 尝试根据样本名称的规律自动分组（这里需要根据你的实际情况调整）
awk '{
    sample_id = $2
    # 尝试根据样本名称规律推断群体
    # 这里只是示例，你需要根据实际的样本命名规则调整
    if (sample_id ~ /Africa/) group = "hpAfrica"
    else if (sample_id ~ /Asia/) group = "hpAsia" 
    else if (sample_id ~ /Europe/) group = "hpEurope"
    else if (sample_id ~ /China/) group = "hpChina"
    else if (sample_id ~ /Japan/) group = "hpJapan"
    else group = "Unknown"
    
    print $1, $2, group
}' "$FAM_FILE" >> "$OUTPUT_FILE"

echo "群体文件模板已创建: $OUTPUT_FILE"
echo ""
echo "请编辑此文件，为每个样本分配正确的群体名称。"
echo "可用的群体名称应该与你的研究设计一致。"
echo ""
echo "编辑完成后，运行主脚本："
echo "bash 0-prepare.sh"
