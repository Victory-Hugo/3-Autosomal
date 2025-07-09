#!/bin/bash

# VCF to EIGENSTRAT format conversion pipeline (Fixed version)
# 修正版，处理细菌基因组和非标准染色体名称
# Author: [Your Name]
# Date: $(date)

set -e  # 遇到错误时停止执行
set -u  # 使用未定义变量时报错

# 定义变量
INPUT_VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.SNP.maf99.vcf.recode.vcf.gz"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp"
PREFIX="7544_filtered"

# 检查输入文件是否存在
if [ ! -f "$INPUT_VCF" ]; then
    echo "错误: 输入VCF文件不存在: $INPUT_VCF"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

echo "开始VCF质量控制和格式转换..."

# 步骤1: VCF质量控制，生成PLINK二进制格式
echo "步骤1: 执行VCF质量控制..."

# 尝试标准方法（需要 --double-id 参数处理复杂样本ID）
echo "尝试标准VCF转换..."
plink --vcf "$INPUT_VCF" \
      --allow-extra-chr \
      --double-id \
      --maf 0.01 \
      --make-bed \
      --out "$OUTPUT_DIR/$PREFIX"

if [ $? -ne 0 ]; then
    echo "标准方法失败，尝试分步处理..."
    
    # 先转换，后过滤
    plink --vcf "$INPUT_VCF" \
          --allow-extra-chr \
          --double-id \
          --make-bed \
          --out "$OUTPUT_DIR/${PREFIX}_raw"
    
    if [ $? -eq 0 ]; then
        echo "基础转换成功，现在进行质量控制..."
        plink --bfile "$OUTPUT_DIR/${PREFIX}_raw" \
              --maf 0.01 \
              --make-bed \
              --out "$OUTPUT_DIR/$PREFIX"
    else
        echo "错误: PLINK VCF转换完全失败"
        exit 1
    fi
fi

if [ $? -ne 0 ]; then
    echo "错误: PLINK VCF质量控制失败"
    exit 1
fi

echo "步骤1完成: PLINK二进制文件已生成"

# 步骤1.5: 修改.bim文件中的染色体名称以兼容EIGENSTRAT
echo "步骤1.5: 修改染色体名称以兼容EIGENSTRAT..."

# 备份原始.bim文件
cp "$OUTPUT_DIR/${PREFIX}.bim" "$OUTPUT_DIR/${PREFIX}.bim.backup"

# 将非标准染色体名称替换为数字（EIGENSTRAT要求）
# NC_000915.1 -> 1 (因为是单个细菌基因组)
sed -i 's/NC_000915\.1/1/g' "$OUTPUT_DIR/${PREFIX}.bim"

echo "染色体名称已修改: NC_000915.1 -> 1"

# 步骤2: 连锁不平衡剪枝
echo "步骤2: 执行连锁不平衡剪枝..."
plink --bfile "$OUTPUT_DIR/$PREFIX" \
      --allow-extra-chr \
      --indep-pairwise 50 10 0.1 \
      --out "$OUTPUT_DIR/${PREFIX}_pruned"

if [ $? -ne 0 ]; then
    echo "错误: 连锁不平衡剪枝失败"
    exit 1
fi

# 步骤3: 提取剪枝后的SNP
echo "步骤3: 提取剪枝后的SNP..."
plink --bfile "$OUTPUT_DIR/$PREFIX" \
      --extract "$OUTPUT_DIR/${PREFIX}_pruned.prune.in" \
      --allow-extra-chr \
      --make-bed \
      --out "$OUTPUT_DIR/${PREFIX}_pruned_data"

if [ $? -ne 0 ]; then
    echo "错误: 提取剪枝后SNP失败"
    exit 1
fi

# 步骤4: 创建CONVERTF参数文件
echo "步骤4: 创建CONVERTF参数文件..."
cat > "$OUTPUT_DIR/parfile.txt" << EOF
genotypename: $OUTPUT_DIR/${PREFIX}_pruned_data.bed
snpname: $OUTPUT_DIR/${PREFIX}_pruned_data.bim
indivname: $OUTPUT_DIR/${PREFIX}_pruned_data.fam
outputformat: EIGENSTRAT
genotypeoutname: $OUTPUT_DIR/${PREFIX}_pruned_data.geno
snpoutname: $OUTPUT_DIR/${PREFIX}_pruned_data.snp
indivoutname: $OUTPUT_DIR/${PREFIX}_pruned_data.ind
EOF

# 步骤5: 转换为EIGENSTRAT格式
echo "步骤5: 转换为EIGENSTRAT格式..."
cd "$OUTPUT_DIR"

# 检查convertf是否可用
if ! command -v convertf &> /dev/null; then
    echo "错误: convertf命令未找到，请确保EIGENSOFT已正确安装"
    exit 1
fi

convertf -p parfile.txt

if [ $? -ne 0 ]; then
    echo "错误: CONVERTF转换失败"
    echo "可能的原因:"
    echo "1. EIGENSOFT未正确安装"
    echo "2. 染色体名称仍不兼容"
    echo "3. 文件格式问题"
    exit 1
fi

echo "CONVERTF转换完成！"

# 步骤6: 清理.ind文件格式
echo "步骤6: 清理.ind文件格式..."

GENO_FILE="$OUTPUT_DIR/${PREFIX}_pruned_data.geno"
SNP_FILE="$OUTPUT_DIR/${PREFIX}_pruned_data.snp"
IND_FILE="$OUTPUT_DIR/${PREFIX}_pruned_data.ind"

# 检查文件是否存在
if [ ! -f "$GENO_FILE" ] || [ ! -f "$SNP_FILE" ] || [ ! -f "$IND_FILE" ]; then
    echo "错误: EIGENSTRAT文件生成不完整"
    echo "检查以下文件是否存在:"
    echo "  $GENO_FILE: $([ -f "$GENO_FILE" ] && echo "存在" || echo "不存在")"
    echo "  $SNP_FILE: $([ -f "$SNP_FILE" ] && echo "存在" || echo "不存在")"
    echo "  $IND_FILE: $([ -f "$IND_FILE" ] && echo "存在" || echo "不存在")"
    exit 1
fi

# 备份原始.ind文件
cp "$IND_FILE" "${IND_FILE}.backup"

# 清理.ind文件格式
echo "清理.ind文件格式问题..."

# 1. 移除重复的样本名（去掉冒号和后面的重复部分）
# 2. 移除前导空格
# 3. 标准化格式为: SampleID U ???
awk '{
    # 移除前导和尾随空格
    gsub(/^[ \t]+|[ \t]+$/, "")
    
    # 提取样本ID（如果有冒号，只取第一部分）
    split($1, parts, ":")
    sample_id = parts[1]
    
    # 输出格式化的行：SampleID U ???
    printf "%-30s U        ???\n", sample_id
}' "${IND_FILE}.backup" > "$IND_FILE"

echo ".ind文件格式已清理完成"

# 步骤7: 文件一致性检查
echo "步骤7: 执行文件一致性检查..."

# 统计文件行数和列数
GENO_ROWS=$(wc -l < "$GENO_FILE")
GENO_COLS=$(head -n 1 "$GENO_FILE" | wc -c)
SNP_ROWS=$(wc -l < "$SNP_FILE")
IND_ROWS=$(wc -l < "$IND_FILE")

echo "文件统计信息:"
echo "  .geno文件: $GENO_ROWS 行, $GENO_COLS 列"
echo "  .snp文件: $SNP_ROWS 行"
echo "  .ind文件: $IND_ROWS 行"

# 验证一致性
if [ "$GENO_ROWS" -eq "$SNP_ROWS" ]; then
    echo "✓ .geno行数与.snp行数一致 ($GENO_ROWS)"
else
    echo "✗ 警告: .geno行数($GENO_ROWS)与.snp行数($SNP_ROWS)不一致"
fi

# 检查样本数一致性
SAMPLES_IN_GENO=$((GENO_COLS - 1))  # 减去换行符
if [ "$SAMPLES_IN_GENO" -eq "$IND_ROWS" ]; then
    echo "✓ .geno列数与.ind行数一致 ($IND_ROWS个样本)"
else
    echo "⚠ 注意: .geno列数($SAMPLES_IN_GENO)与.ind行数($IND_ROWS)可能不一致"
fi

# 步骤8: 显示清理后的.ind文件示例
echo "步骤8: 显示清理后的.ind文件格式..."
echo "清理后的.ind文件前5行内容:"
head -n 5 "$IND_FILE"

echo ""
echo "提示: 现在可以根据实际群体信息修改.ind文件中的群体标签"
echo "当前格式: SampleID U ???"
echo "修改为:   SampleID U PopulationName"
echo ""
echo "示例修改命令:"
echo "  sed -i 's/???$/Helicobacter_pylori/g' $IND_FILE  # 将所有???替换为Helicobacter_pylori"

echo ""
echo "=== 流程完成 ==="
echo "输出文件位置: $OUTPUT_DIR"
echo "EIGENSTRAT格式文件:"
echo "  - ${PREFIX}_pruned_data.geno (基因型数据)"
echo "  - ${PREFIX}_pruned_data.snp  (SNP信息)"
echo "  - ${PREFIX}_pruned_data.ind  (样本信息，格式已清理)"
echo ""
echo "重要说明:"
echo "1. 染色体名称已从 NC_000915.1 修改为 1"
echo "2. 原始文件已备份：${PREFIX}.bim.backup 和 ${PREFIX}_pruned_data.ind.backup"
echo "3. .ind文件格式已清理，移除了重复样本名和前导空格"
echo "4. 下一步: 请根据实际群体信息修改.ind文件中的群体标签"