#!/bin/bash
# 输入文件是bed bim fam文件，通过vcf文件转化而来
# 使用 PLINK 将 VCF 转换为 BED 格式

# 设置基础目录
BASE_DIR="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE"
INPUT_DIR="${BASE_DIR}/data/" # 该文件夹内需要准备vcf.gz和索引文件
OUTPUT_DIR="${BASE_DIR}/output/Chunk_4/" # 空文件夹

# 设置其他必要的路径
SCRIPT_DIR="/home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/"
GENETIC_MAP="${SCRIPT_DIR}/conf/genetic_map_HP.txt" # 这是HP的遗传图谱文件，按照之前的文献计算的
GENETIC_MAP_FINAL="${SCRIPT_DIR}/conf/genetic_map_Finalversion_HP.txt" # 通过genetic_map_HP.txt计算而来
SHAPEIT_BIN="${SCRIPT_DIR}/download/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
FINESTRUCTURE_SCRIPTS="${SCRIPT_DIR}/download/fineSTRUCTURE"
ChromoPainter_SCRIPTS="${SCRIPT_DIR}/download/ChromoPainterv2"
PLINK_BIN="${SCRIPT_DIR}/download/plink_linux_x86_64_20241022/plink"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
# 定义输入VCF文件和输出前缀
VCF_FILE="${INPUT_DIR}/Chunk_4.vcf.gz"
PLINK_PREFIX="$OUTPUT_DIR/Chunk_4"
PLINK_HP_PREFIX="$OUTPUT_DIR/Chunk_4_HP"

Step 1: 使用PLINK过滤VCF文件并生成BED格式
$PLINK_BIN --vcf "$VCF_FILE" \
      --make-bed \
      --double-id \
      --allow-extra-chr \
      --out "$PLINK_PREFIX" \
      --geno 0.1 \
      --mind 0.1 \
      --indep-pairwise 50 10 0.1

# 定义PLINK生成的文件
FILES=("bed" "bim" "fam")
for EXT in "${FILES[@]}"; do
    FILE="$PLINK_PREFIX.$EXT"
    cp "$FILE" "$FILE.bak"
done

# Step 2: 替换BIM文件中的染色体标识符
awk 'BEGIN{OFS="\t"} { $1 = ($1 == "NC_000915.1" ? "1" : $1); print }' "$PLINK_PREFIX.bim" > "$OUTPUT_DIR/temp.bim" && mv "$OUTPUT_DIR/temp.bim" "$PLINK_PREFIX.bim"

echo "替换完成！以下文件已更新："
echo "$PLINK_PREFIX.bim"
echo "原始文件备份在：$PLINK_PREFIX.bim.bak"

# Step 3: 使用PLINK选择染色体1并生成新的BED文件
$PLINK_BIN --bfile "$PLINK_PREFIX" \
      --chr 1 \
      --make-bed \
      --alleleACGT \
      --out "$PLINK_HP_PREFIX"

# Step 4: 使用Shapeit进行相位处理
$SHAPEIT_BIN --input-bed "$PLINK_HP_PREFIX.bed" "$PLINK_HP_PREFIX.bim" "$PLINK_HP_PREFIX.fam" \
             --input-map "$GENETIC_MAP" \
             --output-max "$OUTPUT_DIR/Chunk_4_HP.phased.haps" "$OUTPUT_DIR/Chunk_4_HP.phased.sample" \
             --output-log "$OUTPUT_DIR/Chunk_4_HP.log" \
             --force \
             --burn 30 \
             --prune 30 \
             --main 30 \
             --thread 8

# Step 5: 生成ID文件用于ChromoPainter
#! 特别注意！！这一步生成的ids文件需要进行手动处理，将第二列替换为群体名，第三列替换为1。
awk 'BEGIN{FS=" "; OFS="\t"} {print $2, $1, $6}' "$PLINK_HP_PREFIX.fam" | tr -s '\t ' ' ' > "$OUTPUT_DIR/Chunk_4_HP.ids"

# Step 6: 使用Perl脚本转换Shapeit输出为fineSTRUCTURE格式
echo "使用impute2进行转换，这一步需要耗费很多时间。"
perl "$FINESTRUCTURE_SCRIPTS/impute2chromopainter.pl" -J \
     "$OUTPUT_DIR/Chunk_4_HP.phased.haps" \
     "$OUTPUT_DIR/Chunk_4_HP.phase"
dos2unix "$GENETIC_MAP_FINAL"
echo '保证genetic_map_Finalversion_HP.txt文件是UNIX格式!'

perl "$FINESTRUCTURE_SCRIPTS/convertrecfile.pl" -M \
     hap "$OUTPUT_DIR/Chunk_4_HP.phase" \
     "$GENETIC_MAP_FINAL" \
     "$OUTPUT_DIR/Chunk_4_HP.recombfile"

echo "所有步骤完成！"
# 删除备份的.bak文件
rm "$PLINK_PREFIX.bim.bak"
rm "$PLINK_PREFIX.bed.bak"
rm "$PLINK_PREFIX.fam.bak"
echo "备份文件已经删除!"

# 第三步 
dos2unix "$BASE_DIR/Chunk_4_HP.ids"
echo "保证$BASE_DIR/Chunk_4_HP.ids文件是UNIX格式!"
dos2unix "$BASE_DIR/popDonRec.txt"
echo "保证$BASE_DIR/popDonRec.txt文件是UNIX格式!"

python3 /home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/src/3_replace_ids.py \
      --inf_csv /home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/conf/Anchor.csv \
      /home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/OUTPUT/Chunk_4/Chunk_4_HP.ids

# 执行 ChromoPainterv2 命令
"${ChromoPainter_SCRIPTS}/ChromoPainterv2" \
    -g "${OUTPUT_DIR}/Chunk_4_HP.phase" \
    -r "${OUTPUT_DIR}/Chunk_4_HP.recombfile" \
    -t "${OUTPUT_DIR}/Chunk_4_HP.ids" \
    -f /home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/conf/popDonRec.txt 0 0 -s 0 -n 1089.017 -M 0.038307 \
    -o "${OUTPUT_DIR}/Donor_v_Target"
