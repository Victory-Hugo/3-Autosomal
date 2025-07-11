#!/bin/bash

# ---------------------------------------
# 根据实际环境修改下列路径
# ---------------------------------------
BASE_DIR="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE"
INPUT_DIR="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/data"   # 内含多个待处理 vcf.gz
OUTPUT_DIR="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/output/"  # 作为总的输出目录
GENETIC_MAP="/home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/conf/genetic_map_HP.txt"
GENETIC_MAP_FINAL="/home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/conf/genetic_map_Finalversion_HP.txt"
SHAPEIT_BIN="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"
FINESTRUCTURE_SCRIPTS="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/fineSTRUCTURE"
ChromoPainter_SCRIPTS="/home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/ChromoPainterv2"

mkdir -p "$OUTPUT_DIR"

#############################################
# 先顺序执行：为所有输入文件创建输出文件夹
#############################################
for VCF_FILE in "${INPUT_DIR}"/*.vcf.gz; do
    BASENAME=$(basename "$VCF_FILE" .vcf.gz)
    SAMPLE_OUT_DIR="${OUTPUT_DIR}/${BASENAME}"
    mkdir -p "$SAMPLE_OUT_DIR"
    echo "已创建输出目录: $SAMPLE_OUT_DIR"
done

#############################################
# 定义处理单个 VCF 文件的函数
#############################################
process_vcf() {
    VCF_FILE="$1"
    BASENAME=$(basename "$VCF_FILE" .vcf.gz)
    SAMPLE_OUT_DIR="${OUTPUT_DIR}/${BASENAME}"
    PLINK_PREFIX="${SAMPLE_OUT_DIR}/${BASENAME}"
    PLINK_HP_PREFIX="${SAMPLE_OUT_DIR}/${BASENAME}_plink"

    echo "==================================="
    echo "开始处理文件: $VCF_FILE"
    echo "输出目录: $SAMPLE_OUT_DIR"
    echo "==================================="

    # Step 1: 使用 PLINK 过滤 VCF 文件并生成 BED 格式
    /home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/plink_linux_x86_64_20241022/plink --vcf "$VCF_FILE" \
          --make-bed \
          --double-id \
          --allow-extra-chr \
          --out "$PLINK_PREFIX" \
          --geno 0.1 \
          --mind 0.1 \
          --indep-pairwise 50 10 0.1 

    # 备份三大文件(bed, bim, fam)
    for EXT in bed bim fam; do
        cp "${PLINK_PREFIX}.${EXT}" "${PLINK_PREFIX}.${EXT}.bak"
    done

    # Step 2: 替换 BIM 文件中的染色体标识符
    awk 'BEGIN{OFS="\t"} { $1 = ($1 == "NC_000915.1" ? "1" : $1); print }' \
        "${PLINK_PREFIX}.bim" > "${SAMPLE_OUT_DIR}/temp.bim" && \
        mv "${SAMPLE_OUT_DIR}/temp.bim" "${PLINK_PREFIX}.bim"

    echo "染色体标识替换完成, 新的 .bim 文件: ${PLINK_PREFIX}.bim"

    # Step 3: 使用 PLINK 选择染色体 1 并生成新的 BED 文件（供下游相位使用）
    /home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/plink_linux_x86_64_20241022/plink --bfile "$PLINK_PREFIX" \
          --chr 1 \
          --make-bed \
          --alleleACGT \
          --out "$PLINK_HP_PREFIX"

    # Step 4: 使用 Shapeit 进行相位
    $SHAPEIT_BIN \
        --input-bed  "${PLINK_HP_PREFIX}.bed" "${PLINK_HP_PREFIX}.bim" "${PLINK_HP_PREFIX}.fam" \
        --input-map  "$GENETIC_MAP" \
        --output-max "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phased.haps" \
                     "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phased.sample" \
        --output-log "${SAMPLE_OUT_DIR}/${BASENAME}_plink.log" \
        --force \
        --burn 10 \
        --prune 10 \
        --main 30 \
        --thread 16

    # Step 5: 生成 ChromoPainter 所需的 ids 文件
    awk 'BEGIN{FS=" "; OFS="\t"} {print $2, $1, $6}' "${PLINK_HP_PREFIX}.fam" \
        | tr -s '\t ' ' ' \
        > "${SAMPLE_OUT_DIR}/${BASENAME}_plink.ids"

    # Step 6: 使用 Perl 脚本将 Shapeit 输出转换为 fineSTRUCTURE 格式
    echo "使用 impute2chromopainter.pl 进行转换(耗时较长)..."
    perl "${FINESTRUCTURE_SCRIPTS}/impute2chromopainter.pl" -J \
         "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phased.haps" \
         "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phase"

    # 确保遗传图谱为 UNIX 换行
    dos2unix "$GENETIC_MAP_FINAL"

    perl "${FINESTRUCTURE_SCRIPTS}/convertrecfile.pl" -M hap \
         "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phase" \
         "$GENETIC_MAP_FINAL" \
         "${SAMPLE_OUT_DIR}/${BASENAME}_plink.recombfile"

    # 删除备份文件
    rm "${PLINK_PREFIX}.bim.bak" "${PLINK_PREFIX}.bed.bak" "${PLINK_PREFIX}.fam.bak"
    /home/anaconda3/bin/python3 /home/luolintao/Helicopter/Script/分析结果/fineSTRUCTURE/Script/3_replace_ids.py
    echo "ids文件已经替换完成!"
    # # 以下为估算 ChromoPainter 参数的示例
    # dos2unix "${SAMPLE_OUT_DIR}/${BASENAME}_plink.ids"
    # dos2unix "${BASE_DIR}/popDonRec.txt"

    # "${ChromoPainter_SCRIPTS}/ChromoPainterv2" \
    #     -g "${SAMPLE_OUT_DIR}/${BASENAME}_plink.phase" \
    #     -r "${SAMPLE_OUT_DIR}/${BASENAME}_plink.recombfile" \
    #     -t "${SAMPLE_OUT_DIR}/${BASENAME}_plink.ids" \
    #     -f "${BASE_DIR}/popDonRec.txt" \
    #     0 0 \
    #     -s 0 -i 10 -in -iM \
    #     -o "${SAMPLE_OUT_DIR}/ChromoPainter"

    # echo "ChromoPainter 结果文件前缀为: ${SAMPLE_OUT_DIR}/ChromoPainter"
    echo "===== 完成处理: $VCF_FILE ====="
    echo
}

# 导出函数和需要的变量，以便 GNU parallel 在子进程中能够访问
export -f process_vcf
export BASE_DIR INPUT_DIR OUTPUT_DIR GENETIC_MAP GENETIC_MAP_FINAL SHAPEIT_BIN FINESTRUCTURE_SCRIPTS ChromoPainter_SCRIPTS

#############################################
# 并行执行所有处理步骤
#############################################
parallel -j 16 process_vcf {} ::: "${INPUT_DIR}"/*.vcf.gz

echo "全部文件处理完毕！"
