#!/bin/bash
# 输入文件是bed bim fam文件，通过vcf文件转化而来
# 使用 PLINK 将 VCF 转换为 BED 格式
#!======使用前需要安装如下软件=======
#todo 1. PLINK
#todo 2. Shapeit
#todo 3. Perl
#todo 4. ChromoPainterv2
#todo 5. fineSTRUCTURE
#todo 6. R (用于可视化)
#!=================================
# 设置基础目录
BASE_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/6-fineSTRUCTURE"
INPUT_DIR="${BASE_DIR}/data" # 该文件夹内需要准备vcf.gz和索引文件
OUTPUT_DIR="${BASE_DIR}/OUTPUT/7544_CDS/" # 空文件夹
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

#todo 以下文件参照./conf/目录下的文件
GENETIC_MAP_FINAL="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/conf/genetic_map_Finalversion_HP.txt" # 通过genetic_map_HP.txt计算而来
GENETIC_MAP="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/conf/genetic_map_HP.txt" # 这是HP的遗传图谱文件，按照之前的文献计算的
D_R_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/conf/popDonRec.txt"

# 定义输入VCF文件和输出前缀
VCF_FILE="${INPUT_DIR}/merged_biallelic_7544.NoN.maf99.WGS.recode.SNP.noN.CDS.vcf.gz"
PLINK_PREFIX="${OUTPUT_DIR}/7544_CDS"
PLINK_HP_PREFIX="${PLINK_PREFIX}_HP"


# # Step 1: 使用PLINK过滤VCF文件并生成BED格式
plink --vcf "$VCF_FILE" \
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
    FILE="${PLINK_PREFIX}.${EXT}"
    cp "$FILE" "${FILE}.bak"
done

# Step 2: 替换BIM文件中的染色体标识符
awk 'BEGIN{OFS="\t"} { $1 = ($1 == "NC_000915.1" ? "1" : $1); print }' \
     "${PLINK_PREFIX}.bim" > \
     "${OUTPUT_DIR}/temp.bim" && mv "${OUTPUT_DIR}/temp.bim" "${PLINK_PREFIX}.bim"

echo "替换完成！以下文件已更新："
echo "${PLINK_PREFIX}.bim"
echo "原始文件备份在：${PLINK_PREFIX}.bim.bak"

# Step 3: 使用PLINK选择染色体1并生成新的BED文件
plink --bfile "$PLINK_PREFIX" \
      --chr 1 \
      --make-bed \
      --alleleACGT \
      --out "$PLINK_HP_PREFIX"

# Step 4: 使用Shapeit进行相位处理
shapeit --input-bed "${PLINK_HP_PREFIX}.bed" "${PLINK_HP_PREFIX}.bim" "${PLINK_HP_PREFIX}.fam" \
             --input-map "$GENETIC_MAP" \
             --output-max "${PLINK_PREFIX}_HP.phased.haps" "${PLINK_HP_PREFIX}.phased.sample" \
             --output-log "${PLINK_HP_PREFIX}.log" \
             --force \
             --burn 10 \
             --prune 10 \
             --main 30 \
             --thread 16

# # Step 5: 生成ID文件用于ChromoPainter
# #! 特别注意！！这一步生成的ids文件需要进行手动处理，将第二列替换为群体名，第三列替换为1。
awk 'BEGIN{FS=" "; OFS="\t"} {print $2, $1, $6}' "${PLINK_HP_PREFIX}.fam" | tr -s '\t ' ' ' > \
     "${PLINK_HP_PREFIX}.ids"


# Step 6: 使用Perl脚本转换Shapeit输出为fineSTRUCTURE格式
echo "使用impute2进行转换，这一步需要耗费很多时间。"
impute2chromopainter.pl -J \
     "${PLINK_HP_PREFIX}.phased.haps" \
     "${PLINK_HP_PREFIX}.phase"
dos2unix "$GENETIC_MAP_FINAL"
echo '保证genetic_map_Finalversion_HP.txt文件是UNIX格式!'

impute2chromopainter.pl -J \
#      "/mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/OUTPUT/7544_CDS_HP.phased.haps" \
#      "/mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/OUTPUT/7544_CDS_HP.phase"

convertrecfile.pl -M \
     hap "${PLINK_HP_PREFIX}.phase" \
     "$GENETIC_MAP_FINAL" \
     "${PLINK_HP_PREFIX}.recombfile"

echo "所有步骤完成！"
# 删除备份的.bak文件
rm "${PLINK_PREFIX}.bim.bak"
rm "${PLINK_PREFIX}.bed.bak"
rm "${PLINK_PREFIX}.fam.bak"
echo "备份文件已经删除!"

# 第三步 计算chromopater所需两参数
dos2unix "${PLINK_HP_PREFIX}.ids"
echo "保证${PLINK_HP_PREFIX}.ids文件是UNIX格式!"
dos2unix "${D_R_FILE}"
echo "保证${D_R_FILE}文件是UNIX格式!"

# # 使用ChromoPainter估算参数，这里先选择1000个受体来估算参数！
# "ChromoPainterv2" \
#  -g "${OUTPUT_DIR}/7544_CDS_HP.phase" \
#  -r "${OUTPUT_DIR}/7544_CDS_HP.recombfile" \
#  -t "${BASE_DIR}/7544_CDS_HP.ids" \
#  -f "${BASE_DIR}/popDonRec.txt" 1 1000 -s 0 -i 10 -in -iM \
#  -o "${OUTPUT_DIR}/HP"
# echo "结果文件为：${OUTPUT_DIR}/HP.开头的文件!"

# # 第四步跑代理人群
# echo "为了估算参数，生成22条染色体的伪文件"


# # 直接将${OUTPUT_DIR}/HP.EMprobs.out文件复制22份到"$ChromoPainter_OUTPUT"文件夹，重命名为如下：
# # 循环复制和重命名文件
# for i in {1..22}; do
#     # 构造文件名
#     new_file="output_estimateEM_Chr${i}.EMprobs.out"
#     # 执行复制和重命名
#     cp "${OUTPUT_DIR}/HP.EMprobs.out" "${ChromoPainter_OUTPUT}/${new_file}"
# done

# echo "文件复制和重命名完成。"
# cd ${ChromoPainter_OUTPUT}
# perl "ChromoPainterv2EstimatedNeMutExtractEM.pl" > "${ChromoPainter_OUTPUT}/HP_estimateEM.txt"
# echo "参数估计完成，结果保存在：${ChromoPainter_OUTPUT}/HP_estimateEM.txt"
# echo "-n   -M"
# tail -n 1 "${ChromoPainter_OUTPUT}/HP_estimateEM.txt"

# 再次使用ChromoPainter计算
# # 获取所有 popDonRec_{i}.txt 文件
# input_popDonRecfiles=(${BASE_DIR}/INPUT/popDonRec_*.txt)

# # 循环处理每个文件
# for input_popDonRecfiles in "${input_popDonRecfiles[@]}"
# do
#   # 提取文件名中的数字部分（例如 popDonRec_1.txt -> 1）
#   i=$(basename "$input_popDonRecfiles" .txt | sed 's/popDonRec_//')

#   # 设置输出目录
#   OUTPUT_SUBDIR="${OUTPUT_DIR}/DvR_${i}"
  
#   # 创建输出目录
#   mkdir -p "$OUTPUT_SUBDIR"
  
#   # 执行 ChromoPainterv2 命令
#   "${ChromoPainter_SCRIPTS}/ChromoPainterv2" \
#     -g "${OUTPUT_DIR}/7544_CDS_HP.phase" \
#     -r "${OUTPUT_DIR}/7544_CDS_HP.recombfile" \
#     -t "${BASE_DIR}/7544_CDS_HP.ids" \
#     -f "$input_popDonRecfiles" 0 0 -s 0 -i 10 -in -iM \
#     -o "${OUTPUT_SUBDIR}/Donor_v_Target"
  
#   # 打印结果
#   echo "结果文件为：${OUTPUT_SUBDIR}/Donor_v_Target开头的文件!"
# done


# #删除中间文件
# rm -rf "${ChromoPainter_OUTPUT}"

# 第四步跑代理人群
# #########################Generating the copy vector input file##################################################################
# #To construct the copy vectors, we use the estimated parameters from Section 3.2 and run ChromoPainterv2 again, this time painting all individuals as recipients on chromosome 21:

# ###popfileSurr.txt
# Kalash D
# Basque D
# French D
# Mozabite D
# Dai D
# Tibetan_Lhasa D
# Ulchi D
# Kalash R
# Basque R
# French R
# Mozabite R
# Dai R
# Tibetan_Lhasa R
# Ulchi R

# for x in {1..22}; do
# nohup /home/biosoftware/ChromoPainterv2/ChromoPainterv2 -g 12pops_ibd.chr$x.phase -r genetic_map_Finalversion_GRCh37_chr${x}.recombfile -t 12pops_ibd.ids -f popfileSurr.txt 0 0 -s 0 -n 1089.017 -M 0.038307 -o 12pops_ibd.chr${x}_Donor_v_Donor &
#  done




# 第五步 跑目标人群
# ######text for one chromosome
# ###popfile.txt里面加代理群体D和目标群体R
# Han_HGDP_dul.DG D
# Kalash D
# Basque D
# French D
# Mozabite D
# Dai D
# Tibetan_Lhasa D
# Ulchi D
# Buryat_Ehirit R
# Buryat_Zakamensky R
# Khakass_Koibals R
# Shor_Khakassia R
# Tatar_Zabolotniye R
# Altaian_Chelkans R
# Mogolian_HGDP R
# Shor_Mountain R
# Uzbek_Khorezm R
# Dongxiang_Linxia R
# Khakass_Kachins R
# Mongol_Zasagt R
# Xibo_HGDP R
# Evenk_Transbaikal R
# Hezhen_HGDP R
# Kazakh_Aksay R

# for x in {1..22}; do
# nohup /home/biosoftware/ChromoPainterv2/ChromoPainterv2 -g 12pops_ibd.chr$x.phase -r genetic_map_Finalversion_GRCh37_chr${x}.recombfile -t 12pops_ibd.ids -f popfile.txt 0 0 -s 10 -n 1089.017 -M 0.038307 -o 12pops_ibd.chr${x}_Donor_v_Target &
#  done


# 第六步 合并代第四第五步骤的chunklengthes
# #The command \-f poplistSurr.txt 0 0" specifies that we paint all individuals from the populations listed in popfileSurr.txt as recipients. Repeat the above for each chromosome (i.e. chromosome 22 in this tutorial). We then want to sum the 12pops_ibd.chr21_Donor_v_Donor.chunklengths.out files across all chromosomes, which we can do using a provided perl script:
# perl /home/biosoftware/fastGLOBETROTTER/tutorial/ChromoPainterOutputSum.pl 12pops_ibd.chr _Donor_v_Donor.chunklengths.out
#  ### run fs chromocombine

# fs chromocombine \
#  -u  -o /mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/OUTPUT/chromocombineALLfiles \
#    12pops_ibd.chr1_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr2_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr3_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr4_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr5_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr6_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr7_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr8_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr9_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr10_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr11_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr12_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr13_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr14_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr15_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr16_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr17_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr18_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr19_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr20_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr21_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr22_Donor_v_Target.chunklengths.out \
#    12pops_ibd.chr1_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr2_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr3_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr4_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr5_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr6_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr7_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr8_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr9_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr10_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr11_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr12_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr13_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr14_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr15_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr16_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr17_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr18_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr19_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr20_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr21_Donor_v_Donor.chunklengths.out \
#    12pops_ibd.chr22_Donor_v_Donor.chunklengths.out
#  cd "${OUTPUT_DIR}"
#  cp ../7544_CDS_HP.ids ./
#  #######HPC
#  date > 7544_CDS_HP.time
#  fs 7544_CDS_HP.cp -hpc 1 \
#  -idfile 7544_CDS_HP.ids \
#  -phasefiles 7544_CDS_HP.phase \
#  -recombfiles 7544_CDS_HP.recombfile \
#  -s3iters 100000 \
#  -s4iters 50000 \
#  -s1minsnps 1000 \
#  -s1indfrac 0.1 \
#  -go
#  cat 7544_CDS_HP/commandfiles/commandfile1.txt | parallel
#  fs 7544_CDS_HP.cp -go
#  cat 7544_CDS_HP/commandfiles/commandfile2.txt | parallel
#  fs 7544_CDS_HP.cp -go
#  cat 7544_CDS_HP/commandfiles/commandfile3.txt | parallel
#  fs 7544_CDS_HP.cp -go
#  cat 7544_CDS_HP/commandfiles/commandfile4.txt | parallel
#  fs 7544_CDS_HP.cp -go
#  date >> 7544_CDS_HP.time


#  ##Visulization of fs fineStrcute
# ulimit -s unlimited
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/Fat_Initial_Fs.R $2

# ###delete some files
# #把Fst_mat.csv转化为FstMatrix.csv
# Rscript /home/biosoftware/ppgv1/R_packages/fst.R
#  cd ../
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/Fat_PCA_Admixture_TreeMixVersion1.R  $2
# #/usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/Fat_PCA_Admixture_TreeMix.R  $2
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/FstMatrix.R
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/FstMatrixGroup.R
# /usr/bin/Rscript /home/biosoftware/ppgv1/R_packages/plot_strucAverage.R $2_prunned
# grep -h CV output* > prunned_error.csv
# sort -k4 prunned_error.csv > prunned_error._sort.csv