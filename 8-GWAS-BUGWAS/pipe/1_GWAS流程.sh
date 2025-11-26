#! /bin/bash
#*#######################################################
###* 全基因组关联分析（GWAS），样本不要放太多，在250个左右
#*#######################################################
# .
# ├── input
# │   ├── list.txt
# │   ├── pheno.txt
# └── output
# merge_fasta
# ├── 3697.fasta
# ├── 3738.fasta
# ├── 3745.fasta
#*#######################################################
###* 代码配置区
#*#######################################################
CODE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS" #? 脚本目录
PYTHON_FASTA_SCRIPT="${CODE_DIR}/python/0_从aln提取需要的样本.py" #?该脚本将路径下的多个fasta文件（长度一致）拼接为一个aln文件
PYTHON_BUGWAS_SCRIPT="${CODE_DIR}/python/0_准备bugwas输入文件.py" #?该脚本将输入的txt文件转换为bugwas所需的基因型文件
PYTHON3_PATH="/home/luolintao/miniconda3/envs/BigLin/bin/python3" #? python3的路径
VERYFASTTREE_PATH='VeryFastTree' #? VeryFastTree的路径
#*#######################################################
###* 文件配置区
#*#######################################################
INPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/input" #? 输入目录
FASTA_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/7544_Sample/" #? fasta文件所在目录



${PYTHON3_PATH} ${PYTHON_FASTA_SCRIPT} \
    ${INPUT_DIR}/list.txt \
    ${FASTA_DIR}/ \
    ${INPUT_DIR}/alignmentSub.aln


snp-sites -c -o\
     ${INPUT_DIR}/alignmentSub.temp \
     ${INPUT_DIR}/alignmentSub.aln

nohup ${VERYFASTTREE_PATH} -nt \
    -threads 21 \
    ${INPUT_DIR}/alignmentSub.temp \
    > ${INPUT_DIR}/treeFull.nwk \
    2> ${INPUT_DIR}/treeFull.log & #? 异步运行VeryFastTree生成系统发育树



snp-sites -c -v -o\
     ${INPUT_DIR}/tmp1.vcf \
     ${INPUT_DIR}/alignmentSub.aln



sed "s/#CHROM/CHROM/" \
    ${INPUT_DIR}/tmp1.vcf \
    | grep -v "#" | cut -f 2,4-5,10- | sed "1s/POS/ps/" \
    > ${INPUT_DIR}/tmp2.txt  # 处理VCF文件，提取所需列

#* 使用默认双等位过滤，自动生成输出文件路径
${PYTHON3_PATH} ${PYTHON_BUGWAS_SCRIPT} \
     ${INPUT_DIR}/tmp2.txt \
     --output_file ${INPUT_DIR}/geno_biallelic_SNP.txt

