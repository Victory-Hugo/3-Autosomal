#! /bin/bash
#*#######################################################
###* 全基因组关联分析（GWAS），样本不要放太多，在250个左右
#*#######################################################
## 从全球对齐中提取重点人群的菌株序列
## 输入：strainsSubpop.txt [重点人群菌株列表], DATA/woHpGP_dnds.cvg70.aln [全球对齐]
## 输出：alignmentSub.aln [仅包含重点人群菌株的全球对齐]
PYTHON_FASTA_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/python/0_从aln提取需要的样本.py" #todo该脚本将路径下的多个fasta文件（长度一致）拼接为一个aln文件
PYTHON_BUGWAS_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/python/0_准备bugwas输入文件.py" #todo该脚本将输入的txt文件转换为bugwas所需的基因型文件
PYTHON3_PATH="/home/luolintao/miniconda3/bin/python3" #todo python3的路径
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/8-临床结局相关/data/GWAS/Atrophy-萎缩" #todo 输入目录
FASTA_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/7544_Sample/" #todo fasta文件所在目录
VERYFASTTREE_PATH='/mnt/e/Scientifc_software/veryfasttree-master/veryfasttree-master/VeryFastTree' #todo VeryFastTree的路径

#* 第一步
${PYTHON3_PATH} ${PYTHON_FASTA_SCRIPT} \
    ${INPUT_DIR}/list.txt \
    ${FASTA_DIR}/ \
    ${INPUT_DIR}/alignmentSub.aln

#* 第二步
snp-sites -c -o\
     ${INPUT_DIR}/alignmentSub.temp \
     ${INPUT_DIR}/alignmentSub.aln

nohup ${VERYFASTTREE_PATH} -nt \
    -threads 21 \
    ${INPUT_DIR}/alignmentSub.temp \
    > ${INPUT_DIR}/treeFull.nwk \
    2> ${INPUT_DIR}/treeFull.log & #todo 异步运行VeryFastTree生成系统发育树


#* 生成基因型文件（使用所有SNP - 包括CDS和基因间）
snp-sites -c -v -o\
     ${INPUT_DIR}/tmp1.vcf \
     ${INPUT_DIR}/alignmentSub.aln


#* 处理VCF文件，提取所需列
sed "s/#CHROM/CHROM/" \
    ${INPUT_DIR}/tmp1.vcf \
    | grep -v "#" | cut -f 2,4-5,10- | sed "1s/POS/ps/" \
    > ${INPUT_DIR}/tmp2.txt  # 处理VCF文件，提取所需列

#* 使用默认双等位过滤，自动生成输出文件路径
${PYTHON3_PATH} ${PYTHON_BUGWAS_SCRIPT} \
     ${INPUT_DIR}/tmp2.txt \
     --output_file ${INPUT_DIR}/geno_biallelic_SNP.txt


## 生成输入文件、运行GWAS、后处理和绘制树图
## 输入：
##          1. 生成输入：METADATA/strains_SiberiaIndAm.csv [仅包含来自重点人群的菌株，提取表型；此文件基于不同菌株的元数据手动生成，如其位置、生态种归属等] -> input/pheno.txt [表型文件], strains_pop.txt [所有菌株的人群信息], input/tree [用于GWAS的树]
##          2. GWAS：input/geno_biallelic_SNP.txt [基因型文件], input/pheno.txt, input/tree, gemma可执行文件的路径
##          3. 分析：DATA/26695.gff [参考基因组的注释文件 - 从NCBI下载], bugwas_biallelic_lmmout_allSNPs.txt [bugwas输出]
## 输出：
##          1. 生成输入：input/pheno.txt [表型文件], PLOT/TREE/TREE_GWAS.jpeg [用于GWAS输入的树图]
##          2. GWAS：bugwas输出 + data_GWAS.RData [GWAS输出作为RData对象]
##          3. 分析：PLOT/GWAS/regions/* [曼哈顿图及基因组不同区域的缩放图], significantCDS.csv [显著的CDS列表], significantSNPs.txt [显著的SNP列表]
