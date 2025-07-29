#!/bin/bash
PYTHON3_PATH="/home/luolintao/miniconda3/bin/python3" # todo python3的路径

#* 配置路径
GFF_PATH='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/9-FST-popgenome/conf/NC_000915.gff' #todo gff
FASTA_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/7544_Sample/' # TODO FASTA文件所在目录
#* FST结果目录
BASE_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/temp/' # todo Fst结果目录
TEMP_ALN="${BASE_DIR}/temp.aln" # todo 临时对齐文件
#* 配置代码
PYTHON_T='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/9-FST-popgenome/python/2-1-数据转置.py'
PYTHON_THRESHOLD='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/9-FST-popgenome/python/2-2-合理阈值.py'
PYTHON_VISUALIZATION='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/9-FST-popgenome/python/2-3-FST可视化.py'
R_SCRIPT='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/9-FST-popgenome/pipe/s6_analysis_FST_BATCH.r' # todo R脚本路径

Rscript ${R_SCRIPT} \
  "${BASE_DIR}/META.csv" \
  "${FASTA_DIR}" \
  "${TEMP_ALN}" \
  "${BASE_DIR}/" \
  "${BASE_DIR}/FST.txt"

 ${PYTHON3_PATH} \
   ${PYTHON_T} \
    --base_dir ${BASE_DIR} 


 ${PYTHON3_PATH} \
  ${PYTHON_THRESHOLD} \
    --input_csv "${BASE_DIR}/处理后FST.csv" \
    --output_png "${BASE_DIR}/FST_percentile_distribution.png"

 ${PYTHON3_PATH} ${PYTHON_VISUALIZATION} \
    --base_dir "${BASE_DIR}" \
    --limitation 0.2 \
    --gff_file ${GFF_PATH}
