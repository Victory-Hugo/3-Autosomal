: '
脚本名称: makecluster.sh

用途:
本脚本用于处理群体遗传学分析中的数据预处理和格式转换，主要针对plink生成的数据文件，最终用于TreeMix分析。

主要步骤说明:
1. 设置工作目录和输入文件路径。
2. （已注释）利用plink进行主成分分析（PCA）。
3. （已注释）处理.fam文件，提取个体信息并生成分组文件。
4. （已注释）调用Perl脚本进行配对分析。
5. （已注释）检查并删除已存在的plink.frq.strat.gz文件，然后重新压缩plink.frq.strat文件。
6. 调用Python脚本(plink2treemix.py)将plink.frq.strat.gz文件转换为TreeMix可用的输入文件(TreeMix.gz)。
7. （已注释）调用TreeMix进行群体结构分析。

注意事项:
- 需要提前安装plink、perl、python及相关依赖。
- 路径需根据实际环境进行调整。
- 部分步骤已被注释，如需使用请取消注释。
'
#!/bin/bash

WORK_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/"
file=${WORK_DIR}/7544_filtered_pruned_data


cd ${WORK_DIR}
echo $file

# plink --bfile \
#     ${file} \
#     --pca 10

# cut -d' ' -f1-2 ${file}.fam > ${WORK_DIR}/first.two.csv
# cut -d' ' -f1 ${file}.fam > ${WORK_DIR}/group.name.csv
# paste ${WORK_DIR}/first.two.csv ${WORK_DIR}/group.name.csv >  ${file}.csv
# awk '!seen[$0]++' ${WORK_DIR}/group.name.csv > ${WORK_DIR}/group.csv

# PERL_SRC="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/pairwise.perl"
# perl \
#     ${PERL_SRC}  \
#     ${file} \
#     0 



# # 如果存在 plink.frq.strat.gz 则删除
# if [ -f plink.frq.strat.gz ]; then
#     rm ${WORK_DIR}/plink.frq.strat.gz
# fi

# gzip ${WORK_DIR}/plink.frq.strat

# /home/luolintao/miniconda3/envs/pyg/bin/python3 \
#     /mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/plink2treemix.py \
#     ${WORK_DIR}/plink.frq.strat.gz \
#     ${WORK_DIR}/TreeMix.gz

treemix -i \
    ${WORK_DIR}/TreeMix.gz \
    -m 0 -k 1000 -global 
