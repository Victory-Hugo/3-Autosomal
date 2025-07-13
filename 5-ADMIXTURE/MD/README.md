# 介绍
该脚本用于计算ADMIXTURE值并进行可视化。ADMIXTURE值是指个体在不同种群中的基因组成分比例。

# 使用方法
1. 确保已安装Python和相关依赖库。
2. 使用5-ADMIXTURE/script/0-1-PLINK.sh将`VCF`文件转换为`PLINK`格式。
3. 运行5-ADMIXTURE/script/0-2-ADMIXTURE.sh计算ADMIXTURE值。
   - 注意：HPC可以使用0-2-ADMIXTURE_HPC.sh
4. 使用5-ADMIXTURE/script/0-3-ADMIXTURE_plot.sh进行可视化。
   - 可视化之前请修改`fam`文件第一例为群体名称。