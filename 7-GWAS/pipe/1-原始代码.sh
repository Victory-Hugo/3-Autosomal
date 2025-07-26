#!/usr/bin/env bash
#===============================================================================
# 文件: 7-GWAS/pipe/1-原始代码.sh
# 说明: 一站式 GWAS 分析流水线脚本 —— 完整复现前述步骤
# 作者: 罗林焘
# 日期: 2025年7月26日
#===============================================================================

set -euo pipefail
IFS=$'\n\t'

#----------------------------------------
# 1. GEN/Sample 到 PED/MAP（gtool1）
#----------------------------------------

echo ">>> 1. 生成排除列表：找出多等位基因位点"
awk 'length($4)>1 || length($5)>1 {print $2}' chrH22.impute2 > chrH22.exclude
awk 'length($4)>1 || length($5)>1 {print $2}' chrH1.impute2  > chrH1.exclude

echo ">>> 1.1 用 gtool 将 impute2/SAMPLE 转为二进制 impute2 格式"
gtool -S \
    --g /public/home/BaiHao/impute/chrH22.impute2 \
    --s /public/home/BaiHao/shapeit/bin/chr22.sample \
    --og chrH22.impute2 \
    --exclusion chrH22.exclude

echo ">>> 1.2 导出 PED/MAP 文件（阈值 0.9）"
gtool -G \
    --g chrH22.impute2 \
    --s /public/home/BaiHao/shapeit/bin/chr22.sample \
    --ped chrH22_SNP1.ped \
    --map chrH22_SNP1.map \
    --chr 22 \
    --snp \
    --threshold 0.9

gtool -G \
    --g chrT1.impute2 \
    --s /public/home/BaiHao/shapeit/bin/chrT1.sample \
    --ped chrT1_SNP.ped \
    --map chrT1_SNP.map \
    --chr 1 \
    --snp \
    --threshold 0.9

#----------------------------------------
# 2. 数据合并（merge）
#----------------------------------------

echo ">>> 2. 合并汉族样本 (chrH1)"
plink --file chrH1 --merge-list allfiles.txt --recode --out YBW2
echo "    ※ 共 $(wc -l < YBW2.ped) 个样本，$(wc -l < YBW2.map) 个变异"

echo ">>> 2. 合并藏族样本 (chrT1)"
plink --file chrT1 --merge-list allfiles2.txt --recode --out Tibetan
echo "    ※ 共 $(wc -l < Tibetan.ped) 个样本，$(wc -l < Tibetan.map) 个变异"

echo ">>> 2.1 合并 SNP 子集"
plink --bfile chrH1_SNP --merge-list allfiles.txt   --recode --out Han_SNP
plink --bfile chrT1_SNP --merge-list allfiles1.txt  --recode --out Tibetan_SNP

#----------------------------------------
# 3. INFO ≥ 0.8 提取高质量 SNP
#----------------------------------------

echo ">>> 3. 提取 INFO ≥0.8 的位点"
plink --bfile YBW2             --extract info_extraction.txt --make-bed --out Han
plink --bfile Tibetan          --extract Tibetan_info.txt   --make-bed --out Tibetan1
plink --bfile Han_SNP          --extract info_extraction.txt --make-bed --out Han_SNP1
plink --bfile Tibetan_SNP      --extract info_extraction.txt --make-bed --out Tibetan_SNP1

#----------------------------------------
# 4. SNP 质量控制（QC）
#----------------------------------------

echo ">>> 4.1 缺失率过滤 (geno ≤ 0.05)"
plink --bfile Han        --geno 0.05 --make-bed --out Han1
plink --bfile Tibetan1   --geno 0.05 --make-bed --out Tibetan2

echo ">>> 4.2 次要等位基因频率过滤 (maf ≥ 0.01)"
plink --bfile Han1       --maf 0.01  --make-bed --out Han2
plink --bfile Tibetan2   --maf 0.01  --make-bed --out Tibetan3

echo ">>> 4.3 Hardy–Weinberg 平衡过滤 (hwe P ≥ 1e-5)"
plink --bfile Han2       --hwe 0.00001 --make-bed --out Han3
plink --bfile Tibetan3   --hwe 0.00001 --make-bed --out Tibetan4

echo ">>> 4.4 SNP 子集 QC"
plink --bfile Han_SNP1    --geno 0.05    --make-bed --out Han_SNP2
plink --bfile Han_SNP2    --maf 0.01     --make-bed --out Han_SNP3
plink --bfile Han_SNP3    --hwe 0.00001  --make-bed --out Han_SNP4
plink --bfile Tibetan_SNP1--geno 0.05    --make-bed --out Tibetan_SNP2
plink --bfile Tibetan_SNP2--maf 0.01     --make-bed --out Tibetan_SNP3
plink --bfile Tibetan_SNP3--hwe 0.00001  --make-bed --out Tibetan_SNP4

#----------------------------------------
# 5. 主成分分析（PCA）
#----------------------------------------

echo ">>> 5.1 计算前 20 个主成分"
plink --bfile Han3        --pca 20 --out Han3
plink --bfile Tibetan4    --pca 20 --out Tibetan4
plink --bfile Han_SNP4    --pca 20 --out Han_SNP4
plink --bfile Tibetan_SNP4--pca 20 --out Tibetan_SNP4

echo ">>> 5.2 用 twstats 评估显著主成分"
twstats -t twtable -i Han3.eigenval        -o Han3.out
twstats -t twtable -i Tibetan4.eigenval    -o Tibetan4.out
twstats -t twtable -i Han_SNP4.eigenval    -o Han_SNP4.out
twstats -t twtable -i Tibetan_SNP4.eigenval- o Tibetan_SNP4.out

#----------------------------------------
# 6. Logistic 回归
#----------------------------------------

echo ">>> 6.1 对 Han3 数据进行 logistic 回归（含缺失年龄补均值）"
plink --bfile Han3 \
      --logistic \
      --covar Han3.cov \
      --hide-covar \
      --sex \
      --out Han3

echo ">>> 6.2 仅使用 PCA 协变量"
plink --bfile Han3 \
      --logistic \
      --covar PCA.cov \
      --hide-covar \
      --sex \
      --out Han3_PCAonly

echo ">>> 6.3 对藏族数据 run logistic"
plink --bfile Tibetan4 \
      --logistic \
      --covar Tibetan4.cov \
      --hide-covar \
      --sex \
      --out Tibetan4

echo ">>> 6.4 输出回归系数（--beta）"
plink --bfile Han_SNP4 \
      --logistic \
      --beta \
      --covar Han_SNP4.cov \
      --hide-covar \
      --sex \
      --out Han_SNP5

plink2 --bfile Han_SNP4 \
       --glm hide-covar sex \
       --covar-variance-standardize Han_SNP4.cov \
       --out Han_SNP4_plink2

plink --bfile Tibetan_SNP4 \
      --logistic \
      --beta \
      --covar Tibetan_SNP4.cov \
      --hide-covar \
      --sex \
      --out Tibetan_SNP5

plink2 --bfile Tibetan_SNP4 \
       --glm hide-covar sex \
       --covar-variance-standardize Tibetan_SNP4.cov \
       --out Tibetan_SNP4_plink2

#----------------------------------------
# 7–10. R 脚本：QQ 图、曼哈顿图、SNPID 转换与分列
#----------------------------------------

echo ">>> 7–10. 生成 R 脚本并执行"
cat <<'EOF' > gwas_plots_and_format.R
#!/usr/bin/env Rscript
#===============================================================================
# R 脚本：生成 QQ 图、曼哈顿图；SNPID 转换与列拆分
#===============================================================================

# 7. QQ 图
library(gap)
dat1 <- read.table("Han3.assoc.logistic", header=TRUE)
gc1  <- gcontrol2(dat1$P)
cat("QQ lambda for Han3:", gc1$lambda, "\n")
# 画 QQ 图
pdf("QQ_Han3.pdf")
gap::qq(dat1$P, main="QQ Plot - Han3")
dev.off()

# 8. 曼哈顿图
library(qqman)
# Han3
dat2 <- read.table("Han3.assoc.logistic", header=TRUE)
pdf("Manhattan_Han3.pdf", width=10, height=6)
manhattan(dat2, col=c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62"),
          main="Manhattan Plot - Han3", cex=0.8)
dev.off()
# SNP4
dat3 <- read.table("Han_SNP4_plink2.PHENO1.glm.logistic", header=TRUE)
pdf("Manhattan_Han_SNP4.pdf", width=10, height=6)
manhattan(dat3, col=c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62"),
          main="Manhattan Plot - Han_SNP4", cex=0.8)
qq(dat3$P)
dev.off()

# 9. SNPID 转换
han <- read.table("Han_plink2.PHENO1.glm.logistic", header=TRUE)
han$SNPID <- paste(han$CHR, han$BP, sep="_")
han <- han[ , c("CHR","BP","SNPID", setdiff(names(han), c("CHR","BP"))) ]
write.table(han, file="Han_plink2.logistic", sep="\t", quote=FALSE, row.names=FALSE)

# 10. 分列操作
aa <- read.table("Han_SNP4_plink2.PHENO1.glm.logistic", header=TRUE)
library(dplyr); library(tidyr)
aa$SNPID <- aa$MarkerName
b <- aa %>% separate(SNPID, into=c("CHR","BP"), sep="_")
write.table(b, file="Han_SNP4_splitBP.logistic", sep="\t", quote=FALSE, row.names=FALSE)
EOF

chmod +x gwas_plots_and_format.R
./gwas_plots_and_format.R

echo ">>> 全部步骤完成！"
