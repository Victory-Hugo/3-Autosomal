# 介绍

TreeMix 是一款基于等位基因频率的群体遗传分析工具，可推断群体分化历史并检测基因流事件。通过构建树状图和添加迁移边，TreeMix 直观展示群体进化关系及基因交流，广泛应用于群体遗传学和进化生物学研究。

# 分析流程

1. **VCF 转换**  
    使用 `PLINK` 将 VCF 文件转换为 `bed`、`bim`、`fam` 格式。

2. **编辑 fam 文件**  
    将 `fam` 文件第一列替换为群体名。可用脚本：`1-treemix/python/0-fam_group.py`。

3. **编辑 bim 文件**  
    将 `bim` 文件第二列替换为 SNP 位置（即复制第四列到第二列）。可用脚本：`1-treemix/python/1-bim-site.py`。

4. **计算等位基因频率**  
    用 `PLINK` 计算每个群体的等位基因频率，脚本：`1-treemix/script/1-preprocess-makecluster.sh`。  
    - 需准备与 `bed` 前缀相同的 csv 文件。  
    - csv 第1、3列为群体名，第2列为 ID。  
    - 用空格分隔，且 ID 顺序需与 `bed`、`bim`、`fam` 文件一致。

5. **格式转换**  
    用 `1-treemix/python/plink2treemix.py` 或 `python2-plink2treemix.py` 将 `plink.frq.strat.gz` 转为 `TreeMix.gz`。

6. **运行 TreeMix**  
    执行 `1-treemix/script/1-treemix-global.sh`，根据实际情况调整参数。

7. **结果分析**  
    - 运行 `1-treemix/script/2-Result.R` 分析结果。  
    - 可视化树状图和迁移事件。
    - 残差分析。