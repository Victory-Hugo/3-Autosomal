# ======================================================================
# TreeMix OptM 分析脚本 - 优化版本
# 用于确定最优迁移边数量
# ======================================================================

# 检查并加载必要的包
required_packages <- c("OptM", "tidyverse")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("正在安装包: %s\n", pkg))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("=== TreeMix OptM 分析开始 ===\n")

# 设置目录
setwd('/home/luolintao/Helicopter/Script/分析结果/treemix/output/China/')
# 确保输出目录存在
if (!dir.exists("./analysis")) {
  dir.create("./analysis")
}

# 检查 TreeMix 输出文件是否存在
cat("检查 TreeMix 输出文件...\n")
treemix_files <- list.files("./", pattern = "*.llik$")
if (length(treemix_files) == 0) {
  stop("错误：当前目录中未找到 TreeMix 输出文件 (*.llik)。请确认路径正确。")
}
cat(sprintf("找到 %d 个 TreeMix 输出文件\n", length(treemix_files)))

# 如果文件前缀不是"global"，请修改为实际的前缀
# 分析数据 - 使用linear方法
cat("正在运行 OptM 分析...\n")
linear = optM("./")

# 绘制并保存图形 - 使用 plot_optM 内置的 pdf 参数
cat("正在生成优化图表...\n")
plot_optM(linear, pdf = "./analysis/optM_linear.pdf")


#*保存结果至csv文件
cat("正在保存结果到 CSV 文件...\n")
linear |> as_tibble() -> df_linear

# —— 新增 —— 用 Deltam 最大值来确定最优 m
best_m <- df_linear |>
  filter(Deltam == max(Deltam, na.rm = TRUE)) |>
  slice(1) |>           # 如果有并列，取第一个
  pull(m)

print(paste("根据Deltam 最大值确定最佳的迁移边为:", best_m))
# 保存数据框到CSV文件
write.csv(x = df_linear, file = "./analysis/df_linear.csv", 
          row.names = FALSE, quote = TRUE)

# 打印结果摘要
cat("分析完成！\n")
cat("结果已保存到:\n")
cat("  - 图形文件: ./analysis/optM_linear.pdf\n")
cat("  - 数据文件: ./analysis/df_linear.csv\n")
cat("\n数据摘要:\n")
print(df_linear)

