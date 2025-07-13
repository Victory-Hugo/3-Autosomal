#!/usr/bin/Rscript

#======================================================================
# TreeMix结果分析与可视化脚本
# 使用说明：
# 这个脚本已经设置了固定路径，可以在R环境中逐行运行
# 适用于GlobalTreemix_m{m}_rep{r}命名格式的TreeMix输出文件
#======================================================================

# ======= 路径和参数设置（可根据需要修改） =======

# TreeMix输出文件目录
treemix_dir <- "/home/luolintao/Helicopter/Script/分析结果/treemix/output/global"

# 种群文件路径
pop.uniq <- "/home/luolintao/S00-Github/3-Autosomal/1-treemix/conf/pop.cov" 

# TreeMix输出文件前缀
file_prefix <- "GlobalTreemix" 

# 结果输出目录
results_dir <- file.path(treemix_dir, "analysis")

# 创建结果目录（如果不存在）
if(!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 切换到TreeMix输出目录
setwd(treemix_dir)
cat(sprintf("工作目录设为：%s\n", getwd()))

# 加载TreeMix绘图函数库
source('/home/luolintao/S00-Github/3-Autosomal/1-treemix/func/plotting_funcs.R')

# 检查种群列表文件是否存在
if(!file.exists(pop.uniq)){
  cat("错误：种群文件不存在：", pop.uniq, "\n")
  cat("请提供有效的种群名称列表文件\n")
  stop("程序因上述错误停止运行")
}

# ======= 文件分析 =======
cat("分析输出文件结构...\n")
all_files <- list.files(path=treemix_dir, pattern="covse.gz")

if(length(all_files) == 0) {
  stop(sprintf("在目录 %s 未找到TreeMix输出文件（*.covse.gz）。请检查目录或文件名。", treemix_dir))
}

# 提取不同的迁移边数量
m_values <- unique(gsub(".*_m([0-9]+)_rep.*", "\\1", all_files))
m_values <- sort(as.numeric(m_values))
cat(sprintf("发现迁移边数量: %s\n", paste(m_values, collapse=", ")))

# 确定每个迁移边数量的重复数
rep_counts <- sapply(m_values, function(m) {
  files <- list.files(path=treemix_dir, pattern=paste0("_m", m, "_rep[0-9]+\\.covse\\.gz$"))
  max_rep <- max(as.numeric(gsub(".*_rep([0-9]+)\\.covse\\.gz$", "\\1", files)))
  return(max_rep)
})
names(rep_counts) <- m_values

cat("迁移边数量及其重复次数:\n")
for(m in m_values) {
  cat(sprintf("  m=%s: %d次重复\n", m, rep_counts[m]))
}

max_m <- max(as.numeric(m_values))
cat(sprintf("最大迁移边数量: m=%d\n", max_m))

#======================================================================
# 第一部分：生成各迁移边数量的树结构图
#======================================================================

# 这一部分暂时保留使用第一个重复实验绘图，最佳重复会在后面计算
# 最佳重复的选择将在计算方差解释量的步骤中完成
cat("绘制迁移边树结构图（使用第一个重复实验）...\n")
for(m in m_values){
  rep <- 1  # 默认使用第一次重复的结果，后面会更新为最佳重复
  
  prefix <- sprintf("%s_m%s_rep%d", file_prefix, m, rep)
  pdf_file <- file.path(results_dir, paste0("treemix.m", m, ".pdf"))
  cat(sprintf("  生成图片：%s (使用 %s)\n", basename(pdf_file), prefix))
  
  pdf(file=pdf_file, width=12, height=8)
  plot_tree(prefix)
  dev.off()
}

#======================================================================
# 第二部分：选择最佳重复实验并计算方差解释量
#======================================================================

# 定义函数用于提取.llik文件中的似然值
extract_likelihood <- function(llik_file) {
  if(!file.exists(llik_file)) {
    return(NA)
  }
  
  lines <- readLines(llik_file)
  # 查找"Exiting ln(likelihood)"行，通常是最后的似然值
  exit_line <- grep("Exiting ln\\(likelihood\\)", lines, value=TRUE)
  
  if(length(exit_line) > 0) {
    # 提取似然值（最后一个数字）
    likelihood <- as.numeric(gsub(".*: ([0-9.-]+).*", "\\1", exit_line))
    return(likelihood)
  }
  
  return(NA)
}

# 为每个迁移边数量找到最佳重复实验
cat("为每个迁移边数量查找似然值最高的重复实验...\n")
best_reps <- list()

for(m in m_values) {
  # 收集所有重复实验的似然值
  likelihoods <- numeric(rep_counts[m])
  
  for(r in 1:rep_counts[m]) {
    llik_file <- file.path(treemix_dir, sprintf("%s_m%s_rep%d.llik", file_prefix, m, r))
    likelihoods[r] <- extract_likelihood(llik_file)
  }
  
  # 找出似然值最高的重复实验
  if(all(is.na(likelihoods))) {
    best_rep <- 1  # 如果无法获取似然值，默认使用第一个重复
    cat(sprintf("  警告: m=%s 无法获取任何重复实验的似然值，默认使用 rep=1\n", m))
  } else {
    valid_likelihoods <- likelihoods[!is.na(likelihoods)]
    if(length(valid_likelihoods) > 0) {
      best_rep <- which.max(likelihoods)
      cat(sprintf("  m=%s: 最佳重复实验是 rep=%d (似然值: %.2f)\n", 
                 m, best_rep, likelihoods[best_rep]))
    } else {
      best_rep <- 1
      cat(sprintf("  警告: m=%s 所有似然值都是NA，默认使用 rep=1\n", m))
    }
  }
  
  best_reps[[as.character(m)]] <- best_rep
}

cat("\n计算各迁移边数量的方差解释量...\n")

# 计算每个迁移边数量的方差解释量（使用最佳重复实验）
m_results <- numeric(length(m_values))
names(m_results) <- m_values
best_prefixes <- character(length(m_values))
names(best_prefixes) <- m_values

for(i in 1:length(m_values)){
  m <- m_values[i]
  rep <- best_reps[[as.character(m)]]  # 使用找到的最佳重复实验
  
  prefix <- sprintf("%s_m%s_rep%d", file_prefix, m, rep)
  best_prefixes[i] <- prefix
  m_results[i] <- get_f(prefix)
  cat(sprintf("  m=%s (rep=%d) 方差解释量: %.4f%%\n", m, rep, m_results[i]*100))
}

# 保存方差解释量结果
results_table <- data.frame(
  migration_edges = m_values,
  best_rep = sapply(m_values, function(m) best_reps[[as.character(m)]]),
  best_prefix = best_prefixes,
  variance_explained = m_results * 100
)
csv_file <- file.path(results_dir, "variance_explained.csv")
write.table(results_table, csv_file, 
            quote=FALSE, sep=",", row.names=FALSE)
cat(sprintf("方差解释量结果已保存至 %s\n", csv_file))

# 找出达到99.8%阈值的最小迁移边数量
threshold <- 99.8
threshold_m <- min(m_values[m_results*100 >= threshold], Inf)
if(is.infinite(threshold_m)) {
  cat(sprintf("警告：没有任何迁移边数量达到 %.1f%% 的方差解释量阈值\n", threshold))
  cat(sprintf("最大方差解释量: %.4f%% (m=%s)\n", max(m_results)*100, m_values[which.max(m_results)]))
} else {
  cat(sprintf("达到 %.1f%% 方差解释量阈值的最小迁移边数量: m=%d\n", threshold, threshold_m))
  
  # 创建最佳模型的副本
  best_rep <- best_reps[[as.character(threshold_m)]]
  best_tree_file <- sprintf("%s_m%s_rep%d.treeout.gz", file_prefix, threshold_m, best_rep)
  if(file.exists(file.path(treemix_dir, best_tree_file))) {
    file.copy(
      from = file.path(treemix_dir, best_tree_file),
      to = file.path(results_dir, "best_model.treeout.gz"),
      overwrite = TRUE
    )
    cat(sprintf("已将最佳模型文件复制到: %s\n", file.path(results_dir, "best_model.treeout.gz")))
    
    # 同时重新生成最佳模型的树图
    best_prefix <- sprintf("%s_m%s_rep%d", file_prefix, threshold_m, best_rep)
    pdf_file <- file.path(results_dir, paste0("best_treemix_m", threshold_m, "_rep", best_rep, ".pdf"))
    cat(sprintf("生成最佳模型图片：%s\n", basename(pdf_file)))
    pdf(file=pdf_file, width=12, height=8)
    plot_tree(best_prefix)
    dev.off()
  }
}

# 绘制迁移边数量与方差解释量关系图
cat("绘制迁移边数量选择图...\n")
pdf_file <- file.path(results_dir, "migration.edge.choice.pdf")
pdf(file=pdf_file, width=10, height=7)
plot(as.numeric(m_values), m_results*100, 
     pch=19, cex=1.2, col="blue", type="b",
     xlab="迁移边数量", ylab="方差解释量(%)",
     xaxt="n", # 不显示默认x轴刻度
     main="TreeMix方差解释量分析")

# 添加自定义x轴标签
axis(1, at=as.numeric(m_values), labels=m_values)

# 添加99.8%阈值参考线（根据文献要求）
abline(h=threshold, col="red", lty=2)
text(as.numeric(m_values)[2], threshold+0.5, 
     sprintf("%.1f%% 阈值", threshold), col="red")

# 高亮显示达到阈值的最小迁移边数量
if(!is.infinite(threshold_m)) {
  points(threshold_m, m_results[as.character(threshold_m)]*100, 
         pch=19, col="red", cex=2)
  text(threshold_m, m_results[as.character(threshold_m)]*100-1.5, 
       sprintf("m=%d", threshold_m), col="red")
}

# 添加方差解释量标签
for(i in 1:length(m_values)) {
  text(as.numeric(m_values)[i], m_results[i]*100+0.7, 
       sprintf("%.2f%%", m_results[i]*100), cex=0.8)
}

dev.off()
cat(sprintf("图形已保存到: %s\n", pdf_file))

#======================================================================
# 第三部分：残差分析图
#======================================================================

cat("生成残差热图...\n")
# 为每个迁移边数量绘制残差图（使用最佳重复实验）
for(m in m_values){
  rep <- best_reps[[as.character(m)]]
  
  prefix <- sprintf("%s_m%s_rep%d", file_prefix, m, rep)
  pdf_file <- file.path(results_dir, sprintf("treemix.residuals.m%s.pdf", m))
  cat(sprintf("  生成图片：%s (使用最佳重复 %s)\n", basename(pdf_file), prefix))
  
  pdf(file=pdf_file)
  plot_resid(prefix, pop.uniq)
  dev.off()
  
  # 如果有达到阈值的最佳m值，为其创建额外的高质量残差图
  if(!is.infinite(threshold_m) && m == threshold_m) {
    hi_res_file <- file.path(results_dir, "best_model_residuals.pdf")
    cat(sprintf("  生成最佳模型残差图：%s\n", basename(hi_res_file)))
    pdf(file=hi_res_file, width=12, height=10)
    plot_resid(prefix, pop.uniq)
    dev.off()
  }
}

# 创建汇总报告
report_file <- file.path(results_dir, "treemix_analysis_summary.txt")
cat("创建分析结果汇总报告...\n")

report_con <- file(report_file, "w")
cat("=========================================\n", file=report_con)
cat("TreeMix 方差解释量分析汇总报告\n", file=report_con)
cat("=========================================\n\n", file=report_con)
cat(sprintf("分析日期: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file=report_con)
cat("分析的迁移边数量:\n", file=report_con)
for(i in 1:length(m_values)) {
  cat(sprintf("  m=%s: %.4f%%\n", m_values[i], m_results[i]*100), file=report_con)
}
cat("\n", file=report_con)

if(!is.infinite(threshold_m)) {
  best_rep <- best_reps[[as.character(threshold_m)]]
  cat(sprintf("达到 %.1f%% 方差解释量阈值的最小迁移边数量: m=%d\n", threshold, threshold_m), file=report_con)
  cat(sprintf("最佳模型: %s_m%s_rep%d\n", file_prefix, threshold_m, best_rep), file=report_con)
  cat(sprintf("方差解释量: %.4f%%\n", m_results[as.character(threshold_m)]*100), file=report_con)
} else {
  best_m <- m_values[which.max(m_results)]
  best_rep <- best_reps[[as.character(best_m)]]
  cat(sprintf("警告：没有任何迁移边数量达到 %.1f%% 的方差解释量阈值\n", threshold), file=report_con)
  cat(sprintf("最大方差解释量: %.4f%% (m=%s, rep=%d)\n", 
              max(m_results)*100, best_m, best_rep), file=report_con)
}
close(report_con)

cat("所有分析完成！\n")
cat(sprintf("结果基于 %d 种迁移边数量分析\n", length(m_values)))
cat(sprintf("所有结果已保存到目录: %s\n", results_dir))