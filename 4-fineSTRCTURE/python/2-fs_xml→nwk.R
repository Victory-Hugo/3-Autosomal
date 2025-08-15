#!/usr/bin/env Rscript

## 1) 载入包
suppressMessages({
  if (!requireNamespace("xml2", quietly=TRUE)) install.packages("xml2")
  if (!requireNamespace("ape",   quietly=TRUE)) install.packages("ape")
  if (!requireNamespace("treeio",quietly=TRUE)) install.packages("treeio")
})
library(xml2)
library(ape)
library(treeio)

## 2) 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("用法: convert_tree.R <输入 tree.xml 路径> <输出 .nwk 路径>\n")
  q(status = 1)
}
treefile <- args[1]
outfile  <- args[2]

## 3) 读取 XML 并提取 <Tree> 节点文本
xml_doc   <- read_xml(treefile)
tree_text <- xml_text(xml_find_all(xml_doc, ".//Tree"))

if (length(tree_text) == 0) {
  stop("在 XML 中未找到任何 <Tree> 节点，请检查输入文件。")
}

## 4) 解析为 phylo，并写出 Newick 文件
phy <- read.tree(text = tree_text)
write.tree(phy, file = outfile)

cat("已将", treefile, "中的树保存为 Newick 格式至", outfile, "\n")
