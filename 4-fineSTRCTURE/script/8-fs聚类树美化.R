library(ggtree)
library(treeio)
library(tidyverse)

# 设置工作目录
setwd("/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/")
# 加载树文件
WGS_tree <- read.tree("整体1369_HP_linked_tree.nwk")
#* ================计算从节点到叶子的距离（树的高度 - 节点深度）================
#? BEAST等软件不需要计算，自带
node_heights <- node.depth.edgelength(WGS_tree)
max_height <- max(node_heights)
node_to_tip_distance <- max_height - node_heights
WGS_tree |> as_tibble() |> select(label) |> drop_na() |> write.csv("整体1369_HP_linked_tree_labels.csv", row.names = FALSE)
#* =========================分支着色(如果没有分组请删除下列的aes(color = group))=========================
#* =========================分支着色(如果没有分组请删除下列的scale_color_manual(values = group_colors) =========================
# 文件路径可按需修改
group_file <- "group.txt" #todo 第一列是label(ID)，第二列是group(单倍群)，不包含表头
color_file <- "color.txt" #todo 第一列是group(单倍群)，第二列是color，不包含表头

# 读取分组信息
group_df <- read.table(group_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("label", "group"))
group_list <- split(group_df$label, group_df$group) # 将分组信息转换为列表

# 读取颜色信息
color_df <- read.table(color_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', col.names = c("group", "color"))
group_colors <- setNames(color_df$color, color_df$group) # 将颜色信息转换为命名向量

# 读取树并分组
tree_grouped <- groupOTU(WGS_tree, group_list)

# 找到汉族和藏族的MRCA
ggtree::MRCA(tree_grouped,2650,2555)

#* ================使用ggtree绘制系统发育树================
#* 使用ggtree绘制系统发育树
p0 <- ggtree(tree_grouped, aes(color = group), layout = 'dendrogram',  lwd = 0.2) + #* 如果没有分支着色请删除aes(color = group))
  scale_color_manual(values = group_colors) + #* 如果没有分支着色请删除scale_color_manual(values = group_colors)
  geom_tiplab(size = 0.5, align = TRUE, hjust = -0.5) + #* 添加tiplab，默认显示`label`
  #geom_nodepoint(size = 0.1) +   #* 添加nodepoint
  # geom_nodelab(aes(label = round(branch.length, 3)), hjust = -0.3, size = 0.7) + #* 显示branch.length
  #geom_nodelab(aes(label = node),size = 0.7) + #* 显示branch.length +
  geom_highlight(mapping=aes(subset = node %in% c(2285)),fill = "#56B6C2", #* 高亮显示特定节点,c(节点);
                        type = "gradient", #* 渐变
                        gradient.direction = 'rt', #* 渐变方向
                        alpha = .1) +
  theme(
    axis.text.x = element_text(size = 2),    #* x轴刻度文字大小
    axis.title.x = element_text(size = 2),    #* x轴标题文字大小
    axis.text.y = element_text(size = 2),   #* y轴刻度文字大小
    axis.title.y = element_text(size = 2)     #* y轴标题文字大小
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.03))
  ) + #* 扩展x轴空间，左边2%，右边3%
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.01))) + #* 翻转 y 轴并扩展y轴空间，上下各1%
  coord_flip()  #* 翻转坐标轴，实现从上往下分叉
p0 + layout_dendrogram() -> p0 
p0


#* 保存
ggsave(p0, filename = '系统发育树.pdf', units = "cm", height = 15, width = 30 )

# 方法1: 使用viewClade查看特定节点的子树（推荐）
p2 <- viewClade(p0, 2285) #* 显示节点2285及其后代
p2

# 方法2: 使用tree_subset从原始树提取子树然后重新绘制
library(ape) #* 需要ape包的drop.tip函数
# 获取节点2285的所有后代tip
descendants_2285 <- offspring(WGS_tree, 2285, tiponly = TRUE)
# 提取包含这些tip的子树
sub_tree <- keep.tip(WGS_tree, descendants_2285)
# 重新分组并绘制子树
sub_tree_grouped <- groupOTU(sub_tree, group_list)
p3 <- ggtree(sub_tree_grouped, aes(color = group), layout = 'dendrogram',  lwd = 0.2) + #* 如果没有分支着色请删除aes(color = group))
  scale_color_manual(values = group_colors) + #* 如果没有分支着色请删除scale_color_manual(values = group_colors)
  geom_tiplab(size = 0.5, align = TRUE, hjust = -0.5) + #* 添加tiplab，默认显示`label`
  #geom_nodepoint(size = 0.1) +   #* 添加nodepoint
  # geom_nodelab(aes(label = round(branch.length, 3)), hjust = -0.3, size = 0.7) + #* 显示branch.length
  #geom_nodelab(aes(label = node),size = 0.7) + #* 显示branch.length +
  theme(
    axis.text.x = element_text(size = 2),    #* x轴刻度文字大小
    axis.title.x = element_text(size = 2),    #* x轴标题文字大小
    axis.text.y = element_text(size = 2),   #* y轴刻度文字大小
    axis.title.y = element_text(size = 2)     #* y轴标题文字大小
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.03))
  ) + #* 扩展x轴空间，左边2%，右边3%
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.01))) + #* 翻转 y 轴并扩展y轴空间，上下各1%
  coord_flip()  #* 翻转坐标轴，实现从上往下分叉
p3 + layout_dendrogram() -> p3 
p3 
# 保存子树图
ggsave(p3, filename = '系统发育树_汉藏.pdf', units = "cm", height = 10, width = 20)

# 读取数据
df_meta <- read_csv("Presence_absence.csv")
df_meta

df_meta |> pivot_longer(
    cols = starts_with("hp"),
    names_to = "gene",
    values_to = "presence"
) -> df_meta_long

df_meta_long

library(ggtreeExtra)
library(ggnewscale)
library(MetBrewer)
p4 <- p3 + 
    new_scale_fill() +
    geom_fruit(
      data = df_meta_long,
      geom = geom_tile,
      mapping = aes(y = label, x = gene, fill = presence),
      width = 0.07,
      offset = 0.01,
      axis.params = list(
        axis = "x",
        text.angle = 0,
        hjust = -0.5
      ) #* 把基因名称添加
    ) + 
    scale_fill_met_c(name = "Homer2")  #! 注意fill/color
p4
ggsave(p4, filename = '系统发育树_汉藏_cagA.pdf', units = "cm", height = 10, width = 20)

