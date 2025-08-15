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



library(ggtreeExtra)
library(ggnewscale)
library(MetBrewer)
# 读取数据
df_meta <- read_csv("Presence_absence.csv")
df_meta

df_meta -> df_meta_long

df_meta_long
# 检查CagA_subtype的唯一值数量
unique_subtypes <- length(unique(df_meta_long$CagA_subtype))
cat("CagA_subtype 唯一值数量:", unique_subtypes, "\n")

# 解决方案1: 使用更多颜色的色板
p4 <- p3 + 
    new_scale_fill() +
    geom_fruit(
      data = df_meta_long,
      geom = geom_tile,
      mapping = aes(y = label, fill = CagA_subtype),
      width = 0.2,
      offset = 0.01,
      axis.params = list(axis = "x", text.angle = -45, hjust = 0)
    ) + 
    scale_fill_manual(values = c(
      met.brewer("Troy", n = min(unique_subtypes, 16)),  
      met.brewer("Cassatt2", n = max(0, unique_subtypes - 5)), 
      rainbow(max(0, unique_subtypes - 24))  # 如果还不够用彩虹色
    )) + 
    labs(fill = "CagA Subtype")  #! 自定义图例标题
p4

ggsave(p4, filename = '系统发育树_汉藏_cagA.pdf', units = "cm", height = 10, width = 20)

library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(ggstar)
library(ggnewscale)

# 树文件路径
trfile <- system.file("extdata", "tree.nwk", package = "ggtreeExtra")
# 提示点数据文件路径
tippoint1 <- system.file("extdata", "tree_tippoint_bar.csv", package = "ggtreeExtra")
# 第一层环数据文件路径
ring1 <- system.file("extdata", "first_ring_discrete.csv", package = "ggtreeExtra")
# 第二层环数据文件路径
ring2 <- system.file("extdata", "second_ring_continuous.csv", package = "ggtreeExtra")

# 使用 read.tree 导入树文件
tree <- read.tree(trfile)

# 读取注释数据
dat1 <- read.csv(tippoint1)
knitr::kable(head(dat1))
dat2 <- read.csv(ring1)
knitr::kable(head(dat2))
dat3 <- read.csv(ring2)
knitr::kable(head(dat3))

# 绘制基础树
p <- ggtree(tree, layout = "fan", open.angle = 10, size = 0.5)

# 星状点图层
p2 <- p +
  geom_fruit(
    data = dat1,
    geom = geom_star,
    mapping = aes(
      y = ID,
      fill = Location,
      size = Length,
      starshape = Group
    ),
    position = "identity",
    starstroke = 0.5,
    color = "#FFFFFF"
  ) +
  scale_size_continuous(
    range = c(1, 3),
    guide = guide_legend(
      keywidth = 0.5,
      keyheight = 0.5,
      override.aes = list(starshape = 15),
      order = 2
    )
  ) +
  scale_fill_manual(
    values = c(
      "#F8766D", "#C49A00", "#53B400",
      "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"
    ),
    guide = "none"
  ) +
  scale_starshape_manual(
    values = c(1, 15),
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 1)
  )

# 第一层离散热图
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = dat2,
    geom = geom_tile,
    mapping = aes(y = ID, x = Pos, fill = Type),
    offset = 0.08,
    colour = "white",
    linewidth = 0.2,
    pwidth = 0.25,
    axis.params = list(axis = "x", text.angle = -45, hjust = 0)
  ) +
  scale_fill_manual(
    values = c("#339933", "#dfac03"),
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 3)
  )

# 第二层连续热图
p4 <- p3 +
  new_scale_fill() +
  geom_fruit(
    data = dat3,
    geom = geom_tile,
    mapping = aes(y = ID, x = Type2, alpha = Alpha, fill = Type2),
    pwidth = 0.15,
    axis.params = list(axis = "x", text.angle = -45, hjust = 0)
  ) +
  scale_fill_manual(
    values = c("#b22222", "#005500", "#0000be", "#9f1f9f"),
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 4)
  ) +
  scale_alpha_continuous(
    range = c(0, 0.4),
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 5)
  )

# 条形图层
p5 <- p4 +
  new_scale_fill() +
  geom_fruit(
    data = dat1,
    geom = geom_col,
    mapping = aes(y = ID, x = Abundance, fill = Location),
    pwidth = 0.4,
    axis.params = list(axis = "x", text.angle = -45, hjust = 0)
  ) +
  scale_fill_manual(
    values = c(
      "#F8766D", "#C49A00", "#53B400",
      "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"
    ),
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 6)
  ) +
  theme(
    legend.background = element_rect(fill = NA),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.spacing.y = unit(0.02, "cm")
  )

p5