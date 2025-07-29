# -*- coding: utf-8 -*-
# author: 自动生成
# 功能：读取“Recipient + 多列原始计数”格式的文件，
#       归一化成比例，按照 color.csv 手动映射颜色，绘制默认堆叠图
#* =============================
#* ===========基础功能===========
#* =============================
# 1. 载入依赖包
library(reshape2)
library(tidyverse)
# 2. 参数设置
input_file   <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/6-fineSTRUCTURE/output/东亚受到周边/Chunkcounts.txt"   # 【改成你的数据文件】
color_file   <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/conf/color.csv"
output_file  <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/6-fineSTRUCTURE/output/东亚受到周边/structure_default.pdf"

# 3. 读取并归一化
df <- read.table(input_file, header=TRUE, stringsAsFactors=FALSE)
cluster_cols <- colnames(df)[-1]
df_norm <- df
df_norm[,cluster_cols] <- df[,cluster_cols] / rowSums(df[,cluster_cols])

# 4. 整形数据（长表格）
df_melt <- melt(df_norm,
                id.vars       = "Recipient",
                variable.name = "Cluster",
                value.name    = "Proportion")

# 5. 读取颜色映射文件
colormap_df <- read.csv(color_file, stringsAsFactors = FALSE)
color_map   <- setNames(colormap_df$color, colormap_df$Class)

# 筛选出数据中实际出现的簇
present_clusters   <- unique(df_melt$Cluster)
# 按照 color.csv 的顺序，只保留出现在 present_clusters 中的类别
present_levels     <- intersect(colormap_df$Class, present_clusters)
color_map_present  <- color_map[present_levels]

# 6. 重新定义因子水平，保证画图顺序正确
df_melt$Cluster <- factor(df_melt$Cluster, levels = present_levels)

# 7. 绘图并保存
pdf(output_file, width=12, height=6)
ggplot(df_melt, aes(x=Recipient, y=Proportion, fill=Cluster)) +
  geom_bar(stat="identity", width=1, color=NA) +
  # 只给 present_levels 这几类分配颜色，legend 只会显示它们
  scale_fill_manual(values = color_map_present) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) 
dev.off()

#* =============================
#* ===========分组功能===========
#* =============================
# 读取excel文件
library(readxl)
df_meta <- read_excel("/mnt/f/OneDrive/文档（科研）/P0_幽门螺旋杆菌/基础信息表/HP数据收集2.xlsx", sheet = "HP数据收集")

# 选取ID列和经纬度
df_meta |> 
    filter(`7544个样本` == 'YES') |>
    select(ID,Latitude,Longitude) -> df_7544

# join数据
df_joined <- left_join(df_melt, df_7544, by = c("Recipient" = "ID")) |> as_tibble()

# 移除缺失经纬度的样本
df_joined <- df_joined |> 
    filter(!is.na(Latitude) & !is.na(Longitude))

# 创建经纬度组合标识符
df_joined <- df_joined |> 
    mutate(LatLon_Group = paste(round(Latitude, 4), round(Longitude, 4), sep = "_"))

# 按照每个独特的经纬度组合，在地图上绘制堆叠柱状图
# 统计每个经纬度组合下的样本数量
location_counts <- df_joined |> 
    group_by(LatLon_Group, Latitude, Longitude) |> 
    summarise(n_samples = n_distinct(Recipient), .groups = "drop") |>
    arrange(desc(n_samples))

print("经纬度组合及样本数量：")
print(location_counts)

# 加载地图相关包
library(sf)

# 读取世界地图
world_shp <- read_sf("/mnt/f/OneDrive/文档（科研）/地图/全球行政地图(精细版)/全球国家边界/世界地图国家.shp")

# 准备地图上绘制的数据：计算每个位置的平均比例
df_map_data <- df_joined %>%
    group_by(LatLon_Group, Latitude, Longitude, Cluster) %>%
    summarise(Mean_Proportion = mean(Proportion), .groups = "drop")

# 创建基础地图
base_map <- ggplot() +
    geom_sf(data = world_shp, 
            fill = "#f0f0f0", 
            color = "white", 
            linewidth = 0.2) +
    coord_sf(xlim = c(-180, 180), 
             ylim = c(-60, 80), 
             expand = FALSE) +
    theme_void() +
    theme(
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "bottom",
        legend.title = element_blank()
    )

# 设置柱状图参数
sample_bar_width <- 0.05   # 每个样本柱子的宽度（经度单位）
bar_height_scale <- 8      # 柱状图高度缩放因子
location_spread <- 1.5     # 每个位置样本分布的范围（经度单位）

# 只选择样本数量较多的位置
top_locations <- location_counts %>% 
    filter(n_samples >= 10) %>%  # 显示样本数>=10的位置
    slice_head(n = 30)  # 最多显示前30个位置

print(paste("将在地图上显示", nrow(top_locations), "个主要位置"))

# 为绘制每个样本的堆叠柱状图准备数据
map_bars_data <- data.frame()

for(i in 1:nrow(top_locations)) {
    current_group <- top_locations$LatLon_Group[i]
    current_lat <- top_locations$Latitude[i]
    current_lon <- top_locations$Longitude[i]
    current_n <- top_locations$n_samples[i]
    
    # 获取当前位置的所有样本数据
    current_samples <- df_joined %>%
        filter(LatLon_Group == current_group) %>%
        arrange(Recipient, match(Cluster, present_levels))
    
    # 获取该位置的唯一样本列表
    unique_samples <- unique(current_samples$Recipient)
    
    # 为每个样本分配x位置（在位置周围分布）
    if(length(unique_samples) == 1) {
        sample_x_positions <- current_lon
    } else {
        sample_x_positions <- seq(
            from = current_lon - location_spread/2,
            to = current_lon + location_spread/2,
            length.out = length(unique_samples)
        )
    }
    
    # 为每个样本创建堆叠柱状图数据
    for(j in 1:length(unique_samples)) {
        sample_id <- unique_samples[j]
        sample_x <- sample_x_positions[j]
        
        # 获取该样本的所有簇数据
        sample_data <- current_samples %>%
            filter(Recipient == sample_id) %>%
            arrange(match(Cluster, present_levels))
        
        # 计算累积比例用于堆叠
        sample_data$cumsum_prop <- cumsum(sample_data$Proportion)
        sample_data$cumsum_prop_lag <- c(0, sample_data$cumsum_prop[-nrow(sample_data)])
        
        # 为该样本的每个簇创建矩形数据
        for(k in 1:nrow(sample_data)) {
            bar_data <- data.frame(
                xmin = sample_x - sample_bar_width/2,
                xmax = sample_x + sample_bar_width/2,
                ymin = current_lat + sample_data$cumsum_prop_lag[k] * bar_height_scale,
                ymax = current_lat + sample_data$cumsum_prop[k] * bar_height_scale,
                Cluster = sample_data$Cluster[k],
                Location = current_group,
                Sample = sample_id,
                n_samples = current_n,
                lon = current_lon,
                lat = current_lat
            )
            map_bars_data <- rbind(map_bars_data, bar_data)
        }
    }
}

# 绘制带有每个样本堆叠柱状图的地图
map_with_bars <- base_map +
    geom_rect(data = map_bars_data,
              aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, 
                  fill = Cluster),
              color = "white",
              linewidth = 0.02) +  # 减小边框线宽度，因为柱子很细
    scale_fill_manual(values = color_map_present) +
    # 添加样本数量标签
    geom_text(data = top_locations,
              aes(x = Longitude, y = Latitude - 1),
              label = paste0("n=", top_locations$n_samples),
              size = 2,
              hjust = 0.5,
              vjust = 1) +
    labs(title = "各地理位置的Chromopainter组成分布（每个样本单独显示）",
         subtitle = paste0("显示样本数≥10的前", nrow(top_locations), "个位置，每根细柱代表一个样本"))

# 保存地图
output_map_file <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/6-fineSTRUCTURE/output/东亚受到周边/chromopainter_individual_samples_map.pdf"
ggsave(output_map_file, 
       plot = map_with_bars, 
       width = 20, 
       height = 12, 
       device = "pdf")

cat("地图已保存:", output_map_file, "\n")
cat("共显示了", nrow(top_locations), "个主要位置的Chromopainter组成\n")
cat("总共绘制了", nrow(map_bars_data), "个矩形（每个样本的每个簇一个矩形）\n")

