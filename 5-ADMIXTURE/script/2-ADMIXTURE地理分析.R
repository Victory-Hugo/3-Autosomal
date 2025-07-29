# 安装并加载必要包（如未安装请先 install.packages()）
library(sf)
library(tidyverse)
library(scatterpie)
library(terra)
library(ggspatial)
# 1. 读入世界地图 shapefile
world_shp <- read_sf(
  "/mnt/f/OneDrive/文档（科研）/地图/全球行政地图(精细版)/全球国家边界/世界地图国家.shp"
)

world_raster <- terra::rast(
  "/mnt/f/OneDrive/文档（科研）/地图/全球自然地图包（位图）/HYP_HR_SR_W.tif"
) %>% 
  aggregate(fact = 5) 

# 全球范围
xlim_range <- c(60, 150)    # 经度范围：西到东
ylim_range <- c(-10, 60)     # 纬度范围：南到北

# 创建范围多边形
extent_bbox <- terra::ext(xlim_range[1], xlim_range[2], ylim_range[1], ylim_range[2])
# 裁剪栅格到指定范围
world_raster_cropped <- terra::crop(world_raster, extent_bbox)

# 2. 读入 ADMIXTURE Q 矩阵（含经纬度）
df_Q <- read.csv(
  "/mnt/c/Users/Administrator/Desktop/ADMIXTURE.csv",
  stringsAsFactors = FALSE
) |> as_tibble()

df_Q

# 2.1 自动检测K组分的数量和名称
# 找出所有以"K"开头后跟数字的列名
k_columns <- grep("^K\\d+$", names(df_Q), value = TRUE)
k_count <- length(k_columns)

cat("检测到", k_count, "个K组分:", paste(k_columns, collapse = ", "), "\n")

# 2.2 生成足够的颜色
# 使用RColorBrewer调色板和自定义颜色
base_colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#666666", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
  "#B15928", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5"
)

# 确保有足够的颜色
if (k_count > length(base_colors)) {
  # 如果K数量超过预定义颜色，使用rainbow生成更多颜色
  admix_colors <- rainbow(k_count)
  names(admix_colors) <- k_columns
} else {
  # 使用预定义颜色
  admix_colors <- setNames(base_colors[1:k_count], k_columns)
}

# 筛选出经纬度在亚洲的样本
df_Q_asia <- df_Q |>
  filter(!is.na(Latitude) & !is.na(Longitude)) |>
  filter(Latitude >= -10 & Latitude <= 60, Longitude >= 60 & Longitude <= 150)

df_Q_asia

# 3. 按位点聚合计算平均值（每个地理位置的平均ADMIXTURE比例）
# 使用across()函数动态处理所有K列
df_Q_grouped <- df_Q_asia |>
  group_by(Latitude, Longitude) |>
  summarise(
    across(all_of(k_columns), ~ mean(.x, na.rm = TRUE)),
    n_samples = n(),
    .groups = 'drop'
  ) |>
  # 过滤掉样本数量太少的位点
  filter(n_samples >= 3)

# 4. 动态设置ADMIXTURE组分的颜色（已在上面生成）
cat("使用的颜色配置:\n")
for(i in 1:length(admix_colors)) {
  cat(names(admix_colors)[i], ":", admix_colors[i], "\n")
}

# 5. 创建基础地图
base_map <- ggplot() +
    # 使用 ggspatial 的 layer_spatial() 添加栅格图层
  layer_spatial(world_raster_cropped, alpha = 1) +
  coord_sf(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(title = paste0("ADMIXTURE Analysis Results (K=", k_count, ") - Geographic Distribution"))

# 先定义最小/最大半径
min_r <- 1.0
max_r <- 3.0

# 6. 使用scatterpie绘制饼图在地图上
map_with_pies <- base_map +
  geom_scatterpie(
    data = df_Q_grouped,
    aes(
      x = Longitude,
      y = Latitude,
      # 先按原来 sqrt(n_samples)*0.28 算半径，再用 pmax/pmin 钳制
      r = pmin(
            pmax(sqrt(n_samples) * 0.28, min_r),
            max_r
          )
    ),
    cols  = k_columns,  # 使用动态检测的K列名
    color     = "#A0C7E3",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  scale_fill_manual(values = admix_colors) +
  guides(fill = guide_legend(title = "ADMIXTURE Components"))

# 显示地图
print(map_with_pies)

# 7. 保存地图
ggsave(paste0("/mnt/c/Users/Administrator/Desktop/ADMIXTURE_K", k_count, "_geographic_map.pdf"), 
       plot = map_with_pies, width = 12, height = 8, dpi = 300)

# 8. 创建柱状图版本的地理分布图
# 首先准备数据 - 将宽格式转换为长格式
df_Q_long <- df_Q_grouped |>
  pivot_longer(
    cols = all_of(k_columns),
    names_to = "Component", 
    values_to = "Proportion"
  ) |>
  mutate(Component = factor(Component, levels = k_columns))

# 创建自定义函数来绘制迷你柱状图
create_mini_barplot <- function(data, lon, lat, width = 1.5, height = 2) {
  # 为每个位置创建柱状图的坐标
  n_components <- length(unique(data$Component))
  bar_width <- width / n_components
  
  # 计算每个柱子的位置
  x_positions <- seq(lon - width/2, lon + width/2, length.out = n_components + 1)
  x_positions <- x_positions[-length(x_positions)] + bar_width/2
  
  # 创建柱状图数据框
  bar_data <- data.frame(
    xmin = x_positions - bar_width/2,
    xmax = x_positions + bar_width/2,
    ymin = lat,
    ymax = lat + data$Proportion * height,
    Component = data$Component,
    Proportion = data$Proportion
  )
  
  return(bar_data)
}

# 为每个位置创建柱状图数据
all_bar_data <- data.frame()
for(i in 1:nrow(df_Q_grouped)) {
  location_data <- df_Q_long |> 
    filter(Latitude == df_Q_grouped$Latitude[i], 
           Longitude == df_Q_grouped$Longitude[i])
  
  bar_data <- create_mini_barplot(
    location_data, 
    df_Q_grouped$Longitude[i], 
    df_Q_grouped$Latitude[i],
    width = 2,  # 柱状图宽度
    height = 3  # 柱状图高度
  )
  
  all_bar_data <- rbind(all_bar_data, bar_data)
}

# 9. 创建带柱状图的地图
map_with_bars <- base_map +
  geom_rect(
    data = all_bar_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Component),
    color = "black", 
    size = 0.1,
    alpha = 0.8
  ) +
  scale_fill_manual(values = admix_colors) +
  guides(fill = guide_legend(title = "ADMIXTURE Components")) +
  labs(title = paste0("ADMIXTURE Analysis Results (K=", k_count, ") - Geographic Distribution (Bar Charts)"))

# 显示柱状图地图
print(map_with_bars)

# 保存柱状图地图
ggsave(paste0("/mnt/c/Users/Administrator/Desktop/ADMIXTURE_K", k_count, "_geographic_barmap.pdf"), 
       plot = map_with_bars, width = 12, height = 8, dpi = 300)

