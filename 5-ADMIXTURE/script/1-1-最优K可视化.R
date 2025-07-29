library(tidyverse)

# 读取CSV文件
K_CV <- read.csv('/mnt/d/幽门螺旋杆菌/Script/分析结果/ADMIXTURE/data/cv_results.csv')

# 方法1: 保持K为数值型，以便线条正确连接
K_CV |> 
  ggplot(aes(x = K, y = CV, group = 1)) +  # 添加group=1使所有点连接在一起
  geom_line() +
  geom_point(size = 3) +  # 增大点的大小以便更清晰
  labs(title = 'Cross-validation Results for K', 
       x = 'K', 
       y = 'Cross-Validation Error') +
  scale_x_continuous(breaks = K_CV$K) +  # 确保X轴只显示数据中的K值
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) 

