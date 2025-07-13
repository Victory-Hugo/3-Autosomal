library(OptM)
# 如果文件前缀不是"global"，请修改为实际的前缀
# 分析rep1的数据
linear_rep1 = optM("/home/luolintao/Helicopter/Script/分析结果/treemix/output")
plot_optM(linear_rep1, main="Rep1")

# 如果有rep2，可以类似分析
