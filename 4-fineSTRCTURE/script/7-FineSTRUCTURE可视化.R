###########################################
#### FineSTRUCTURE 分析脚本 (优化示例) ####
###########################################

## 前置说明：
#todo 1. 需要提前安装 FineSTRUCTURE v2 软件，并将其可执行文件(fs)加入到系统路径或在代码中使用绝对路径。
#todo 2. R 环境中需安装 "gplots", "XML" 等依赖包。
#todo 3. 本脚本中的路径、文件名等请根据实际数据进行修改。

#################################
#### 1. 加载所需的 R 包和函数 ####
#################################
rm(list=ls())
# conda install -c conda-forge r-gplots
# install.packages("ape")
# install.packages("XML")
library("gplots") 
source("/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/python/Matrix_FinestructureLibrary.R")
## 上述 FinestructureLibrary.R 中会加载本分析需要的函数和依赖包


##########################################
#### 2. 定义 FineSTRUCTURE 输出文件路径 ####
##########################################
## 以下三个文件为 FineSTRUCTURE 主流程输出文件
chunkfile <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked.chunkcounts.out"
mcmcfile  <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked_mcmc.xml"
treefile  <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked_tree.xml"

########################################################
#### 3. 提取 FineSTRUCTURE 额外生成的结果文件并保存 ####
########################################################
## (1) population-by-population chunkcount 文件
mappopchunkfile <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked.mutationprobs.out"
## 使用 fs 命令行执行 -X -Y -e X2

command_X2 <- paste("fs fs -X -Y -e X2", chunkfile, treefile, mappopchunkfile)
system(command_X2)

## (2) pairwise coincidence 文件 (即 MCMC 中个体共聚类的平均比例)
meancoincidencefile <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/整体1369_HP_linked.meancoincidence.csv"
## 使用 fs 命令行执行 -X -Y -e meancoincidence
command_meancoin <- paste("fs fs -X -Y -e meancoincidence", chunkfile, mcmcfile, meancoincidencefile)
system(command_meancoin)

##########################################
#### 4. 读取数据：chunkcounts、MCMC 和树 ####
##########################################
## 4.1 读取 chunkcounts (coancestry matrix)
dataraw <- as.matrix(read.table(chunkfile, 
                                row.names = 1, 
                                header = TRUE, 
                                skip = 1))

## 4.2 读取并解析 MCMC 文件（XML 格式）
library("XML")
mcmcxml  <- xmlTreeParse(mcmcfile)
mcmcdata <- as.data.frame.myres(mcmcxml)

## 4.3 读取并提取树文件（XML 格式）
treexml <- xmlTreeParse(treefile)
ttree   <- extractTree(treexml)     # 将树转换为 ape 包的 phylo 格式
ttree$node.label <- NULL            # 移除内部节点标签，避免之后绘图时绘制多余信息
tdend   <- myapetodend(ttree, 
                       factor = 1)  # 将 phylo 转换为 dendrogram 格式

###########################################
#### 5. 提取并处理 FineSTRUCTURE 聚类结果 ####
###########################################
## 获取 MAP (Maximum A Posteriori) 状态下的群体划分
mapstate     <- extractValue(treexml, "Pop")  # finestructure 的聚类分配结果
mapstatelist <- popAsList(mapstate)           # 转换为各群体所含个体列表

## 生成群体名称
popnames      <- lapply(mapstatelist, NameSummary)      # 可逆格式（无信息损失）
popnamesplot  <- lapply(mapstatelist, NameMoreSummary)  # 更易读的可视化名称
names(popnames)     <- popnamesplot
names(popnamesplot) <- popnamesplot

## 生成 dendrogram 并修正其中的中点分支信息
popdend      <- makemydend(tdend, mapstatelist)
popdend      <- fixMidpointMembers(popdend)
popdendclear <- makemydend(tdend, mapstatelist, "NameMoreSummary")
popdendclear <- fixMidpointMembers(popdendclear)

##################################################
#### 6. 读取并处理 Pairwise coincidence 矩阵 ####
##################################################
## 按树的顺序对 coincidence 矩阵重新排序
fullorder       <- labels(tdend) 
mcmcmatrixraw   <- as.matrix(read.csv(meancoincidencefile, 
                                      row.names = 1))
mcmcmatrix      <- mcmcmatrixraw[fullorder, fullorder]

## 获取 MAP state 的二元矩阵形式 (groupingAsMatrix)
mapstatematrix  <- groupingAsMatrix(mapstatelist)[fullorder, fullorder]

######################################################
#### 7. 处理和保存 coancestry 矩阵（chunkcounts） ####
######################################################
datamatrix <- dataraw[fullorder, fullorder]

## 将 coancestry 矩阵写入 CSV 文件，以便后续在外部软件中替换 ID 
write.csv(datamatrix, 
          "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/datamatrix.csv")

## 注意：可在后续使用外部程序（如 Excel 的 VLOOKUP）将该矩阵的 INDXXX 替换为实际样本 ID

############################################################
#!### 8. 从替换过 ID 的矩阵文件中再次读取并进行后续操作 ####
#! 在终端输入 awk -v FS=',' '{gsub("\"", "", $1); print $1}' "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/datamatrix.csv" > /mnt/f/1_唐小琼项目/5_fineSTRUCTURE/script/整体1369_list.csv
############################################################
InputFileData  <- "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/datamatrix.csv"
datamatrix_R   <- read.table(file = InputFileData, 
                             header = TRUE, 
                             sep = ",", 
                             check.names = FALSE)

## 将第一列取回行名
names_pop            <- datamatrix_R[, 1]
datamatrix_R[, 1]    <- NULL
rownames(datamatrix_R) <- names_pop
datamatrix_R         <- as.matrix(datamatrix_R)

##################################################
#### 9. 定义新的颜色函数，并示例化获取颜色序列 ####
##################################################
MakeCustomColors <- function(numColors = 300, useDarkEnd = FALSE) {
  ## 定义一组基础配色
  baseColors <- c("#003936", "#31646C", "#4E9280", "#96B89B", 
                  "#DCDFD2", "#D49C87", "#B86265", "#50184E")
  
  ## 如果需要在末尾使用深色
  if (useDarkEnd) {
    baseColors[length(baseColors)] <- "#333333"  # 将最后一个颜色替换为深灰色
  }

  ## 使用 colorRampPalette 创建一个颜色渐变函数
  colorFunction <- colorRampPalette(baseColors)
  
  ## 生成并返回指定数量的颜色
  return(colorFunction(numColors))
}

## 获取两种颜色序列示例
some.colors    <- MakeCustomColors()
some.colorsEnd <- MakeCustomColors(numColors = 100, TRUE)  # 获取带深灰末尾的颜色序列

## 将矩阵中过大的值限定在 tmatmax 以内，以便更好地显示热图
tmatmax       <- 200
tmpmat        <- datamatrix_R
tmpmat[tmpmat > tmatmax] <- tmatmax



######################################################
#### 10. 读取样本群体信息，按需为每个群体定义颜色 ####
######################################################
## 若要对行和列进行着色，需要先根据群体/分群定义颜色
cols      <- rep('#828282', ncol(datamatrix_R))  # 初始默认颜色
cols_rows <- rep('#828282', nrow(datamatrix_R))
#! 请记得修改如下！
list_color <- read.table(file = "/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/script/整体1369_list.csv", 
                         header = TRUE, 
                         sep = ",")

## 根据群体名称(POP)分配颜色
Population <- rep('#CEF9FF', nrow(list_color))

Population[list_color$POP %in% "Han"]                     <- '#FDD262'
Population[list_color$POP %in% "hpAfrica1"]               <- '#E2D200'
Population[list_color$POP %in% "hpAfrica2"]               <- '#CB2314'
Population[list_color$POP %in% "hpAsia2"]                 <- '#79402E'
Population[list_color$POP %in% "hspEAsia"]                <- '#0B775E'
Population[list_color$POP %in% "hpEurope"]                <- '#DC863B'
Population[list_color$POP %in% "hpEuropeSahul"]           <- '#FDDDA0'
Population[list_color$POP %in% "hpNEAfrica"]              <- '#C27D38'
Population[list_color$POP %in% "hpNorthAsia"]             <- '#9B110E'
Population[list_color$POP %in% "hpSahul"]                 <- '#C52E19'
Population[list_color$POP %in% "hspAfrica1MiscAmerica"]   <- '#E79805'
Population[list_color$POP %in% "hspAfrica1NAmerica"]      <- '#D3DDDC'
Population[list_color$POP %in% "hspAfrica1SAfrica"]       <- '#E6A0C4'
Population[list_color$POP %in% "hspAfrica1WAfrica"]       <- '#ABDDDE'
Population[list_color$POP %in% "hspEurasia"]               <- '#ECCBAE'
Population[list_color$POP %in% "hspIndigenousAmerica"]    <- '#E1AF00'
Population[list_color$POP %in% "hspNEurope"]               <- '#F8AFA8'
Population[list_color$POP %in% "hspSWEurope"]              <- '#CCBA72'
Population[list_color$POP %in% "hspUral"]                  <- '#0E2C68'
Population[list_color$POP %in% "T"]                        <- '#000000'

## 根据 DAPC 分群信息分配颜色
colours_list_DAPC <- rep('#CEF9FF', nrow(list_color))
colours_list_DAPC[list_color$DAPC %in% "East_Asia"] <- '#008041'
colours_list_DAPC[list_color$DAPC %in% "Southeast_Asia"] <- '#0081cf'
colours_list_DAPC[list_color$DAPC %in% "South_Asia"] <- '#1951af'
colours_list_DAPC[list_color$DAPC %in% "Other"] <- '#000000'
colours_list_DAPC[list_color$DAPC %in% "Oceania"] <- '#aB3D4F'
colours_list_DAPC[list_color$DAPC %in% "Europe"] <- '#d11fae'
colours_list_DAPC[list_color$DAPC %in% "Central_Asia"] <- '#C1D2D2'
colours_list_DAPC[list_color$DAPC %in% "Papua_New_Guinea"] <- '#CB3D4F'
colours_list_DAPC[list_color$DAPC %in% "America"] <- '#947145'
colours_list_DAPC[list_color$DAPC %in% "South_Africa"] <- '#7F7F7F'
colours_list_DAPC[list_color$DAPC %in% "North_America"] <- '#DD8D29'
colours_list_DAPC[list_color$DAPC %in% "South_America"] <- '#FAD510'
colours_list_DAPC[list_color$DAPC %in% "West_Africa"] <- '#FAEFD1'
colours_list_DAPC[list_color$DAPC %in% "Central_Africa"]  <- '#9AA83A'
colours_list_DAPC[list_color$DAPC %in% "West_Asia"]  <- '#00668c'
colours_list_DAPC[list_color$DAPC %in% "Europe-Asia"]  <- '#81A88D'
colours_list_DAPC[list_color$DAPC %in% "North_Africa"]  <- '#000000'
# # 根据分群信息分配颜色
# colours_list_Classification <- rep('#7F7F7F', nrow(list_color))
# colours_list_Classification[list_color$Classification %in% "Other"] <- '#7F7F7F'
# colours_list_Classification[list_color$Classification %in% "Russia"] <- '#81A88D'
# colours_list_Classification[list_color$Classification %in% "Korea"] <- '#FAD510'
# colours_list_Classification[list_color$Classification %in% "Mongolia"] <- '#DD8D29'
# colours_list_Classification[list_color$Classification %in% "China"] <- '#de283b'
# colours_list_Classification[list_color$Classification %in% "Australia"] <- '#00668c'
# colours_list_Classification[list_color$Classification %in% "Indonesia"] <- '#1951af'
# colours_list_Classification[list_color$Classification %in% "India"] <- '#71c4ef'
# colours_list_Classification[list_color$Classification %in% "Thailand"] <- '#afffff'
# colours_list_Classification[list_color$Classification %in% "Myanmar"] <- '#25b1bf'
# colours_list_Classification[list_color$Classification %in% "Malaysia"] <- '#008041'
# colours_list_Classification[list_color$Classification %in% "Japan"] <- '#ffaea0'
# colours_list_Classification[list_color$Classification %in% "Nepal"] <- '#228B22'
# colours_list_Classification[list_color$Classification %in% "Vietnam"] <- '#FAEFD1'
# colours_list_Classification[list_color$Classification %in% "Bhutan"] <- '#3A9AB2'

# 根据分群信息分配颜色
colours_list_Classification <- rep('#7F7F7F', nrow(list_color))
colours_list_Classification[list_color$Classification %in% "South"]  <- '#1951af'
colours_list_Classification[list_color$Classification %in% "North"] <- '#de283b'


colours_list_Classification2 <- rep('#7F7F7F', nrow(list_color))
colours_list_Classification2[list_color$Classification2 %in% "Other"] <- '#000000'
colours_list_Classification2[list_color$Classification2 %in% "Inner_Mongolia"] <- '#CB3D4F'
colours_list_Classification2[list_color$Classification2 %in% "Heilongjiang"] <- '#d11fae'
colours_list_Classification2[list_color$Classification2 %in% "Ningxia"] <- '#FF3D3D'
colours_list_Classification2[list_color$Classification2 %in% "Beijing"] <- '#de283b'
colours_list_Classification2[list_color$Classification2 %in% "Sichuan"] <- '#00668c'
colours_list_Classification2[list_color$Classification2 %in% "Xizang"] <- '#1951af'
colours_list_Classification2[list_color$Classification2 %in% "Guangxi"] <- '#71c4ef'
colours_list_Classification2[list_color$Classification2 %in% "Yunnan"] <- '#afffff'
colours_list_Classification2[list_color$Classification2 %in% "Guizhou"] <- '#25b1bf'
colours_list_Classification2[list_color$Classification2 %in% "Fujian"] <- '#5B2C6F'
colours_list_Classification2[list_color$Classification2 %in% "Zhejiang"] <- '#3F51B5'
colours_list_Classification2[list_color$Classification2 %in% "Taiwan"] <- '#757de8'
colours_list_Classification2[list_color$Classification2 %in% "Shenzhen"] <- '#0E549B'
colours_list_Classification2[list_color$Classification2 %in% "Hongkong"] <- '#9AA83A'
colours_list_Classification2[list_color$Classification2 %in% "Shandong"] <- '#6E211F'
colours_list_Classification2[list_color$Classification2 %in% "Shanghai"] <- '#52A400'
colours_list_Classification2[list_color$Classification2 %in% "Hunan"] <- '#0E73CF'
colours_list_Classification2[list_color$Classification2 %in% "Chengdu"] <- '#00668c'
colours_list_Classification2[list_color$Classification2 %in% "Gansu"] <- '#DC5C00'
colours_list_Classification2[list_color$Classification2 %in% "Shaanxi"] <- '#FF4E38'
colours_list_Classification2[list_color$Classification2 %in% "Taipei"] <- '#757de8'


# 构建行、列注释矩阵，以在热图中实现侧边颜色标记
#! 自由选择，如果需要添加多余的标签，请取消注释
rlab=as.matrix(t(cbind(colours_list_Classification)))
clab=as.matrix((cbind(colours_list_Classification)))
rlab=as.matrix(t(cbind(colours_list_Classification,colours_list_Classification2)))
clab=as.matrix((cbind(colours_list_Classification,colours_list_Classification2)))
rlab=as.matrix(t(cbind(Population,colours_list_DAPC)))
clab=as.matrix((cbind(Population,colours_list_DAPC)))

# rlab=as.matrix(t(cbind(Population,colours_list_DAPC,colours_list_Classification)))
# clab=as.matrix((cbind(Population,colours_list_DAPC,colours_list_Classification)))

# rlab=as.matrix(t(cbind(Population,colours_list_DAPC,colours_list_Classification,colours_list_Classification2)))
# clab=as.matrix((cbind(Population,colours_list_DAPC,colours_list_Classification,colours_list_Classification2)))
###############################
#### 11. 绘制热图并保存为 PDF ####
###############################
# 使用自定义的 heatmap.3 函数绘制
dev.off()
pdf("/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/Fs_coancestry_heatmap.pdf",
    height = 100, 
    width  = 100)


source("/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/python/heatmap.3.R")

p <- heatmap.3(tmpmat,
               scale            = "none",
               Rowv             = tdend,
               Colv             = tdend,
               dendrogram       = "both",
               trace            = "none",
               key              = TRUE,
               keysize          = 0.5,
               cexRow           = 0.3,
               cexCol           = 0.3,
               col              = some.colorsEnd,  # 使用带深灰结尾的配色
               RowSideColorsSize= 1,
               ColSideColorsSize= 1,
               ColSideColors    = clab,
               RowSideColors    = rlab
)
dev.off()

#! 如果太大，尝试生成位图
dev.off()
# 使用 tiff() 函数代替 pdf() 函数，生成 TIFF 格式的图像
# 设置TIFF图形设备的尺寸；此处使用10英寸 x 10英寸，分辨率300dpi
tiff("/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/output/Fs_coancestry_heatmap.tif",
     width = 20,      
     height = 20,     
     units = "in",    
     res = 1000,       
     compression = "lzw")

# 可选：自定义更小的边距（必要时）
par(mar = c(4,4,2,2))

# 加载自定义函数
source("/mnt/f/1_唐小琼项目/5_fineSTRUCTURE/script/heatmap.3.R")

# 绘图
p <- heatmap.3(tmpmat,
               scale            = "none",
               Rowv             = tdend,
               Colv             = tdend,
               dendrogram       = "both",
               trace            = "none",
               key              = TRUE,
               keysize          = 0.5,
               cexRow           = 0.3,
               cexCol           = 0.3,
               col              = some.colorsEnd,
               RowSideColorsSize= 1,
               ColSideColorsSize= 1,
               ColSideColors    = clab,
               RowSideColors    = rlab
)
dev.off()
