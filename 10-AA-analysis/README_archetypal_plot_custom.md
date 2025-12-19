# 自定义Archetypal Plot工具使用说明

## 概述

这是一个自定义版本的`archetypal-plot`工具，主要功能与原版工具一致，但增加了按照FAM文件中群体信息进行着色的功能。

## 主要特点

1. **群体着色**: 根据FAM文件中的群体信息为每个样本分配不同颜色
2. **多种图形类型**: 支持条形图、单纯形图及其交互式版本
3. **多种颜色方案**: 提供tab20、Set3、viridis等颜色方案
4. **高质量输出**: 支持自定义DPI设置
5. **完全兼容**: 参数接口与原版archetypal-plot保持兼容

## 文件说明

- `archetypal_plot_custom.py`: 主要脚本文件
- `run_custom_plot.sh`: 使用示例脚本

## 依赖环境

```bash
# 必需的Python包
pip install pandas numpy matplotlib plotly
```

## 使用方法

### 基本语法

```bash
python3 archetypal_plot_custom.py -i INPUT_Q_FILE -f FAM_FILE -p PLOT_TYPE [OPTIONS]
```

### 参数说明

#### 必需参数
- `-i, --input_file`: Q文件路径（archetypal analysis输出的.Q文件）
- `-f, --fam_file`: FAM文件路径（包含群体信息的.fam文件）
- `-p, --plot_type`: 图形类型，选择：
  - `bar_simple`: 静态条形图
  - `bar_html`: 交互式HTML条形图
  - `plot_simplex`: 静态单纯形图
  - `plot_simplex_html`: 交互式HTML单纯形图
  - `all`: 生成所有类型图形

#### 可选参数
- `-so, --sorted`: 按最大原型组分排序（默认：False）
- `-dpi, --dpi_number`: 图像分辨率（默认：200）
- `-title, --data_title`: 图形标题（默认：空）
- `-c, --color_scheme`: 颜色方案，选择：
  - `tab20`: 20色调色板（默认）
  - `Set3`: 12色柔和调色板
  - `viridis`: 连续色彩映射
- `-o, --output_prefix`: 输出文件前缀（默认：aa_custom）

### 使用示例

#### 1. 生成基本条形图
```bash
python3 archetypal_plot_custom.py \
    -i ADMIXTURE_800_top5_AA.5.Q \
    -f ADMIXTURE_800_filtered_top5.fam \
    -p bar_simple \
    -title "Population Analysis" \
    -o population_analysis
```

#### 2. 生成高分辨率单纯形图
```bash
python3 archetypal_plot_custom.py \
    -i ADMIXTURE_800_top5_AA.5.Q \
    -f ADMIXTURE_800_filtered_top5.fam \
    -p plot_simplex \
    -title "Simplex Plot" \
    -c Set3 \
    -dpi 500 \
    -o high_res_simplex
```

#### 3. 生成所有类型图形
```bash
python3 archetypal_plot_custom.py \
    -i ADMIXTURE_800_top5_AA.5.Q \
    -f ADMIXTURE_800_filtered_top5.fam \
    -p all \
    -title "Complete Analysis" \
    -c viridis \
    -o complete_analysis
```

#### 4. 生成交互式图形
```bash
python3 archetypal_plot_custom.py \
    -i ADMIXTURE_800_top5_AA.5.Q \
    -f ADMIXTURE_800_filtered_top5.fam \
    -p bar_html \
    -title "Interactive Bar Plot"
```

## 输出文件

根据选择的图形类型，脚本会生成以下文件：

### 静态图形（.jpg）
- `{prefix}_bar_k{K}.jpg`: 条形图
- `{prefix}_simplex_k{K}.jpg`: 单纯形图

### 交互式图形（.html）
- `{prefix}_bar_k{K}_interactive.html`: 交互式条形图
- `{prefix}_simplex_k{K}_interactive.html`: 交互式单纯形图

其中`{prefix}`为指定的输出前缀，`{K}`为原型数量。

## 输入文件格式

### Q文件格式
- 空格分隔的数值文件
- 每行代表一个样本
- 每列代表一个原型的组分
- 行数应与FAM文件的样本数一致

示例：
```
0.0 0.028 0.0 0.0 0.972
0.0 0.0 0.0 0.0 1.000
...
```

### FAM文件格式
- 制表符分隔文件，包含6列：
  1. 群体标识
  2. 样本ID
  3. 父系ID（通常为0）
  4. 母系ID（通常为0）
  5. 性别（通常为0）
  6. 表型（通常为-9）

示例：
```
hspUral	HEL_AA7341AA_AS	0	0	0	-9
hpEurope	HEL_AA5958AA_AS	0	0	0	-9
...
```

## 功能特点详解

### 1. 群体着色系统
- 自动为每个群体分配唯一颜色
- 支持多种科学色彩方案
- 样本按群体排序显示，便于观察模式

### 2. 条形图功能
- 堆叠条形图显示每个样本的原型组分
- 群体区域用半透明色带标记
- 自动调整标签密度避免重叠
- 支持群体标签和样本标签

### 3. 单纯形图功能
- K维原型空间投影到2D平面
- 每个群体用不同颜色和符号表示
- 包含原型顶点标签
- 支持图例和注释

### 4. 交互式功能
- HTML输出支持缩放和悬停
- 悬停显示详细样本信息
- 群体信息集成到交互界面

## 与原版工具对比

| 功能 | 原版archetypal-plot | 自定义版本 |
|------|-------------------|-----------|
| 条形图 | ✓ | ✓ |
| 单纯形图 | ✓ | ✓ |
| 交互式图形 | ✓ | ✓ |
| 群体着色 | ✗ | ✓ |
| 多种颜色方案 | ✗ | ✓ |
| 群体标签 | ✗ | ✓ |
| 自动排序 | 有限 | 增强 |

## 故障排除

### 常见问题

1. **"文件不匹配"错误**
   - 检查Q文件行数是否与FAM文件样本数一致
   - 确认文件路径正确

2. **"单纯形图至少需要3个原型"警告**
   - 对于K<3的情况，只能生成条形图
   - 考虑增加原型数量

3. **图形显示不完整**
   - 调整DPI设置
   - 检查输出目录权限

4. **颜色显示异常**
   - 尝试不同的颜色方案
   - 检查群体数量是否超过颜色方案限制

### 调试技巧

1. 使用`-title`参数添加描述性标题
2. 先生成小样本数据测试
3. 检查Python依赖包版本
4. 确认工作目录有写入权限

## 更新日志

### v1.0 (2025-12-19)
- 初始版本发布
- 实现群体着色功能
- 支持多种图形类型
- 添加交互式HTML输出
- 兼容原版参数接口

## 作者与许可

基于原版archetypal-plot工具开发，增加了群体着色和可视化增强功能。

使用中如有问题，请检查输入文件格式和参数设置。