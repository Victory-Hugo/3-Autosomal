# popstats
种群遗传学汇总统计

POPSTATS 是一个用于计算 f 统计量、种群对称性检验和其他种群遗传学量的 Python 程序。它还允许估计 h4 统计量，该统计量首次在 Skoglund 等人 (2015, Nature) 中使用。这是一个初步发布，后续将提供更多文档和一个精炼的版本。

如果您使用 POPSTATS，请引用：

P Skoglund, S Mallick, MC Bortolini, N Chennagiri, T Hünemeier, ML Petzl-Erler, FM Salzano, N Patterson, D Reich (2015) Genetic evidence for two founding populations of the Americas, Nature, 525:104-108

### 基本用法

POPSTATS 使用 PLINK 转置文件，我们提供脚本 vcf2tped.py 将 VCF 文件转换为此格式。

下面的示例用于计算 D 统计量 (Green 等人 2010, Science; Patterson 等人 2012, Genetics)，以检测非非洲人群中尼安德特人混入的证据：

```python
# 我们在 myfile.tped 和 myfile.tfam 中拥有来自 1000 基因组计划和阿尔泰尼尼安德特人基因组 (Prufer 等人 2014, Nature) 的数据
python popstats.py --file myfile --pops chimpanzee,Neandertal,Yoruba,Japanese --informative

# 结果：
chimpanzee 	Neandertal 	Yoruba Japanese 	0.0566273725392 	0.00376607791239 	15.036165968 	1002084 	530 	2 	2 	208 	214
```

输出结果中的列代表：
```
1. 群体 A
2. 群体 B
3. 群体 X
4. 群体 Y
5. D(A, B; X, Y)
6. 加权块自助法标准误 (SE)
7. Z 分数 (D/SE)
8. 用于分析的位点数
9. 自助法块数
10. 群体 A 的染色体数
11. 群体 B 的染色体数
12. 群体 X 的染色体数
13. 群体 Y 的染色体数
```

### 统计量

以下所有统计量的标准误均使用块自助法计算（默认块大小为 5mb，并按位点数加权，详见下文选项）。这些统计量使用 **--pops POP1,POP2,POP3,POP4** 定义的四个群体，分别对应等位基因频率 p1, p2, p3, p4。

**--pi**  
通过从 POP1 中随机抽取两个个体的一个等位基因，计算不匹配概率来估计杂合度

**--FST**  
使用 Hudson 估计量 (Bhatia 等人) 估计 POP1 与 POP2 之间的 F_ST。使用 **--FSTWC** 可使用 Weir 和 Cockerham 估计量。

**--f2**  
估计 f2 统计量（平均平方等位基因频率差），用于 POP1 与 POP2 (Reich 等人 2009, Nature)

**--f3**  
估计 f3 统计量 (Reich 等人 2009)。默认使用 Patterson 等人 2012, Genetics 中描述的杂合度校正。若使用简单的 f3 统计量 f3=(p3-p1)(p3-p2)，请加上 **--f3vanilla**。

**--f4**  
估计 f4 统计量 (p1-p2)*(p3-p4) (Reich 等人 2009)。

**--D**  
估计 D 统计量 (Green 等人 2010, Science; Patterson 等人 2012, Genetics)。如果未指定其他选项，默认使用此统计量。

**--symmetry**  
通过在 POP1 和 POP2 中分别随机抽取一个基因拷贝，并在它们等位基因不同的位点上条件抽样，计算一个全基因组平均对称性统计量，以检测 POP1 在位点上衍生等位基因过剩 (正值) 或不足 (负值)。祖先等位基因通过 **--outgroup** 指定的外群确定。该统计量类似于 Do 等人 2014, Nature Genetics 使用的方法。

**--LD** [距离]  
估计指定距离处的 h4 统计量 (Skoglund 等人 2015, Nature)。必须同时使用 **--LDwindow** [距离] 和 **--withinfreq**。

通过指定 **--testpop**，可以计算 **f4 比率**，该比率将估计第五群体 POPT 的等位基因频率 pt。计算的统计量为两个 f4 统计量的比值 ((p1-p2)*(pt-p4))/((p1-p2)*(p3-p4))，在特定系统发育假设下可作为混合群体中祖先比例的无偏估计量。详情见 Patterson 等人 2012, Genetics。

**--FAB**  
估计在随机抽取的两个 POP1 拷贝等位基因不同的位点上，POP2 携带衍生等位基因的概率。祖先等位基因通过 **--ancestor** 指定的外群确定。该统计量可在对 POP1 遗传漂变做出假设后用于估计群体分化时间。详见 Green 等人 2010, Science 中对尼安德特人分化时间的估计。

### 可选参数

**--informative**  
在块自助法加权时，仅使用在 POP1+POP2 与 POP3+POP4 中均多态的 SNP，可在某些情况下略微减小标准误

**--morgan**  
使用遗传距离（默认 5 cM）而非物理距离来定义块大小

**--noweighting**  
执行非加权块自助法，不考虑各块的位点数

**--chromblocks**  
将整条染色体作为块执行自助法

**--nojackknife**  
不估计标准误

**--not23**  
使用输入文件中提供的所有染色体。默认仅使用 1-22 号染色体。该选项完全支持非人类生物。

**--haploidize**  
从每个群体随机抽取一个单倍型基因型进行统计。该选项与 **--D** 一起使用时，将计算 Green 等人 2010, Science 中经典的 ABBA-BABA 统计量。
