
# pip install pandas matplotlib seaborn aquarel
# aquarel 可选；缺失则自动跳过主题
#? 输入文件有2个，一个是数据，一个是颜色
#? 数据可以是csv或者txt，修改--sep即可
#? 数据文件应当符合如下格式：
#* Recipient,供体1,供体2,供体3,供体4,供体5,供体6,供体7,供体8,供体9,供体10,供体11,供体12……,属性
#* ID1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2……,hpAfrica1
#* ID2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2……,hpAfrica1
#* ID3,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2……,hpAfrica2
#* ……
#? 颜色文件应当符合如下格式：
#* hpAfrica1,#123029
#*# hpAfrica2,#000000
#! 数据文件中的供体名在--value-vars指定
#! 希望绘制的受体在--valid-recipients指定
#! 颜色在--colors-csv指定
PYTHON3="python3"
BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/"


"$PYTHON3" ${BASE_DIR}/python/1-CP-箱线图.py \
  --input ${BASE_DIR}/conf/1-CP-箱线图.csv \
  --sep "," \
  --value-vars "['hpAfrica1','hpAfrica2','hpAsia2','hpEurope','hpNEAfrica','hpNorthAsia','hpSahul','hpSiberia','hspHAEAsiaN','hspHAEAsiaS','hspLASEAsia','hspLAEAJapan']" \
  --colors-csv ${BASE_DIR}/conf/1-CP-箱线图color.csv \
  --valid-recipients "['hspHAEAsiaN','hspHAEAsiaS','hspLAEAChinaN','hspLAEAChinaS','hspLAEAChinaNS']" \
  --n-cols 3 \
  --out /mnt/c/Users/Administrator/Desktop/Chromopainter_2.pdf
