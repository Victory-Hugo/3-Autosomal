#*======== 从ADMIXTURE的LOG文件中获取最优K =======
#*==============================================
LOG_DIR='/mnt/f/1_唐小琼项目/13-ADMIXTURE/output/'
PREFIX='TXQ_filtered'
OUTPUT_FILE="${LOG_DIR}/cv_results.csv"

# 创建CSV文件并添加标题行
echo "K,CV" > $OUTPUT_FILE

# 进入日志目录
cd $LOG_DIR
# 获取所有以PREFIX开头的日志文件
LOG_FILES=$(ls ${PREFIX}*.log)
# 遍历每个日志文件
for LOG_FILE in $LOG_FILES; do
    # 提取K值和CV error值并格式化为CSV
    grep 'CV error' "$LOG_FILE" | sed -E 's/CV error \(K=([0-9]+)\): ([0-9.]+)/\1,\2/' >> $OUTPUT_FILE
done
# 将结果按K值排序
sort -n -t ',' -k1 $OUTPUT_FILE -o $OUTPUT_FILE
echo "CSV结果已保存到 $OUTPUT_FILE 文件中"
