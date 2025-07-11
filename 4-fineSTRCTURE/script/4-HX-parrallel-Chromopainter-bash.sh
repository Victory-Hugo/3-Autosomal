#!/bin/bash

# 设置工作目录
WORK_DIR="/home/luolintao/S00-Github/3-Autosomal/4-fineSTRCTURE/download/Batch_chromopainter"

# 检查目录是否存在
if [ ! -d "$WORK_DIR" ]; then
    echo "错误：目录 $WORK_DIR 不存在"
    exit 1
fi

# 切换到工作目录
cd "$WORK_DIR"

# 循环处理所有.sh文件
for script in *.sh; do
    # 检查文件是否存在且是常规文件
    if [ -f "$script" ]; then
        # 获取脚本名称（不包含.sh后缀）
        base_name="${script%.sh}"
        
        # 输出正在处理的文件名
        echo "正在提交任务: $script"
        
        # 使用nohup运行脚本，并将输出重定向到对应的log文件
        nohup bash "$script" > "${base_name}.log" 2>&1 &
        
        # 等待1秒，避免同时提交太多任务
        sleep 1
    fi
done

echo "所有任务已提交完成"