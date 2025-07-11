#!/bin/bash
set -e  # 如果任何命令返回非零状态，则退出脚本
set -o pipefail  # 如果管道中的任何命令失败，则整个管道失败

# 定义输出目录
OUTPUTDIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/OUTPUT/718中国信息样本/'
cd "${OUTPUTDIR}" || { echo "无法切换到目录 ${OUTPUTDIR}"; exit 1; }

# 记录开始时间
date > 718中国信息样本_HP.time
echo "分析开始于 $(date)" >> 718中国信息样本_HP.time

# 运行 fineSTRUCTURE 的初始命令
echo "运行 fineSTRUCTURE 初始命令..."
fs 718中国信息样本_HP.cp -hpc 1 \
    -idfile 718中国信息样本_HP.ids \
    -phasefiles 718中国信息样本_HP.phase \
    -recombfiles 718中国信息样本_HP.recombfile \
    -s3iters 200000 \
    -s4iters 50000 \
    -s1minsnps 1000 \
    -s1indfrac 0.1 \
    -go

# 定义一个函数来执行命令文件
run_commandfile() {
    local cmdfile=$1
    local stage=$2

    if [ -f "${cmdfile}" ]; then
        echo "运行 Stage${stage} 命令文件: ${cmdfile}"
        # 使用 GNU Parallel 执行命令文件中的命令，并等待完成
        parallel -j 16 < "${cmdfile}" || { echo "Stage${stage} 的并行执行失败。"; exit 1; }
        echo "Stage${stage} 完成，继续分析..."
        fs 718中国信息样本_HP.cp -go || { echo "Stage${stage} 后的 fs 命令失败。"; exit 1; }
    else
        echo "警告: ${cmdfile} 不存在，跳过 Stage${stage}。"
    fi
}

# 依次运行 Stage1 到 Stage4
for stage in {1..4}; do
    cmdfile="718中国信息样本_HP/commandfiles/commandfile${stage}.txt"
    run_commandfile "${cmdfile}" "${stage}"
done

# 记录结束时间
date >> 718中国信息样本_HP.time
echo "分析结束于 $(date)" >> 718中国信息样本_HP.time

echo "fineSTRUCTURE 分析完成。"
