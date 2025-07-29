#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="读取 FST.txt，将其转置后提取第一列（重命名为 Fst），并输出到 CSV 文件。"
    )
    parser.add_argument(
        "--base_dir",
        required=True,
        help="基础目录路径，例如 /mnt/d/幽门螺旋杆菌/Script/分析结果/FST/output/临床指标/AML-阿莫西林"
    )
    args = parser.parse_args()
    BASE_DIR = args.base_dir

    # 以下代码保持与原来一模一样
    df_结果 = pd.read_csv(f'{BASE_DIR}/FST.txt', sep='\t')
    df_结果 = df_结果.T
    df_结果 = df_结果.rename(columns={0: 'Fst'})
    df_结果 = df_结果.reset_index()
    df_结果 = df_结果.loc[:, ['Fst']]
    df_结果.to_csv(f'{BASE_DIR}/处理后FST.csv', header=True, index_label='Location')
    # 如果需要在脚本里对 data 继续操作，可在此处添加

if __name__ == "__main__":
    main()
