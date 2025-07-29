#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from aquarel import load_theme

def main():
    parser = argparse.ArgumentParser(
        description="对 F_ST 值做经验分位数分析并绘制直方图，标注 99th/99.5th percentile 以及阈值线。"
    )
    parser.add_argument(
        "--input_csv", "-i",
        required=True,
        help="输入的处理后 FST CSV 文件路径（必须包含 'Fst' 列），例如 /path/处理后FST.csv"
    )
    parser.add_argument(
        "--output_png", "-o",
        required=True,
        help="输出 PNG 文件路径，例如 /path/FST_percentile_distribution.png"
    )

    args = parser.parse_args()
    input_csv = args.input_csv
    output_png = args.output_png

    # 加载主题
    theme = load_theme('arctic_light')
    theme.apply()
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # 1. 读取 CSV（假设含有 'Location' 和 'Fst' 两列）
    df = pd.read_csv(input_csv, header=0, encoding='utf-8')
    if "Fst" not in df.columns:
        raise KeyError("在处理后的 CSV 中找不到 'Fst' 列，请检查文件格式。")

    # 2. 提取 Fst 值向量
    fst_vals = df["Fst"].values

    # 3. 计算全基因组位点总数
    total_sites = len(fst_vals)

    # 4. 计算经验分位数
    p99  = np.percentile(fst_vals, 99)   # 99th percentile
    p995 = np.percentile(fst_vals, 99.5) # 99.5th percentile

    # 5. 设定“强分化”阈值（以 99th percentile 作为阈值）
    threshold = p99

    # 6. 打印说明文字（保留两位小数）
    print(
        f"全基因组共 {total_sites} 个位点，其中 99th percentile 的 F_ST 值约为 {p99:.2f}"
        f"（99.5th percentile 为 {p995:.2f}），本研究将 F_ST > {threshold:.2f} 视为强分化区段。"
    )

    # 7. 可视化分布：绘制直方图并标出 99th / 99.5th 以及阈值线
    plt.figure(figsize=(8, 5))
    plt.hist(fst_vals, bins=100, color='#71A48D', edgecolor='black', alpha=0.7)
    plt.xlabel("Fst value")
    plt.ylabel("Count")
    plt.title("Genome-wide Fst distribution (empirical quantile notation)")

    plt.axvline(x=p99,  color='red',   linestyle='--', linewidth=1.5, label=f"99th percentile = {p99:.2f}")
    plt.axvline(x=p995, color='orange',linestyle='--', linewidth=1.5, label=f"99.5th percentile = {p995:.2f}")
    plt.axvline(x=threshold, color='blue', linestyle='-',  linewidth=1.5, label=f"Threshold = {threshold:.2f}")

    plt.legend(loc="upper right", fontsize=10)
    plt.tight_layout()

    # 保存图片到指定位置
    plt.savefig(output_png, dpi=300)
    plt.close()

    print(f"直方图已保存到：{output_png}")

if __name__ == "__main__":
    main()
