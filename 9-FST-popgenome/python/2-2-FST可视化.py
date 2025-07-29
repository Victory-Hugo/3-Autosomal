#!/usr/bin/env python3
# 文件名：2-2-FST可视化.py
# 用法示例：
# /home/luolintao/miniconda3/bin/python 2-2-FST可视化.py \
#     --base_dir "${BASE_DIR}" \
#     --limitation 0.2 \
#     --gff_file "/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/script/R方法/NC_000915.gff"

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from aquarel import load_theme

def parse_args():
    parser = argparse.ArgumentParser(
        description="FST 可视化及注释脚本：\n"
                    "1) 根据 Fst 值绘制上下两个子图（上方和下方散点均为位图，坐标轴与注释保持矢量）。\n"
                    "2) 对 Fst > limitation 的位点进行 CDS 注释，并输出到注释文件。\n"
                    "3) 保存最终带注释的 PDF 图。\n"
                    "示例：\n"
                    "  python 2-2-FST可视化.py \\\n"
                    "    --base_dir \"/path/to/your/BASE_DIR\" \\\n"
                    "    --limitation 0.2 \\\n"
                    "    --gff_file \"/path/to/NC_000915.gff\""
    )
    parser.add_argument(
        "--base_dir", "-b",
        required=True,
        help="输入输出的基目录，脚本会从 {base_dir}/处理后FST.csv 读取数据，"
             "并将结果保存到 base_dir 下。"
    )
    parser.add_argument(
        "--limitation", "-l",
        type=float,
        required=True,
        help="Fst 分界值，范围在 0~1 之间，所有 Fst > limitation 的位点会做注释、"
             "上下子图的分割线也设为该值。"
    )
    parser.add_argument(
        "--gff_file", "-g",
        required=True,
        help="NC_000915.gff 文件路径，用于提取 CDS 区间及 product 注释信息。"
    )
    return parser.parse_args()

def read_cds_list(gff_path):
    """
    从 GFF 文件中提取所有 CDS 项目的信息，返回列表：
    [(start1, end1, product1), (start2, end2, product2), ...]
    """
    cds_list = []
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            # GFF 文件至少有 9 列，第三列为 "CDS" 的行才是我们要的
            if len(cols) >= 9 and cols[2] == "CDS":
                try:
                    start = int(cols[3])
                    end   = int(cols[4])
                except ValueError:
                    # 如果解析失败，跳过这一行
                    continue
                annotation_field = cols[8]
                # annotation_field 里各个子域用分号分隔，最后一段通常包含 product 信息
                fields = annotation_field.split(";")
                last_field = fields[-1].strip()
                if last_field.startswith("product="):
                    product_info = last_field[len("product="):]
                else:
                    product_info = last_field
                cds_list.append((start, end, product_info))
    return cds_list

def annotate_sites(data_df, cds_list, threshold, ax_top, txt_output_path):
    """
    对 data_df 中 Fst > threshold 的行进行注释：
    1) 如果某个位置 pos 在 cds_list 中的某个区间 [start, end] 内，则在 ax_top 上添加注释文字与箭头；
    2) 同时将注释结果写入 txt_output_path。
    """
    # 筛选 Fst > threshold 的行
    annotate_df = data_df[data_df["Fst"] > threshold].copy()
    results = []

    # 遍历每个需要注释的位点
    for idx, row in annotate_df.iterrows():
        pos = int(row["Location"])
        fst_value = float(row["Fst"])
        # 查找这个 pos 是否落在某个 CDS 区间内
        for (start, end, product_info) in cds_list:
            if start <= pos <= end:
                # 在图中做注释（矢量文字 + 箭头）
                ax_top.annotate(
                    product_info,
                    xy=(pos, fst_value),
                    xytext=(pos, min(fst_value + 0.03, 1.0)),
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.5),
                    fontsize=8
                )
                # 写到结果列表里：位置、Fst 值、产物说明、注释字段
                results.append(f"位置: {pos}, Fst: {fst_value:.4f}, CDS产物: {product_info}")
                break

    # 将结果写入文件
    with open(txt_output_path, "w", encoding="utf-8") as out_fh:
        for line in results:
            out_fh.write(line + "\n")

    return len(results)

def main():
    args = parse_args()
    BASE_DIR = args.base_dir.rstrip("/")  # 去掉末尾可能的斜杠
    limitation = args.limitation
    gff_file = args.gff_file

    # 检查输入路径与文件
    fst_csv_path = os.path.join(BASE_DIR, "处理后FST.csv")
    if not os.path.exists(fst_csv_path):
        raise FileNotFoundError(f"未找到 Fst 数据文件：{fst_csv_path}")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"未找到 GFF 文件：{gff_file}")

    # 读取 Fst 数据
    data = pd.read_csv(fst_csv_path, header=0, encoding="utf-8")
    # 必须包含 'Location' 和 'Fst' 两列
    if "Location" not in data.columns or "Fst" not in data.columns:
        raise KeyError(f"CSV 文件必须包含 'Location' 和 'Fst' 两列，请检查 {fst_csv_path}")

    # ---------- 第 1 步：绘制上下子图（上下散点都为位图，坐标轴、注释保持矢量） ----------
    # 加载主题
    theme = load_theme("arctic_light")
    theme.apply()
    plt.rcParams["font.family"]     = "Arial"
    plt.rcParams["pdf.fonttype"]    = 42
    plt.rcParams["ps.fonttype"]     = 42

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1,
        figsize=(12, 6),
        sharex=True,
        gridspec_kw={"height_ratios": [6, 1]}
    )
    plt.subplots_adjust(hspace=0)

    # 下方子图：所有散点栅格化
    scatter_bottom = ax_bottom.scatter(
        data["Location"],
        data["Fst"],
        c="#71A48D",
        alpha=1,
        s=7,
        rasterized=True
    )
    ax_bottom.set_ylim(0, limitation)
    ax_bottom.get_yaxis().set_visible(False)
    ax_bottom.set_ylabel("Fst", fontsize=10)
    ax_bottom.grid(False)

    # 上方子图：散点同样栅格化，其余元素保持矢量
    scatter_top = ax_top.scatter(
        data["Location"],
        data["Fst"],
        c="#AB5962",
        alpha=1,
        s=7,
        rasterized=True
    )
    ax_top.set_ylim(limitation, 1)
    ax_top.set_ylabel("Fst", fontsize=10)
    ax_top.grid(False)

    # 隐藏上方子图的 x 轴刻度和轴线
    ax_top.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax_top.spines["bottom"].set_visible(False)

    # 在 y=limitation 处画一条黑色矢量线
    ax_top.axhline(
        y=limitation,
        color="black",
        linewidth=4,
        alpha=1,
        linestyle="-"
    )

    # 公共设置
    ax_bottom.set_xlabel("Location (Position)", fontsize=10)
    ax_bottom.set_title("")
    ax_bottom.set_xlim(left=0, right=data["Location"].max() + 10000)
    ax_bottom.ticklabel_format(style="plain", axis="x")

    # ---------- 第 2 步：对 Fst > limitation 的位点进行 CDS 注释，并输出到 TXT ----------
    cds_list = read_cds_list(gff_file)
    txt_out = os.path.join(BASE_DIR, "注释.txt")
    n_hits = annotate_sites(data, cds_list, limitation, ax_top, txt_out)
    print(f"共找到 {n_hits} 个 Fst > {limitation} 的位点并完成注释。注释结果已保存到：{txt_out}")

    # ---------- 第 3 步：保存最终带注释的 PDF 图 ----------
    # 先应用主题的后处理
    theme.apply_transforms()
    pdf_out = os.path.join(BASE_DIR, "FST_注释.pdf")
    plt.savefig(pdf_out, dpi=500)
    plt.close(fig)
    print(f"带注释的 PDF 图已保存到：{pdf_out}")

if __name__ == "__main__":
    main()
