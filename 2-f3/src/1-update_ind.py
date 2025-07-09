#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(
        description="将分组 TXT 文件中的分组信息替换到 .ind 文件的第三列，并打印出未匹配的 ID"
    )
    parser.add_argument("-i", "--ind", required=True,
                        help="输入的 .ind 文件路径（空格分隔，三列：ID, 原标签, 备注）")
    parser.add_argument("-g", "--groups", required=True,
                        help="分组 TXT 文件路径（无表头，列 1 为 ID，列 2 为分组标签，用空格或制表符分隔）")
    parser.add_argument("-o", "--output", required=True,
                        help="输出 .ind 文件路径")
    args = parser.parse_args()

    # 读取 .ind 文件
    try:
        ind_df = pd.read_csv(
            args.ind,
            sep=r"\s+",
            header=None,
            names=["ID", "OriginalLabel", "Remark"],
            dtype=str
        )
    except Exception as e:
        sys.exit(f"读取 .ind 文件失败：{e}")

    # 读取分组文件
    try:
        groups_df = pd.read_csv(
            args.groups,
            sep=r"\s+",
            header=None,
            names=["ID", "Group"],
            dtype=str
        )
    except Exception as e:
        sys.exit(f"读取分组文件失败：{e}")

    # 合并
    merged = ind_df.merge(groups_df, on="ID", how="left")

    # 找出缺失 Group 的样本
    missing = merged.loc[merged["Group"].isnull(), "ID"].tolist()
    if missing:
        sys.stderr.write(
            f"警告：{len(missing)} 个样本在分组文件中未找到对应 Group，已保留原 Remark。\n"
            "未匹配的 ID 列表如下：\n"
        )
        for idx in missing:
            sys.stderr.write(f"{idx}\n")

    # 用 assign 一步到位，避免 chained assignment
    merged = merged.assign(
        Remark=merged["Group"].fillna(merged["Remark"])
    )

    # 输出
    try:
        merged[["ID", "OriginalLabel", "Remark"]].to_csv(
            args.output,
            sep=" ",
            header=False,
            index=False
        )
        print(f"成功生成：{args.output}")
    except Exception as e:
        sys.exit(f"写出文件失败：{e}")

if __name__ == "__main__":
    main()
