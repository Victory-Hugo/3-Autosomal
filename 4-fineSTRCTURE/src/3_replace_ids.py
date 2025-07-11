#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import argparse

def load_id_anchor_map(csv_file):
    """
    从 CSV 文件中加载 ID→Anchor 的映射关系。
    假设第一行为表头，ID 在第一列，Anchor 在第二列，文件以逗号分隔。
    """
    id_anchor_dict = {}
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None)  # 跳过表头
        for row in reader:
            if not row:
                continue
            the_id = row[0].strip()
            anchor = row[1].strip()
            id_anchor_dict[the_id] = anchor
    return id_anchor_dict

def process_ids_file(ids_file, id_anchor_dict):
    """
    逐行读取 .ids 文件，根据第一列查找对应的 Anchor 值，
    将第二列替换为对应的 Anchor（若未找到则保持原值）。
    """
    temp_file = ids_file + ".tmp"
    with open(ids_file, 'r', encoding='utf-8') as fin, \
         open(temp_file, 'w', encoding='utf-8') as fout:
        
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # 假定格式：第一列为 ID，第二列为默认值（通常与ID相同），第三列为其他信息
            if len(parts) < 3:
                # 格式异常时直接跳过该行
                continue
            
            the_id = parts[0]
            original_val = parts[1]
            col3 = parts[2]
            
            if the_id in id_anchor_dict:
                fout.write(f"{the_id} {id_anchor_dict[the_id]} {col3}\n")
            else:
                fout.write(f"{the_id} {original_val} {col3}\n")
    
    # 用处理后的临时文件覆盖原文件
    os.replace(temp_file, ids_file)

def main():
    parser = argparse.ArgumentParser(
        description="根据传入 CSV 映射文件处理 .ids 文件，替换第二列为对应 Anchor 值"
    )
    parser.add_argument("--inf_csv", required=True,
                        help="包含 ID 和 Anchor 映射关系的 CSV 文件路径")
    parser.add_argument("ids_files", nargs="+",
                        help="需要处理的一个或多个 .ids 文件路径")
    args = parser.parse_args()

    id_anchor_dict = load_id_anchor_map(args.inf_csv)
    
    for ids_file in args.ids_files:
        if os.path.isfile(ids_file):
            process_ids_file(ids_file, id_anchor_dict)
        else:
            print(f"警告: {ids_file} 不是一个有效的文件")

if __name__ == "__main__":
    main()
