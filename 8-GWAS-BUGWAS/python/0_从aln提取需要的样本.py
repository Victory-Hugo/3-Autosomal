########
#!方案一#
########
# # 1. 使用集合存储目标菌株，提升匹配速度
# with open("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/script/strainsSubpop.txt", 'r') as f:
#     target_strains = {line.strip().split()[0] for line in f if line.strip()}

# # 2. 逐行读取 FASTA 文件，实时写入匹配记录（支持多行序列）
# with open("/mnt/d/幽门螺旋杆菌/Script/分析结果/merged_fasta/8622_noInDel.aln", 'r') as fin, \
#      open("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/alignmentSub.aln", 'w') as fout:
    
#     write_flag = False  # 标记是否写入当前记录
#     for line in fin:
#         if line.startswith('>'):
#             # 遇到新的 header，判断菌株名称是否在目标集合中
#             strain_name = line[1:].split()[0]
#             write_flag = strain_name in target_strains
#             if write_flag:
#                 fout.write(line)
#                 fout.flush()  # 实时写入 header
#         else:
#             # 非 header 行，如果上一个 header 匹配，则写入当前行（可能为多行序列）
#             if write_flag:
#                 fout.write(line)
#                 fout.flush()  # 实时写入序列行
########
#!方案二#
########
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="合并多个 FASTA 文件到一个输出文件")
    parser.add_argument("list_file", help="包含目标ID的列表文件")
    parser.add_argument("fasta_dir", help="存放 FASTA 文件的目录")
    parser.add_argument("output_file", help="输出合并后的文件路径")
    args = parser.parse_args()

    list_file = args.list_file
    fasta_dir = args.fasta_dir
    output_file = args.output_file

    # 读取目标ID列表
    try:
        with open(list_file, 'r') as f:
            target_ids = [line.strip() for line in f if line.strip()]
    except IOError as e:
        print(f"Error: 无法读取列表文件 {list_file}: {e}")
        sys.exit(1)

    # 打开输出文件（覆盖模式）
    try:
        with open(output_file, 'w') as fout:
            # 遍历所有目标ID
            for tid in target_ids:
                fasta_file = os.path.join(fasta_dir, tid + ".fasta")
                # 判断对应的 FASTA 文件是否存在
                if os.path.exists(fasta_file):
                    with open(fasta_file, 'r') as fin:
                        # 逐行读取，并写入输出文件，实现实时写入
                        for line in fin:
                            fout.write(line)
                    fout.flush()  # 每个文件处理完后刷新写入缓冲区
                else:
                    print(f"Warning: 文件 {fasta_file} 不存在")
    except IOError as e:
        print(f"Error: 无法写入输出文件 {output_file}: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
