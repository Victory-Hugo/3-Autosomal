#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="生成用于Bugwas分析的基因型输入文件，将0/1/2/3转换为A/T/C/G"
    )
    parser.add_argument(
        "input_file", 
        help="输入文件路径（例如：/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/tmp2.txt）"
    )
    parser.add_argument(
        "--output_file", 
        help="输出文件路径。如果不提供，则根据 biallelic 参数自动生成默认路径",
        default=None
    )
    # 默认只保留双等位变异；使用 --no-biallelic 关闭此过滤（即保留所有变异）
    parser.add_argument(
        "--no-biallelic", 
        dest="biallelic", 
        action="store_false", 
        help="关闭双等位过滤，保留所有变异"
    )
    parser.set_defaults(biallelic=True)

    args = parser.parse_args()
    input_file = args.input_file
    biallelic = args.biallelic

    # 如果未指定输出文件，则根据 biallelic 选项自动生成默认输出文件路径
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = (
            "/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/geno_biallelic_SNP.txt"
            if biallelic else
            "/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/geno_multiallelic_SNP.txt"
        )

    header = []          # 存储表头信息：第一个元素为位点标识名称，第二个为样本名称列表
    positions = []       # 存储各个SNP位点的标识
    genotype_data = []   # 存储每个位点经过转换后的基因型数据

    try:
        with open(input_file, 'r') as f:
            for line_number, line in enumerate(f, start=1):
                tokens = line.strip().split()
                if not tokens:
                    continue  # 跳过空行

                # 第一行作为表头
                if line_number == 1:
                    if len(tokens) < 4:
                        raise ValueError("表头至少需要包含4列（位点标识、Ref、Alt以及样本名）")
                    header = [tokens[0], tokens[3:]]  # 第1列为位点标识，其余样本名从第4列开始
                    continue

                # 数据行至少需要4列
                if len(tokens) < 4:
                    print(f"Warning: 第 {line_number} 行数据不足，已跳过")
                    continue

                pos_id = tokens[0]
                ref_allele = tokens[1]
                alt_alleles = tokens[2].split(",")
                genotypes = tokens[3:]

                # 如果要求双等位但该行替代等位基因不为1个，则跳过此行
                if biallelic and len(alt_alleles) != 1:
                    continue

                # 将每个样本的基因型编码转换为对应的碱基
                try:
                    converted_genotypes = [
                        ref_allele if g == '0' else alt_alleles[int(g) - 1]
                        for g in genotypes
                    ]
                except (ValueError, IndexError) as e:
                    print(f"Error processing line {line_number}: {e}. 行内容: {tokens}")
                    continue

                positions.append(pos_id)
                genotype_data.append(converted_genotypes)
                # 输出进度信息
                print(f"Processed line {line_number}")
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}")
        sys.exit(1)

    try:
        with open(output_file, 'w') as f_out:
            # 写入表头
            f_out.write(header[0] + '\t' + '\t'.join(header[1]) + '\n')
            # 写入每一行数据
            for pos, geno in zip(positions, genotype_data):
                f_out.write(pos + '\t' + '\t'.join(geno) + '\n')
    except Exception as e:
        print(f"Error writing output file {output_file}: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
