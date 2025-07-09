#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import pandas as pd

def main():
    # ---------------------------
    # 1. 解析命令行参数
    # ---------------------------
    parser = argparse.ArgumentParser(description="将 plink 分群的等位基因频率转换为 treemix 输入格式。")
    parser.add_argument("plink_freq_file", help="plink 输出的 frq.strat 文件路径")
    args = parser.parse_args()

    # ---------------------------
    # 2. 读取输入文件
    # ---------------------------
    try:
        df = pd.read_csv(args.plink_freq_file, delim_whitespace=True)
    except Exception as e:
        sys.stderr.write(f"读取文件失败：{e}\n")
        sys.exit(1)

    # ---------------------------
    # 3. 检查重复的 SNP ID
    # ---------------------------
    # 检查是否存在重复的 SNP ID
    duplicated_snps = df["SNP"][df["SNP"].duplicated()].unique()
    if len(duplicated_snps) > 0:
        sys.stderr.write(f"警告: 发现 {len(duplicated_snps)} 个重复的 SNP ID\n")
        sys.stderr.write(f"将为重复的 SNP ID 添加唯一标识符...\n")
        
        # 创建一个计数字典，用于记录每个 SNP ID 出现的次数
        snp_counter = {}
        
        # 为重复的 SNP ID 添加后缀
        new_snp_ids = []
        for snp_id in df["SNP"]:
            if snp_id in snp_counter:
                snp_counter[snp_id] += 1
                new_snp_ids.append(f"{snp_id}_{snp_counter[snp_id]}")
            else:
                snp_counter[snp_id] = 0
                new_snp_ids.append(snp_id)
        
        # 更新 DataFrame 中的 SNP ID
        df["SNP"] = new_snp_ids
    
    # ---------------------------
    # 4. 获取所有 SNP 和群体 ID，保留原始出现顺序
    # ---------------------------
    snp_ids = df["SNP"].unique()
    pop_ids = df["CLST"].unique()
    nsnps = len(snp_ids)
    npops = len(pop_ids)

    # 在 stderr 输出基本信息
    sys.stderr.write(f"SNPs in file: {nsnps}\n")
    sys.stderr.write(f"Populations in input (n = {npops}):\n")
    for pop in pop_ids:
        sys.stderr.write(f"  {pop}\n")
    sys.stderr.write("\n")

    # ---------------------------
    # 5. 构造透视表
    #    - mac_table: index=SNP, columns=群体, 值=MAC（次要等位基因计数）
    #    - nchrobs_table: index=SNP, columns=群体, 值=NCHROBS（染色体观察数）
    # ---------------------------
    try:
        mac_table = df.pivot(index="SNP", columns="CLST", values="MAC")
        nchrobs_table = df.pivot(index="SNP", columns="CLST", values="NCHROBS")
    except ValueError as e:
        sys.stderr.write(f"透视表创建失败: {e}\n")
        sys.stderr.write("尝试改用 groupby 方法...\n")
        
        # 使用 groupby 作为备选方案
        mac_table = df.groupby(["SNP", "CLST"])["MAC"].first().unstack()
        nchrobs_table = df.groupby(["SNP", "CLST"])["NCHROBS"].first().unstack()

    # ---------------------------
    # 6. 过滤：只保留在所有群体都有记录的 SNP（无缺失）
    # ---------------------------
    # 检查是否有缺失值
    na_count_before = mac_table.isna().sum().sum()
    if na_count_before > 0:
        sys.stderr.write(f"警告：发现 {na_count_before} 个缺失值，将过滤掉含有缺失值的 SNP...\n")
    
    valid_snps = mac_table.dropna(axis=0, how="any").index
    
    # 如果没有有效的 SNP，输出错误并退出
    if len(valid_snps) == 0:
        sys.stderr.write("错误：过滤后没有有效的 SNP！请检查输入文件。\n")
        sys.exit(1)
    
    sys.stderr.write(f"过滤后保留了 {len(valid_snps)}/{len(mac_table)} 个 SNP\n")
    
    # 确保所有群体都存在于数据中
    missing_pops = set(pop_ids) - set(mac_table.columns)
    if missing_pops:
        sys.stderr.write(f"警告：以下群体在数据中缺失：{missing_pops}\n")
        # 只保留存在的群体
        pop_ids = [p for p in pop_ids if p in mac_table.columns]
    
    mac_table = mac_table.loc[valid_snps, pop_ids]
    nchrobs_table = nchrobs_table.loc[valid_snps, pop_ids]

    # ---------------------------
    # 6. 输出 treemix 格式
    #    - 首行：群体 ID，用空格分隔
    #    - 每行：每个 SNP 的 major,minor 计数对，用空格分组，不同群体之间用空格
    # ---------------------------
    # 写出群体 ID 行
    header = " ".join(pop_ids)
    sys.stdout.write(header + "\n")

    # 检查是否有有效的SNP可以输出
    if len(valid_snps) == 0:
        sys.stderr.write("错误：没有有效的SNP可以输出！\n")
        sys.exit(1)
        
    nsnps_out = 0
    # 对每个有效的 SNP 逐行处理并输出
    for snp in valid_snps:
        nsnps_out += 1
        if nsnps_out % 1000 == 0:  # 每1000个SNP更新一次进度
            sys.stderr.write(f"SNPs written: {nsnps_out}/{len(valid_snps)}\r")

        try:
            # 次要等位基因计数
            minor_counts = mac_table.loc[snp]
            # 主要等位基因计数 = NCHROBS - MAC
            major_counts = nchrobs_table.loc[snp] - minor_counts
            
            # 检查是否有负值（错误）
            if (major_counts < 0).any() or (minor_counts < 0).any():
                sys.stderr.write(f"警告：SNP {snp} 存在负数计数，将被跳过\n")
                continue
        except Exception as e:
            sys.stderr.write(f"处理 SNP {snp} 时出错: {e}，将被跳过\n")
            continue

        # 将每个群体的 "major,minor" 拼接成字符串
        # 处理可能的浮点数和缺失值
        allele_strs = []
        for pop in pop_ids:
            try:
                major = major_counts[pop]
                minor = minor_counts[pop]
                # 检查是否为 NaN
                if pd.isna(major) or pd.isna(minor):
                    allele_strs.append("?,?")  # TreeMix 中的缺失值表示
                else:
                    # 确保是整数
                    allele_strs.append(f"{int(major)},{int(minor)}")
            except Exception as e:
                sys.stderr.write(f"处理 SNP {snp} 群体 {pop} 时出错: {e}\n")
                allele_strs.append("?,?")
        line = " ".join(allele_strs)
        sys.stdout.write(line + "\n")

    # 换行以清晰结束 stderr 输出
    sys.stderr.write("\n")

if __name__ == "__main__":
    main()
