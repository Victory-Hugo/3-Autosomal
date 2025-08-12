#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import ast
import os
import sys
import math
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 强制使用非交互式后端，避免 Qt 依赖
import matplotlib.pyplot as plt


# seaborn / aquarel 可选（若缺失，脚本仍能运行）
try:
    import seaborn as sns
    _HAS_SEABORN = True
except Exception:
    _HAS_SEABORN = False

def _try_apply_aquarel_theme():
    try:
        from aquarel import load_theme
        theme = load_theme('boxy_light')
        theme.apply()
        return theme
    except Exception:
        return None

def parse_list(s: str):
    """
    解析列表参数：支持
    1) JSON/py-literal 形式：'["a","b"]' 或 "['a','b']"
    2) 逗号分隔：'a,b,c'
    """
    if s is None:
        return []
    s = s.strip()
    if not s:
        return []
    # 尝试用 literal_eval
    if (s.startswith('[') and s.endswith(']')) or (s.startswith('(') and s.endswith(')')):
        try:
            val = ast.literal_eval(s)
            if isinstance(val, (list, tuple)):
                return [str(x).strip() for x in val]
        except Exception:
            pass
    # 逗号分隔兜底
    return [x.strip() for x in s.split(',') if x.strip()]

def read_colors_map(path_or_str: str):
    """
    按需求从“无表头两列 CSV”读取颜色映射：
    第1列：供体名
    第2列：HEX 颜色
    若传入的不是路径，尝试按 literal_eval 解析成 dict 以增强兼容性（非必需）。
    """
    if path_or_str and os.path.exists(path_or_str):
        dfc = pd.read_csv(path_or_str, header=None)
        if dfc.shape[1] < 2:
            raise ValueError(f"colors CSV 至少需要两列，无表头。当前列数={dfc.shape[1]}")
        dfc = dfc.iloc[:, :2]
        dfc.columns = ['donor', 'color']
        # 去重并保留首次出现
        dfc = dfc.dropna().drop_duplicates(subset=['donor'], keep='first')
        return dict(zip(dfc['donor'].astype(str), dfc['color'].astype(str)))
    # 兜底：如果不是路径，尝试解析成字典（可选）
    try:
        maybe = ast.literal_eval(path_or_str)
        if isinstance(maybe, dict):
            return {str(k): str(v) for k, v in maybe.items()}
    except Exception:
        pass
    raise FileNotFoundError("未找到 colors CSV 文件，且无法从字面量解析为字典。")

def main():
    parser = argparse.ArgumentParser(
        description="Make donor-wise horizontal boxplots with global normalization."
    )
    parser.add_argument("--input", required=True, help="输入表的绝对路径（如 1.txt）")
    parser.add_argument("--sep", default="\\t", help="输入分隔符，默认 '\\t'")
    parser.add_argument("--value-vars", required=True,
                        help="供体列名列表；支持 JSON/py-literal 或逗号分隔")
    parser.add_argument("--colors-csv", required=True,
                        help="颜色映射：无表头两列 CSV 路径（donor, hex）。")
    parser.add_argument("--valid-recipients", required=True,
                        help="有效受体列表；支持 JSON/py-literal 或逗号分隔")
    parser.add_argument("--n-cols", type=int, default=3, help="每行子图数，默认 3")
    parser.add_argument("--out", default="Chromopainter_2.pdf", help="输出 PDF 路径")
    args = parser.parse_args()

    # 解析参数
    value_vars = parse_list(args.value_vars)
    valid_recipients = parse_list(args.valid_recipients)
    if not value_vars:
        raise SystemExit("value_vars 为空。请提供至少一个供体列名。")
    if not valid_recipients:
        raise SystemExit("valid_recipients 为空。请至少提供一个受体名称。")

    colors = read_colors_map(args.colors_csv)

    # 读入数据
    df = pd.read_csv(args.input, sep=args.sep, dtype=str)
    # 转数值列：只对 value_vars 做数值转换，保留非数值为 NaN
    for c in value_vars:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # 校验必须列
    required_cols = {"Recipient", "属性"}
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"输入表缺少必需列：{missing}。需要包含列：{sorted(required_cols)}")

    # 校验 value_vars 存在性
    missing_vals = [v for v in value_vars if v not in df.columns]
    if missing_vals:
        print(f"警告：以下 value_vars 列在输入数据中不存在，将被忽略：{missing_vals}", file=sys.stderr)
        value_vars = [v for v in value_vars if v in df.columns]
        if not value_vars:
            raise SystemExit("在数据中找不到任何有效的 value_vars 列。")

    # melt
    df_long = df.melt(
        id_vars=["Recipient", "属性"],
        value_vars=value_vars,
        var_name="供体",
        value_name="贡献值"
    )

    # 全局归一化（仅在有效受体子集上计算 min/max）
    global_data = df_long[df_long["属性"].isin(valid_recipients)].copy()

    # 若子集为空，无法归一化——提前终止
    if global_data.empty:
        raise SystemExit("有效受体子集为空：请检查 valid_recipients 与数据中的“属性”列是否匹配。")

    global_min = pd.to_numeric(global_data["贡献值"], errors="coerce").min()
    global_max = pd.to_numeric(global_data["贡献值"], errors="coerce").max()

    # 处理 max == min 的退化情况，避免除零
    denom = (global_max - global_min)
    if pd.isna(global_min) or pd.isna(global_max):
        raise SystemExit("贡献值存在全 NaN，无法归一化。请检查输入数据。")

    if denom == 0:
        # 全部相等时，归一化结果统一为 0.5（或 0 任一固定值都可）
        df_long["global_normalized"] = 0.5
    else:
        df_long["global_normalized"] = (pd.to_numeric(df_long["贡献值"], errors="coerce") - global_min) / denom

    # 拟绘制的供体：既在数据中出现、又在颜色表中定义
    donors_in_data = set(df_long["供体"].dropna().unique().tolist())
    donor_list = [d for d in colors.keys() if d in donors_in_data]
    if not donor_list:
        raise SystemExit("没有可绘制的供体：请检查 colors 中的 donor 名称与数据中的列是否一致。")

    # 主题与全局字体设置
    theme = _try_apply_aquarel_theme()
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] = False

    # 计算子图行列
    n_donors = len(donor_list)
    n_cols = max(1, int(args.n_cols))
    n_rows = math.ceil(n_donors / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4 * n_cols, 2.5 * n_rows),
        squeeze=False
    )

    # 作图
    for i, donor in enumerate(donor_list):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row][col]

        donor_data = df_long[df_long["供体"] == donor].copy()
        donor_data = donor_data[donor_data["属性"].isin(valid_recipients)]

        if donor_data.empty:
            ax.set_visible(False)
            continue

        # 若无 seaborn，则用 matplotlib 的 boxplot 兜底
        if _HAS_SEABORN:
            import seaborn as sns
            sns.boxplot(
                data=donor_data,
                x="global_normalized",
                y="属性",
                color=colors.get(donor, None),
                width=0.5,
                showfliers=True,
                flierprops=dict(
                    marker='o',
                    markersize=3,
                    markerfacecolor=colors.get(donor, None),
                    markeredgecolor=colors.get(donor, None),
                    alpha=0.5,
                ),
                boxprops={'edgecolor': colors.get(donor, None)},
                ax=ax
            )
        else:
            # matplotlib 兜底：按“属性”分类分别绘制
            ycats = list(donor_data["属性"].dropna().unique())
            y_index = {v: i+1 for i, v in enumerate(ycats)}
            positions = [y_index[v] for v in donor_data["属性"]]
            grouped = donor_data.groupby("属性")["global_normalized"].apply(list)
            bp = ax.boxplot(
                grouped.values,
                vert=False,
                patch_artist=True,
                showfliers=True
            )
            # 着色
            for patch in bp['boxes']:
                patch.set_facecolor(colors.get(donor, "#BBBBBB"))
                patch.set_edgecolor(colors.get(donor, "#666666"))
            for fl in bp.get('fliers', []):
                fl.set_marker('o')
                fl.set_markersize(3)
                fl.set_alpha(0.5)
                fl.set_markerfacecolor(colors.get(donor, "#888888"))
                fl.set_markeredgecolor(colors.get(donor, "#888888"))
            ax.set_yticks(range(1, len(ycats)+1))
            ax.set_yticklabels(ycats)

        ax.set_title(f"Donor: {donor}")
        ax.grid(False)
        ax.set_xlabel("Global normalized contribution")
        ax.set_ylabel("Recipient")

    # 删除多余子图
    total_subplots = n_rows * n_cols
    if n_donors < total_subplots:
        for j in range(n_donors, total_subplots):
            row = j // n_cols
            col = j % n_cols
            fig.delaxes(axes[row][col])

    if theme is not None:
        try:
            theme.apply_transforms()
        except Exception:
            pass

    plt.tight_layout()
    out_path = args.out
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path)
    print(f"Saved figure to: {out_path}")
    # 如需显示请手动 -- 但命令行脚本默认不 show
    # plt.show()

if __name__ == "__main__":
    main()
