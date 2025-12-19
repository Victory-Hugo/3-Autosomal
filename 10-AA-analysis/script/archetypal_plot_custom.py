#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
自定义版本的archetypal-plot工具
功能与原版一致但按照FAM文件中的群体进行着色

Created on: December 19, 2025
Author: Custom implementation based on archetypal-plot
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from math import pi, cos, sin


def read_fam_file(fam_file_path):
    """读取FAM文件获取群体信息"""
    fam_data = pd.read_csv(fam_file_path, sep='\t', header=None, 
                          names=['population', 'sample_id', 'father', 'mother', 'sex', 'phenotype'])
    return fam_data


def create_population_colors(populations, color_scheme='tab20'):
    """为每个群体创建颜色映射"""
    unique_pops = list(pd.Series(populations).unique())
    n_pops = len(unique_pops)
    
    if color_scheme == 'tab20':
        if n_pops <= 20:
            colors = plt.cm.tab20(np.linspace(0, 1, 20))[:n_pops]
        else:
            colors = plt.cm.nipy_spectral(np.linspace(0, 1, n_pops))
    elif color_scheme == 'Set3':
        colors = plt.cm.Set3(np.linspace(0, 1, n_pops))
    elif color_scheme == 'viridis':
        colors = plt.cm.viridis(np.linspace(0, 1, n_pops))
    else:
        colors = plt.cm.tab20(np.linspace(0, 1, n_pops))
    
    color_map = dict(zip(unique_pops, colors))
    return color_map, unique_pops


def read_q_file(q_file_path, sorted_by_max=False):
    """读取Q文件"""
    data = pd.read_csv(q_file_path, header=None, sep=" ")
    
    if sorted_by_max:
        data['Max'] = data.idxmax(axis=1)
        data = data.sort_values(by='Max')
        data.drop('Max', axis=1, inplace=True)
    
    return data


def bar_simple_population(q_data, populations, fam_data, title="", dpi=200, 
                         color_scheme='tab20', output_prefix="aa_custom"):
    """创建按群体着色的条形图"""
    fig, ax = plt.subplots(figsize=(15, 8))
    k = q_data.shape[1]
    
    # 创建颜色映射
    color_map, unique_pops = create_population_colors(populations, color_scheme)
    
    # 按群体排序数据
    sort_df = pd.DataFrame({'pop': populations, 'index': range(len(populations))})
    sort_df = sort_df.sort_values('pop')
    sorted_indices = sort_df['index'].values
    
    sorted_q_data = q_data.iloc[sorted_indices]
    sorted_populations = [populations[i] for i in sorted_indices]
    sorted_sample_ids = [fam_data['sample_id'].iloc[i] for i in sorted_indices]
    
    # 创建x轴位置
    x_pos = np.arange(len(sorted_q_data))
    
    # 绘制堆叠条形图
    bottom = np.zeros(len(sorted_q_data))
    for i in range(k):
        ax.bar(x_pos, sorted_q_data.iloc[:, i], bottom=bottom, 
               width=1.0, label=f'Archetype {i+1}')
        bottom += sorted_q_data.iloc[:, i]
    
    # 设置x轴标签颜色
    if len(sorted_q_data) > 100:
        ax.set_xticks([])
        # 在底部添加群体标签
        pop_positions = {}
        for i, pop in enumerate(sorted_populations):
            if pop not in pop_positions:
                pop_positions[pop] = []
            pop_positions[pop].append(i)
        
        # 绘制群体区域标记
        y_offset = -0.1
        for pop, positions in pop_positions.items():
            start_pos = min(positions)
            end_pos = max(positions)
            center_pos = (start_pos + end_pos) / 2
            ax.axvspan(start_pos-0.5, end_pos+0.5, alpha=0.3, 
                      color=color_map[pop], ymin=0, ymax=0.05)
            if end_pos - start_pos > 5:  # 只标记较大的群体
                ax.text(center_pos, y_offset, pop, ha='center', va='top',
                       rotation=45, fontsize=8, color=color_map[pop])
    else:
        ax.set_xticks(x_pos[::max(1, len(x_pos)//50)])
        ax.set_xticklabels([sorted_sample_ids[i] for i in range(0, len(x_pos), max(1, len(x_pos)//50))],
                          rotation=90, fontsize=8)
    
    ax.set_ylabel('Archetypal Composition')
    ax.set_title(f'{title} K={k}', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_bar_k{k}.jpg', dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"保存条形图: {output_prefix}_bar_k{k}.jpg")


def plot_simplex_population(q_data, populations, fam_data, title="", 
                           color_scheme='tab20', output_prefix="aa_custom"):
    """创建按群体着色的单纯形图"""
    k = q_data.shape[1]
    
    if k < 3:
        print("警告: 单纯形图至少需要3个原型")
        return
    
    # 创建颜色映射
    color_map, unique_pops = create_population_colors(populations, color_scheme)
    
    # 转置数据
    alfa = q_data.T
    archetype_num = alfa.shape[0]
    
    # 创建单纯形基
    sides = archetype_num
    basis = np.array([
        [cos(2*i*pi/sides + 90*pi/180), sin(2*i*pi/sides + 90*pi/180)]
        for i in range(sides)
    ])
    
    # 投影数据到2D
    data = alfa.T
    newdata = np.dot(data, basis)
    
    # 创建图形
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    
    # 绘制单纯形边界
    boundary_x = [basis[i, 0] for i in range(sides)] + [basis[0, 0]]
    boundary_y = [basis[i, 1] for i in range(sides)] + [basis[0, 1]]
    ax.plot(boundary_x, boundary_y, 'k-', linewidth=2, zorder=1)
    
    # 添加原型标签
    labels = [f'A{i+1}' for i in range(archetype_num)]
    label_offset = 0.1
    for i, label in enumerate(labels):
        x = basis[i, 0]
        y = basis[i, 1]
        ax.text(x*(1 + label_offset), y*(1 + label_offset), label,
                ha='center', va='center', fontsize=12, fontweight='bold')
    
    # 按群体绘制散点
    for pop in unique_pops:
        mask = [p == pop for p in populations]
        if any(mask):
            pop_data = newdata[mask]
            ax.scatter(pop_data[:, 0], pop_data[:, 1], 
                      c=[color_map[pop]], label=pop, alpha=0.7, s=50, zorder=3)
    
    # 设置图形属性
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.set_aspect('equal')
    ax.set_title(f'{title} K={archetype_num}', fontsize=16)
    
    # 添加图例
    legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
                      frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_alpha(0.9)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_simplex_k{archetype_num}.jpg', 
                dpi=500, bbox_inches='tight')
    plt.close()
    
    print(f"保存单纯形图: {output_prefix}_simplex_k{archetype_num}.jpg")


def bar_html_population(q_data, populations, fam_data, title="", output_prefix="aa_custom"):
    """创建交互式HTML条形图"""
    k = q_data.shape[1]
    
    # 准备数据
    plot_data = []
    for i in range(len(q_data)):
        sample_id = fam_data['sample_id'].iloc[i]
        population = populations[i]
        for j in range(k):
            plot_data.append({
                'sample_index': i,
                'sample_id': sample_id,
                'population': population,
                'archetype': f'A{j+1}',
                'value': q_data.iloc[i, j],
                'archetype_num': j
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # 创建堆叠条形图
    fig = px.bar(plot_df, x='sample_index', y='value', color='archetype',
                hover_data=['sample_id', 'population'],
                color_discrete_sequence=px.colors.qualitative.Set3,
                title=f"{title} K={k}")
    
    fig.update_layout(xaxis_title="Sample Index", yaxis_title="Archetypal Composition",
                     xaxis_type='category', height=600)
    
    fig.write_html(f'{output_prefix}_bar_k{k}_interactive.html')
    print(f"保存交互式条形图: {output_prefix}_bar_k{k}_interactive.html")


def plot_simplex_html_population(q_data, populations, fam_data, title="", output_prefix="aa_custom"):
    """创建交互式HTML单纯形图"""
    k = q_data.shape[1]
    
    if k < 3:
        print("警告: 单纯形图至少需要3个原型")
        return
    
    # 转置数据并投影
    alfa = q_data.T
    archetype_num = alfa.shape[0]
    sides = archetype_num
    
    basis = np.array([
        [cos(2*i*pi/sides + 90*pi/180), sin(2*i*pi/sides + 90*pi/180)]
        for i in range(sides)
    ])
    
    data = alfa.T
    newdata = np.dot(data, basis)
    
    # 创建边界
    boundary_x = [basis[i, 0] for i in range(sides)] + [basis[0, 0]]
    boundary_y = [basis[i, 1] for i in range(sides)] + [basis[0, 1]]
    
    # 创建散点图
    fig = px.scatter(x=newdata[:, 0], y=newdata[:, 1], color=populations,
                    hover_data={'sample_id': fam_data['sample_id']},
                    title=f'{title} K={archetype_num}')
    
    # 添加边界
    fig.add_trace(go.Scatter(x=boundary_x, y=boundary_y, mode='lines',
                           line=dict(color='black', width=2), name='Boundary'))
    
    # 添加原型标签
    labels = [f'A{i+1}' for i in range(archetype_num)]
    label_offset = 0.1
    for i, label in enumerate(labels):
        x = basis[i, 0] * (1 + label_offset)
        y = basis[i, 1] * (1 + label_offset)
        fig.add_annotation(x=x, y=y, text=label, showarrow=False,
                          font=dict(size=14, color="black"))
    
    fig.update_layout(yaxis_scaleanchor="x", yaxis_scaleratio=1, height=700)
    fig.write_html(f'{output_prefix}_simplex_k{archetype_num}_interactive.html')
    print(f"保存交互式单纯形图: {output_prefix}_simplex_k{archetype_num}_interactive.html")


def main():
    parser = argparse.ArgumentParser(description="自定义archetypal-plot工具，按群体着色")
    parser.add_argument('-i', '--input_file', required=True, type=str, 
                       help='Q文件路径')
    parser.add_argument('-f', '--fam_file', required=True, type=str,
                       help='FAM文件路径')
    parser.add_argument('-p', '--plot_type', required=True, type=str,
                       choices=['bar_simple', 'bar_html', 'plot_simplex', 'plot_simplex_html', 'all'],
                       help='图形类型')
    parser.add_argument('-so', '--sorted', action='store_true', default=False,
                       help='按最大原型组分排序')
    parser.add_argument('-dpi', '--dpi_number', type=int, default=200,
                       help='图像分辨率')
    parser.add_argument('-title', '--data_title', type=str, default="",
                       help='图形标题')
    parser.add_argument('-c', '--color_scheme', type=str, default='tab20',
                       choices=['tab20', 'Set3', 'viridis'],
                       help='颜色方案')
    parser.add_argument('-o', '--output_prefix', type=str, default="aa_custom",
                       help='输出文件前缀')
    
    args = parser.parse_args()
    
    # 读取数据
    print("读取数据文件...")
    q_data = read_q_file(args.input_file, args.sorted)
    fam_data = read_fam_file(args.fam_file)
    populations = fam_data['population'].tolist()
    
    print(f"样本数量: {len(q_data)}")
    print(f"原型数量: {q_data.shape[1]}")
    print(f"群体数量: {len(set(populations))}")
    
    # 生成图形
    if args.plot_type == 'bar_simple' or args.plot_type == 'all':
        print("生成条形图...")
        bar_simple_population(q_data, populations, fam_data, args.data_title, 
                            args.dpi_number, args.color_scheme, args.output_prefix)
    
    if args.plot_type == 'bar_html' or args.plot_type == 'all':
        print("生成交互式条形图...")
        bar_html_population(q_data, populations, fam_data, args.data_title, args.output_prefix)
    
    if args.plot_type == 'plot_simplex' or args.plot_type == 'all':
        print("生成单纯形图...")
        plot_simplex_population(q_data, populations, fam_data, args.data_title, 
                               args.color_scheme, args.output_prefix)
    
    if args.plot_type == 'plot_simplex_html' or args.plot_type == 'all':
        print("生成交互式单纯形图...")
        plot_simplex_html_population(q_data, populations, fam_data, args.data_title, args.output_prefix)
    
    print("所有图形已生成完成！")


if __name__ == '__main__':
    main()