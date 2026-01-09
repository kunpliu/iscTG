import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import os
import argparse

def extract_title_from_filename(filename):
    """从文件名中提取标题部分"""
    base_name = os.path.basename(filename)
    
    # 去掉 .csv 扩展名
    if base_name.endswith('.csv'):
        base_name_no_ext = base_name[:-4]
    else:
        base_name_no_ext = base_name
    
    # 去掉固定的前缀（如果存在）
    prefix_to_remove = "freq61_hbf2000_"
    if base_name_no_ext.startswith(prefix_to_remove):
        title_part = base_name_no_ext[len(prefix_to_remove):]
    else:
        title_part = base_name_no_ext
    
    return title_part

def plot_cell_type_proportions(file_path, thresholds=[0.01, 0.05, 0.10, 0.15], alpha=0.05):
    """绘制细胞类型比例热图"""
    title_part = extract_title_from_filename(file_path)
    print(f"处理文件: {file_path}")
    print(f"图表标题: {title_part}")
    
    # 读取数据
    df = pd.read_csv(file_path)
    
    # 计算每种细胞类型的总个数
    total_counts = df['type'].value_counts()
    count_matrix = {thresh: {} for thresh in thresholds}
    
    # 对每个细胞类型进行处理
    for cell_type in df['type'].unique():
        subset = df[df['type'] == cell_type]
        for thresh in thresholds:
            count = (subset['p_norm'] <= thresh).sum()
            count_matrix[thresh][cell_type] = count

    # 转换为 DataFrame 并计算比例
    count_df = pd.DataFrame(count_matrix)
    proportion_df = count_df.div(total_counts, axis=0)
    
    # 卡方检验
    chi2_results = {}
    for thresh in thresholds:
        p_values = {}
        for cell_type in count_df.index:
            total_count = total_counts[cell_type]
            less_equal = count_df.at[cell_type, thresh]
            greater = total_count - less_equal
            
            other_less_equal = count_df.loc[count_df.index != cell_type, thresh].sum()
            other_greater = count_df.loc[count_df.index != cell_type, :].sum().sum() - other_less_equal
            
            observed = [less_equal, greater]
            expected = [other_less_equal, other_greater]
            
            current_ratio = less_equal / total_count if total_count > 0 else 0
            background_ratio = other_less_equal / (other_less_equal + other_greater) if (other_less_equal + other_greater) > 0 else 0
            
            if current_ratio > background_ratio:
                chi2, p_value, _, _ = chi2_contingency([observed, expected])
                p_values[cell_type] = p_value
            else:
                p_values[cell_type] = 1.0
        
        chi2_results[thresh] = p_values
    
    chi2_results_df = pd.DataFrame(chi2_results)

    # 创建显著性标注
    annotations = proportion_df.copy()
    annotations[:] = ''
    for thresh in thresholds:
        for cell_type in count_df.index:
            if chi2_results_df.at[cell_type, thresh] < alpha:
                annotations.at[cell_type, thresh] = '*'

    # 绘制热图
    plt.figure(figsize=(8, 5))
    
    thresholds_reversed = thresholds[::-1]
    ax = sns.heatmap(
        proportion_df.T.reindex(thresholds_reversed),
        annot=annotations.T.reindex(thresholds_reversed),
        fmt='s',
        cmap="Blues",
        cbar=True,
        linewidths=0.5,
        linecolor='gray',
        vmin=0,
        vmax=proportion_df.max().max(),
        square=True,
        annot_kws={"fontsize": 12},
        cbar_kws={"shrink": 0.5}
    )
    
    plt.title(title_part, fontsize=14)
    plt.xlabel('Cell Type', fontsize=12)
    plt.ylabel('Threshold', fontsize=12)
    
    #plt.xticks(rotation=45, ha='right')
    plt.xticks()
    plt.tight_layout()
    
    # 保存图像
    output_filename = f'heatmap_{title_part}.png'
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    print(f"图表已保存: {output_filename}")
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='绘制细胞类型比例热图')
    
    # 必需参数
    parser.add_argument('input_file', help='输入CSV文件路径')
    
    # 可选参数
    parser.add_argument('--thresholds', nargs='+', type=float, 
                       default=[0.01, 0.05, 0.10, 0.15],
                       help='p值阈值（默认: 0.01 0.05 0.10 0.15）')
    parser.add_argument('--alpha', type=float, default=0.05,
                       help='显著性水平（默认: 0.05）')
    
    args = parser.parse_args()
    
    # 验证文件是否存在
    if not os.path.exists(args.input_file):
        print(f"错误：文件 '{args.input_file}' 不存在！")
        return
    
    # 执行绘图
    plot_cell_type_proportions(
        file_path=args.input_file,
        thresholds=args.thresholds,
        alpha=args.alpha
    )

if __name__ == "__main__":
    main()

#python celltype_heat.py freq61_hbf2000_mono_example.csv