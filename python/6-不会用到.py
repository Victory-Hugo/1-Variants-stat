#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-韦恩图数据生成_Indel_Common.py

从变异统计结果中提取Freq为Common且Type为Indel的变异，
生成韦恩图所需的数据文件，用于比较东亚和全球数据集。

功能：
1. 读取东亚和全球的变异统计CSV文件
2. 筛选出Freq=Common且Type=Indel的变异
3. 对数据进行去重处理（解决重复的(POS,REF,ALT)组合）
4. 合并两个数据集并分配唯一的Variant_ID
5. 生成韦恩图所需的CSV文件

输出文件格式：
- 第一列：EA_Variant_ID（东亚特有的变异ID）
- 第二列：Global_Variant_ID（全球特有的变异ID）
"""

import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description="生成韦恩图数据文件（Freq=Common, Type=Indel）"
    )
    parser.add_argument("--ea-file", required=True, 
                       help="东亚变异统计CSV文件路径")
    parser.add_argument("--global-file", required=True, 
                       help="全球变异统计CSV文件路径")
    parser.add_argument("--output", required=True, 
                       help="输出韦恩图数据CSV文件路径")
    parser.add_argument("--dedup-method", choices=['max', 'min', 'first'], 
                       default='max',
                       help="去重方法：max(保留AC最大), min(保留AC最小), first(保留第一个) [默认: max]")
    args = parser.parse_args()
    
    print("=== 韦恩图数据生成工具 ===")
    print(f"东亚数据文件: {args.ea_file}")
    print(f"全球数据文件: {args.global_file}")
    print(f"输出文件: {args.output}")
    print(f"去重方法: {args.dedup_method}")
    
    # 检查输入文件是否存在
    if not os.path.exists(args.ea_file):
        sys.exit(f"错误：东亚数据文件不存在: {args.ea_file}")
    if not os.path.exists(args.global_file):
        sys.exit(f"错误：全球数据文件不存在: {args.global_file}")
    
    try:
        # 1. 读取数据
        print("\n步骤1：读取数据...")
        df_EA = pd.read_csv(args.ea_file)
        df_Global = pd.read_csv(args.global_file)
        print(f"东亚数据: {len(df_EA)} 行")
        print(f"全球数据: {len(df_Global)} 行")
        
        # 2. 筛选Common且Indel的变异
        print("\n步骤2：筛选Freq=Common且Type=Indel的变异...")
        df_EA_filtered = df_EA[(df_EA['Freq'] == 'Common') & (df_EA['Type'] == 'Indel')]
        df_Global_filtered = df_Global[(df_Global['Freq'] == 'Common') & (df_Global['Type'] == 'Indel')]
        
        print(f"东亚筛选后: {len(df_EA_filtered)} 行")
        print(f"全球筛选后: {len(df_Global_filtered)} 行")
        
        # 3. 检查并处理重复
        print("\n步骤3：检查重复情况...")
        
        # 检查东亚数据重复
        ea_duplicated = df_EA_filtered.duplicated(subset=['POS', 'REF', 'ALT']).sum()
        if ea_duplicated > 0:
            print(f"东亚数据中发现 {ea_duplicated} 个重复的(POS,REF,ALT)组合，正在去重...")
            if args.dedup_method == 'max':
                df_EA_clean = df_EA_filtered.loc[df_EA_filtered.groupby(['POS', 'REF', 'ALT'])['AC'].idxmax()]
            elif args.dedup_method == 'min':
                df_EA_clean = df_EA_filtered.loc[df_EA_filtered.groupby(['POS', 'REF', 'ALT'])['AC'].idxmin()]
            else:  # first
                df_EA_clean = df_EA_filtered.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first')
            print(f"东亚去重后: {len(df_EA_clean)} 行")
        else:
            df_EA_clean = df_EA_filtered
            print("东亚数据无重复")
        
        # 检查全球数据重复
        global_duplicated = df_Global_filtered.duplicated(subset=['POS', 'REF', 'ALT']).sum()
        if global_duplicated > 0:
            print(f"全球数据中发现 {global_duplicated} 个重复的(POS,REF,ALT)组合，正在去重...")
            if args.dedup_method == 'max':
                df_Global_clean = df_Global_filtered.loc[df_Global_filtered.groupby(['POS', 'REF', 'ALT'])['AC'].idxmax()]
            elif args.dedup_method == 'min':
                df_Global_clean = df_Global_filtered.loc[df_Global_filtered.groupby(['POS', 'REF', 'ALT'])['AC'].idxmin()]
            else:  # first
                df_Global_clean = df_Global_filtered.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first')
            print(f"全球去重后: {len(df_Global_clean)} 行")
        else:
            df_Global_clean = df_Global_filtered
            print("全球数据无重复")
        
        # 4. 合并数据并分配Variant_ID
        print("\n步骤4：合并数据并分配变异ID...")
        df_merged = pd.merge(df_EA_clean, df_Global_clean, 
                            on=['POS', 'REF', 'ALT'], 
                            suffixes=('_EA', '_Global'), 
                            how='outer')
        
        # 为每个唯一的(POS,REF,ALT)组合分配ID
        df_merged['Variant_ID'] = df_merged.groupby(['POS', 'REF', 'ALT']).ngroup()
        
        print(f"合并后总行数: {len(df_merged)}")
        print(f"唯一变异数: {df_merged['Variant_ID'].nunique()}")
        
        # 5. 生成韦恩图数据
        print("\n步骤5：生成韦恩图数据...")
        
        # 提取东亚和全球的Variant_ID
        df_EA_venn = df_merged.loc[df_merged['Source_EA'] == 'East_Asia', ['Variant_ID']]
        df_Global_venn = df_merged.loc[df_merged['Source_Global'] == 'Global', ['Variant_ID']]
        
        # 重命名列
        df_EA_venn = df_EA_venn.rename(columns={'Variant_ID': 'EA_Variant_ID'})
        df_Global_venn = df_Global_venn.rename(columns={'Variant_ID': 'Global_Variant_ID'})
        
        # 合并为最终的韦恩图数据
        df_venn = pd.concat([df_EA_venn, df_Global_venn], axis=1)
        
        print(f"东亚独有变异: {df_venn['EA_Variant_ID'].dropna().nunique()}")
        print(f"全球独有变异: {df_venn['Global_Variant_ID'].dropna().nunique()}")
        
        # 计算交集
        ea_variants = set(df_EA_venn['EA_Variant_ID'].dropna())
        global_variants = set(df_Global_venn['Global_Variant_ID'].dropna())
        intersection = len(ea_variants & global_variants)
        print(f"共同变异: {intersection}")
        
        # 6. 保存结果
        print(f"\n步骤6：保存结果到 {args.output}...")
        
        # 创建输出目录（如果不存在）
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 保存韦恩图数据
        df_venn.to_csv(args.output, index=False, float_format='%.0f')
        
        # 生成统计报告
        report_file = args.output.replace('.csv', '_report.txt')
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("韦恩图数据生成报告\n")
            f.write("=" * 30 + "\n\n")
            f.write(f"筛选条件: Freq=Common, Type=Indel\n")
            f.write(f"去重方法: {args.dedup_method}\n\n")
            f.write("数据统计:\n")
            f.write(f"- 东亚原始数据: {len(df_EA)} 行\n")
            f.write(f"- 全球原始数据: {len(df_Global)} 行\n")
            f.write(f"- 东亚筛选后: {len(df_EA_filtered)} 行\n")
            f.write(f"- 全球筛选后: {len(df_Global_filtered)} 行\n")
            f.write(f"- 东亚去重后: {len(df_EA_clean)} 行\n")
            f.write(f"- 全球去重后: {len(df_Global_clean)} 行\n\n")
            f.write("韦恩图统计:\n")
            f.write(f"- 东亚独有变异: {df_venn['EA_Variant_ID'].dropna().nunique()}\n")
            f.write(f"- 全球独有变异: {df_venn['Global_Variant_ID'].dropna().nunique()}\n")
            f.write(f"- 共同变异: {intersection}\n")
            f.write(f"- 总唯一变异: {df_merged['Variant_ID'].nunique()}\n")
        
        print(f"✅ 成功！")
        print(f"   韦恩图数据: {args.output}")
        print(f"   统计报告: {report_file}")
        
    except Exception as e:
        sys.exit(f"处理过程中出现错误: {e}")

if __name__ == "__main__":
    main()
