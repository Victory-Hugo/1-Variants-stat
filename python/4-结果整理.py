#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_variants.py

将多个地区的 variants_with_source CSV 文件合并成一个整合表格，
每个变异（CHROM, POS, REF, ALT）对应各地区的 AC 计数。

用法：
    python merge_variants.py \
        --var-dir /path/to/variants_csvs \
        --out /path/to/merged_all_sources.csv
"""

import argparse
import glob
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="合并多个地区的变异详情 CSV，输出按地区 AC 分列的整合表"
    )
    parser.add_argument(
        "--var-dir", "-d", required=True,
        help="包含各地区 var CSV 文件的目录（例如 Africa.var.csv、Central_Asia.var.csv 等）"
    )
    parser.add_argument(
        "--out", "-o", required=True,
        help="整合后的输出 CSV 文件路径"
    )
    args = parser.parse_args()

    # 支持两种命名模式：*.var.csv 和 *.var_*.csv
    pattern1 = os.path.join(args.var_dir, "*.var.csv")
    pattern2 = os.path.join(args.var_dir, "*.var_*.csv")
    files = glob.glob(pattern1) + glob.glob(pattern2)
    files = sorted(set(files))
    if not files:
        sys.exit(f"在目录 {args.var_dir} 中未找到任何 '.var.csv' 或 '.var_*.csv' 文件")

    # 读取并拼接所有文件的关键信息
    dfs = []
    for f in files:
        df = pd.read_csv(f, dtype={'POS': str})
        if not {'CHROM','POS','REF','ALT','AC','Source'}.issubset(df.columns):
            sys.exit(f"文件 {f} 中缺少必需的列：CHROM, POS, REF, ALT, AC, Source")
        dfs.append(df[['CHROM','POS','REF','ALT','AC','Source']])

    all_df = pd.concat(dfs, ignore_index=True)

    # 透视：不同 Source 的 AC 值分列
    pivot = (
        all_df
        .pivot_table(
            index=['CHROM','POS','REF','ALT'],
            columns='Source',
            values='AC',
            aggfunc='sum',
            fill_value=0
        )
        .reset_index()
    )
    pivot.columns.name = None  # 去除列名层次

    # 写出整合结果
    pivot.to_csv(args.out, index=False, encoding='utf-8')
    print(f"已输出整合文件：{args.out}，共包含 {len(pivot)} 条变异记录")

if __name__ == "__main__":
    main()
