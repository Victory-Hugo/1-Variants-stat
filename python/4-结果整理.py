#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_variants.py

功能一：将多个地区的 variants_with_source CSV 文件合并成一个按地区 AC 分列的整合表格。
功能二：将指定目录下所有不以 .var.csv 结尾的 .csv 文件合并成一个 CSV，只保留一次表头。

用法示例：
    # 仅做 .var.csv 透视整合
    python merge_variants.py \
        --var-dir /path/to/variants_csvs \
        --out pivoted.csv

    # 仅做普通 CSV 合并
    python merge_variants.py \
        --merge-dir /path/to/csvs \
        --merge-out merged.csv

    # 两项都做
    python merge_variants.py \
        --var-dir /path/to/variants_csvs \
        --out pivoted.csv \
        --merge-dir /path/to/csvs \
        --merge-out merged.csv
"""

import argparse
import glob
import os
import sys
import pandas as pd

def merge_plain_csv(input_dir, output_file):
    """合并普通 CSV（不包括以 .var.csv 或 .var_*.csv 结尾的），保留一次表头"""
    import glob, os, sys

    pattern = os.path.join(input_dir, "*.csv")
    # 排除所有 *.var.csv 和 *.var_*.csv
    all_files = glob.glob(pattern)
    files = sorted(f for f in all_files
                   if not (f.endswith(".var.csv") or os.path.basename(f).startswith("*.var_")))
    if not files:
        sys.exit(f"在目录 {input_dir} 中未找到符合条件的普通 CSV 文件。")

    with open(output_file, 'w', encoding='utf-8', newline='') as fo:
        for idx, f in enumerate(files):
            with open(f, 'r', encoding='utf-8') as fi:
                lines = fi.readlines()
                if idx == 0:
                    # 第一份文件，连同表头一起写入
                    fo.writelines(lines)
                else:
                    # 后续文件，若仅有表头（len<=1）则跳过，否则写入除表头外的行
                    if len(lines) > 1:
                        fo.writelines(lines[1:])
    print(f"[普通 CSV 合并] 已将 {len(files)} 个文件合并，并保存为：{output_file}")


def merge_and_pivot_vars(var_dir, pivot_out):
    """合并 .var.csv / .var_*.csv 并按 Source 透视 AC"""
    # 支持两种命名模式：*.var.csv 和 *.var_*.csv
    pattern1 = os.path.join(var_dir, "*.var.csv")
    pattern2 = os.path.join(var_dir, "*.var_*.csv")
    files = sorted(set(glob.glob(pattern1) + glob.glob(pattern2)))
    if not files:
        sys.exit(f"在目录 {var_dir} 中未找到 '.var.csv' 或 '.var_*.csv' 文件。")

    dfs = []
    for f in files:
        df = pd.read_csv(f, dtype={'POS': str})
        needed = {'CHROM','POS','REF','ALT','AC','Source'}
        if not needed.issubset(df.columns):
            sys.exit(f"文件 {f} 中缺少必需的列：{','.join(needed)}")
        dfs.append(df[['CHROM','POS','REF','ALT','AC','Source']])

    all_df = pd.concat(dfs, ignore_index=True)
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
    pivot.to_csv(pivot_out, index=False, encoding='utf-8')
    print(f"[.var.csv 透视] 已输出整合文件：{pivot_out}，共 {len(pivot)} 条记录")

def main():
    parser = argparse.ArgumentParser(
        description="合并普通 CSV 或者按地区透视 .var.csv 变异表"
    )
    parser.add_argument(
        "--var-dir", "-d",
        help="包含各地区 var CSV 文件的目录（例如 Africa.var.csv、Central_Asia.var.csv 等）"
    )
    parser.add_argument(
        "--out", "-o",
        help="透视整合后的 .var.csv 输出文件路径"
    )
    parser.add_argument(
        "--merge-dir", "-m",
        help="要合并的普通 CSV 文件目录（排除 *.var.csv）"
    )
    parser.add_argument(
        "--merge-out", "-r",
        help="普通 CSV 合并后的输出文件路径"
    )
    args = parser.parse_args()

    # 校验参数
    if args.var_dir and not args.out:
        sys.exit("指定了 --var-dir，请同时指定 --out 输出路径。")
    if args.out and not args.var_dir:
        sys.exit("指定了 --out，请同时指定 --var-dir 输入目录。")
    if args.merge_dir and not args.merge_out:
        sys.exit("指定了 --merge-dir，请同时指定 --merge-out 输出路径。")
    if args.merge_out and not args.merge_dir:
        sys.exit("指定了 --merge-out，请同时指定 --merge-dir 输入目录。")
    if not ( (args.var_dir and args.out) or (args.merge_dir and args.merge_out) ):
        parser.print_help()
        sys.exit(1)

    # 执行功能
    if args.merge_dir and args.merge_out:
        merge_plain_csv(args.merge_dir, args.merge_out)
    if args.var_dir and args.out:
        merge_and_pivot_vars(args.var_dir, args.out)

if __name__ == "__main__":
    main()
