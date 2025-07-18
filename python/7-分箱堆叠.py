#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
maf_bin_analysis.py

对 Global.var.csv 和 East_Asia.var.csv 两份文件：
1. 基于全局 MAF 做分箱，标记东亚是否出现，输出 Bin_MAF_Comparison.csv
2. 分别对全局/东亚 MAF 做分箱统计，输出 Bin_MAF_Comparison_with_eas.csv
"""

import os
import argparse
import pandas as pd


def load_and_prepare(path: str) -> pd.DataFrame:
    """
    读取 CSV 并把 '0.06%' 形式的 MAF 转为浮点数（单位 %）。
    """
    df = pd.read_csv(path, dtype=str)
    # 去掉 '%' 并转换为 float
    df['MAF_pct'] = df['MAF'].str.rstrip('%').astype(float)
    return df


def define_bins():
    """
    返回分箱阈值和对应的标签。
    """
    bins = [0, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
    labels = [
        '<0.1', '0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-1.0',
        '1-2', '2-5', '5-10', '10-20', '20-30', '30-40', '40-50'
    ]
    return bins, labels


def compute_presence(df_global: pd.DataFrame,
                     df_eas: pd.DataFrame,
                     bins, labels) -> pd.DataFrame:
    """
    基于全局 MAF 分箱，并标记该变异在东亚中是否出现。
    """
    # 1. 全局分箱
    df = df_global.copy()
    df['maf_bin'] = pd.cut(
        df['MAF_pct'],
        bins=bins,
        right=False,
        labels=labels
    )

    # 2. 构建东亚出现标记
    eas_keys = (
        df_eas
        [['CHROM', 'POS', 'REF', 'ALT']]
        .drop_duplicates()
    )
    eas_keys['present_eas'] = 1

    # 3. 左合并
    df = df.merge(
        eas_keys,
        on=['CHROM', 'POS', 'REF', 'ALT'],
        how='left'
    )
    df['present_eas'] = df['present_eas'].fillna(0).astype(int)

    # 4. 分箱统计
    grp = df.groupby('maf_bin', observed=False)
    result = pd.DataFrame({
        'global_count'    : grp.size(),
        'east_asia_count' : grp['present_eas'].sum()
    })
    result['non_east_asia_count'] = (
        result['global_count'] - result['east_asia_count']
    )

    return result


def compute_bin_only(df_global: pd.DataFrame,
                     df_eas: pd.DataFrame,
                     bins, labels) -> pd.DataFrame:
    """
    分别对全局和东亚 MAF 做分箱统计，比较两者数量。
    """
    df_g = df_global.copy()
    df_e = df_eas.copy()

    df_g['maf_bin_global'] = pd.cut(
        df_g['MAF_pct'], bins=bins, right=False, labels=labels
    )
    df_e['maf_bin_eas'] = pd.cut(
        df_e['MAF_pct'], bins=bins, right=False, labels=labels
    )

    global_counts = (
        df_g['maf_bin_global']
        .value_counts()
        .sort_index()
        .rename('global_count')
    )
    eas_counts = (
        df_e['maf_bin_eas']
        .value_counts()
        .sort_index()
        .rename('east_asia_count')
    )

    result = pd.concat([global_counts, eas_counts], axis=1).fillna(0).astype(int)
    result['non_east_asia_count'] = (
        result['global_count'] - result['east_asia_count']
    )

    return result


def main():
    parser = argparse.ArgumentParser(
        description="MAF 分箱对比脚本"
    )
    parser.add_argument(
        "--global_csv",
        required=True,
        help="Global.var.csv 路径"
    )
    parser.add_argument(
        "--eas_csv",
        required=True,
        help="East_Asia.var.csv 路径"
    )
    parser.add_argument(
        "--out_dir",
        default="output",
        help="输出目录（默认 ./output）"
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 加载并预处理
    df_global = load_and_prepare(args.global_csv)
    df_eas    = load_and_prepare(args.eas_csv)
    bins, labels = define_bins()

    # 1) 全局分箱 + 出现标记
    res1 = compute_presence(df_global, df_eas, bins, labels)
    out1 = os.path.join(args.out_dir, "Bin_MAF_Comparison.csv")
    res1.to_csv(out1, index=True, index_label="MAF")
    print(f"已保存：{out1}")

    # 2) 分别分箱统计
    res2 = compute_bin_only(df_global, df_eas, bins, labels)
    out2 = os.path.join(args.out_dir, "Bin_MAF_Comparison_with_eas.csv")
    res2.to_csv(out2, index=True, index_label="MAF")
    print(f"已保存：{out2}")


if __name__ == "__main__":
    main()
