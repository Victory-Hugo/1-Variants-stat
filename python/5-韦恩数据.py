#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_venn_sets.py

将 merged_all_sources.csv 转换为 Venn 图所需的“列是集合”格式：
- 每列对应一个地区（Africa, Central_Asia, ...）
- 列中的值是唯一的变异 ID，带有文字前缀（例如 var1, var2, ...）
- 列长度相同，不足部分填空（""）

用法：
    python generate_venn_sets.py \
        --input  /path/to/merged_all_sources.csv \
        --output /path/to/venn_sets.csv
"""

import argparse
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="生成 Venn 图集合列格式（带文字前缀）")
    parser.add_argument("--input", "-i", required=True, help="merged_all_sources.csv 路径")
    parser.add_argument("--output", "-o", required=True, help="输出的 venn_sets.csv 路径")
    args = parser.parse_args()

    # 读取 merged_all_sources.csv
    df = pd.read_csv(args.input, dtype={"POS": str})
    # 地区列
    regions = [c for c in df.columns if c not in ["CHROM", "POS", "REF", "ALT"]]
    if not regions:
        sys.exit("未发现地区列，请检查输入文件格式")

    # 为每行变异分配带前缀的唯一 ID
    # 前缀使用 "var"，序号从1开始
    df["var_id"] = df.index.to_series().add(1).apply(lambda x: f"var{x}")

    # 为每个地区生成 var_id 列表
    sets = {}
    for region in regions:
        ids = df.loc[df[region] > 0, "var_id"].tolist()
        sets[region] = ids

    # 计算最长列表长度
    max_len = max(len(v) for v in sets.values())

    # 填充每个列表至相同长度
    for region in regions:
        lst = sets[region]
        padded = lst + ["" for _ in range(max_len - len(lst))]
        sets[region] = padded

    # 构建 DataFrame
    out_df = pd.DataFrame(sets, columns=regions)

    # 保存为 CSV
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out_df.to_csv(args.output, index=False, encoding="utf-8")
    print(f"已生成集合格式文件：{args.output}（共 {max_len} 行，ID前缀为 'var'）")

if __name__ == "__main__":
    main()
