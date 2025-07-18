#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
count_variants_per_sample.py

对超大 VCF 文件流式统计每个样本的 variants per genome（任何非 0/0 都计作一次变异），
输出 CSV，并自动从 VCF 基本名填写 Source 列。
"""

import os
import sys
import argparse
import csv
from cyvcf2 import VCF

def count_variants(vcf_path):
    """
    遍历 VCF，返回 (samples, counts)：
      - samples: 样本列表
      - counts: 每个样本的变异计数（任何非 0/0 的基因型都算一次变异）
    """
    vcf = VCF(vcf_path)
    samples = vcf.samples
    counts = [0] * len(samples)

    for var in vcf:
        # 只统计 FILTER=PASS 的记录
        if var.FILTER not in (None, [], 'PASS'):
            continue
        # 遍历每个样本的基因型
        for i, gt in enumerate(var.genotypes):
            a0, a1, _ = gt
            # 跳过缺失或无效
            if a0 is None or a1 is None or a0 < 0 or a1 < 0:
                continue
            # 只要任一等位基因 != 0，就算一次变异
            if a0 != 0 or a1 != 0:
                counts[i] += 1

    return samples, counts

def derive_source(vcf_path):
    """
    从文件名中提取基本名作为 Source，去掉 .vcf 或 .vcf.gz 后缀
    """
    base = os.path.basename(vcf_path)
    for ext in ('.vcf.gz', '.vcf'):
        if base.endswith(ext):
            return base[:-len(ext)]
    return base

def main():
    parser = argparse.ArgumentParser(
        description="流式统计每个样本的 variants per genome 并输出 CSV"
    )
    parser.add_argument(
        "--vcf", required=True,
        help="输入 VCF(.gz) 文件路径"
    )
    parser.add_argument(
        "--out", required=True,
        help="输出 CSV 文件路径"
    )
    args = parser.parse_args()

    vcf_path = args.vcf
    out_csv  = args.out

    if not os.path.exists(vcf_path):
        sys.exit(f"Error: 找不到 VCF 文件 {vcf_path}")

    source = derive_source(vcf_path)
    print(f"[INFO] 开始统计：{vcf_path} （Source={source}）")

    samples, counts = count_variants(vcf_path)

    # 流式写出 CSV
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'variant_count', 'Source'])
        for sample, cnt in zip(samples, counts):
            writer.writerow([sample, cnt, source])

    print(f"[INFO] 完成，结果已保存到：{out_csv}")

if __name__ == "__main__":
    main()
