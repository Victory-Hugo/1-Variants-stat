#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-真二倍体文件统计_带MAF.py

在真二倍体 VCF.gz 基础上，逐等位基因统计并输出 MAF：
- 按频率分类：Common/LowFreq/Rare/UltraRare
- 类型分类：SNV/Indel
- singleton/doubleton 统计

输出两份 CSV：
1) 汇总统计：Frequency/Type/Special 分类计数（含 Source 列）
2) 变异详情：CHROM, POS, REF, ALT, AC, Source, Freq, Type, Special, MAF
"""

import argparse
import sys
import csv
import os
from cyvcf2 import VCF

def main():
    parser = argparse.ArgumentParser(
        description="真二倍体 VCF 统计脚本（带 MAF 列）"
    )
    parser.add_argument("-i", "--vcf", required=True, help="输入 VCF.gz 文件")
    parser.add_argument("-o", "--out", required=True, help="输出汇总统计 CSV")
    parser.add_argument("-v", "--var-out", required=True, help="输出变异详情 CSV")
    args = parser.parse_args()

    # 计算 Source 名称
    base = os.path.basename(args.vcf)
    for ext in (".vcf.gz", ".vcf"):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    # 打开 VCF
    try:
        vcf = VCF(args.vcf)
    except Exception as e:
        sys.exit(f"无法打开 VCF：{e}")

    # 写入变异详情 CSV（含 MAF 列）
    try:
        var_f = open(args.var_out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建 {args.var_out}：{e}")
    var_writer = csv.writer(var_f)
    var_writer.writerow([
        "CHROM", "POS", "REF", "ALT", "AC", "Source",
        "Freq", "Type", "Special", "MAF"
    ])

    # 遍历每个位点
    for var in vcf:
        # 从 INFO 拿 AN 和 AC
        an = var.INFO.get("AN")
        ac_info = var.INFO.get("AC")
        if an is None or an == 0 or ac_info is None:
            continue

        # 可能多等位
        ac_list = ac_info if isinstance(ac_info, list) else [ac_info]

        # 变异类型
        type_label = "SNV" if var.is_snp else "Indel"

        # 逐等位统计
        for alt, ac_val in zip(var.ALT, ac_list):
            if ac_val == 0:
                continue

            # 计算 AF 和 MAF
            af = ac_val / an
            maf = af if af <= 0.5 else 1 - af

            # 格式化 MAF（百分比，两位小数）
            maf_str = f"{maf * 100:.2f}%"

            # 频率分类
            if maf >= 0.05:
                freq_label = "Common"
            elif maf >= 0.01:
                freq_label = "LowFreq"
            elif maf >= 0.001:
                freq_label = "Rare"
            else:
                freq_label = "UltraRare"

            # special 分类
            if ac_val == 1:
                special_label = "Singleton"
            elif ac_val == 2:
                special_label = "Doubleton"
            else:
                special_label = ""

            # 写入详情行
            var_writer.writerow([
                var.CHROM,
                var.POS,
                var.REF,
                alt,
                ac_val,
                base,
                freq_label,
                type_label,
                special_label,
                maf_str
            ])

    var_f.close()

    # 从详情 CSV 读回，统计汇总
    freq_counts    = {}
    type_counts    = {}
    special_counts = {}
    with open(args.var_out, newline="", encoding="utf-8") as vf:
        reader = csv.DictReader(vf)
        for row in reader:
            f = row["Freq"]
            freq_counts[f] = freq_counts.get(f, 0) + 1
            t = row["Type"]
            type_counts[t] = type_counts.get(t, 0) + 1
            s = row["Special"]
            special_counts[s] = special_counts.get(s, 0) + 1

    # 写入汇总统计 CSV
    try:
        with open(args.out, "w", newline="", encoding="utf-8") as out_f:
            writer = csv.writer(out_f)
            writer.writerow(["Category", "Class", "Count", "Source"])
            for cls, cnt in freq_counts.items():
                writer.writerow(["Frequency", cls, cnt, base])
            for cls, cnt in type_counts.items():
                writer.writerow(["Type", cls, cnt, base])
            for cls, cnt in special_counts.items():
                writer.writerow(["Special", cls, cnt, base])
    except Exception as e:
        sys.exit(f"无法写入 {args.out}：{e}")

    print(
        f"Done. 汇总统计：{args.out}；"
        f"变异详情（含 MAF 列）：{args.var_out}"
    )

if __name__ == "__main__":
    main()
