#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-纯单倍体文件统计.py

统计单倍体 VCF.gz 中每个位点的变异：
- 按频率分类：Common/LowFreq/Rare/UltraRare
- 类型分类：SNV/Indel
- singleton/doubleton 统计

输出两份 CSV：
1) 汇总统计：Frequency/Type/Special 分类计数（含 Source 列）
2) 变异详情：CHROM, POS, REF, ALT, AC, Source, Freq, Type, Special
"""

import argparse
import sys
import csv
import os
from cyvcf2 import VCF

def main():
    parser = argparse.ArgumentParser(
        description="纯单倍体 VCF 统计脚本（伪二倍体逻辑）"
    )
    parser.add_argument("-i", "--vcf",     required=True, help="输入 VCF.gz 文件")
    parser.add_argument("-o", "--out",     required=True, help="输出统计结果 CSV")
    parser.add_argument("-v", "--var-out", required=True, help="输出变异详情 CSV")
    args = parser.parse_args()

    # 计算来源基础名
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

    # ----- 写入变异详情（含 Freq, Type, Special 列） -----
    try:
        var_f = open(args.var-out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建 {args.var_out}：{e}")
    var_writer = csv.writer(var_f)
    var_writer.writerow([
        "CHROM", "POS", "REF", "ALT", "AC", "Source",
        "Freq", "Type", "Special"
    ])

    # 遍历每个位点
    for var in vcf:
        n_alt = len(var.ALT)
        ac_list = [0] * n_alt
        an = 0

        # 伪二倍体逻辑：只统计同型非缺失基因型
        for g in var.genotypes:
            a0, a1, _ = g
            if a0 is None or a1 is None or a0 != a1:
                continue
            an += 1
            if a0 == 0:
                continue
            idx = a0 - 1
            if 0 <= idx < n_alt:
                ac_list[idx] += 1

        total_ac = sum(ac_list)
        if an == 0 or total_ac == 0:
            continue

        # 类型分类（位点层面）
        type_label = "SNV" if var.is_snp else "Indel"

        # 写入每个 ALT 的一行
        for alt, ac_val in zip(var.ALT, ac_list):
            if ac_val == 0:
                continue

            # 计算等位基因频率 AF 和最小等位基因频率 MAF
            # ```python
            # af_ind  = ac_val / an
            # maf_ind = af_ind if af_ind <= 0.5 else 1 - af_ind
            # ```
            af_ind  = ac_val / an
            maf_ind = af_ind if af_ind <= 0.5 else 1 - af_ind

            # 频率分类
            if maf_ind >= 0.05:
                freq_label = "Common"
            elif maf_ind >= 0.01:
                freq_label = "LowFreq"
            elif maf_ind >= 0.001:
                freq_label = "Rare"
            else:
                freq_label = "UltraRare"

            # special 标签按每个 ALT 的 AC
            if ac_val == 1:
                special_label = "Singleton"
            elif ac_val == 2:
                special_label = "Doubleton"
            else:
                special_label = ""

            var_writer.writerow([
                var.CHROM,
                var.POS,
                var.REF,
                alt,
                ac_val,
                base,
                freq_label,
                type_label,
                special_label
            ])

    var_f.close()

    # ----- 读取详情文件，重新统计汇总 -----
    freq_counts    = {}
    type_counts    = {}
    special_counts = {}
    try:
        with open(args.var_out, newline="", encoding="utf-8") as vf:
            reader = csv.DictReader(vf)
            for row in reader:
                freq_counts[row["Freq"]]    = freq_counts.get(row["Freq"], 0) + 1
                type_counts[row["Type"]]    = type_counts.get(row["Type"], 0) + 1
                special_counts[row["Special"]] = special_counts.get(row["Special"], 0) + 1
    except Exception as e:
        sys.exit(f"无法读取详情文件 {args.var_out}：{e}")

    # ----- 写入汇总统计 CSV -----
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
        f"Done. 统计文件：{args.out}（基于详情文件计算）；"
        f"变异详情：{args.var_out}（含 Freq、Type、Special 列）"
    )

if __name__ == "__main__":
    main()
