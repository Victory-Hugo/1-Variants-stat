#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-二倍体文件统计.py

统计真二倍体 VCF.gz 中每个位点的变异（按等位基因）：
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
        description="真二倍体 VCF 统计脚本"
    )
    parser.add_argument("-i", "--vcf", required=True, help="输入 VCF.gz 文件")
    parser.add_argument("-o", "--out", required=True, help="输出统计结果 CSV")
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

    # ----- 写入变异详情（按等位基因） -----
    try:
        var_f = open(args.var_out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建 {args.var_out}：{e}")
    var_writer = csv.writer(var_f)
    var_writer.writerow([
        "CHROM", "POS", "REF", "ALT", "AC", "Source",
        "Freq", "Type", "Special"
    ])

    for var in vcf:
        # 从 INFO 里取 AC 和 AN
        ac_info = var.INFO.get("AC")
        an      = var.INFO.get("AN")
        # 跳过无效记录
        if an is None or an == 0 or ac_info is None:
            continue

        # 将 AC 转为列表（多等位基因时）
        ac_list = ac_info if isinstance(ac_info, list) else [ac_info]

        # 变异类型：同所有 ALT 共享
        type_label = "SNV" if var.is_snp else "Indel"

        # 对每个 ALT 单独计算
        for alt, ac_val in zip(var.ALT, ac_list):
            # 跳过未观测到的等位基因
            if ac_val == 0:
                continue

            # 计算等位基因频率 AF 和最小等位基因频率 MAF
            # ```python
            # af_ind  = ac_val / an
            # maf_ind = min(af_ind, 1 - af_ind)
            # ```
            af_ind  = ac_val / an
            maf_ind = af_ind if af_ind <= 0.5 else 1 - af_ind

            # 基于个体 MAF 进行频率分类
            if maf_ind >= 0.05:
                freq_label = "Common"
            elif maf_ind >= 0.01:
                freq_label = "LowFreq"
            elif maf_ind >= 0.001:
                freq_label = "Rare"
            else:
                freq_label = "UltraRare"

            # 基于 AC 计算 special 标签
            if ac_val == 1:
                special_label = "Singleton"
            elif ac_val == 2:
                special_label = "Doubleton"
            else:
                special_label = ""

            # 写入详情
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

    # ----- 从详情文件读回，重新统计保证一致 -----
    freq_counts    = {}
    type_counts    = {}
    special_counts = {}
    with open(args.var_out, newline="", encoding="utf-8") as vf:
        reader = csv.DictReader(vf)
        for row in reader:
            # 频率
            f = row["Freq"]
            freq_counts[f] = freq_counts.get(f, 0) + 1
            # 类型
            t = row["Type"]
            type_counts[t] = type_counts.get(t, 0) + 1
            # 特殊
            s = row["Special"]
            special_counts[s] = special_counts.get(s, 0) + 1

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
