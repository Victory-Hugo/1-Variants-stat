#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1-伪二倍体文件统计_带MAF.py

在原有功能基础上，新增 MAF 列，输出每个 allele 的具体 MAF 值（如 0.25%），暂不做分箱。
"""

import argparse
import sys
import csv
import os
from cyvcf2 import VCF

def main():
    parser = argparse.ArgumentParser(
        description="伪二倍体 VCF 统计脚本（带 MAF 列）"
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

    try:
        vcf = VCF(args.vcf)
    except Exception as e:
        sys.exit(f"无法打开 VCF：{e}")

    # ----- 写入变异详情（含 Freq, Type, Special, MAF 列） -----
    try:
        var_f = open(args.var_out, "w", newline="", encoding="utf-8")
    except Exception as e:
        sys.exit(f"无法创建 {args.var_out}：{e}")
    var_writer = csv.writer(var_f)
    var_writer.writerow([
        "CHROM", "POS", "REF", "ALT", "AC", "Source",
        "Freq", "Type", "Special", "MAF"
    ])

    for var in vcf:
        genos = var.genotypes
        n_alt = len(var.ALT)
        ac_list = [0] * n_alt
        an = 0

        # 统计等位基因个数（只算同型非缺失的基因型）
        for g in genos:
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

        # 类型分类（基于位点，对所有ALT相同）
        type_label = "SNV" if var.is_snp else "Indel"

        for alt, ac_val in zip(var.ALT, ac_list):
            if ac_val > 0:
                # 计算等位基因频率 AF
                af_individual = ac_val / an

                # 计算 MAF：
                # MAF = min(af_individual, 1 - af_individual)
                maf_individual = af_individual if af_individual <= 0.5 else 1 - af_individual

                # 格式化为百分比字符串，保留两位小数
                maf_pct_str = f"{maf_individual * 100:.2f}%"

                # 基于个体 MAF 进行频率分类（原有逻辑）
                if maf_individual >= 0.05:
                    freq_label_individual = "Common"
                elif maf_individual >= 0.01:
                    freq_label_individual = "LowFreq"
                elif maf_individual >= 0.001:
                    freq_label_individual = "Rare"
                else:
                    freq_label_individual = "UltraRare"

                # 基于个体 AC 计算 special 标签
                if ac_val == 1:
                    special_label_individual = "Singleton"
                elif ac_val == 2:
                    special_label_individual = "Doubleton"
                else:
                    special_label_individual = ""

                var_writer.writerow([
                    var.CHROM,
                    var.POS,
                    var.REF,
                    alt,
                    ac_val,
                    base,
                    freq_label_individual,
                    type_label,
                    special_label_individual,
                    maf_pct_str
                ])

    var_f.close()

    # ----- 剩余部分保持不变：从详情文件读回统计汇总 -----
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
        f"Done. 统计文件：{args.out}；"
        f"变异详情（含 MAF 列）：{args.var_out}"
    )

if __name__ == "__main__":
    main()
