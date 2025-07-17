#!/bin/bash
# -*- coding: utf-8 -*-
"""
使用示例：3-韦恩图数据生成_Indel_Common.py

这个脚本展示了如何使用韦恩图数据生成工具
"""

# 设置文件路径
BASE_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/比较东亚和全球/var"
EA_FILE="${BASE_DIR}/East_Asia.var.csv"
GLOBAL_FILE="${BASE_DIR}/Global.var.csv"
OUTPUT_FILE="${BASE_DIR}/df_Venn_Indel_Common.csv"

# 运行韦恩图数据生成脚本
echo "=== 生成韦恩图数据（Indel + Common）==="
/home/luolintao/miniconda3/envs/pyg/bin/python3 \
    /mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/6-不会用到.py \
    --ea-file "$EA_FILE" \
    --global-file "$GLOBAL_FILE" \
    --output "$OUTPUT_FILE" \
    --dedup-method max

echo "完成！"
