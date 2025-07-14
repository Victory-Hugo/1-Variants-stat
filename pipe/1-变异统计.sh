#! /bin/bash
# -----------------------------------------------------------------------------
# 脚本功能说明：
# 本脚本调用三个 Python 脚本，对指定的 VCF 文件（Africa.vcf.gz）进行不同类型的变异统计分析，并将结果输出到指定目录。
#
# 1. 1-二倍体文件统计.py
#    - 对 VCF 文件进行二倍体统计分析。
#    - 输出总体统计结果（Africa_伪二倍体.csv）和变异位点统计（Africa.var_伪二倍体.csv）。
#
# 2. 2-伪二倍体文件统计.py
#    - 对 VCF 文件进行伪二倍体统计分析。
#    - 输出总体统计结果（Africa_二倍体.csv）和变异位点统计（Africa.var_二倍体.csv）。
#
# 3. 3-单倍体文件统计.py
#    - 对 VCF 文件进行单倍体统计分析。
#    - 输出总体统计结果（Africa_单倍体.csv）和变异位点统计（Africa.var_单倍体.csv）。
#
# 每个 Python 脚本均通过命令行参数指定输入 VCF 文件、输出统计文件和变异位点统计文件。
# -----------------------------------------------------------------------------
#! 请根据自己的文件类型使用下列任意一个脚本
# 统计
/home/luolintao/miniconda3/envs/pyg/bin/python3 \
 "/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/2-伪二倍体文件统计.py" \
 --vcf "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/data/Africa.vcf.gz" \
 --out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.csv" \
 --var-out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.var.csv" \

# /home/luolintao/miniconda3/envs/pyg/bin/python3 \
#  "/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/1-二倍体文件统计.py" \
#  --vcf "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/data/Africa.vcf.gz" \
#  --out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.csv" \
#  --var-out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.var.csv" \



# /home/luolintao/miniconda3/envs/pyg/bin/python3 \
#  "/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/3-单倍体文件统计.py" \
#  --vcf "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/data/Africa.vcf.gz" \
#  --out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.csv" \
#  --var-out "/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/Africa.var.csv" \