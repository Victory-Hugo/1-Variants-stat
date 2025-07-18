# !/bin/bash

VAR_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/conf/细分/'

/home/luolintao/miniconda3/envs/pyg/bin/python3 \
     /mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/python/7-分箱堆叠.py \
  --global_csv ${VAR_DIR}/Global.var.csv \
  --eas_csv    ${VAR_DIR}/East_Asia.var.csv \
  --out_dir    "/mnt/f/OneDrive/文档（科研）/脚本/Download/1-Variants-stat/output/"

#TODO 下一步使用pipe/4-可视化MAF.R 进行可视化