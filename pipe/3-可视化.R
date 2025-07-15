library(tidyverse)
library(tidyplots)

# 读取数据文件
df_variants <- read.csv('/mnt/d/幽门螺旋杆菌/Script/分析结果/2-变异统计/output/粗分/csv/merged.csv', sep = ',')


df_variants |>
   filter(Category == 'Frequency' | Category == 'Special') |>
   tidyplot(x = Class ,y = Count, color = Class) |>
   add_barstack_absolute()
