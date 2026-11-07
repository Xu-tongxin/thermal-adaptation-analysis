# 载入需要的包
library(pheatmap)
library(grDevices) 

# 1. 导入由 snp-dists 生成的距离矩阵
snp_dist <- read.csv("snp_matrix_no_outgroup.csv", row.names = 1, check.names = FALSE)
snp_mat <- as.matrix(snp_dist)

# 2. 过渡色（深红 -> 浅色 -> 深蓝）
# 使用提取的 HEX 色值构建 100 种渐变色
my_colors <- colorRampPalette(c("#6A0624", "#FDD2BC", "#053061"))(100)

# 3. 绘制聚类热图
pheatmap(snp_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = my_colors,             
         border_color = "grey60",        
         display_numbers = FALSE)