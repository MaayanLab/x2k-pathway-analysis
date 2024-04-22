library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# User parameters
data = "x2k" # x2k, chea, kea
tracks <- c('baseline.tumor_stage_pathological', 'binned_age', 'Sex') # sample (column) annotation tracks
N = 50 # number of genes (rows) to show in heatmap
column_clust_method = 'average'
row_clust_method = 'average'
low_color = "blue"
high_color = "red"

# Data to plot
top <- read.csv(sprintf("data/up_%s_zscores.csv", data), header=T, row.names=1, check.names = F)
bottom <- read.csv(sprintf("data/dn_%s_zscores.csv", data), header=T, row.names=1, check.names = F)

metadata <- read.csv("data/metadata.csv", header=T, row.names=1, check.names = F)

# Align sample order
col_order <- rownames(metadata)[rownames(metadata) %in% names(top)]
top <- top %>% select(all_of(col_order))
bottom <- bottom %>% select(all_of(col_order))
metadata <- metadata[col_order, ]

# Heatmap annotation dataframe
anno_df <- metadata[tracks]
ha = HeatmapAnnotation(df = anno_df,
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_gp = gpar(fontsize = 8))

# Rows to plot
row_vars <- apply(top, 1, sd)
sorted_indices <- order(row_vars, decreasing = TRUE)
top_rows <- row.names(top)[sorted_indices[1:N]]
row_vars <- apply(bottom, 1, sd)
sorted_indices <- order(row_vars, decreasing = TRUE)
bottom_rows <- row.names(bottom)[sorted_indices[1:N]]

mat1 = data.matrix(top[rownames(top) %in% top_rows, ])
mat2 = data.matrix(bottom[rownames(bottom) %in% bottom_rows, ])

if (data == "kea") {
  split_df = c(rep('Enriched in up-phosphosites', length(top_rows)), rep('Enriched in down-phosphosites', length(bottom_rows)))
} else {
  split_df = c(rep('Enriched in up-genes', length(top_rows)), rep('Enriched in down-genes', length(bottom_rows)))
}

colors = colorRamp2(seq(min(min(mat1), min(mat2)), max(max(mat1), max(mat2)), length = 3), c(low_color, "#EEEEEE", high_color))

hm = Heatmap(rbind(mat1, mat2), 
                 name = "Enrichment score",
                 clustering_method_columns = column_clust_method,
                 clustering_method_rows = row_clust_method,
                 top_annotation = ha, 
                 row_split = split_df,
                 show_column_names = F,
                 column_title = NULL,
                 col = colors,
                 row_names_gp = grid::gpar(fontsize = 8),
                 heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop", legend_width = unit(5, "cm")),
                 row_gap = unit(3.2, "mm"))
png(file = sprintf("results/%s_heatmap.png", data), width = 225, height = 205, units='mm', res = 300)
draw(hm, heatmap_legend_side = "bottom")
dev.off()
