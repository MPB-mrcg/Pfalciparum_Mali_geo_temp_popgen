
library(fs)
library(tidyverse)
library(ComplexHeatmap)
library(magick)
library(reshape2)
library(circlize)
library(gplots)



Path <- "results/tables/Fst/Geographical"

get_country <- function(filename){
    str_match(filename, paste0(Path, "/(.*)_wcFst.txt"))[,2]
}

# files <- dir_ls("results/tables/Fst/Geographical")
# 
# df_list <- map(files, read_tsv, col_names = F)
# names(df_list)
# 
# bind_rows(df_list, .id = 'Population')
# 
# bind_rows(df_list, .id = 'Population') %>% 
#     mutate(Population = str_match(Population, paste0(Path, "/(.*)_wcFst.txt"))[,2])

data <- dir_ls("results/tables/Fst/Geographical") %>% 
    map(read_tsv, col_names = F) %>% 
    bind_rows(.id = 'Population') %>% 
    mutate(Population = get_country(Population)) %>% 
    mutate(chr = recode(X1, 
                        "Pf3D7_01_v3" = 1, "Pf3D7_02_v3" = 2, "Pf3D7_03_v3" = 3, "Pf3D7_04_v3" = 4, 
                        "Pf3D7_05_v3" = 5, "Pf3D7_06_v3" = 6, "Pf3D7_07_v3" = 7, "Pf3D7_08_v3" = 8, 
                        "Pf3D7_09_v3" = 9, "Pf3D7_10_v3" = 10, "Pf3D7_11_v3" = 11, "Pf3D7_12_v3" = 12, 
                        "Pf3D7_13_v3" = 13, "Pf3D7_14_v3" = 14))

data.format <- as_tibble(reshape2::dcast(data, X1+chr+X2 ~ Population, value.var="X5"))

data.format1 <- as.matrix(t(data.format[-c(1:3)]))

#================================
## Annotated Heatmap with boxplot
#================================
myBreaks <- seq(0.0 , 0.5,  length.out =10)
myCol <- colorRampPalette(c("lightgreen","red"))(10)

ha.box <- HeatmapAnnotation(FST = anno_boxplot(data.format1, height = unit(1, "cm"), 
                                               pch = 20, size = unit(1, "mm"), axis = TRUE,
                                               gp = gpar(fill = "#CC0000",fontface = "bold"), 
                                               axis_param = list(gp=gpar(fontsize=7,fontface = "bold"))))

Heatmap(data.format1, name = "FST", column_split=data.format$chr, 
        col=colorRamp2(myBreaks, myCol), heatmap_width = unit(25, "cm"), #col = bluered(64),
        heatmap_height = unit(20, "cm"), cluster_columns = F,
        show_column_dend = F, cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 9, fontface = "bold"),
        show_column_names =F, row_names_g = gpar(fontsize = 7,fontface = "bold"),
        row_names_rot = 0, column_title_rot = 0, top_annotation = ha.box,
        heatmap_legend_param = list(title_gp=gpar(fontsize=7, fontface="bold"),legend_width=unit(8,"cm"),  
                                    legend_height=unit(10,"cm"), labels_gp=gpar(fontsize=8, fontface="bold")))

# top <- data %>% 
#     group_by(Population) %>% 
#     summarise(Outlier = quantile(X5, c(0.99)), q = 0.99) %>%
#     left_join(data, ., by="Population") %>% 
#     filter(X5 > Outlier)

as_tibble(reshape2::dcast(top, X1+chr+X2 ~ Population, value.var="X5")) %>% 
    write.table("results/tables/Fst/geographical.xlsx", 
                col.names = T, row.names = F, quote = F, sep = "\t")
