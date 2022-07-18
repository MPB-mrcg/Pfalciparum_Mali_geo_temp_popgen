

rm(list = ls())
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

setwd("~/Documents/Aoua_PhD/data/data_filt/Tajima_results/")

Bamako <- readxl::read_xlsx("Bamako.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()

Bougoula_Hameau <- readxl::read_xlsx("Bougoula-Hameau.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()
Dangassa <- readxl::read_xlsx("Dangassa.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()

Faladje <- readxl::read_xlsx("Faladje.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()
Kenieroba <- readxl::read_xlsx("Kenieroba.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()
Kolle <- readxl::read_xlsx("Kolle.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()
Nioro <- readxl::read_xlsx("nioro.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()

## COMBINE DATAFRAMES
df1 <- left_join(Bamako, Bougoula_Hameau, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Bougoula_Hameau", sep = "_"))

df2 <- left_join(Bamako, Dangassa, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Dangassa", sep = "_"))

df3 <- left_join(Bamako, Faladje, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Faladje", sep = "_"))

df4 <- left_join(Bamako, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Kenieroba", sep = "_"))

df5 <- left_join(Bamako, Kolle, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Kolle", sep = "_"))

df6 <- left_join(Bamako, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bamako", "Nioro", sep = "_"))

df7 <- left_join(Bougoula_Hameau, Dangassa, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Dangassa", sep = "_"))

df8 <- left_join(Bougoula_Hameau,Faladje, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Faladje", sep = "_"))

df9 <- left_join(Bougoula_Hameau, Kenieroba, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Kenieroba", sep = "_"))

df10 <- left_join(Bougoula_Hameau, Kolle, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Kolle", sep = "_"))

df11 <- left_join(Bougoula_Hameau, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Nioro", sep = "_"))

df12 <- left_join(Dangassa, Faladje, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Dangassa", "Faladje", sep = "_"))

df13 <- left_join(Dangassa, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Dangassa", "Kenieroba", sep = "_"))

df14 <- left_join(Dangassa, Kolle, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Dangassa", "Kolle", sep = "_"))

df15 <- left_join(Dangassa, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Dangassa", "Nioro", sep = "_"))

df16 <- left_join(Faladje, Kenieroba, by = c("Gene", "Chrom", "Start", "End")) %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Faladje", "Kenieroba", sep = "_"))

df17 <- left_join(Faladje, Kolle, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Faladje", "Kolle", sep = "_"))

df18 <- left_join(Faladje, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Faladje", "Nioro", sep = "_"))

df19 <- left_join(Kenieroba, Kolle, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Kenieroba", "Kolle", sep = "_"))

df20 <- left_join(Kenieroba, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Kenieroba", "Nioro", sep = "_"))

df21 <- left_join(Kolle, Nioro, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    select(Gene, Chrom, Start, End, TajimasD.x, TajimasD.y) %>% 
    add_column(Pairwise = paste("Kolle", "Nioro", sep = "_"))



## IDENTIFY TOP GENES
top.genes.Bam <- Bamako %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.BH <- Bougoula_Hameau %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.DG <- Dangassa %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.FL <- Faladje %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.KEN <- Kenieroba %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.KOL <- Kolle %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.NR <- Nioro %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

## COMMON GENES BETWEEN POP
outliers1 <- inner_join(top.genes.Bam, top.genes.BH, by = c("Gene", "Chrom", "Start", "End")) %>% 
    add_column(Pairwise = paste("Bamako", "Bougoula_Hameau", sep = "_"))

outliers2 <- inner_join(top.genes.Bam, top.genes.DG, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bamako", "Dangassa", sep = "_"))

outliers3 <- inner_join(top.genes.Bam, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bamako", "Faladje", sep = "_"))

outliers4 <- inner_join(top.genes.Bam, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bamako", "Kenieroba", sep = "_"))

outliers5 <- inner_join(top.genes.Bam, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bamako", "Kolle", sep = "_"))

outliers6 <- inner_join(top.genes.Bam, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bamako", "Nioro", sep = "_"))

outliers7 <- inner_join(top.genes.BH, top.genes.DG, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Dangassa", sep = "_"))

outliers8 <- inner_join(top.genes.BH, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Faladje", sep = "_"))

outliers9 <- inner_join(top.genes.BH, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Kenieroba", sep = "_"))

outliers10 <- inner_join(top.genes.BH, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Kolle", sep = "_"))

outliers11 <- inner_join(top.genes.BH, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Bougoula_Hameau", "Nioro", sep = "_"))

outliers12 <- inner_join(top.genes.DG, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Dangassa", "Faladje", sep = "_"))

outliers13 <- inner_join(top.genes.DG, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Dangassa", "Kenieroba", sep = "_"))

outliers14 <- inner_join(top.genes.DG, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Dangassa", "Kolle", sep = "_"))

outliers15 <- inner_join(top.genes.DG, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Dangassa", "Nioro", sep = "_"))

outliers16 <- inner_join(top.genes.FL, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Faladje", "Kenieroba", sep = "_"))

outliers17 <- inner_join(top.genes.FL, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Faladje", "Kolle", sep = "_"))

outliers18 <- inner_join(top.genes.FL, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Faladje", "Nioro", sep = "_"))

outliers19 <- inner_join(top.genes.KEN, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Kenieroba", "Kolle", sep = "_"))

outliers20 <- inner_join(top.genes.KEN, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))  %>% 
    add_column(Pairwise = paste("Kenieroba", "Nioro", sep = "_"))

outliers21 <- inner_join(top.genes.KOL, top.genes.NR, by = c("Gene", "Chrom", "Start", "End")) %>% 
    add_column(Pairwise = paste("Kolle", "Nioro", sep = "_"))


## COMBINE DATAFRAMES
df <- rbind.data.frame(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, 
                       df13, df14, df15, df16, df17, df18, df19, df20, df21)

## COMBINE OUTLIERS
outliers <- rbind.data.frame(outliers1, outliers2, outliers3, outliers4, outliers5, outliers6, 
                             outliers7, outliers8, outliers9, outliers10, outliers11, outliers12, 
                             outliers13, outliers14, outliers15, outliers16, outliers17, 
                             outliers18, outliers19, outliers20, outliers21)

## CREATE PLOT

ggplot(df, aes(x = TajimasD.x, y = TajimasD.y)) + 
    geom_point() + 
    geom_abline(color = "red", slope = 1) + 
    xlim(-2, 5) + ylim(-3, 4) +
    
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    geom_point(data = outliers, 
               aes(x = TajimasD.x, y = TajimasD.y), 
               color ='red') +
    geom_label_repel(data = outliers, aes(x = TajimasD.x, y = TajimasD.y, fill = Gene, label = Gene), 
                     force = 5, size = 1.5, force_pull = 5, max.overlaps = Inf,
                     segment.colour = "black", fontface = 'bold', color = 'black', 
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines")) +
    facet_wrap(~Pairwise)

## SAVE PLOT
ggsave(filename = "correlation_plots.jpeg", width = 20, height = 10)
