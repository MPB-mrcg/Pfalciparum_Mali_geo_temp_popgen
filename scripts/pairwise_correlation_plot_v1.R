
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
df1 <- left_join(Bamako, Bougoula_Hameau, by = c("Gene", "Chrom", "Start", "End"))
df2 <- left_join(Bamako, Dangassa, by = c("Gene", "Chrom", "Start", "End"))
df3 <- left_join(Bamako, Faladje, by = c("Gene", "Chrom", "Start", "End"))
df4 <- left_join(Bamako, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))
df5 <- left_join(Bamako, Kolle, by = c("Gene", "Chrom", "Start", "End"))
df6 <- left_join(Bamako, Nioro, by = c("Gene", "Chrom", "Start", "End"))
df7 <- left_join(Bougoula_Hameau, Dangassa, by = c("Gene", "Chrom", "Start", "End"))
df8 <- left_join(Bougoula_Hameau,Faladje, by = c("Gene", "Chrom", "Start", "End"))
df9 <- left_join(Bougoula_Hameau, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))
df10 <- left_join(Bougoula_Hameau, Kolle, by = c("Gene", "Chrom", "Start", "End"))
df11 <- left_join(Bougoula_Hameau, Nioro, by = c("Gene", "Chrom", "Start", "End"))
df12 <- left_join(Dangassa, Faladje, by = c("Gene", "Chrom", "Start", "End"))
df13 <- left_join(Dangassa, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))
df14 <- left_join(Dangassa, Kolle, by = c("Gene", "Chrom", "Start", "End"))
df15 <- left_join(Dangassa, Nioro, by = c("Gene", "Chrom", "Start", "End"))
df16 <- left_join(Faladje, Kenieroba, by = c("Gene", "Chrom", "Start", "End"))
df17 <- left_join(Faladje, Kolle, by = c("Gene", "Chrom", "Start", "End"))
df18 <- left_join(Faladje, Nioro, by = c("Gene", "Chrom", "Start", "End"))
df19 <- left_join(Kenieroba, Kolle, by = c("Gene", "Chrom", "Start", "End"))
df20 <- left_join(Kenieroba, Nioro, by = c("Gene", "Chrom", "Start", "End"))
df21 <- left_join(Kolle, Nioro, by = c("Gene", "Chrom", "Start", "End"))



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
outliers1 <- inner_join(top.genes.Bam, top.genes.BH, by = c("Gene", "Chrom", "Start", "End"))
outliers2 <- inner_join(top.genes.Bam, top.genes.DG, by = c("Gene", "Chrom", "Start", "End"))
outliers3 <- inner_join(top.genes.Bam, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))
outliers4 <- inner_join(top.genes.Bam, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))
outliers5 <- inner_join(top.genes.Bam, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))
outliers6 <- inner_join(top.genes.Bam, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))
outliers7 <- inner_join(top.genes.BH, top.genes.DG, by = c("Gene", "Chrom", "Start", "End"))
outliers8 <- inner_join(top.genes.BH, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))
outliers9 <- inner_join(top.genes.BH, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))
outliers10 <- inner_join(top.genes.BH, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))
outliers11 <- inner_join(top.genes.BH, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))
outliers12 <- inner_join(top.genes.DG, top.genes.FL, by = c("Gene", "Chrom", "Start", "End"))
outliers13 <- inner_join(top.genes.DG, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))
outliers14 <- inner_join(top.genes.DG, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))
outliers15 <- inner_join(top.genes.DG, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))
outliers16 <- inner_join(top.genes.FL, top.genes.KEN, by = c("Gene", "Chrom", "Start", "End"))
outliers17 <- inner_join(top.genes.FL, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))
outliers18 <- inner_join(top.genes.FL, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))
outliers19 <- inner_join(top.genes.KEN, top.genes.KOL, by = c("Gene", "Chrom", "Start", "End"))
outliers20 <- inner_join(top.genes.KEN, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))
outliers21 <- inner_join(top.genes.KOL, top.genes.NR, by = c("Gene", "Chrom", "Start", "End"))

## CALCULATE CORRELATION COEF
cor.val1 <- round(cor(df1$TajimasD.x, df1$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val2 <- round(cor(df2$TajimasD.x, df2$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val3 <- round(cor(df3$TajimasD.x, df3$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val4 <- round(cor(df4$TajimasD.x, df4$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val5 <- round(cor(df5$TajimasD.x, df5$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val6 <- round(cor(df6$TajimasD.x, df6$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val7 <- round(cor(df7$TajimasD.x, df7$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val8 <- round(cor(df8$TajimasD.x, df8$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val9 <- round(cor(df9$TajimasD.x, df9$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val10 <- round(cor(df10$TajimasD.x, df10$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val11 <- round(cor(df11$TajimasD.x, df11$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val12 <- round(cor(df12$TajimasD.x, df12$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val13 <- round(cor(df13$TajimasD.x, df13$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val14 <- round(cor(df14$TajimasD.x, df14$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val15 <- round(cor(df15$TajimasD.x, df15$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val16 <- round(cor(df16$TajimasD.x, df16$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val17 <- round(cor(df17$TajimasD.x, df17$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val18 <- round(cor(df18$TajimasD.x, df18$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val19 <- round(cor(df19$TajimasD.x, df19$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val20 <- round(cor(df20$TajimasD.x, df20$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)
cor.val21 <- round(cor(df21$TajimasD.x, df21$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)

## LABEL
cor.label1 <- paste0("Correlation: ", cor.val1)
cor.label2 <- paste0("Correlation: ", cor.val2)
cor.label3 <- paste0("Correlation: ", cor.val3)
cor.label4 <- paste0("Correlation: ", cor.val4)
cor.label5 <- paste0("Correlation: ", cor.val5)
cor.label6 <- paste0("Correlation: ", cor.val6)
cor.label7 <- paste0("Correlation: ", cor.val7)
cor.label8 <- paste0("Correlation: ", cor.val8)
cor.label9 <- paste0("Correlation: ", cor.val9)
cor.label10 <- paste0("Correlation: ", cor.val10)
cor.label11 <- paste0("Correlation: ", cor.val11)
cor.label12 <- paste0("Correlation: ", cor.val12)
cor.label13 <- paste0("Correlation: ", cor.val13)
cor.label14 <- paste0("Correlation: ", cor.val14)
cor.label15 <- paste0("Correlation: ", cor.val15)
cor.label16 <- paste0("Correlation: ", cor.val16)
cor.label17 <- paste0("Correlation: ", cor.val17)
cor.label18 <- paste0("Correlation: ", cor.val18)
cor.label19 <- paste0("Correlation: ", cor.val19)
cor.label20 <- paste0("Correlation: ", cor.val20)
cor.label21 <- paste0("Correlation: ", cor.val21)

Mypalette = c("gold4 ", "red3", "darkmagenta ", "blue2") #,"darkorange", "darkgreen", "yellow", "pink")

## CREATE PLOT
 plot1 <- ggplot(df1, aes(x = TajimasD.x, y = TajimasD.y)) + 
    geom_point() + 
    geom_abline(color = "red", slope = 1) + 
    xlim(-2, 5) + ylim(-3, 4) +
    annotate(x = -1.3, y = 3.5,  geom = "text", label = cor.label1, size = 4) +
    labs( x = "Bamako", y = "Bougoula_Hameau") +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    geom_point(data = outliers1, 
               aes(x = TajimasD.x, y = TajimasD.y), 
               color ='red',
               size = 2) +
    geom_label_repel(data = outliers1, 
                     aes(x = TajimasD.x, y = TajimasD.y, fill = Gene, label = Gene), 
                    segment.colour = "black", fontface = 'bold', color = 'black', 
                    box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines")) +
     scale_fill_manual(values = setNames(RColorBrewer::brewer.pal(nrow(outliers1), 'Dark2'), levels(outliers1$Gene))) +
     scale_color_manual(values = setNames(RColorBrewer::brewer.pal(nrow(outliers1), 'Dark2'), levels(outliers1$Gene)))

plot2 <- ggplot(df2, aes(x = TajimasD.x, y = TajimasD.y)) + 
  geom_point() + 
  geom_abline(color = "red", slope = 1) + 
  xlim(-2, 5) + ylim(-3, 4) +
  annotate(x = -1.3, y = 3.5,  geom = "text", label = cor.label2, size = 4) +
  labs( x = "Bamako", y = "Dangassa") +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=8.5, face="bold"),
        axis.title=element_text(size=10,face="bold")) +
  geom_point(data = outliers2, 
             aes(x = TajimasD.x, y = TajimasD.y), 
             color ='red',
             size = 2) +
    geom_label_repel(data = outliers2, 
                     aes(x = TajimasD.x, y = TajimasD.y, fill = Gene, label = Gene), 
                     segment.colour = "black", fontface = 'bold', color = 'black', 
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines")) +
    scale_fill_manual(values = setNames(RColorBrewer::brewer.pal(nrow(outliers2), 'Dark2'), levels(outliers2$Gene))) +
    scale_color_manual(values = setNames(RColorBrewer::brewer.pal(nrow(outliers2), 'Dark2'), levels(outliers2$Gene)))


#ggarrange(plot1, plot2, heights = 0.25, widths = 0.5)
plot1 + plot2 
ggsave(filename = "correlation_plots.jpeg", width = 20, height = 10)
