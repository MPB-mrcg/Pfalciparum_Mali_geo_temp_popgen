
rm(list = ls())
library(tidyverse)
library(ggrepel)

Bamako <- readxl::read_xlsx("Bamako.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()

Bougoula_Hameau <- readxl::read_xlsx("Bougoula-Hameau.filtered_fuli.xlsx") %>% 
    separate(Pos, c("Start", "End")) %>% 
    type_convert()

## COMBINE DATAFRAMES
df <- left_join(Bamako, Bougoula_Hameau, by = c("Gene", "Chrom", "Start", "End"))

## IDENTIFY TOP GENES
top.genes.Bam <- Bamako %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

top.genes.BH <- Bougoula_Hameau %>% 
    arrange(desc(TajimasD)) %>%
    head(n = 10) %>% 
    select(c(1:5))

## COMMON GENES BETWEEN POP
outliers <- inner_join(top.genes.Bam, top.genes.BH, by = c("Gene", "Chrom", "Start", "End"))

## CALCULATE CORRELATION COEF
cor.val <- round(cor(df$TajimasD.x, df$TajimasD.y,  method = "pearson", use = "pairwise.complete.obs"), 2)

## LABEL
cor.label <- paste0("Correlation: ", cor.val)

## CREATE PLOT
p1 <- ggplot(df, aes(x = TajimasD.x, y = TajimasD.y)) + 
    geom_point() + 
    geom_abline(color = "red", slope = 1) + 
    xlim(-2, 5) + ylim(-3, 4) +
    annotate(x = -1.3, y = 3.5,  geom = "text", label = cor.label, size = 4) +
    labs( x = "Bamako", y = "Bougoula_Hameau") +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    geom_point(data = outliers, 
               aes(x = TajimasD.x, y = TajimasD.y), 
               color ='red',
               size = 2) +
    geom_label_repel(data = outliers, aes(x = TajimasD.x, y = TajimasD.y, 
                    label = Gene, fill = factor(Gene)), force = 5,
                    size = 2.5, force_pull = 5, max.overlaps = Inf)

p2 <- ggplot(df, aes(x = TajimasD.x, y = TajimasD.y)) + 
    geom_point() + 
    geom_abline(color = "red", slope = 1) + 
    xlim(-2, 5) + ylim(-3, 4) +
    annotate(x = -1.3, y = 3.5,  geom = "text", label = cor.label, size = 4) +
    labs( x = "Bamako", y = "Bougoula_Hameau") +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    geom_point(data = outliers, 
               aes(x = TajimasD.x, y = TajimasD.y), 
               color ='red',
               size = 2) +
    geom_label_repel(data = outliers, aes(x = TajimasD.x, y = TajimasD.y, 
                                          label = Gene, fill = factor(Gene)), force = 5,
                     size = 2.5, force_pull = 5, max.overlaps = Inf)

p1 + p2

ggsave(filename = "correlation.jpeg", width = 20, height = 10)
