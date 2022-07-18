
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()
extrafont::loadfonts(device = "all", quiet = T)
extrafont::fonts()
library(tidyverse)

metadata <- read_tsv("raw_data/Mali_Matadata.txt")

pca <- read_tsv("results/tables/Mali_PCA.txt")
mds <- read_tsv("results/tables/Mali_MDS.txt")

pca <- pca %>% inner_join(metadata, by = c("Taxa" = "Sample"))
mds <- mds %>% inner_join(metadata, by = c("Taxa" = "Sample"))

legend.color <- c("#000000", "#eb1f10", "#30ff08", "#05f0fc", "#0509fc", 
                  "#FEFD03", "#9003FE") # , "#900C3F"

ggplot(pca, aes(PC1, PC2)) +
    geom_point(colour="Black", shape=21, size = 5, aes(fill = factor(Location))) +
    scale_y_continuous(breaks=seq(-8, 5, by = 1))+
    scale_x_continuous(breaks=seq(-15, 7, by = 1)) +
    xlab("PC2") + ylab("PC1") +
    scale_fill_manual("Location", values= legend.color) +
    theme(text = element_text(family = "Arial", size = 15, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))



