
rm(list=ls())
library(tidyverse)
library(readxl)
library(reshape2)
library(ggtext)

## Load Gene Name file
Gene_name <- "../PF_GeneName.tsv"
Gene_name <- read_tsv(Gene_name, show_col_types = FALSE)

## List files
tjd_files <- list.files(path = ".", pattern = ".xlsx", full.names = TRUE)

names(tjd_files) <- str_replace(string = tjd_files,
                                pattern = paste0(".", "/(.*).filtered_fuli.xlsx"),
                                replacement = "\\1")

#===================================
## Combine files by adding filenames
#===================================
tjd <- map_dfr(.x = tjd_files, .f = read_xlsx, .id = "Population")

tajima <- tjd %>% 
    separate(Pos, c("Start","End"), sep = "-") %>% 
    type.convert() 

tajima <- tajima %>% 
    
    # Compute chromosome size
    group_by(Chrom) %>% 
    summarise(chr_len = max(Start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(tajima, ., by=c("Chrom"="Chrom")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chrom, Start) %>%
    mutate( BPcum = Start + tot)

axis_set <- tajima %>%
    group_by(Chrom) %>%
    summarize(center = mean(BPcum))

# ylim <- abs(floor(min(taj$TajimasD))) + 1
ylimits <- c(floor(min(tajima$TajimasD)), abs(floor(min(tajima$TajimasD))) + 2)

# Ready to make the plot using ggplot2:
ggplot(tajima, aes(x=BPcum, y=TajimasD, 
               color = as_factor(Chrom))) +
    # Show all points
    geom_point(alpha = 0.75, size = 1.2) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chrom, breaks = axis_set$center, expand = c(0.02,0.05)) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +    # remove space between plot area and x axis
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$Chrom)))) +
    scale_size_continuous(range = c(0.5,3)) +

    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold.italic"),
        plot.subtitle = element_markdown(face = "bold", color = "gray60", size = 11, family = "Playfair"),
        axis.text.x = element_text(),
        axis.title.x = element_text(vjust = 0),
        axis.title.y = element_text(vjust = 2),
        axis.title = element_text(face = "bold.italic")
    ) + 
    facet_wrap(~Population, ncol = 4)
