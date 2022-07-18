
library(patchwork)
library(tidyverse)

aoua_data <- read_tsv("data/reference/Genotypes.txt", col_names = F) %>% 
    dplyr::select(c(1:4))

fst <- read_csv("results/tables/Fst/maligeopops_ebfst.csv") %>% 
    mutate(chrom = recode(CHR, 
                          "1" = "Pf3D7_01_v3", "2" = "Pf3D7_02_v3", "3" = "Pf3D7_03_v3", "4" = "Pf3D7_04_v3", 
                          "5" = "Pf3D7_05_v3", "6" = "Pf3D7_06_v3", "7" = "Pf3D7_07_v3", "8" = "Pf3D7_08_v3", 
                          "9" = "Pf3D7_09_v3", "10" = "Pf3D7_10_v3", "11" = "Pf3D7_11_v3", "12" = "Pf3D7_12_v3", 
                          "13" = "Pf3D7_13_v3", "14" = "Pf3D7_14_v3"))

fst <- fst %>% 
    inner_join(aoua_data, by = c("chrom"="X1", "POS"="X2")) %>% 
    relocate(chrom, .before = SNP) %>% 
    dplyr::select(-SNP)


fst.format <- fst %>% 
    
    # Compute chromosome size
    group_by(chrom) %>% 
    summarise(chr_len = max(POS)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(fst, ., by=c("chrom"="chrom")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chrom, POS) %>%
    mutate( BPcum = POS + tot) %>% 
    ungroup()

axis <- fst.format %>%
    group_by(CHR) %>%
    summarize(center = mean(BPcum))

limits <- c(min(fst.format$Fst), max(fst.format$Fst)+0.1)

# Ready to make the plot using ggplot2:
ebfst.plot <- ggplot(fst.format, aes(x = BPcum, y = Fst, color = as_factor(CHR))) +
    # Show all points
    geom_point(size = .3) +
    geom_hline(yintercept = quantile(fst.format$Fst, 0.99), color = "darkred", size = 0.3) +
    
    # custom X axis:
    scale_x_continuous(label = axis$CHR, breaks = axis$center) +
    scale_y_continuous(expand = c(0,0), limits = limits) +
    scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Fst") +
    ggtitle('Geographical') +
    
    # Custom the theme:
    theme_minimal() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 6),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(face = "bold.italic", size = 8),
        axis.text.y = element_text(colour = "black", face = "bold", size = 6))


Fst <- read_tsv("results/tables/Fst/Heterogeneity/wcFst_07_13.txt", col_names = F)
Fst1 <- read_tsv("results/tables/Fst/Heterogeneity/wcFst_07_1517.txt", col_names = F)
Fst2 <- read_tsv("results/tables/Fst/Heterogeneity/wcFst_13_1517.txt", col_names = F)

Fst.format <- Fst %>% 
    
    # Compute chromosome size
    group_by(X1) %>% 
    summarise(chr_len = max(X2)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(Fst, ., by=c("X1"="X1")) %>%
    
    # Add a cumulative position of each SNP
    arrange(X1, X2) %>%
    mutate( BPcum = X2 + tot) %>% 
    ungroup()


Fst.format <- Fst.format %>% 
    mutate(Chr = as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", X1)))) %>%
    mutate(X5 = ifelse(X5 <0, 0, X5))

axis_set <- Fst.format %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylimits <- c(min(Fst.format$X5), max(Fst.format$X5)+0.1)

# Ready to make the plot using ggplot2:
plot1 <- ggplot(Fst.format, aes(x = BPcum, y = X5, color = as_factor(Chr))) +
    # Show all points
    geom_point(size = .3) +
    geom_hline(yintercept = quantile(Fst.format$X5, 0.99), color = "darkred", size = 0.3) +
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis_set$Chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Fst") +
    ggtitle('Faladje 2007-2013') +

    # Custom the theme:
    theme_minimal() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 6), #, hjust = 0, vjust = 0
        # plot.title.position = "plot",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(face = "bold.italic", size = 8),
        axis.text.y = element_text(colour = "black", face = "bold", size = 6))

#=============================================================
Fst1.format <- Fst1 %>% 
    
    # Compute chromosome size
    group_by(X1) %>% 
    summarise(chr_len = max(X2)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(Fst1, ., by=c("X1"="X1")) %>%
    
    # Add a cumulative position of each SNP
    arrange(X1, X2) %>%
    mutate( BPcum = X2 + tot) %>% 
    ungroup()


Fst1.format <- Fst1.format %>% 
    mutate(Chr = as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", X1)))) %>%
    mutate(X5 = ifelse(X5 <0, 0, X5))

axis_set1 <- Fst1.format %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))


# Ready to make the plot using ggplot2:
plot2 <- ggplot(Fst1.format, aes(x = BPcum, y = X5, color = as_factor(Chr))) +
    # Show all points
    geom_point(size = .3) +
    geom_hline(yintercept = quantile(Fst1.format$X5, 0.99), color = "darkred", size = 0.3) +
    
    # custom X axis:
    scale_x_continuous(label = axis_set1$Chr, breaks = axis_set1$center) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis_set1$Chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Fst") +
    ggtitle('Faladje 2007-2015_17') +
    
    # Custom the theme:
    theme_minimal() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 6),
        axis.line = element_line(colour = "black"),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(colour = "black", face = "bold", size = 6))

#=======================================================
# Fst.format2 <- Fst2 %>% 
#     
#     # Compute chromosome size
#     group_by(X1) %>% 
#     summarise(chr_len = max(X2)) %>% 
#     
#     # Calculate cumulative position of each chromosome
#     mutate(tot = cumsum(chr_len) - chr_len) %>%
#     dplyr::select(-chr_len) %>%
#     
#     # Add this info to the initial dataset
#     left_join(Fst2, ., by=c("X1"="X1")) %>%
#     
#     # Add a cumulative position of each SNP
#     arrange(X1, X2) %>%
#     mutate( BPcum = X2 + tot) %>% 
#     ungroup()
# 
# 
# Fst.format2 <- Fst.format2 %>% 
#     mutate(Chr = as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", X1)))) %>%
#     mutate(X5 = ifelse(X5 <0, 0, X5))
# 
# axis_set2 <- Fst.format2 %>%
#     group_by(Chr) %>%
#     summarize(center = mean(BPcum))
# 
# # Ready to make the plot using ggplot2:
# plot3 <- ggplot(Fst.format2, aes(x = BPcum, y = X5, color = as_factor(Chr))) +
#     # Show all points
#     geom_point(size = .3) +
#     geom_hline(yintercept = quantile(Fst.format2$X5, 0.99), color = "darkred", size = 0.3) +
#     
#     # custom X axis:
#     scale_x_continuous(label = axis_set2$Chr, breaks = axis_set2$center) +
#     scale_y_continuous(expand = c(0,0), limits = ylimits) +
#     scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis_set2$Chr)))) +
#     scale_size_continuous(range = c(0.5,3)) +
#     labs(x = "Chromosomes", y = "Fst") +
#     ggtitle('Faladje 2013-2015_17') +
#     
#     # Custom the theme:
#     theme_minimal() +
#     theme( 
#         legend.position = "none",
#         plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 6),
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(face = "bold.italic", size = 8),
#         axis.text = element_text(colour = "black", face = "bold", size = 6))
# 
# 
# plot1 / plot2 / plot3
# ggsave("results/figures/Fst_heterogeneity.pdf", width = 180, height = 150, units = "mm", dpi = 600)

ebfst.plot / plot1 / plot2

ggsave("results/figures/Fst.pdf", 
       width = 180, height = 150, units = "mm", dpi = 600)


#=======================
# SUMMARY TABLES
#=======================

outliers1 <- fst %>% 
    filter(Fst > quantile(fst.format$Fst, 0.99)) %>% 
    dplyr::select(-CHR,-X3,-X4)

outliers2 <-  Fst.format %>% 
    filter(X5 > quantile(Fst.format$X5, 0.99)) %>% 
    dplyr::select(-Chr,-X3,-X4,-tot , -BPcum,   -Chr)

names(outliers2) <- c("chrom","POS","Fst")

outliers3 <- Fst1.format %>% 
    filter(X5 > quantile(Fst1.format$X5, 0.99)) %>% 
    dplyr::select(-Chr,-X3,-X4,-tot , -BPcum,   -Chr)
names(outliers3) <- c("chrom","POS","Fst")

#put all data frames into list
df_list <- list(outliers1, outliers2, outliers3)

# Common genes across all indices
top_outliers <- df_list %>% reduce(full_join, by=c('chrom', "POS")) %>% 
    arrange(chrom, POS)

names(top_outliers)[3:5] <- c("Geo", "Faladje_07_1517", "Faladje_07_13")

write.table(top_outliers, "results/tables/Fst/Heterogeneity.xlsx", 
            col.names = T, row.names = F, quote = F, sep = "\t")
