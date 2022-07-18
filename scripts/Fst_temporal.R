
Fst <- read_tsv("results/tables/Fst/Temporal/2007_2013_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2007_2013")
Fst1 <- read_tsv("results/tables/Fst/Temporal/2007_2015_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2007_2015")
Fst2 <- read_tsv("results/tables/Fst/Temporal/2007_2016_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2007_2016")
Fst3 <- read_tsv("results/tables/Fst/Temporal/2007_2017_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2007_2017")
Fst4 <- read_tsv("results/tables/Fst/Temporal/2013_2015_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2013_2015")
Fst5 <- read_tsv("results/tables/Fst/Temporal/2013_2016_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2013_2016")
Fst6 <- read_tsv("results/tables/Fst/Temporal/2013_2017_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2013_2017")
Fst7 <- read_tsv("results/tables/Fst/Temporal/2015_2016_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2015_2016")
Fst8 <- read_tsv("results/tables/Fst/Temporal/2015_2017_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2015_2017")
Fst9 <- read_tsv("results/tables/Fst/Temporal/2016_2017_wcFst.txt", col_names = F) %>% 
    add_column(Population = "2016_2017")

combined <- rbind.data.frame(Fst, Fst1, Fst2, Fst3, Fst4, 
                 Fst5, Fst6, Fst7, Fst8, Fst9)

combined.format <- combined %>% 
    
    # Compute chromosome size
    group_by(X1) %>% 
    summarise(chr_len = max(X2)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(combined, ., by=c("X1"="X1")) %>%
    
    # Add a cumulative position of each SNP
    arrange(X1, X2) %>%
    mutate( BPcum = X2 + tot) %>% 
    ungroup()


combined.format <- combined.format %>% 
    mutate(Chr = as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", X1)))) %>%
    mutate(X5 = ifelse(X5 <0, 0, X5))

axis_set <- combined.format %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylimits <- c(min(combined.format$X5), max(combined.format$X5))

# Ready to make the plot using ggplot2:
combined.plot <- ggplot(combined.format, aes(x = BPcum, y = X5, color = factor(Chr), size = X5)) +
    # Show all points
    geom_point(size = .15) + #, colour="Black"
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0)) + # , limits = ylimits
    scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis_set$Chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Fst") +

    # Custom the theme:
    theme_light() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 5), #, hjust = 0, vjust = 0
        # # panel.grid.major = element_line(size = 0.5),
        # panel.grid.major.y = element_line(size = 1.5), 
        # panel.grid.minor.y = element_line(size = 0.5),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(face = "bold.italic", size = 8),
        axis.text = element_text(colour = "black", face = "bold", size = 6),
        strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
        # strip.background = element_rect(color="black", fill="#92C5DE", size=1.5, linetype="solid")) + #"#FC4E07"
    facet_wrap(~factor(Population), scales = "free_y", ncol = 2)

print(combined.plot)

ggsave("results/figures/Fst_temporal.pdf", width = 180, height = 150, units = "mm",)
