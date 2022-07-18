rm(list=ls())

library(tidyverse)
library(rlang)
library(ggrepel)
library(stringr)
theme_set(theme_bw())

# load New Gene Name file
Gene_name <- "PF_GeneName.tsv" # File with Chrom Start End and GeneName
Gene_name <- read_tsv(Gene_name)

#************************************
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

combined.tajima <- combined.tajima %>% 
    inner_join(Gene_name) %>% 
    select(-c(10, 12, 13))

#================================================
# Identifying outliers genes for each population
#================================================
top.genes.pop1 <- pop1.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))

top.genes.pop2 <- pop2.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))

top.genes.pop3 <- pop3.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))
top.genes.pop4 <- pop4.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))
top.genes.pop5 <- pop5.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))
top.genes.pop6 <- pop6.tajima %>% 
    arrange(desc(Tajima)) %>%
    head(n = 5) %>% 
    inner_join(Gene_name) %>% 
    select(-c(10:13))

outliers <- bind_rows(
    top.genes.pop1, 
    top.genes.pop2,
    top.genes.pop3,
    top.genes.pop4,
    top.genes.pop5,
    top.genes.pop6
)

#=================================================
## Correlation between Tajima's D and Fu & Li F*
#=================================================
combined.tajima %>% 
    ggplot(aes(x=Tajima, y=Fu.Li_F)) + 
    geom_point(size=1, alpha = 3/5) +
    geom_label_repel(data=outliers, 
                     aes(x=Tajima, 
                         y = Fu.Li_F, 
                         label = GeneName, 
                         fill = factor(GeneName)), force = 5,
                     size = 1.5, force_pull = 5, max.overlaps = Inf) + 
    geom_point(size=1, data = combined.tajima[combined.tajima$Tajima  >1 & combined.tajima$Fu.Li_F > 1,], color = "red") +
    labs( x = "Tajima's D", y = "Fu & Li F*") + theme_bw() + 
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    facet_wrap(~Population, nrow = 3)


ggsave("../results/Result1.jpeg")   
 


