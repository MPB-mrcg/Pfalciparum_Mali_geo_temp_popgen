
rm(list = ls())

library(tidyverse)
library(rehh)

VCF <- "data/reference/Mali_snps_data.recode.vcf"
system(paste0("bgzip -f ", VCF))
system(paste0("tabix -f ", VCF, ".gz"))
VCF <- paste0(VCF, ".gz")

OUTPUT <- "results/tables/"

data <- read.table(VCF, header = FALSE, stringsAsFactors = FALSE)

metadata <- read_tsv("data/metadata/Mali_metadata.txt") %>% 
    filter(Location == "Faladje") %>% 
    dplyr::select(Sample, Year, Location) 

metadata %>% 
    filter(Year != 2014) %>% 
    count(Year)

years <- metadata %>% 
    filter(Year != 2014) %>% 
    pull(Year) %>% unique()

ihs_file <- "results/tables/Faladje_iHs.xlsx"

chromosomes <- data %>% pull(V1) %>% unique()

for (i in years) {
    cat("Performing iHs for", i)
    
    metadata %>% 
        filter(Year == i) %>%
        pull(Sample) %>% 
        write.table(paste0(OUTPUT, i, ".txt"), col.names = F, row.names = F, quote = F)
    
    for(c in 1:length(chromosomes)) {
        # Split VCF by chromosome
        system(paste0("bcftools view -r ", chromosomes[c],
                      " -S ", paste0(OUTPUT, i, ".txt"),
                      " -Oz -o ", paste0(OUTPUT, i, "_", chromosomes[c], ".vcf.gz"),
                      " ", VCF))
        
        # vcf file name for each chromosome
        vcf_file <- paste0(OUTPUT, i, "_", chromosomes[c], ".vcf.gz")
        
        # create internal representation
        hh <- data2haplohh(hap_file = vcf_file,
                           polarize_vcf = FALSE,
                           # vcf_reader = "data.table")
                           vcf_reader = "vcfR")
        
        # perform scan on a single chromosome (calculate iHH values)
        scan <- scan_hh(hh)
        
        # concatenate chromosome-wise data frames to
        # a data frame for the whole genome
        # (more efficient ways certainly exist...)
        
        if (c == 1) {
            wgscan <- scan
        } else {
            wgscan <- rbind(wgscan, scan)
        }
        
        file.remove(paste0(OUTPUT, i, "_", chromosomes[c], ".vcf.gz"))
    }
    
    # calculate genome-wide iHS values
    wgscan.ihs <- ihh2ihs(wgscan)
    ihs <- cbind.data.frame(i, wgscan.ihs$ihs)
    
    ## SAVE MAF FILE
    if(!file.exists(ihs_file)){
        file.create(ihs_file)
        write.table(ihs, ihs_file, col.names = TRUE,
                    row.names = FALSE, quote = FALSE, sep = '\t')
    }
    else write.table(ihs, ihs_file, col.names = FALSE, 
                     append = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
    
    file.remove(paste0(OUTPUT, i, ".txt"))
}


ihs_data <- read_tsv(ihs_file) %>% 
    drop_na(IHS)

ihs_data$Chr <- as.integer(factor(ihs_data$CHR), levels = unique(ihs_data$CHR))
ihs_data$P = 10**(-ihs_data$LOGPVALUE)

# Adjust pvalue using qvalues
ihs_data$P = 10**(-ihs_data$LOGPVALUE)

adjusted.p <- p.adjust(ihs_data$P, method = "fdr",  n = length(ihs_data$P))
fdr = fdrtool(ihs_data$P, statistic = "pvalue")
ihs_data$Qval <- fdr$qval

qvalues <- qvalue(ihs_data$P, fdr.level = 0.01)
ihs_data$Qvalues <- qvalues$qvalues

test <- ihs_data %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len = max(POSITION)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(ihs_data, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    # arrange(X1, X2) %>%
    mutate( BPcum = POSITION + tot) %>% 
    ungroup()

axis_set <- test %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

test$P = 10**(-test$LOGPVALUE)
ylab <- bquote("-" ~ log[10] ~ scriptstyle(italic(.("pvalue"))) )

# Ready to make the plot using ggplot2:
ggplot(test, aes(x = BPcum, y = -log10(P), color = factor(Chr), size = LOGPVALUE)) + # -log10(LOGPVALUE)
    # Show all points
    geom_point(size = .15) + 
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = rep(c("gray30", "firebrick"),
                                    unique(length(axis_set$Chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = ylab) +
    
    # Custom the theme:
    theme_light() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 5),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(face = "bold.italic", size = 9),
        axis.title.y = element_text(face = "bold.italic", size = 12),
        axis.text = element_text(colour = "black", face = "bold", size = 6),
        strip.text.x = element_text(size = 8, color = "black", face = "bold.italic")) +
    facet_wrap(~factor(i), scales = "free_y", ncol = 1)


ggsave("results/figures/Faladje_iHs.pdf", width = 180, height = 150, units = "mm",)
