
rm(list=ls())

library(tidyverse)

VCF <- "/media/Data/Data/Documents_Karim/Fadel/Aoua/Mali_bi_final_samples_maf01.vcf.gz"
metadata <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/Mali_metadata.txt")
    
pop <- metadata %>% select(Location) %>% 
    pull() %>% unique() %>% sort()

for (i in pop) {
    cat("====================\n")
    cat("Population ", i, "\n")
    cat("====================\n")
    samples <- paste0("/media/Data/Data/Documents_Karim/Fadel/Aoua/", i, ".txt")
    metadata %>% 
        filter(Location == i) %>% 
        select(Sample_id) %>% 
        pull() %>% 
        write.table(samples, col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    OUT <- paste0("/media/Data/Data/Documents_Karim/Fadel/Aoua/", i)
    system(paste0("vcftools --gzvcf ", VCF,
                  " --keep ", samples, 
                  " --recode --recode-INFO-all --out ", OUT))
    
    file.remove(samples)
    
    system(paste0("mv ", OUT, ".recode.vcf ", OUT, ".vcf"))
    system(paste0("bgzip ", OUT, ".vcf"))
    system(paste0("tabix ", OUT, ".vcf.gz"))
}
