
rm(list = ls())

FILE <- "raw_data/Mali"

plink <-"~/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"

common_isolates <- read_tsv("raw_data/Mali_bi_final_samples_maf01.filtered.nosex", col_names = F) %>% 
    inner_join(metadata, by = c("X1" = "Sample"))

common_isolates %>% 
    dplyr::select(X1, X2) %>% 
    write.table("samples_in_metadata.txt", col.names = F, row.names = F, quote = F)

common_isolates %>% 
    rename(Samples = X1) %>% 
    dplyr::select(Samples, Location) %>% 
    write.table("metadata.txt", col.names = T, row.names = F, quote = F, sep = '\t')

# Generate the input file in plink format
system(paste0(plink, " --vcf ", FILE, 
              "_bi_final_samples_maf01.filtered.vcf.gz --make-bed --out ", FILE,
              " --keep samples_in_metadata.txt",
              " --geno 0.1 --mind 0.1 --maf 0.05 ",
              " --allow-extra-chr"))

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
system(paste0("awk '{$1=\"0\";print $0}' ", FILE, ".bim > ", FILE, ".bim.tmp"))

system(paste0("mv ", FILE, ".bim.tmp ", FILE, ".bim"))