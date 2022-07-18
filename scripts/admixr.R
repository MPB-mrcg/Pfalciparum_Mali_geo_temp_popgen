

rm(list = ls())

install.packages("admixr")

# To install the development version from Github (which might be slightly ahead in terms of 
# new features and bugfixes compared to the stable release on CRAN), you need the package devtools. 
# You can run:
    
# install.packages("devtools")
# devtools::install_github("bodkan/admixr")

usethis::edit_r_environ()

library(admixr)
library(tidyverse)
library(data.table)

metadata <- read_tsv("raw_data/Mali_Matadata.txt")

vcf <- "raw_data/Mali_bi_final_samples_maf01.filtered.vcf.gz"
Genotypes <- 'Genotypes.txt'
expression <- '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

bcftools <-"/usr/local/bin/bcftools"
system(paste0(bcftools, " query -f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ", vcf, " > ",Genotypes))

#------- Get the sample names from the VCF file

sampleList <- 'list_of_isolates.txt'
system(paste0(bcftools, " query -l ", vcf, " > ",sampleList))
list_of_isolates <- fread(sampleList, header = FALSE)

#---- put the first four columns in a variable
Genotype <- fread(Genotypes, header = FALSE)
names(Genotype) <- c("CHROM", "POS", "REF", "ALT", list_of_isolates$V1)

first4Column <- subset(Genotype, select=c(1:4))
Genotype <- subset(Genotype, select=-c(1:4))

index <- which(colnames(Genotype) %in% metadata$Sample)

Genotype <- as.data.frame(Genotype)[, index]

#==============================
# Create Input file for admixr
#==============================

# 1. Individual file
 
metadata %>% 
    filter(Sample %in% colnames(Genotype)) %>% 
    dplyr::select(Sample, Location) %>% 
    mutate(Sex = "U") %>% 
    relocate(Sex, .before = Location) %>% 
    write.table("admixure/mali.ind", col.names = F, row.names = F, quote = F, sep = "\t")


# 2. SNP file
first4Column %>% 
    as_tibble() %>% 
    mutate(snp_name = paste(CHROM, POS, sep = "_")) %>% 
    mutate(CHROM = as.numeric(gsub("[Pf3D7_, _v3]", "", CHROM))) %>% 
    mutate(genetic_position = 0.0) %>% 
    relocate(snp_name, .before = CHROM) %>% 
    relocate(genetic_position, .before = POS) %>% 
    write.table("admixure/mali.snp", col.names = F, row.names = F, quote = F, sep = "\t")


# 3. Genotype file
Genotype <- as.matrix(Genotype)

Genotype[Genotype == "1/1"] <- 0
Genotype[Genotype == "0/1"] <- 1
Genotype[Genotype == "0/0"] <- 2
Genotype[Genotype == "./."] <- 9

Genotype <- Genotype %>% as_tibble()

Genotype$Genotypes <- apply(Genotype, 1, function(x) paste(x, collapse = ""))

write.table(Genotype$Genotypes, "admixure/mali.geno", 
            col.names = F, row.names = F, quote = F)

file.remove(sampleList, Genotypes)

snps <- eigenstrat("admixure/mali")

pops <- metadata %>% pull(Location) %>% unique()

# Using the admixr package we can then calculate our D statistic simply by running:
    
result <- d(W = pops, X = "Kolle", Y = "Faladje", Z = "Bamako", 
            outdir = "admixure/", data = snps)


result

result_f4ratio <- f4ratio(
    X = pops, A = "Kolle", B = "Faladje", C = "Bamako", O = "Kenieroba",
    data = snps
)
