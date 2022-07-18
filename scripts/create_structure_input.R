

# Data check and preparation to examine crossbreds with the Structure software
# See video on the Genomics Boot Camp YouTube channel

# Prerequisites:
# 1) Download and install Structure https://web.stanford.edu/group/pritchardlab/structure.html
# 2) Get data https://datadryad.org/stash/dataset/doi:10.5061/dryad.v8g21pt


# Clear workspace and load packages
rm(list = ls())
library(tidyverse)


# perform Quality control

samples <- colnames(Genotype)[-ncol(Genotype)]

pop.info <- metadata %>% 
    filter(Sample %in% samples) %>% 
    dplyr::select(Sample, Location)

plink <- "/home/karim/Documents/Mes_Programmes/plink_linux_x86_64_20181202/plink"
vcfname <- "raw_data/Mali_bi_final_samples_maf01.filtered.vcf.gz"

system(paste0(plink , " --vcf ",vcfname, " --allow-extra-chr ",
              " --geno 0.1 --mind 0.01 --maf 0.05 ",
              " --make-bed --out results/afterQC"))

#################
# Check PCA plot
#################

system(paste0(plink , " --bfile results/afterQC --allow-extra-chr --pca --out results/plinkPCA"))

###
# Visualize PCA results
###

# read in result files
eigenValues <- read_delim("results/plinkPCA.eigenval", delim = " ", col_names = F)
eigenVectors <- read_delim("results/plinkPCA.eigenvec", delim = " ", col_names = F) %>% 
    inner_join(pop.info, by = c("X1" = "Sample"))

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

# prepare file for the Structure software
# system("plink --bfile afterQC --chr-set 29 --recode structure --out forStructure")
system(paste0(plink, " --bfile results/afterQC --allow-extra-chr --recode structure --out results/forStructure"))













rm(list = ls())
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

Genotype <- as.matrix(Genotype)

Genotype[Genotype == "1/1"] <- 3
Genotype[Genotype == "0/1"] <- 2
Genotype[Genotype == "0/0"] <- 1
Genotype[Genotype == "./."] <- -9

population <- metadata[index, ] %>% 
    select(Sample, Location) %>% 
    mutate(Pop = case_when(Location == "Kolle" ~ 1, 
                           Location == "Faladje" ~ 2,
                           Location == "Bamako" ~ 3,
                           Location == "Bougoula-Hameau" ~ 4,
                           Location == "Kenieroba" ~ 5,
                           Location == "Nioro" ~ 6,
                           Location == "Dangassa" ~ 7))

Genotype <- Genotype %>% as_tibble()
Genotype <- as_tibble(t(Genotype))

Genotype <- cbind(rownames(Genotype), Genotype) %>% as_tibble() %>% 
    inner_join(population, by = c("V1" = "Sample")) %>% 
    relocate(Pop, .after = V1) %>% 
    relocate(Location, .after = Pop)

colnames(Genotype)[-(1:3)] <- paste0("Locus", 1:(ncol(Genotype)-3))

write.table(Genotype, "Structure/MaliforStructure.txt", 
            col.names = T, row.names = F, quote = F)
