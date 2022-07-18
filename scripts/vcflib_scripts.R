
rm(list = ls())

library(tidyverse)

vcf <- "data/reference/Mali_snps_data.recode.vcf"

metadata <- read_tsv("data/metadata/Mali_metadata.txt") %>% 
    dplyr::select(Sample, Location, Year)

sampleIDs <- scan("sampleIDs.txt", what = character())

metadata <- metadata[metadata$Sample %in% sampleIDs,]

#===================
## Geographical Fst
#===================

metadata$index <- seq(0, nrow(metadata)-1)

population <- metadata %>% pull(Location) %>% unique()

wcFst <- "/home/karim/Documents/Mes_Programmes/vcflib/bin/wcFst "
vcflib_file <- "scripts/vcflib_scripts.sh"
vcflib <- file(vcflib_file, open='w')
cat("#!/bin/bash", file = vcflib, sep = '\n')
cat(" ", file = vcflib, sep = '\n')

for (i in 1:(length(population)-1)) {
    
    for (j in (i+1):length(population)) {

        output <- paste0("results/tables/Fst/Geographical/", 
                         population[i], "_", population[j], "_wcFst.txt")
        
        target <- metadata %>% filter(Location == population[i]) %>% pull(index)
        target <- paste(target, collapse=",")
        
        background <- metadata %>% filter(Location == population[j]) %>% pull(index)
        background <- paste(background, collapse=",")
        
        cat(" ", file = vcflib, sep = '\n')
        cat("#=========================", file = vcflib, sep = '\n')
        cat(paste0("# ", population[i], "_", population[j]), file = vcflib, sep = '\n')
        cat("#=========================", file = vcflib, sep = '\n')
        cat(" ", file = vcflib, sep = '\n')
        cat(paste0(wcFst, "--target ", target,
                   " --background ", background, 
                   " --file ", vcf,
                   " --deltaaf 0.05 --type GT > ",
                   output),
            file = vcflib, sep = '\n')
    }
}

close(vcflib)

system(paste0("bash ", vcflib_file))

#===============
# Temporal Fst
#===============

year <- metadata %>% pull(Year) %>% unique()


for (i in 1:(length(year)-1)) {
    
    target <- metadata %>% filter(Year == year[i]) %>% pull(index)
    target <- paste(target, collapse=",")
    
    for (j in (i+1):length(year)) {
        
        output <- paste0("results/tables/Fst/Temporal/", 
                         year[i], "_", year[j], "_wcFst.txt")
        
    
        background <- metadata %>% filter(Year == year[j]) %>% pull(index)
        background <- paste(background, collapse=",")
        
        system(paste0(wcFst, "--target ", target,
                   " --background ", background, 
                   " --file ", vcf,
                   " --deltaaf 0.05 --type GT > ",
                   output))
    }
}

#=======================
# Heterogeneity Faladje
#======================
# NB: We combined samples from 2015 and 2017 into one

#============
# 2007 -2013
#============
output <- paste0("results/tables/Fst/Heterogeneity/wcFst_07_13.txt")

target <- metadata %>% filter(Location == "Faladje" & Year == 2007) %>% pull(index)
target <- paste(target, collapse=",")

background <- metadata %>% filter(Location == "Faladje" & Year == 2013) %>% pull(index)
background <- paste(background, collapse=",")

system(paste0(wcFst, "--target ", target,
              " --background ", background, 
              " --file ", vcf,
              " --deltaaf 0.05 --type GT > ",
              output))

#==================
# 2007 - 2015/2017
#==================
output <- paste0("results/tables/Fst/Heterogeneity/wcFst_07_1517.txt")

background <- metadata %>% filter(Location == "Faladje" & (Year == 2015 | Year == 2017)) %>% pull(index)
background <- paste(background, collapse=",")

system(paste0(wcFst, "--target ", target,
              " --background ", background, 
              " --file ", vcf,
              " --deltaaf 0.05 --type GT > ",
              output))


#============
# 2013 -2015/2017
#============
output <- paste0("results/tables/Fst/Heterogeneity/wcFst_13_1517.txt")

target <- metadata %>% filter(Location == "Faladje" & Year == 2013) %>% pull(index)
target <- paste(target, collapse=",")

# background <- metadata %>% filter(Location == "Faladje" & (Year == 2015 | Year == 2017)) %>% pull(index)
# background <- paste(background, collapse=",")

system(paste0(wcFst, "--target ", target,
              " --background ", background, 
              " --file ", vcf,
              " --deltaaf 0.05 --type GT > ",
              output))
