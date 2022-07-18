

rm(list = ls())

# install.packages(c("data.table", "tictoc", "statip"))

library(data.table)
library(tictoc)
library(statip)
library(rlang)
options(scipen=999)

#-------- computing missingness on individuals
computeMissingness = function(dataFrame)
{
    dataFrame = t(dataFrame)
    missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
    for(i in 1:dim(dataFrame)[1])
    {
        m=0
        for(j in 1:dim(dataFrame)[2])
        {
            if(dataFrame[i,j] == './.')
                m=m+1
        }
        missingness[i] = m/(dim(dataFrame)[2])
    }
    return(missingness)
}

#-------- computing missingness on SNPs
computeMissingnessOnSNPs = function(dataFrame)
{
    dataFrame = as.matrix(dataFrame)
    missingness = vector(mode = "numeric", length = dim(dataFrame)[1])
    for(i in 1:dim(dataFrame)[1])
    {
        m=0
        for(j in 1:dim(dataFrame)[2])
        {
            if(dataFrame[i,j] == './.')
                m=m+1
        }
        missingness[i] = m/(dim(dataFrame)[2]) #-1
    }
    return(missingness)
}

## TRIQUAD ##
# Calculates the number of different homozygous calls present in the data (1 - 4). 
# Take heterozygous calls into account. 
# Also does not take into account the frequency of each allele.

triquad <- function(x)
{
    xx <- x[x != "./."];
    res <- 1*(sum(1*(xx=="0/0"))>0) + 1*(sum(1*(xx=="1/1"))>0) + 1*(sum(1*(xx=="0/1"))>0);
    res
}

#------------ Extracting genotypes from the raw data:
vcf <- "raw_data/Mali_bi_final_samples_maf01.vcf.gz"

Genotypes <- 'Genotypes.txt'

expression <- '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, Genotypes))

#------- Computing the missingness on SNPs and samples
#------- Get the sample names from the VCF file

sampleList <- 'list_of_isolates.txt'
system(sprintf("bcftools query -l %s > %s", vcf, sampleList))

list_of_isolates <- fread(sampleList, header = FALSE)

#---- put the first four columns in a variable
Genotype <- fread(Genotypes, header = FALSE)

first4Column <- subset(Genotype, select=c(1:4))
Genotype <- subset(Genotype, select=-c(1:4))

#===================
## Remove invariants
#==================
varsnp <- apply(Genotype, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, Genotype)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#======================================
#--- computing missingness on SNPs
#======================================

tic();snpMissingness = computeMissingnessOnSNPs(geno);toc()

geno <- cbind(first4Column, geno, snpMissingness)

geno <- geno[which(geno$snpMissingness <= 0.2), ]
first4Column <- subset(geno, select = c(1:4))
geno <- subset(geno, select = -c(1:4, ncol(geno)))

#================================
#=== Compute Missingness on samples
#================================
tic();sampleMissingness = computeMissingness(geno);toc()

geno <- t(geno)
rownames(geno) <- list_of_isolates$V1
geno <- as.data.frame(geno)
geno$Missingness <- sampleMissingness

index <- which(geno$Missingness > 0.2)

geno <- subset(geno, select = -c(ncol(geno)))
geno <- as.data.frame(t(geno))

if(!is_empty(index)) geno <- geno[, -index]

#===================
## Remove invariants
#==================
varsnp <- apply(geno, 1, triquad)
keepvar <- which(varsnp != 1)

data2 <- cbind(first4Column, geno)
data2 <- data2[keepvar,]
geno <- data2[,5:ncol(data2)]
first4Column <- data2[,1:4]
nrow(geno)
ncol(geno)

#======================================
#--- 2nd Run: computing missingness on SNPs
#======================================

tic();snpMissingness = computeMissingnessOnSNPs(geno);toc()

geno <- cbind(first4Column, geno, snpMissingness)

geno <- geno[which(geno$snpMissingness < 0.01), ]
first4Column <- subset(geno, select = c(1:4))
geno <- subset(geno, select = -c(1:4, ncol(geno)))

#================================
#=== 2nd Run: Compute Missingness on samples
#================================
tic();sampleMissingness = computeMissingness(geno);toc()

geno <- t(geno)
rownames(geno) <- list_of_isolates$V1[-index]
geno <- as.data.frame(geno)
geno$Missingness <- sampleMissingness

index <- which(geno$Missingness > 0.01)

geno <- subset(geno, select = -c(ncol(geno)))
geno <- as.data.frame(t(geno))

if(!is_empty(index)) geno <- geno[, -index]

## Save positions and isolates to keep for downstream analysis
write.table(colnames(geno), 'samplesTokeep.txt', 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(first4Column[,1:2], "snpsTokeep.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

#=========================================================
#---- Removing SNPs and Isolates to discard from vcf file
#========================================================
samplesToKeep = "samplesTokeep.txt"
snpsToKeep = "snpsTokeep.txt"

filtered.vcf <- gsub(".vcf.gz", ".filtered", vcf)
system(paste0("vcftools --gzvcf ", vcf,
              " --keep ", samplesToKeep, 
              " --positions ", snpsToKeep,
              " --not-chr Pf3D7_API_v3", 
              " --recode --recode-INFO-all --out ", filtered.vcf))

file.remove(Genotypes, samplesToKeep, snpsToKeep, paste0(filtered.vcf, ".log"))

system(paste0("mv ", filtered.vcf, ".recode.vcf ", filtered.vcf, ".vcf"))
system(paste0("bgzip -f ", filtered.vcf, ".vcf"))
system(paste0("tabix -f ", filtered.vcf, ".vcf.gz"))
