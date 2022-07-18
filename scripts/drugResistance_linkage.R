

BiocManager::install("snpStats")
devtools::install_github("SFUStatgen/LDheatmap")

library(data.table)
library(LDheatmap)
library(genetics)
library(vcfR)
library(gaston) # Can be also used to estimate LD
# https://developpaper.com/r-language-draw-a-heat-map-of-a-triangle/

rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")

#=====================
# 1. 1st method: VCF
#=====================
gene_list  <- list.files("results/Allele_frequency", full.names = T, recursive = T)

for (file in gene_list) {
    #read the VCF file and store the vcfR object
    vcf <- read.vcfR("results/Allele_frequency/drugResistance_loci.vcf")
    
    #extract genetic distances, subject IDs and genotypes from the vcfR object
    list_mali <- vcfR2SnpMatrix(vcf)
    
    SNP_name <- vcf@fix %>% 
        as_tibble() %>% 
        mutate(ID = paste(CHROM, POS, sep = "_")) %>% 
        pull(ID)
    
    #draw the heatmap
    MyHeatmap <- LDheatmap(list_mali$data, list_mali$genetic.distance, add.map=T, 
                           add.key=T, color= rgb.palette(18), flip=T,text=T)
    
    LDheatmap(MyHeatmap, SNP.name = SNP_name, flip=T,text=T)
}

#=====================
# 1. 2nd method: VCF
#=====================

vcf <- "results/Allele_frequency/drugResistance_loci.vcf"
Genotypes <- 'Genotypes.txt'
expression <- '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'

system(paste0("bcftools query -f'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ", vcf, " > ",Genotypes))

#---- put the first four columns in a variable
Genotype <- fread(Genotypes, header = FALSE)

first4Column <- subset(Genotype, select=c(1:4))
Genotype <- as.data.frame(subset(Genotype, select=-c(1:4)))

Genotypes <- matrix(NA, dim(Genotype)[1], dim(Genotype)[2])
# Recode Genotype
for (i in 1:nrow(Genotype)) {
    
    for (j in 1:ncol(Genotype)) {
        if(Genotype[i,j] == "0/0") Genotypes[i,j] <- paste(first4Column$V3[i],first4Column$V3[i], sep = '/')
        if(Genotype[i,j] == "1/1") Genotypes[i,j] <- paste(first4Column$V4[i],first4Column$V4[i], sep = '/')
        if(Genotype[i,j] == "0/1") Genotypes[i,j] <- paste(first4Column$V3[i],first4Column$V4[i], sep = '/')
        if(Genotype[i,j] == "./.") Genotypes[i,j] <- NA
    }
    
}

Genotypes <- as.data.frame(Genotypes)


genetic.distance <- first4Column$V2
SNP_name <- paste(first4Column$V1,first4Column$V2, sep = '_')

Genotypes <- t(Genotypes) %>% as_tibble()
names(Genotypes) <- as.character(SNP_name)

Genotypes <- Genotypes %>% as.data.frame() %>% 
    mutate_if(is.character,as.factor)

# convert the columns of the original data frame into  genotype objects
for(i in 1:20){
    Genotypes[,i]<-as.genotype(Genotypes[,i]) 
}

class(Genotypes[,1])

LDheatmap(Genotypes, genetic.distances = genetic.distance,add.map=F,
          add.key=F, color= rgb.palette(18), flip=F, text=T, SNP.name = SNP_name,
          title=NULL)

MyHeatmap <- LDheatmap(Genotypes, genetic.distances = genetic.distance,add.map=F,
                       add.key=F, color= rgb.palette(18), flip=F, text=T, SNP.name = SNP_name,
                       title=NULL)
#================
# 3rd Method
#================
# Use of gaston package to estimate pairwise LD

LDmatrix <- MyHeatmap$LDmatrix
    
LD.plot( LDmatrix, snp.positions = genetic.distance, draw.chr = FALSE, cex.snp=.8,cex.ld = 0.7, 
         color.scheme = function(ld) rgb(1,1-abs(ld),1-abs(ld)),
         pdf.file="results/figures/drugResistance_linkage.pdf", finalize.pdf = TRUE)

LDmatrix <- LDmatrix[ order(row.names(LDmatrix)), order(colnames(LDmatrix))]
write.table(LDmatrix, "results/tables/LD_matrix.xlsx", col.names = T, row.names = T, quote = F, sep = '\t')

testmatrix <- readxl::read_xlsx("results/tables/LD_matrix.xlsx") %>% 
    mutate_at(vars(-SNPs), as.numeric) 

test<-as.matrix(testmatrix[-1])
rownames(test)=testmatrix$SNPs   #function(ld) rgb(1,1-max(ld),1-max(ld))

LD.plot( test, snp.positions = genetic.distance, draw.chr = FALSE, cex.snp=.6,cex.ld = 0.6, 
         color.scheme = function(ld) rgb(1,1-abs(ld),1-abs(ld)),
         polygon.par = list(border = "white"),
         pdf.file="results/figures/drugResistance_linkage.pdf", finalize.pdf = TRUE)
    