
library(tidyverse)
vcfName <- "../../Aoua/Mali_bi_Snps20_indv20_dp5_MQ30_maf2.vcf.recode.vcf"

skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), TRUE))
vcf  = read.table( vcfName, skip=skipNum, header=TRUE, 
                   comment.char="", stringsAsFactors = FALSE, check.names=FALSE)


genotype <- "../../Aoua/Genotypes_Mali_Data.txt"

genotype <- read_tsv(genotype)
samples <- genotype %>% select(Samples) %>% pull()

POS <- as.numeric(sapply(strsplit(colnames(genotype)[-c(1:2)], "_"), "[", 4))

genotype <- genotype %>% select(-c(1,2)) %>%
    tibble::rownames_to_column() %>%  
    pivot_longer(-rowname) %>% 
    pivot_wider(names_from=rowname, values_from=value) %>% 
    select(-1)

names(genotype) <- samples

## Replace 0 per 0/0
genotype <- as_tibble(apply(genotype, MARGIN = 2, function(x){ gsub("0", "0/0", x)}))

## Replace 1 per 1/1
genotype <- as_tibble(apply(genotype, MARGIN = 2, function(x){ gsub("1", "1/1", x)}))

genotype[is.na(genotype)] <- "./."


firstColumns <- as_tibble(vcf[,1:9])

missing_position <- setdiff(firstColumns$POS, as.integer(POS))

firstColumns <- firstColumns[!firstColumns$POS %in% missing_position, ]

newVcfFilename = "../../Aoua/Mali_bi_Snps20_indv20_dp5_MQ30_maf2.vcf"

tmpVcf <- cbind(firstColumns, genotype)

system ( paste("grep \"##\"", vcfName, ">", newVcfFilename) )
names(tmpVcf)[1] = "#CHROM"


write.table(tmpVcf, file = newVcfFilename, append = T, sep = "\t", quote = F, row.names = F)




