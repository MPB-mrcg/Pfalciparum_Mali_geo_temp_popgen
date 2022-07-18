library(data.table)
library(tidyverse)

vcf = "/media/Data/Data/Documents_Karim/Fadel/Aoua/Mali_bi_Snps20_indv20_dp5_MQ30_maf2.vcf"

prefix = "/media/Data/Data/Documents_Karim/Fadel/Aoua/Mali"

source('../Codes/craig_functions.r')
## tajima's D function

tajima <- function(haps,lege,Gene,msnp,minn)
{
    keep1 <- lege[,"gene"] == Gene	
    haps1 <- haps[keep1,]
    lege1 <- lege[keep1,]
    
    rr <- sum(1*( keep1));
    
    if (rr > 1)
    {
        
        ## Filter isolates with missing data
        isomiss <- apply(haps1, 2, function(x) length(which((x != 1 & x != 0) | is.na(x))))
        ## Calculate the number remaining and number total, need to add this to results
        isototal <- length(isomiss)
        isokeeptotal <- length(which(isomiss == 0))
        haps1 <- haps1[,isomiss==0]
        isopercent <- isokeeptotal / isototal
        ## Remove monomorphics
        # SNPs with two alleles will have a length(table()) of 2
        hapskeep <- apply(haps1, 1, function(x) length(table(x)))
        haps1 <- haps1[hapskeep == 2,]
        lege1 <- lege1[hapskeep == 2,]
        
    }
    
    keep1 <- lege1[,"gene"] == Gene
    
    results <- matrix(NA,nrow=1,ncol=9) #generate results matrix
    colnames(results) <- c("Gene","H","S","theta","khat","tajimasd","#Iso_used","Total#Iso", "Iso%") #Add colnames
    results[1,1] <- Gene #add the gene id to the results table
    ss <- sum(1*( keep1));
    results[1,3] <- ss #put this info into the table
    
    if(ss > 1) #if more than one segregating site
    {
        temp <- haps1

        temps <- apply(temp,2,collap); count.temp <- table(temps) #collapses the snps for each isolate into a single entry (ie the haplotype for the gene) then uses the table function to determine the number of unique haplotypes and number of isolates with each unique haplotype
        
        n.ind <- sum(count.temp); freqs.temp <- count.temp / n.ind #determines the number of samples represented then determines the frequency of each individual haplotype within the population
        results[1,2] <- length(count.temp); #calculates the number of haplotypes and puts it into the table
        results[1,7] <- isokeeptotal
        results[1,8] <- isototal
        results[1,9] <- isopercent
        if(as.numeric(results[1,3]) >= msnp) #if there are at least 3 segregating sites within the gene...
        {
            all  <- as.matrix( dist(t(temp), "manhattan",diag=T,upper=T)) #calculates the distance between haplotypes
            khat <- mean(c((all)[lower.tri(all)]),na.rm=T)
            a1 <- sum( 1/(1:(n.ind-1))) # calculates the sum of 1 / (1:(number of samples-1))
            a2 <- sum( 1/(1:(n.ind-1))^2); #calculates sum of 1 / (1:(num samples-1))^2
            b1 <- (n.ind + 1)/3/(n.ind-1) #calculates (samples + 1) / 3 / (samples - 1)
            b2 <- 2*(n.ind^2 + n.ind +3)/9/n.ind/(n.ind-1); c1 <- b1 - 1/a1
            c2 <- b2 - (n.ind + 2)/ a1 / n.ind + a2/a1/a1; e1 <- c1/a1
            e2 <- c2 / (a1^2 + a2)
            results[1,4] <- ss / a1
            results[1,5] <- khat
            results[1,6] <- (khat - ss/a1) / sqrt ( e1 * ss + e2 * ss * (ss-1))
        }
    }
    results[!is.na(results[1,"tajimasd"]),]
}

## end of function

vcf <- read.table(vcf)
firstColumns <- vcf %>% select(c(1:9))
vcf <- vcf %>% select(-c(1:9))

genome <- as.matrix(vcf)
genome[genome == "0/0"]= 0
genome[genome == "1/1"]= 1
genome[genome == "./."]= "."

genome <- as.data.frame(cbind(firstColumns, genome))

Triquad <- function(x)
{
    xx <- x[x != "-" & x != "N"];
    res <- 1*(sum(1*(xx=="0"))>0) + 1*(sum(1*(xx=="1"))>0);
    res
}
## Remove invariants
geno <- genome[,10:ncol(genome)]
varsnp <- apply(geno, 1, Triquad)
keepvar <- which(varsnp == 2)

genome <- genome[keepvar,]

## Leaves 97971 SNPs

genic <- genome[-grep("INTERGENIC", genome$V8),]
exonic <- genome[grep("CODING", genome$V8), ]

gids <- strsplit(as.vector(exonic$V8), ";")

geneIDs <- NULL
for (i in 1:length(gids))
{
    print(paste0('i= ', i))
    
    pp <- grep("SNPEFF_TRANSCRIPT_ID", gids[[i]])
    geneIDs <- c(geneIDs, gids[[i]][pp])
    geneIDs <- gsub("SNPEFF_TRANSCRIPT_ID=\"","", geneIDs)
    geneIDs <- gsub("\"","", geneIDs)
}

exonic$gene <- geneIDs

exonic <- as.data.frame(exonic)
exonic <- exonic[ , c(1:9, ncol(exonic), 10:(ncol(exonic)-1))]

datalegend <- exonic[,1:10]
datagenotypes <- exonic[,11:ncol(exonic)]
genos <- datagenotypes
genos <- apply(genos, 2, as.numeric)

genenames <- names(table(datalegend[,"gene"]))
genenames <- genenames[genenames != "-"]

tajimascores3 <- NULL
for (i in genenames)
{
    tajimascores3 <- rbind(tajimascores3, tajima(genos, datalegend, i, 3, 1))
}

write.table(tajimascores3, paste0(prefix,'_TajimaD.txt'), col.names = T, row.names = F, quote = F)
