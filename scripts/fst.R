library(radiator)
library(diveRsity)
library(PGDspideR)
library("devtools")
library(adegenet)
library("hierfstat")
library(vcfR)
library(MASS) 
library(dplyr)
library(data.table)
install.packages("AnalyzeFMRI")
install.packages("GeneNet")
library("GeneNet")
install.packages("FDRsampsize")
library(FDRsampsize)
install.packages("sommer")
library(qvalue)
library("sommer")
library(ggplot2)
if (!("qvalue" %in% installed.packages())){TODO}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 
devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)
setwd("~/Documents/Aoua_PhD/data/FST/")
#######read data #########
Maj_read_data <- fread("~/Documents/Aoua_PhD/data/MetaData/Maj_Mali_Data.txt", sep = "\t", header = T)
Maj_genind <-  df2genind(as.data.frame(Maj_read_data[, -c(1,2)]), ploidy = 1, ind.names= Maj_read_data$Samples)
Maj_genind$pop <- factor(Maj_read_data$Locations)
allele_freq <- makefreq(Maj_genind, missing = NA, truenames = TRUE)
####split Data by population#########
bko_boug <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "Bougoula-Hameau")
bko_dang <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "Dangassa")
bko_ken <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "Kenieroba")
bko_fala <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "Faladje")
bko_kol <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "Kolle")
bko_nioro <- Maj_read_data %>% filter( Locations == "Bamako" |  Locations == "nioro")
boug_dang <- Maj_read_data %>% filter( Locations == "Dangassa" |  Locations == "Bougoula-Hameau")
boug_fala <- Maj_read_data %>% filter( Locations == "Faladje" |  Locations == "Bougoula-Hameau")
boug_ken <- Maj_read_data %>% filter( Locations == "Kenieroba" |  Locations == "Bougoula-Hameau")
boug_kol <- Maj_read_data %>% filter( Locations == "Kolle" |  Locations == "Bougoula-Hameau")
boug_nioro <- Maj_read_data %>% filter( Locations == "nioro" |  Locations == "Bougoula-Hameau")
dang_fala <- Maj_read_data %>% filter( Locations == "Dangassa" |  Locations == "Faladje")
dang_ken <- Maj_read_data %>% filter( Locations == "Dangassa" |  Locations == "Kenieroba")
dang_kol <- Maj_read_data %>% filter( Locations == "Dangassa" |  Locations == "Kolle")
dang_nioro <- Maj_read_data %>% filter( Locations == "Dangassa" |  Locations == "nioro")
fala_ken <- Maj_read_data %>% filter( Locations == "Faladje" |  Locations == "Kenieroba")
fala_kol <- Maj_read_data %>% filter( Locations == "Faladje" |  Locations == "Kolle")
fala_nioro <- Maj_read_data %>% filter( Locations == "Faladje" |  Locations == "nioro")
ken_kol <- Maj_read_data %>% filter( Locations == "Kenieroba" |  Locations == "Kolle")
ken_nioro <- Maj_read_data %>% filter( Locations == "Kenieroba" |  Locations == "nioro")
kol_nioro <- Maj_read_data %>% filter( Locations == "Kolle" |  Locations == "nioro")

##########convert df to genind ###############
bko_boug_genind <- df2genind(as.data.frame(bko_boug[, -c(1,2)]), ploidy = 1, ind.names= bko_boug$Samples)
bko_boug_genind$pop <- factor(bko_boug$Locations)
bko_dang_genind <- df2genind(bko_dang[, -c(1,2)], ploidy = 1, ind.names= bko_dang$Samples)
bko_dang_genind$pop <- factor(bko_dang$Locations)
bko_ken_genind <- df2genind(bko_ken[, -c(1,2)], ploidy = 1, ind.names= bko_ken $Samples)
bko_ken_genind$pop <- factor(bko_ken$Locations)
bko_fala_genind <- df2genind(bko_fala[, -c(1,2)], ploidy = 1, ind.names= bko_fala$Samples)
bko_fala_genind$pop <- factor(bko_fala$Locations)
bko_kol_genind <- df2genind(bko_kol[, -c(1,2)], ploidy = 1, ind.names= bko_kol$Samples)
bko_kol_genind$pop <- factor(bko_kol$Locations)
bko_nioro_genind <- df2genind(bko_nioro[, -c(1,2)], ploidy = 1, ind.names=  bko_nioro$Samples)
bko_nioro_genind$pop <-  factor(bko_nioro$Locations)
boug_dang_genind <- df2genind(boug_dang[, -c(1,2)], ploidy = 1, ind.names= boug_dang$Samples)
boug_dang_genind$pop <- factor(boug_dang$Locations)
boug_fala_genind <- df2genind(boug_fala[, -c(1,2)], ploidy = 1, ind.names= boug_fala$Samples)
boug_fala_genind$pop <- factor(boug_fala$Locations)
boug_ken_genind <- df2genind( boug_ken[, -c(1,2)], ploidy = 1, ind.names=  boug_ken$Samples)
boug_ken_genind$pop <-  factor(boug_ken$Locations)
boug_kol_genind <- df2genind(boug_kol[, -c(1,2)], ploidy = 1, ind.names= boug_kol$Samples)
boug_kol_genind$pop <- factor(boug_kol$Locations)
boug_nioro_genind <- df2genind(boug_nioro[, -c(1,2)], ploidy = 1, ind.names= boug_nioro$Samples)
boug_nioro_genind$pop <- factor(boug_nioro$Locations)

dang_fala_genind <- df2genind(dang_fala [, -c(1,2)], ploidy = 1, ind.names= dang_fala$Samples)
dang_fala_genind$pop <-factor(dang_fala$Locations)
dang_ken_genind <- df2genind(dang_ken[, -c(1,2)], ploidy = 1, ind.names= dang_ken$Samples)
dang_ken_genind$pop <- factor(dang_ken$Locations)
dang_kol_genind <- df2genind(dang_kol[, -c(1,2)], ploidy = 1, ind.names= dang_kol$Samples)
dang_kol_genind$pop <- factor(dang_kol$Locations)
dang_nioro_genind <- df2genind(dang_nioro[, -c(1,2)], ploidy = 1, ind.names= dang_nioro$Samples)
dang_nioro_genind$pop <- factor(dang_nioro$Locations)
fala_ken_genind <- df2genind(fala_ken[, -c(1,2)], ploidy = 1, ind.names= fala_ken$Samples)
fala_ken_genind$pop <- factor(fala_ken$Locations)
fala_kol_genind <- df2genind(fala_kol[, -c(1,2)], ploidy = 1, ind.names= fala_kol$Samples)
fala_kol_genind$pop <- factor(fala_kol$Locations)
fala_nioro_genind <- df2genind(fala_nioro[, -c(1,2)], ploidy = 1, ind.names= fala_nioro$Samples)
fala_nioro_genind$pop <- factor(fala_nioro$Locations)
ken_kol_genind <- df2genind(ken_kol[, -c(1,2)], ploidy = 1, ind.names= ken_kol$Samples)
ken_kol_genind$pop <- factor(ken_kol$Locations)
ken_nioro_genind <- df2genind(ken_nioro[, -c(1,2)], ploidy = 1, ind.names= ken_nioro$Samples)
ken_nioro_genind$pop <- factor(ken_nioro$Locations)
kol_nioro_genind <- df2genind(kol_nioro[, -c(1,2)], ploidy = 1, ind.names= kol_nioro$Samples)
kol_nioro_genind$pop <- factor(kol_nioro$Locations)
####convert genind object to hierfstat data
data1 <- genind2hierfstat(bko_boug_genind)
data1$pop <- bko_boug_genind$pop
data1[,1] <- as.integer(data1[,1])
bko_boug_stat_data <- genind2hierfstat(bko_boug_genind)
bko_boug_stat_data$pop <- bko_boug_genind$pop

bko_dang_stat_data <- genind2hierfstat(bko_dang_genind)
bko_dang_stat_data$pop <- bko_dang_genind$pop
bko_fala_stat_data <- genind2hierfstat(bko_fala_genind)
bko_fala_stat_data$pop <- bko_fala_genind$pop
bko_ken_stat_data <- genind2hierfstat(bko_ken_genind)
bko_ken_stat_data$pop <- bko_ken_genind$pop
bko_kol_stat_data <- genind2hierfstat(bko_kol_genind)
bko_kol_stat_data$pop <- bko_kol_genind$pop
bko_nioro_stat_data <- genind2hierfstat(bko_nioro_genind)
bko_nioro_stat_data$pop <- bko_nioro_genind$pop
boug_dang_stat_data <- genind2hierfstat(boug_dang_genind)
boug_dang_stat_data$pop <- boug_dang_genind$pop
boug_fala_stat_data <- genind2hierfstat(boug_fala_genind)
boug_fala_stat_data$pop <- boug_fala_genind$pop
boug_kol_stat_data <- genind2hierfstat(boug_kol_genind)
boug_kol_stat_data$pop <- boug_kol_genind$pop
boug_ken_stat_data <- genind2hierfstat(boug_ken_genind)
boug_ken_stat_data$pop <- boug_ken_genind$pop
boug_nioro_stat_data <- genind2hierfstat(boug_nioro_genind)
boug_nioro_stat_data$pop <- boug_nioro_genind$pop
dang_fala_stat_data <- genind2hierfstat(dang_fala_genind)
dang_fala_stat_data$pop <- dang_fala_genind$pop
dang_ken_stat_data <- genind2hierfstat(dang_ken_genind)
dang_ken_stat_data$pop <- dang_ken_genind$pop
dang_kol_stat_data <- genind2hierfstat(dang_kol_genind)
dang_kol_stat_data$pop <- dang_kol_genind$pop
dang_nioro_stat_data <- genind2hierfstat(dang_nioro_genind)
dang_nioro_stat_data$pop <- dang_nioro_genind$pop
fala_ken_stat_data <- genind2hierfstat(fala_ken_genind)
fala_ken_stat_data$pop <- fala_ken_genind$pop
fala_kol_stat_data <- genind2hierfstat(fala_kol_genind)
fala_kol_stat_data$pop <- fala_kol_genind$pop
fala_nioro_stat_data <- genind2hierfstat(fala_nioro_genind)
fala_nioro_stat_data$pop <- fala_nioro_genind$pop
ken_kol_stat_data <- genind2hierfstat(ken_kol_genind)
ken_kol_stat_data$pop <- ken_kol_genind$pop
ken_nioro_stat_data <- genind2hierfstat(ken_nioro_genind)
ken_nioro_stat_data$pop <- ken_nioro_genind$pop
kol_nioro_stat_data <- genind2hierfstat(kol_nioro_genind)
kol_nioro_stat_data$pop <- kol_nioro_genind$pop

#############bootstraping##############
data1_btt <- boot.ppfst(dat = data1,nboot=100,quant=c(0.01,0.99),diploid=FALSE )
###############fst calcularion ###########
data1_fst <- basic.stats(data1)
perloci_data1_fst <- data1_fst$perloc
bko_boug_fst <- basic.stats(bko_boug_stat_data)
perloci_bko_boug_fst <- bko_boug_fst$perloc
bko_dang_fst <- basic.stats(bko_dang_stat_data)
perloci_bko_dang_fst <- bko_dang_fst$perloc
bko_fala_fst <- basic.stats(bko_fala_stat_data)
perloci_bko_fala_fst <- bko_fala_fst$perloc
bko_ken_fst <- basic.stats(bko_ken_stat_data)
perloci_bko_ken_fst <- bko_ken_fst$perloc
bko_kol_fst <- basic.stats(bko_kol_stat_data)
perloci_bko_kol_fst <- bko_kol_fst$perloc
bko_nioro_fst <- basic.stats( bko_nioro_stat_data)
perloci_bko_nioro_fst <- bko_nioro_fst$perloc
boug_dang_fst <- basic.stats( boug_dang_stat_data)
perloci_boug_dang_fst <- boug_dang_fst$perloc
boug_fala_fst <- basic.stats( boug_fala_stat_data)
perloci_boug_fala_fst <- boug_fala_fst$perloc
boug_kol_fst <- basic.stats( boug_kol_stat_data)
perloci_boug_kol_fst <- boug_kol_fst$perloc
boug_ken_fst <- basic.stats(boug_ken_stat_data)
perloci_boug_ken_fst <- boug_ken_fst$perloc
boug_nioro_fst <- basic.stats( boug_nioro_stat_data)
perloci_boug_nioro_fst <- boug_nioro_fst$perloc
dang_fala_fst <- basic.stats( dang_fala_stat_data)
perloci_dang_fala_fst <- dang_fala_fst$perloc
dang_ken_fst <- basic.stats( dang_ken_stat_data)
perloci_dang_ken_fst <- dang_ken_fst$perloc
dang_kol_fst <- basic.stats( dang_kol_stat_data)
perloci_dang_kol_fst <- dang_kol_fst$perloc
dang_nioro_fst <- basic.stats(dang_nioro_stat_data)
perloci_dang_nioro_fst <- dang_nioro_fst$perloc
fala_ken_fst <- basic.stats(fala_ken_stat_data)
perloci_fala_ken_fst <- fala_ken_fst$perloc
fala_kol_fst <- basic.stats(fala_kol_stat_data)
perloci_fala_kol_fst <- fala_kol_fst$perloc
fala_nioro_fst <- basic.stats(fala_nioro_stat_data)
perloci_fala_nioro_fst <- fala_nioro_fst$perloc
ken_kol_fst <- basic.stats( ken_kol_stat_data)
perloci_ken_kol_fst <- ken_kol_fst$perloc
ken_nioro_fst <- basic.stats(ken_nioro_stat_data)
perloci_ken_nioro_fst <- ken_nioro_fst$perloc
kol_nioro_fst <- basic.stats(kol_nioro_stat_data)
perloci_kol_nioro_fst <- kol_nioro_fst$perloc

#######################merging all perloci ##########

all_perloci_fst <- rbind(perloci_bko_boug_fst, perloci_bko_dang_fst, perloci_bko_fala_fst, perloci_bko_ken_fst, perloci_bko_kol_fst, perloci_bko_nioro_fst, perloci_boug_dang_fst, perloci_boug_fala_fst, perloci_boug_ken_fst,
                         perloci_boug_kol_fst, perloci_boug_nioro_fst, perloci_dang_fala_fst, perloci_dang_ken_fst, perloci_dang_kol_fst, perloci_dang_nioro_fst,
                         perloci_fala_ken_fst, perloci_fala_kol_fst, perloci_fala_nioro_fst, perloci_ken_kol_fst, perloci_ken_nioro_fst, perloci_kol_nioro_fst)
###########write perloci fst to file ##########
write.table(perloci_bko_boug_fst, file = "bko_boug_fst.txt", sep = "\t", quote = F )
write.table(perloci_bko_dang_fst, file = "bko_dang_fst.txt", sep = "\t", quote = F )
write.table(perloci_bko_fala_fst, file = "bko_fala_fst.txt", sep = "\t", quote = F)
write.table(perloci_bko_ken_fst, file = "bko_ken_fst.txt", sep = "\t", quote = F)
write.table(perloci_bko_kol_fst, file = "bko_kol_fst.txt", sep = "\t", quote = F)
write.table(perloci_bko_nioro_fst, file = "bko_nioro_fst.txt", sep = "\t", quote = F)
write.table(perloci_boug_dang_fst, file = "boug_dang_fst.txt", sep = "\t", quote = F)
write.table(perloci_boug_fala_fst, file = "boug_dang_fst.txt", sep = "\t", quote = F)
write.table(perloci_boug_ken_fst, file = "boug_ken_fst.txt", sep = "\t", quote = F)
write.table(perloci_boug_kol_fst, file = "boug_kol_fst.txt", sep = "\t", quote = F)
write.table(perloci_boug_nioro_fst, file = "boug_nioro_fst.txt", sep = "\t", quote = F)
write.table(perloci_dang_fala_fst, file = "dang_fala_fst.txt", sep = "\t", quote = F)
write.table(perloci_dang_ken_fst, file = "dang_ken_fst.txt", sep = "\t", quote = F)
write.table(perloci_dang_kol_fst, file = "dang_kol_fst.txt", sep = "\t", quote = F)
write.table(perloci_dang_nioro_fst, file = "dang_nioro_fst.txt", sep = "\t", quote = F)
write.table(perloci_fala_ken_fst, file = "fala_ken_fst.txt", sep = "\t", quote = F)
write.table(perloci_fala_kol_fst, file = "fala_kol_fst.txt", sep = "\t", quote = F)
write.table(perloci_fala_nioro_fst, file = "fala_nioro_fst.txt", sep = "\t", quote = F)
write.table(perloci_ken_nioro_fst, file = "ken_nioro_fst.txt", sep = "\t", quote = F)
write.table(perloci_ken_kol_fst, file = "ken_kol_fst.txt", sep = "\t", quote = F)
write.table(perloci_kol_nioro_fst, file = "kol_nioro_fst.txt", sep = "\t", quote = F)

####
############calculate pairwise fst ####################
####fst for alt_data ######
Maj_stat_data <- genind2hierfstat(Maj_genind)
Maj_fst <- pairwise.neifst(Maj_stat_data)
Maj_basic_stat <- basic.stats(Maj_stat_data)
Maj_per_loci <- Maj_basic_stat$perloc
write.table(Maj_fst, file = "mali_pairwise_fst.txt", quote = F, sep = "\t")
write.table()
#################allele freq per population############
makefreq() pop.freq()
allele_freq2 <- pop.freq(Maj_read_data[, -1], diploid = FALSE)
write.table(allele_freq2, file = "allele_freq.txt", col.names = T, row.names = F,quote = F, sep = "\t" )
############read fst files ############
bko_boug_fst <- read.table("bko_boug_fst.txt", header = T, sep ="\t" )
bko_dang_fst <- read.table("bko_dang_fst.txt", header = T, sep ="\t" )
bko_fala_fst <- read.table("bko_fala_fst.txt", header = T, sep ="\t" )
bko_ken_fst <- read.table("bko_ken_fst.txt", header = T, sep ="\t" )
bko_kol_fst <- read.table("bko_kol_fst.txt", header = T, sep ="\t" )
bko_nioro_fst <- read.table("bko_nioro_fst.txt", header = T, sep ="\t" )
boug_dang_fst <- read.table("boug_dang_fst.txt", header = T, sep ="\t" )
boug_fala_fst <- read.table("boug_fala_fst.txt", header = T, sep ="\t" )
boug_ken_fst <- read.table("boug_ken_fst.txt", header = T, sep ="\t" )
boug_kol_fst <- read.table("boug_kol_fst.txt", header = T, sep ="\t" )
boug_nioro_fst <- read.table("boug_nioro_fst.txt", header = T, sep ="\t" )
dang_fala_fst<- read.table("dang_fala_fst.txt", header = T, sep ="\t" ) 
dang_ken_fst <- read.table("dang_ken_fst.txt", header = T, sep ="\t" )
dang_kol_fst <- read.table("dang_kol_fst.txt", header = T, sep ="\t" )
dang_nioro_fst <- read.table("dang_nioro_fst.txt", header = T, sep ="\t" )
fala_ken_fst <- read.table("fala_ken_fst.txt", header = T, sep ="\t" )
fala_kol_fst <- read.table("fala_kol_fst.txt", header = T, sep ="\t" )
fala_nioro_fst <- read.table("fala_nioro_fst.txt", header = T, sep ="\t" )
ken_kol_fst <- read.table("ken_kol_fst.txt", header = T, sep ="\t" )
ken_nioro_fst <- read.table("ken_nioro_fst.txt", header = T, sep ="\t" )
kol_nioro_fst<- read.table("kol_nioro_fst.txt", header = T, sep ="\t" )
##########find outlilers############
outliers=function(){
  # looks at Fst out;ier
  library(data.table)
  library("devtools")
  # first time - must install qvalue first
 # BiocManager::install("qvalue")
#install_github("whitlock/OutFLANK")
  
  
  C=fread("Fst_values.txt",sep="\t")
  hist(C$Fst,50)
  qqnorm(C$Fst);qqline(C$Fst,col=2)
  
  


################fst to zscore #############
perloci_bko_boug_fst <- na.omit(perloci_bko_boug_fst)
perloci_bko_boug_fst$Fst[perloci_bko_boug_fst$Fst<0] = 0
fst.values = perloci_bko_boug_fst$Fst
n = length(perloci_bko_boug_fst$Fst)
#fst.values[fst.values<0] = 0.000001
K = 3
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
z <- z.transform(perloci_bko_boug_fst$Fst)
p_bko_boug  = 2*pnorm(-abs(z.scores))
fdr = fdrtool(z.scores)
q_value <- qvalue(p_bko_boug, fdr.level = 0.01)
####################################
#####visualize fst###############
fst_all_data <- read.table("final_fst_maj.txt", sep = "\t", header = T)
 library(qqman)
pdf("fst_maj_manhattan.pdf",paper = "a4r", width = 10)
par(mfrow=c(4,2),oma=c(3, 3, 2, 2), mar=c(2,1,2,1))
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_boug_fst", snp = "SNP",
    col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst", genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_dang_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_fala_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_ken_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_kol_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="bko_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="boug_dang_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="boug_fala_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="boug_ken_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="boug_kol_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="boug_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="dang_fala_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="dang_ken_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="dang_kol_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="dang_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="fala_ken_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="fala_kol_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="fala_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="ken_kol_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="ken_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
qqman::manhattan(fst_all_data, chr = "CHR", bp = "BP", p ="kol_nioro_fst", snp = "SNP",
                 col = c("blue", "red"), chrlabs = NULL, cex=1.0,cex.axis = 0.8, logp = F,  ylab = "Fst",genomewideline = F, suggestiveline = F, ylim=c(0,0.10))       
dev.off()

fst_all_data <- read.table("final_fst_maj.txt", header = T, sep = "\t")
########HEATMAP#################
##plotting heatmaps
library(heatmap.plus)
library(ComplexHeatmap)
library(circlize)
mat5<-fst_all_data[,4:24]
tmat<-t(mat5)
tmat<-as.matrix(tmat)
myBreaks <- seq(0.0, 0.1, length.out=100)
myCol <- colorRampPalette(c("lightgreen", "red"))(100)
pdf('~/Documents/Aoua_PhD/data/Manuscript1/Figures3/mali_pop_fst.pdf',width = 20, height = 20)
ha = HeatmapAnnotation(FST=anno_boxplot(tmat, height = unit(2, "cm"),
                                        pch = 20,size = unit(2, "mm")))
Heatmap(tmat, name = "FST", column_split=fst_all_data$CHR, heatmap_width = unit(40, "cm"),
        heatmap_height = unit(25, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),show_column_names =F,
        row_names_g = gpar(fontsize = 16,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
        top_annotation = ha,col=colorRamp2(myBreaks, myCol))
dev.off()

 ##############get the data with fst > 0.05 ###############
write.table( sign_1 <- bko_boug_fst %>% filter( P > 0.05), file = "bko_boug_fst0025.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_2 <- bko_dang_fst %>% filter( P > 0.05), file = "bko_dang_fst0025.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_3 <- bko_fala_fst %>% filter( P > 0.05), file = "bko_fala_fst0025.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_4 <- bko_ken_fst %>% filter( P > 0.05), file = "bko_ken_fst0025.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_5 <- bko_kol_fst %>% filter( P > 0.05), file = "bko_kol_fst0025.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_6 <- bko_nioro_fst %>% filter( P > 0.05), file = "bko_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_7 <- boug_dang_fst %>% filter( P > 0.05), file = "boug_dang_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_8 <- boug_fala_fst %>% filter( P > 0.05), file = "boug_fala_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_9 <- boug_ken_fst %>% filter( P > 0.05), file = "boug_ken_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_10 <- boug_kol_fst %>% filter( P > 0.05), file = "boug_kol_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_11 <- boug_nioro_fst %>% filter( P > 0.05), file = "boug_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_12 <- dang_fala_fst %>% filter( P > 0.05), file = "dang_fala_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_13 <- dang_ken_fst %>% filter( P > 0.05), file = "dang_ken_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_14 <- dang_kol_fst %>% filter( P > 0.05), file = "dang_kol_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_15 <- dang_nioro_fst %>% filter( P > 0.05), file = "dang_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_16 <- fala_ken_fst %>% filter( P > 0.05), file = "fala_ken_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_17 <- fala_kol_fst %>% filter( P > 0.05), file = "fala_kol_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_18 <- fala_nioro_fst %>% filter( P > 0.05), file = "fala_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_19 <- ken_kol_fst %>% filter( P > 0.05), file = "ken_kol_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_20 <- ken_nioro_fst %>% filter( P > 0.05), file = "ken_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")
write.table(sign_21 <- kol_nioro_fst %>% filter( P > 0.05), file = "kol_nioro_fst005.txt", row.names = F, quote = F, sep = "\t")

#all_fst_0025 <-
all_fst_005 <-merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(sign_1[, c(1,3,4)], sign_2[, 3:4], by ="SNP", all = TRUE), sign_3[,3:4], by = "SNP", all = TRUE), sign_4[,3:4], by = "SNP", all = TRUE), sign_5[, 3:4], by = "SNP", all = TRUE), sign_6[, 3:4], by = "SNP", all = TRUE),
                                    sign_7[, 3:4], by = "SNP", all = TRUE), sign_8[,3:4], by = "SNP", all = TRUE), sign_9[,3:4], by = "SNP", all = TRUE), sign_10[,3:4], by = "SNP", all = TRUE), sign_11[,3:4], by = "SNP", all = TRUE), sign_12[,3:4], by = "SNP", all = TRUE),
                                    sign_13[,3:4], by = "SNP", all = TRUE), sign_14[,3:4], by = "SNP", all = TRUE),sign_15[,3:4], by = "SNP", all = TRUE),sign_16[,3:4], by = "SNP", all = TRUE),
                                    sign_17[,3:4], by = "SNP", all = TRUE), sign_18[,3:4], by = "SNP", all = TRUE), sign_19[,3:4], by = "SNP", all = TRUE),
                                    sign_20[,3:4], by = "SNP", all = TRUE), sign_21[,3:4], by = "SNP", all = TRUE)




colnames(all_fst_0025)=c("SNP","CHR","bko vs boug",
                        "bko vs dang","bko vs fala","bko vs ken",
                        "bko vs kol","bko vs nioro","boug vs dang",
                        "boug vs fala","boug vs ken","boug vs kol",
                        "boug vs nioro", "dang vs fala","dang vs ken",
                        "dang vs kol","dang vs nioro","fala vs ken",
                        "fala vs kol","fala vs nioro","ken vs kol",
                        "ken vs nioro","kol vs nioro")
write.table(all_fst_005, file = "fst_005.txt", row.names = F, quote = F, sep = "\t")

###add an extra column to identify each pair########
bko_boug_fst$PAIRS <- 1 
bko_dang_fst$PAIRS <- 2
bko_fala_fst$PAIRS <- 3
bko_ken_fst$PAIRS <- 4
bko_kol_fst$PAIRS <- 5
bko_nioro_fst$PAIRS <- 6
boug_dang_fst$PAIRS <- 7
boug_fala_fst$PAIRS <- 8
boug_ken_fst$PAIRS <- 9
boug_kol_fst$PAIRS <- 10
boug_nioro_fst$PAIRS <- 11
dang_fala_fst$PAIRS <- 12
dang_ken_fst$PAIRS <- 13
dang_kol_fst$PAIRS <- 14
dang_nioro_fst$PAIRS <- 15
ken_nioro_fst$PAIRS <- 20 
kol_nioro_fst$PAIRS<- 21
###########adjust p value ############
all_fst_values <- rbind(bko_boug_fst, bko_dang_fst, bko_fala_fst, bko_ken_fst, bko_kol_fst, bko_nioro_fst, boug_dang_fst,boug_fala_fst, boug_ken_fst,
                        boug_kol_fst, boug_nioro_fst, dang_fala_fst, dang_ken_fst, dang_kol_fst, dang_nioro_fst,
                        fala_ken_fst, fala_kol_fst, fala_nioro_fst, ken_kol_fst, ken_nioro_fst, kol_nioro_fst)
#############selesct snp with fst value > 0.05
fst_more_0.05 <- subset(all_fst_values, all_fst_values$P > 0.05)
write.table(fst_more_0.05, file = "fst_more_0.05.txt", col.names = T, row.names = F, quote = F, sep = "\t")
library(data.table)
setDT(fst_more_0.05)[, .N, "SNP" ]
#############adjusted pvalue ###################
############determine significance ########################a
all_fst_values <- na.omit(all_fst_values) 

score_z <- z.transform(all_fst_values$P)
#adjusted.p <- p.adjust(p, method = "fdr",  n = length(p))
fdr = fdrtool(score_z)
all_fst_values$QVAL = fdr$qval
write.table(all_fst_values, file = "all_fst_values.txt", row.names = F, quote = F, sep = "\t")



all_fst_values.qqman <- data.frame(
  CHR = as.integer(factor(all_fst_values$CHR, 
                          levels = unique(all_fst_values$CHR))),
  # chromosomes as integers
  BP = all_fst_values$BP,         # base pairs
  P = all_fst_values$P , # transform back to p-values
  SNP = all_fst_values$SNP     # SNP names
)
library(qqman)
png("Fst_plot.png", width=800, height=600)

manhattan(all_fst_values.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(all_fst_values.qqman$CHR)),
          ylab = "Fst",
          suggestiveline = F,
           #annotatePval = 0.08,
          genomewideline = F,
          cex.axis = 1 ,
          cex.lab=1,
       # highlight = c("Pf3D7_07_v3_404407", "Pf3D7_05_v3_323247"),
          logp = FALSE ,
          ylim =c(0,0.12)
       
       
      
)
dev.off()

manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray10", "gray60"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # The old way to initialize the plot
  # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
  #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  
  # The new way to initialize the plot.
  ## See http://stackoverflow.com/q/23922130/654296
  ## First, define your default arguments
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))),
                   xlab=xlabel, ylab=expression(-log[10](italic(p))))
  ## Next, get a list of ... arguments
  #dotargs <- as.list(match.call())[-1L]
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs)==length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, col=col[1], ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
      icol=icol+1
    }
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
  
  # Highlight snps from a character vector
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col="red", pch=20, ...)) 
  }
  
  # Highlight top SNPs
  if (!is.null(annotatePval)) {
    # extract top SNPs at given p-val
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    # annotate these SNPs
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), 
           textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45), ...)
    }
    else {
      # could try alternative, annotate top SNP of each sig chr
      topHits <- topHits[order(topHits$P),]
      topSNPs <- NULL
      
      for (i in unique(topHits$CHR)) {
        
        chrSNPs <- topHits[topHits$CHR == i,]
        topSNPs <- rbind(topSNPs, chrSNPs[1,])
        
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs =topSNPs$SNP , cex = 0.5, ...)
    }
  }  
  par(xpd = FALSE)
}
