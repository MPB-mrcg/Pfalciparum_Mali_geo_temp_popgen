setwd("~/Documents/Aoua_PhD/data/Selection/")
library(devtools)
devtools::install_github("bahlolab/isoRelate")
BiocManager::install("rehh")
library(rehh)
library(dplyr)
library(qqman)
library(fdrtool)
#####maj data to be splited by population
Maj_genotype_data <- read.table("Maj_hap_Data.txt", header = T, sep = "\t")
map_file <- read.table("map_file.txt", header = F, sep = "\t")
######spliting data by population#########
bko_hap <- Maj_genotype_data %>% filter( Locations == "Bamako")
#bko_hap <- Maj_genotype_data[grepl("Bamako", Maj_genotype_data$Locations),]
bko_hap <- bko_hap[, -1]
boug_hap <- Maj_genotype_data %>% filter( Locations == "Bougoula-Hameau")
boug_hap <- boug_hap[,-1]
dang_hap <- Maj_genotype_data %>% filter( Locations == "Dangassa")
dang_hap <- dang_hap[, -1]
fala_hap <- Maj_genotype_data %>% filter( Locations == "Faladje")
fala_hap <- fala_hap[, -1]
ken_hap <- Maj_genotype_data %>% filter( Locations == "Kenieroba")
ken_hap <- ken_hap[, -1]
kol_hap <- Maj_genotype_data %>% filter( Locations == "Kolle")
kol_hap <- kol_hap[, -1]
nioro_hap <- Maj_genotype_data %>% filter( Locations == "nioro")
nioro_hap <- nioro_hap[, -1]

##############split by chromosome###################
write.table(nioro_hap[, grepl("Pf3D7_01_v3", colnames(nioro_hap))], "Nioro/nioro_chr1.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_02_v3", colnames(nioro_hap))], "Nioro/nioro_chr2.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_03_v3", colnames(nioro_hap))], "Nioro/nioro_chr3.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_04_v3", colnames(nioro_hap))], "Nioro/nioro_chr4.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_05_v3", colnames(nioro_hap))], "Nioro/nioro_chr5.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_06_v3", colnames(nioro_hap))], "Nioro/nioro_chr6.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_07_v3", colnames(nioro_hap))], "Nioro/nioro_chr7.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_08_v3", colnames(nioro_hap))], "Nioro/nioro_chr8.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_09_v3", colnames(nioro_hap))], "Nioro/nioro_chr9.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_10_v3", colnames(nioro_hap))], "Nioro/nioro_chr10.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_11_v3", colnames(nioro_hap))], "Nioro/nioro_chr11.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_12_v3", colnames(nioro_hap))], "Nioro/nioro_chr12.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_13_v3", colnames(nioro_hap))], "Nioro/nioro_chr13.hap",row.names = T, col.names=F, sep='\t',quote = F)
write.table(nioro_hap[, grepl("Pf3D7_14_v3", colnames(nioro_hap))], "Nioro/nioro_chr14.hap",row.names = T, col.names=F, sep='\t',quote = F)

###################computing iHS, Rsb and XP-EHH within each population.################
for(i in 1:14){ 
  hap_file=paste("Bamako/bko_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80 )
  res<-scan_hh(data)
  if(i==1){wg.res.bko<-res}else{wg.res.bko<-rbind(wg.res.bko,res)}
  bko.ihs<-ihh2ihs(wg.res.bko)
}
manhattanplot(bko.ihs,pval = T)

for(i in 1:14){ 
  hap_file=paste("Bougoula/boug_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80 )
  res<-scan_hh(data)
  if(i==1){wg.res.boug<-res}else{wg.res.boug<-rbind(wg.res.boug,res)}
  boug.ihs<-ihh2ihs(wg.res.boug)
}
manhattanplot(boug.ihs, pval = T)
for(i in 1:14){
  hap_file=paste("Dangassa/dang_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80)
  res<-scan_hh(data)
  if(i==1){wg.res.dang<-res}else{wg.res.dang<-rbind(wg.res.dang,res)}
  dang.ihs<-ihh2ihs(wg.res.dang)
}
manhattanplot(dang.ihs,pval = T)

for(i in 1:14){ 
  hap_file=paste("Faladje/fala_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80)
  res<-scan_hh(data)
  if(i==1){wg.res.fala<-res}else{wg.res.fala<-rbind(wg.res.fala,res)}
  fala.ihs<-ihh2ihs(wg.res.fala)
}
manhattanplot(fala.ihs, pval = T)
for(i in 1:14){
  hap_file=paste("Kenieroba/ken_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80)
  res<-scan_hh(data)
  if(i==1){wg.res.ken<-res}else{wg.res.ken<-rbind(wg.res.ken,res)}
  ken.ihs<-ihh2ihs(wg.res.ken)
}
manhattanplot(ken.ihs,pval = T)

for(i in 1:14){ 
  hap_file=paste("Kolle/kol_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80)
  res<-scan_hh(data)
  if(i==1){wg.res.kol<-res}else{wg.res.kol<-rbind(wg.res.kol,res)}
  kol.ihs<-ihh2ihs(wg.res.kol)
}
manhattanplot(kol.ihs, pval = T)

for(i in 1:14){ 
  hap_file=paste("Nioro/nioro_chr",i,".hap",sep="")
  data<-data2haplohh(hap_file,"map_file.txt",chr.name=i, min_perc_geno.hap=80, min_perc_geno.mrk=80)
  res<-scan_hh(data)
  if(i==1){wg.res.nioro<-res}else{wg.res.nioro<-rbind(wg.res.nioro,res)}
  nioro.ihs<-ihh2ihs(wg.res.nioro)
}
manhattanplot(nioro.ihs, pval = T)
##########addlocations#############
bko.ihs$ihs$Location<- "Bamako"
boug.ihs$ihs$Location <- "Bougoula-Hameau"
dang.ihs$ihs$Location <- "Dangassa"
fala.ihs$ihs$Location <- "Faladje"
ken.ihs$ihs$Location <- "Kenioroba"
kol.ihs$ihs$Location <- "Kolle"
nioro.ihs$ihs$Location <- "Nioro"
##########remove Nas ##############
bko.ihs$ihs <- na.omit(bko.ihs$ihs)
boug.ihs$ihs <- na.omit(boug.ihs$ihs)
fala.ihs$ihs <- na.omit(fala.ihs$ihs)
dang.ihs$ihs <- na.omit(dang.ihs$ihs)
ken.ihs$ihs <- na.omit(ken.ihs$ihs)
kol.ihs$ihs <- na.omit(kol.ihs$ihs)
nioro.ihs$ihs <- na.omit(nioro.ihs$ihs)
###########adjusted pvalue for each location #########
nioro.ihs$ihs$P <- 10^(-nioro.ihs$ihs$LOGPVALUE)
fdr = fdrtool(nioro.ihs$ihs$P, statistic = "pvalue")
nioro.ihs$ihs$P <- fdr$qval
###########adjust pvalue using qvalues##########
nioro.ihs$ihs$P <- 10^(-nioro.ihs$ihs$LOGPVALUE)
p <- nioro.ihs$ihs$P
nioro_qval <- qvalue(p, fdr.level = 0.01)
nioro.ihs$ihs$qval <- nioro_qval$qvalues
ihs_nioro <-nioro.ihs$ihs
ihs_nioro$SNP <- rownames(ihs_nioro)
write.table(ihs_nioro, file = "~/Documents/Aoua_PhD/data/Selection/IHS/ihs_NIORO.txt", quote = F, row.names = F, sep ="\t")
write.table(nioro_ihs_sign <- subset(ihs_nioro, ihs_nioro$qval < 0.00001), file = "~/Documents/Aoua_PhD/data/Selection/IHS/snp_of_interest_ihs_NIORO.txt", quote = F, row.names = F, sep = "\t")
all_ihs_sign <- subset(all_ihs, all_ihs$qval < 0.00001)
 ######################manhattan plots ###############
nioro_ihs.qqman <- data.frame(
  CHR = as.integer(factor(ihs_nioro$CHR),
                          levels = unique(ihs_nioro$CHR)),
  # chromosomes as integers
  BP = ihs_nioro$POSITION,         # base pairs
  P = ihs_nioro$qval,  # transform back to p-values
  SNP = ihs_nioro$SNP      # SNP names
)
library(qqman) 
pdf("~/Documents/Aoua_PhD/data/Manuscript1/all_pop_ihs.pdf",paper = "a4r", width = 10)
png("~/Documents/Aoua_PhD/data/Manuscript1/pop_ihs.png", width=1000, height=700)
par(mfrow=c(7,2),oma=c(4, 4, 3, 3), mar=c(1,1,1,1))
manhattan(bko_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_ihs.qqman$CHR)),
          suggestiveline = F,
         #annotatePval = 0.00001,
         genomewideline = F,
         highlight = bko_ihs_sign$SNP 
         #annotateTop = T
      )
manhattan(boug_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight =  boug_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(dang_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = dang_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(fala_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(ken_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(ken_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = ken_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(kol_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(kol_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = kol_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(nioro_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(nioro_ihs.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = nioro_ihs_sign$SNP,
          ylim =c(0,15))
manhattan(all_ihs.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(all_ihs.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = all_ihs_sign$SNP,
          ylim =c(0,15))
dev.off()

####################combine all ihs values #######
all_ihs <- rbind(ihs_bko,ihs_boug, ihs_dang, ihs_fala, ihs_ken, ihs_kol, ihs_nioro)
#all_ihs <- na.omit(all_ihs)
all_ihs$P <- 10^(-all_ihs$LOGPVALUE)
p <- all_ihs$P
adjusted.p <- p.adjust(p, method = "fdr",  n = length(p))
fdr = fdrtool(p, statistic = "pvalue")
all_ihs$QVAL <- fdr$qval
write.table(all_ihs, file = "all_ihs.txt", row.names = F, quote = F, sep = "\t")
ihs <- wgscan.ihs.cgu$ihs
# create new data frame
#all_ihs$P <- NULL 
all_ihs.qqman <- data.frame(
  CHR = as.integer(factor(all_ihs$CHR, 
                          levels = unique(all_ihs$CHR))),
  # chromosomes as integers
  BP = all_ihs$POSITION,         # base pairs
  P = all_ihs$qval , # transform back to p-values
  SNP = all_ihs$SNP       # SNP names
  )
library(qqman)
pdf("./manhattanplot_ihs_new")
manhattan(all_ihs.qqman,
          col = c("lightgreen", "red"),
          chrlabs = as.character(unique(all_ihs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          #highlight = "Pf3D7_7_v3_465508",
          ylab="-log10(q)",
          genomewideline = F,
          ylim =c(0,15))
dev.off()
bko_ihs <- bko_ihs.qqman
boug_ihs <- boug_ihs.qqman
dang_ihs <- dang_ihs.qqman
fala_ihs <- fala_ihs.qqman
ken_ihs <- ken_ihs.qqman
kol_ihs <- kol_ihs.qqman
nioro_ihs <- nioro_ihs.qqman
 ########heatmap###########
#######   
ihs_merged <-merge(merge(merge(merge(merge(merge(bko_ihs[, c(1,3,4)], boug_ihs[, c(3,4)], by = "SNP", all = TRUE), dang_ihs[,c(3,4)], by ="SNP",all = TRUE ), fala_ihs[,c(3,4)], by = "SNP",all = TRUE ), ken_ihs[, c(3,4)], by = "SNP", all = TRUE), kol_ihs[, c(3,4)], by = "SNP", all = TRUE), nioro_ihs[, c(3,4)], by="SNP", all = TRUE)
  colnames(ihs_merged) <- c("SNP", "CHR", "bko", "boug","dang", "fala", "ken", "kol","nioro")       
  write.table(ihs_merged, file = "ihs_merged.txt", row.names = F, quote = F, sep = "\t")  
  ##plotting heatmaps
  library(heatmap.plus)
  library(ComplexHeatmap)
  library(circlize)
  mat5<-ihs_merged_new[,3:9]
  tmat<-t(mat5)
  tmat<-as.matrix(tmat)
  myBreaks <- seq(0.0 , 5,  length.out =100)
  myCol <- colorRampPalette(c("lightgreen","red"))(100)
  pdf('ihs_mali_new.pdf',width = 20, height = 20)
  ha = HeatmapAnnotation(log10qval=anno_boxplot(tmat, height = unit(2, "cm"),
                                            pch = 20, size = unit(2, "mm")))
  Heatmap(tmat, name = "-log10(qval)", column_split=ihs_merged_new$CHR, heatmap_width = unit(40, "cm"),
          heatmap_height = unit(30, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
          column_title_gp = gpar(fontsize = 16, fontface = "bold"),show_column_names =F,
          row_names_g = gpar(fontsize = 16,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
          top_annotation = ha,col=colorRamp2(myBreaks, myCol))
  dev.off() 
  ##############significant snps############
write.table(bko_ihs_sign <- subset(bko_ihs, bko_ihs$QVAL < 0.00001) ,file = "~/Documents/Aoua_PhD/data/Selection/IHS/bko_ihs_snp.txt", sep= "\t", row.names = F )
  write.table(boug_ihs_sign <- subset(boug_ihs, boug_ihs$QVAL < 0.00001), file = "~/Documents/Aoua_PhD/data/Selection/IHS/boug_ihs_snp.txt", sep= "\t", row.names = F )
  
  ####manhattan plot ihs, xpehh, rsb #############
pdf("./ihs_mali_pop.pdf",paper = "a4r", width = 10)
par(mfrow=c(4,2),oma=c(3, 3, 2, 2), mar=c(2,1,2,1))
manhattanplot(bko.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(boug.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(dang.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(fala.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(ken.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(kol.ihs,  pval = T, ylim=c(0, 10))
manhattanplot(nioro.ihs,  pval = T, ylim=c(0, 10))
manha
dev.off()
nioro_ihs_can_regions <- calc_candidate_regions(nioro.ihs,
                                 threshold = 4,
                                 pval = TRUE,
                                 #window_size = 1E6,
                                 #overlap = 1E5,
                                 min_n_extr_mrk = 2)
all_ihs_can_regions <- calc_candidate_regions(all_ihs,
                                                threshold = 5,
                                                pval = TRUE,
                                                #window_size = 1E6,
                                                #overlap = 1E5,
                                                min_n_extr_mrk = 2)

###### Rsb pairwise population statistic ###################
ken_kol.rsb<-ines2rsb(wg.res.ken,wg.res.kol)
ken_kol.rsb$SNP <- rownames(ken_kol.rsb)
ken_kol.rsb <- na.omit(ken_kol.rsb)
#manhattanplot(fala_ken.rsb, pval = T)
write.table(ken_kol.rsb,"RSB/ken_kol.rsb.txt",sep = "\t", row.names = F)
ken_nioro.rsb<-ines2rsb(wg.res.ken,wg.res.nioro)
ken_nioro.rsb$SNP <- rownames(ken_nioro.rsb)
ken_nioro.rsb <- na.omit(ken_nioro.rsb)
#manhattanplot(fala_kol.rsb, pval = T)
write.table(ken_nioro.rsb,"RSB/ken_nioro.rsb.txt",sep = "\t", row.names = F)
bko_dang.rsb<-ines2rsb(wg.res.bko,wg.res.dang)
bko_dang.rsb$SNP <- rownames(bko_dang.rsb)
bko_dang.rsb <- na.omit(bko_dang.rsb)
#manhattanplot(fala_nioro.rsb, pval = T)
write.table(kol_nioro.rsb,"RSB/kol_nioro.rsb.txt",sep = "\t", row.names = F)
dang_nioro.rsb<-ines2rsb(wg.res.dang,wg.res.nioro)
dang_nioro.rsb$SNP <- rownames(dang_nioro.rsb)
dang_nioro.rsb <- na.omit(dang_nioro.rsb)
#manhattanplot(ken_kol.rsb, pval = T)
write.table(dang_nioro.rsb,"RSB/dang_nioro.rsb.txt",sep = "\t", row.names = F)
boug_nioro.rsb<-ines2rsb(wg.res.boug,wg.res.nioro)
boug_nioro.rsb$SNP <- rownames(boug_nioro.rsb)
boug_nioro.rsb <- na.omit(boug_nioro.rsb)
#manhattanplot(ken_nioro.rsb, pval = T)
write.table(boug_nioro.rsb,"RSB/boug_nioro.rsb.txt",sep = "\t", row.names = F)
bko_nioro.rsb<-ines2rsb(wg.res.bko,wg.res.nioro)
bko_nioro.rsb$SNP <- rownames(bko_nioro.rsb)
bko_nioro.rsb <- na.omit(bko_nioro.rsb)
#manhattanplot(kol_nioro.rsb, pval = T)
write.table(bko_nioro.rsb,"RSB/bko_nioro.rsb.txt",sep = "\t", row.names = F)

fala_nioro.rsb<-ines2rsb(wg.res.fala,wg.res.nioro)
manhattanplot(fala_nioro.rsb, pval = T)
write.table(fala_nioro.rsb,"RSB/fala_nioro.rsb.txt",sep = "\t", row.names = F)
ken_kol.rsb<-ines2rsb(wg.res.ken,wg.res.kol)
manhattanplot(ken_kol.rsb, pval = T)
write.table(ken_kol.rsb,"RSB/ken_kol.rsb.txt",sep = "\t", row.names = F)
ken_nioro.rsb<-ines2rsb(wg.res.ken,wg.res.nioro)
manhattanplot(ken_nioro.rsb, pval = T)
write.table(ken_nioro.rsb,"RSB/ken_nioro.rsb.txt",sep = "\t", row.names = F)
kol_nioro.rsb<-ines2rsb(wg.res.kol,wg.res.nioro)
manhattanplot(kol_nioro.rsb, pval = T)
write.table(kol_nioro.rsb,"RSB/kol_nioro.rsb.txt",sep = "\t", row.names = F)
###########adjust rsb pvalue using qvalues##########
bko_dang.rsb$P <- 10^(-bko_dang.rsb$LOGPVALUE)
p <- bko_dang.rsb$P
bko_dang_qval <- qvalue(p, fdr.level = 0.01)
bko_dang.rsb$qval <- bko_dang_qval$qvalues
write.table(bko_dang.rsb, file = "~/Documents/Aoua_PhD/data/Selection/RSB/bko_dang.rsb.txt", quote = F, row.names = F, sep ="\t")
write.table(bko_dang.rsb_sign <- subset(bko_dang.rsb, bko_dang.rsb$qval < 0.00001), file = "~/Documents/Aoua_PhD/data/Selection/RSB/snp_of_interest_bko_dang.txt", quote = F, row.names = F, sep = "\t")

######################manhattan plots ###############
bko_dang.rsb.qqman <- data.frame(
  CHR = as.integer(factor(bko_dang.rsb$CHR),
                   levels = unique(bko_dang.rsb$CHR)),
  # chromosomes as integers
  BP = bko_dang.rsb$POSITION,         # base pairs
  P =bko_dang.rsb$qval,  # transform back to p-values
  SNP =bko_dang.rsb$SNP      # SNP names
)

######################rsb manhattan plot #######################################
all_rsb_values <- rbind(bko_boug.rsb, bko_dang.rsb, bko_fala.rsb, bko_ken.rsb, bko_kol.rsb, bko_nioro.rsb, boug_dang.rsb, boug_fala.rsb,
                        boug_ken.rsb, boug_kol.rsb, boug_nioro.rsb, dang_fala.rsb, dang_ken.rsb, dang_kol.rsb, dang_nioro.rsb, fala_ken.rsb, fala_kol.rsb, fala_nioro.rsb, ken_kol.rsb, ken_nioro.rsb,kol_nioro.rsb)
all_rsb_values.qqman <- data.frame(
CHR = as.integer(factor(all_rsb_values$CHR, 
                        levels = unique(all_rsb_values$CHR))),
# chromosomes as integers
BP = all_rsb_values$POSITION,         # base pairs
P = all_rsb_values$RSB, # rsb values
SNP = all_rsb_values$SNP     # SNP names
)
library(qqman)
png("rsb_plot.png", width=800, height=600)

manhattan(all_rsb_values.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(all_rsb_values.qqman$CHR)),
          ylab ="Mean RSB",
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          cex.axis = 1 ,
          cex.lab=1,
          # highlight = c("Pf3D7_07_v3_404407", "Pf3D7_05_v3_323247"),
          logp = FALSE ,
          ylim =c(0,20)
          
          
          
)
dev.off()
write.table(all_rsb_values, file = "~/Documents/Aoua_PhD/data/Selection/RSB/all_rsb_values.txt", quote = F, row.names = F, sep ="\t")
top_all_rsb_values <-  filter(allRSB > 5) %>% distinct(all_rsb_values, SNP, .keep_all = TRUE)
Top_rsb <- all_rsb_values %>% filter(RSB >5) %>% distinct(SNP, .keep_all = TRUE)
########Get gene Names####
gene_names <- fread(file = "~/Documents/Aoua_PhD/data/vcf-6.0/filtering/All_sample_snpeff_new.txt")
gene_data <- left_join(gene_names, Top_rsb, by = c("CHR", "POSITION"))

write.table(Top_rsb, file = "~/Documents/Aoua_PhD/data/Selection/RSB/top_rsb_values.txt", quote = F, row.names = F, sep ="\t")
write.table(gene_data, file = "~/Documents/Aoua_PhD/data/Selection/RSB/top_rsb_gene_info.txt", quote = F, row.names = F, sep ="\t")

########### rsb Manhattan plots #########################
pdf("~/Documents/Aoua_PhD/data/Manuscript1/all_pop_rsb.pdf",paper = "a4r", width = 10)
par(mfrow=c(8,2),oma=c(3, 3, 2, 2), mar=c(1,1,1,1))
manhattan(bko_boug.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_boug.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
        highlight = bko_boug_rsb_sign$SNP,
          ylim =c(0,20)
)
manhattan(bko_dang.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_dang.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_dang.rsb_sign$SNP,
          ylim =c(0,20))
manhattan(bko_fala.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_fala.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_fala_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(bko_ken.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_ken.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_ken_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(bko_kol.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_kol.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
        highlight = bko_kol_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(bko_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = bko_nioro_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(boug_dang.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_dang.rsb.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = boug_dang_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(boug_fala.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_fala.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_fala_rsb_sign$SNP,
          ylim =c(0,20)
)
manhattan(boug_ken.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_ken.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_ken_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(boug_kol.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_kol.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_kol_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(boug_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_nioro_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(dang_fala.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_fala.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
        highlight =   dang_fala_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(dang_ken.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_ken.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = dang_ken.rsb_sign$SNP,
          ylim =c(0,20))
manhattan(dang_kol.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_kol.rsb.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = dang_kol.rsb_sign$SNP,
          ylim =c(0,20))
manhattan(dang_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = dang_nioro.rsb_sign$SNP,
          ylim =c(0,20)
)
manhattan(fala_ken.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_ken.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_ken_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(fala_kol.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_kol.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_kol_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(fala_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_nioro_rsb_sign$SNP,
          ylim =c(0,20))
manhattan(ken_kol.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(ken_kol.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = ken_kol.rsb_sign$SNP,
          ylim =c(0,20))
manhattan(ken_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(ken_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = ken_nioro.rsb_sign$SNP,
          ylim =c(0,20))
manhattan(kol_nioro.rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(kol_nioro.rsb.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = kol_nioro.rsb_sign$SNP,
          ylim =c(0,20))


dev.off()
bko_boug.xpehh
###########merge all population###########
ken_kol.rsb$LOGQVAL <- -log10(ken_kol.rsb$qval)
kol_nioro.rsb$LOGQVAL <- -log10(kol_nioro.rsb$qval)
fala_nioro.rsb$LOGQVAL <- -log10(fala_nioro.rsb$qval)
#dang_nioro.rsb$LOGQVAL <- -log10(dang_nioro.rsb$qval)
#boug_nioro.rsb$LOGQVAL <- -log10(boug_nioro.rsb$qval)
#bko_nioro.rsb$LOGQVAL <- -log10(bko_nioro.rsb$qval)
bko_boug.rsb$
merge1 <- merge(bko_boug.rsb[, c(1,5,8)], bko_dang.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge2 <- merge(merge1, bko_fala.rsb[,c(5,8)], by ="SNP", all = TRUE)
merge3 <- merge(merge2, bko_ken.rsb[,c(5,8)], by = "SNP", all = TRUE)
merge4 <- merge(merge3, bko_kol.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge5 <- merge(merge4, bko_nioro.rsb[, c(5,8)], by = "SNP",all = TRUE)
merge6 <- merge(merge5, boug_dang.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge7 <- merge(merge6,boug_fala.rsb[, c(5,8)], by= "SNP", all = TRUE)
merge8 <- merge(merge7,boug_ken.rsb[, c(5,8)], by = "SNP",all = TRUE)
merge9 <- merge(merge8, boug_kol.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge10 <- merge(merge9, boug_nioro.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge11 <- merge(merge10, dang_fala.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge12 <- merge(merge11, dang_ken.rsb[, c(5,8)], by ="SNP", all = TRUE)
merge13 <- merge(merge12,dang_kol.rsb[, c(5,8)], by="SNP", all = TRUE)
merge14 <- merge(merge13, dang_nioro.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge15 <- merge(merge14, fala_ken.rsb[, c(5,8)], by = "SNP", all = TRUE )
merge16 <- merge(merge15, fala_kol.rsb[, c(5,8)],by = "SNP", all = TRUE)
merge17 <- merge(merge16, fala_nioro.rsb[, c(5,8)],by = "SNP", all = TRUE)
merge18 <- merge(merge17, ken_kol.rsb[, c(5,8)] ,by = "SNP", all = TRUE)
merge19 <- merge(merge18, ken_nioro.rsb[, c(5,8)], by = "SNP", all = TRUE)
merge20 <- merge(merge19,kol_nioro.rsb[, c(5,8)], by = "SNP", all = TRUE) 

colnames(merge20)=c("SNP","CHR","bko vs boug", 
                    "bko vs dang","bko vs fala","bko vs ken",
                    "bko vs kol","bko vs nioro","boug vs dang",
                    "boug vs fala","boug vs ken","boug vs kol",
                    "boug vs nioro", "dang vs fala","dang vs ken",
                    "dang vs kol","dang vs nioro","fala vs ken",
                    "fala vs kol","fala vs nioro","ken vs kol",
                    "ken vs nioro","kol vs nioro")
rsb_pop_mali <- na.omit(merge20 )
write.table(rsb_pop_mali, file = "~/Documents/Aoua_PhD/data/Selection/RSB/rsb_pop_mali.txt", quote = F, row.names = F, sep ="\t")
###############
##plotting heatmaps
library(heatmap.plus)
library(ComplexHeatmap)
library(circlize)
mat5<-rsb_pop_mali[ ,3:23]
tmat<-t(mat5)
tmat<-as.matrix(tmat)
myBreaks <- seq(0.0 , 15,  length.out =100)
myCol <- colorRampPalette(c("lightgreen","red"))(100)
pdf('~/Documents/Aoua_PhD/data/Manuscript1/Figures3/mali_pop_rsb_new.pdf',width = 20, height = 20)
ha = HeatmapAnnotation("-log10(q-value)"=anno_boxplot(tmat, height = unit(3, "cm"),
                               pch = 20,size = unit(1.5, "mm"), axis = TRUE,
                       axis_param = list(gp=gpar(fontsize=15))))

Heatmap(tmat, name = "-log10(q-value)", column_split=rsb_pop_mali$CHR, heatmap_width = unit(40, "cm"),
        heatmap_height = unit(25, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),show_column_names =F,
        row_names_g = gpar(fontsize = 16,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
        heatmap_legend_param = list( title_gp=gpar(fontsize=16, fontface="bold"),legend_width=unit(8,"cm"),  legend_height=unit(10,"cm"), labels_gp=gpar(fontsize=15, fontface="bold") ),
        top_annotation = ha,col=colorRamp2(myBreaks, myCol))
dev.off()

################ computing XP-EHH pairwise population statistic#############

bko_dang.xpehh<-ies2xpehh(wg.res.bko,wg.res.dang)
#ayout(matrix(1:2,2,1))
#manhattanplot(boug_dang.xpehh, pval = T)
#write.table(ken_kol.xpehh,"xpehh/ken_kol.xpehh.txt",sep = "\t", row.names = F)
fala_kol.xpehh<-ies2xpehh(wg.res.fala,wg.res.kol)
#layout(matrix(1:2,2,1))
#xpehhplot(bko_dang.xpehh)
#write.table(boug_dang.xpehh,"xpehh/boug_dang.xpehh.txt",sep = "\t")
fala_nioro.xpehh<-ies2xpehh(wg.res.fala,wg.res.nioro)
#layout(matrix(1:2,2,1))
#xpehhplot(bko_fala.xpehh)
#write.table(fala_nioro.xpehh,"xpehh/fala_nioro.xpehh.txt",sep = "\t")
ken_kol.xpehh<-ies2xpehh(wg.res.ken,wg.res.kol)
#layout(matrix(1:2,2,1))
#xpehhplot(bko_ken.xpehh)
#write.table(dang_nioro.xpehh,"xpehh/dang_nioro.xpehh.txt",sep = "\t")
ken_nioro.xpehh<-ies2xpehh(wg.res.ken,wg.res.nioro)
#layout(matrix(1:2,2,1))
#xpehhplot(bko_kol.xpehh)
#write.table(boug_nioro.xpehh,"xpehh/boug_nioro.xpehh.txt",sep = "\t")
kol_nioro.xpehh<-ies2xpehh(wg.res.kol,wg.res.nioro)
#layout(matrix(1:2,2,1))
#xpehhplot(bko_nioro.xpehh)
#write.table(bko_nioro.xpehh,"bko_nioro.xpehh.txt",sep = "\t")
##########################pvalue to qvalue ####################
bko_nioro.xpehh$SNP <- rownames(bko_nioro.xpehh)
bko_nioro.xpehh <- na.omit(bko_nioro.xpehh)
bko_dang.xpehh$P <- 10^(-bko_dang.xpehh$LOGPVALUE)
p <-bko_nioro.xpehh$P
bko_nioroxpehh_qval <- qvalue(p, fdr.level = 0.01)
bko_nioro.xpehh$qval <- bko_nioroxpehh_qval$qvalues
write.table(bko_nioro.xpehh, file = "~/Documents/Aoua_PhD/data/Selection/xpehh/bko_nioro.xpehh.txt", quote = F, row.names = F, sep ="\t")
write.table(bko_nioro.xpehh_sign <- subset(bko_nioro.xpehh, bko_nioro.xpehh$qval < 0.00001), file = "~/Documents/Aoua_PhD/data/Selection/xpehh/snp_of_interest_bko_bioro.txt", quote = F, row.names = F, sep = "\t")

######################manhattan plots ###############
bko_dang.xpehh.qqman <- data.frame(
  CHR = as.integer(factor(bko_dang.xpehh$CHR),
                   levels = unique(bko_dang.xpehh$CHR)),
  # chromosomes as integers
  BP = bko_dang.xpehh$POSITION,         # base pairs
  P =bko_dang.xpehh$qval,  # transform back to p-values
  SNP =bko_dang.xpehh$SNP      # SNP names
)
pdf("~/Documents/Aoua_PhD/data/Manuscript1/all_pop_xpehh.pdf",paper = "a4r", width = 10)
png("~/Documents/Aoua_PhD/data/Manuscript1/pop_xpehh.png", width=700, height=500)
par(mfrow=c(8,2),oma=c(3, 3, 2, 2), mar=c(1,1,1,1))
manhattan(bko_boug.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_boug.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_boug.xpehh_sign$SNP,
          ylim =c(0,20)
)
manhattan(bko_dang.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_dang.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_dang.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(bko_fala.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_fala.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_fala.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(bko_ken.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_ken.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_ken.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(bko_kol.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_kol.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = bko_kol.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(bko_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(bko_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = bko_nioro.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(boug_dang.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_dang.xpehh.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = boug_dang.xpehh$SNP,
          ylim =c(0,20))
manhattan(boug_fala.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_fala.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_fala.xpehh_sign$SNP,
          ylim =c(0,20)
)
manhattan(boug_ken.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_ken.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_ken.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(boug_kol.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_kol.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_kol.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(boug_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(boug_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = boug_nioro.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(dang_fala.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_fala.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = dang_fala.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(dang_ken.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_ken.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = dang_ken.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(dang_kol.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_kol.xpehh.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = dang_kol.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(dang_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(dang_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = dang_nioro.xpehh_sign$SNP,
          ylim =c(0,20)
)
manhattan(fala_ken.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_ken.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_ken.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(fala_kol.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_kol.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_kol.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(fala_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(fala_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = fala_nioro.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(ken_kol.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(ken_kol.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.01,
          genomewideline = F,
          highlight = ken_kol.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(ken_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(ken_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.00001,
          genomewideline = F,
          highlight = ken_nioro.xpehh_sign$SNP,
          ylim =c(0,20))
manhattan(kol_nioro.xpehh.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(kol_nioro.xpehh.qqman$CHR)),
          suggestiveline = F,
          genomewideline = F,
          highlight = kol_nioro.xpehh_sign$SNP,
          ylim =c(0,20))


dev.off()
###########merge all population###########
bko_dang.xpehh$LOGQVAL <- -log10(bko_dang.xpehh$qval)
ken_nioro.xpehh$LOGQVAL <- -log10(ken_nioro.xpehh$qval)
kol_nioro.xpehh$LOGQVAL <- -log10(kol_nioro.xpehh$qval)
dang_nioro.xpehh$LOGQVAL <- -log10(dang_nioro.xpehh$qval)
boug_nioro.xpehh$LOGQVAL <- -log10(boug_nioro.xpehh$qval)
bko_nioro.xpehh$LOGQVAL <- -log10(bko_nioro.xpehh$qval)
bko_boug.xpehh$
merge1 <- merge(bko_boug.xpehh[, c(1,5,8)], bko_dang.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge2 <- merge(merge1, bko_fala.xpehh[,c(5,8)], by ="SNP", all = TRUE)
merge3 <- merge(merge2, bko_ken.xpehh[,c(5,8)], by = "SNP", all = TRUE)
merge4 <- merge(merge3, bko_kol.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge5 <- merge(merge4, bko_nioro.xpehh[, c(5,8)], by = "SNP",all = TRUE)
merge6 <- merge(merge5, boug_dang.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge7 <- merge(merge6,boug_fala.xpehh[, c(5,8)], by= "SNP", all = TRUE)
merge8 <- merge(merge7,boug_ken.xpehh[, c(5,8)], by = "SNP",all = TRUE)
merge9 <- merge(merge8, boug_kol.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge10 <- merge(merge9, boug_nioro.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge11 <- merge(merge10, dang_fala.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge12 <- merge(merge11, dang_ken.xpehh[, c(5,8)], by ="SNP", all = TRUE)
merge13 <- merge(merge12,dang_kol.xpehh[, c(5,8)], by="SNP", all = TRUE)
merge14 <- merge(merge13, dang_nioro.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge15 <- merge(merge14, fala_ken.xpehh[, c(5,8)], by = "SNP", all = TRUE )
merge16 <- merge(merge15, fala_kol.xpehh[, c(5,8)],by = "SNP", all = TRUE)
merge17 <- merge(merge16, fala_nioro.xpehh[, c(5,8)],by = "SNP", all = TRUE)
merge18 <- merge(merge17, ken_kol.xpehh[, c(5,8)] ,by = "SNP", all = TRUE)
merge19 <- merge(merge18, ken_nioro.xpehh[, c(5,8)], by = "SNP", all = TRUE)
merge20 <- merge(merge19,kol_nioro.xpehh[, c(5,8)], by = "SNP", all = TRUE) 

colnames(merge20)=c("SNP","CHR","bko vs boug", 
                    "bko vs dang","bko vs fala","bko vs ken",
                    "bko vs kol","bko vs nioro","boug vs dang",
                    "boug vs fala","boug vs ken","boug vs kol",
                    "boug vs nioro", "dang vs fala","dang vs ken",
                    "dang vs kol","dang vs nioro","fala vs ken",
                    "fala vs kol","fala vs nioro","ken vs kol",
                    "ken vs nioro","kol vs nioro")
xpehh_pop_mali <- na.omit(merge20 )
###############
##plotting heatmaps
library(heatmap.plus)
library(ComplexHeatmap)
library(circlize)
mat5<-xpehh_pop_mali[ ,3:23]
tmat<-t(mat5)
tmat<-as.matrix(tmat)
myBreaks <- seq(0.0 , 20,  length.out =100)
myCol <- colorRampPalette(c("lightgreen","red"))(100)
pdf('~/Documents/Aoua_PhD/data/Manuscript1/Figures3/mali_pop_xpehh.pdf',width = 20, height = 20)
ha = HeatmapAnnotation("-log10(q-value)"=anno_boxplot(tmat, height = unit(3, "cm"),
                                                      pch = 20,size = unit(1.5, "mm")))
Heatmap(tmat, name = "-log10(q-value)", column_split=xpehh_pop_mali$CHR, heatmap_width = unit(40, "cm"),
        heatmap_height = unit(25, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),show_column_names =F,
        row_names_g = gpar(fontsize = 16,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
        top_annotation = ha,col=colorRamp2(myBreaks, myCol))
dev.off()

heat
####read rsb and xpehh
bko_boug.rsb <- read.table("RSB/bko_boug.rsb.txt", sep = "\t", header = T)
bko_boug.xpehh <- read.table("xpehh/bko_boug.xpehh.txt", sep = "\t", header = T)
bko_dang.rsb <- read.table("RSB/bko_dang.rsb.txt", sep = "\t", header = T)
bko_dang.xpehh <- read.table("xpehh/bko_dang.xpehh.txt", sep = "\t", header = T)
bko_fala.rsb <- read.table("RSB/bko_fala.rsb.txt", sep = "\t", header = T)
bko_fala.xpehh <- read.table("xpehh/bko_fala.xpehh.txt", sep = "\t", header = T)
bko_ken.rsb <- read.table("RSB/bko_ken.rsb.txt", sep = "\t", header = T)
bko_ken.xpehh <- read.table("xpehh/bko_ken.xpehh.txt", sep = "\t", header = T)
bko_kol.rsb <- read.table("RSB/bko_kol.rsb.txt", sep = "\t", header = T)
bko_kol.xpehh <- read.table("xpehh/bko_kol.xpehh.txt", sep = "\t", header = T)
bko_nioro.rsb <- read.table("RSB/bko_nioro.rsb.txt", sep = "\t", header = T)
bko_nioro.xpehh <- read.table("xpehh/bko_nioro.xpehh.txt", sep = "\t", header = T)

boug_dang.rsb <- read.table("RSB/boug_dang.rsb.txt", sep = "\t", header = T)
boug_dang.xpehh <- read.table("xpehh/boug_dang.xpehh.txt", sep = "\t", header = T)
boug_fala.rsb <- read.table("RSB/boug_fala.rsb.txt", sep = "\t", header = T)
boug_fala.xpehh <- read.table("xpehh/boug_fala.xpehh.txt", sep = "\t", header = T)
boug_ken.rsb <- read.table("RSB/boug_ken.rsb.txt", sep = "\t", header = T)
boug_ken.xpehh <- read.table("xpehh/boug_ken.xpehh.txt", sep = "\t", header = T)
boug_kol.rsb <- read.table("RSB/boug_kol.rsb.txt", sep = "\t", header = T)
boug_kol.xpehh <- read.table("xpehh/boug_kol.xpehh.txt", sep = "\t", header = T)
boug_nioro.rsb <- read.table("RSB/boug_nioro.rsb.txt", sep = "\t", header = T)
boug_nioro.xpehh <- read.table("xpehh/boug_nioro.xpehh.txt", sep = "\t", header = T)

dang_fala.rsb <- read.table("RSB/dang_fala.rsb.txt", sep = "\t", header = T)
dang_fala.xpehh <- read.table("xpehh/dang_fala.xpehh.txt", sep = "\t", header = T)
dang_ken.rsb <- read.table("RSB/dang_ken.rsb.txt", sep = "\t", header = T)
dang_ken.xpehh <- read.table("xpehh/dang_ken.xpehh.txt", sep = "\t", header = T)
dang_kol.rsb <- read.table("RSB/dang_kol.rsb.txt", sep = "\t", header = T)
dang_kol.xpehh <- read.table("xpehh/dang_kol.xpehh.txt", sep = "\t", header = T)
dang_nioro.rsb <- read.table("RSB/dang_nioro.rsb.txt", sep = "\t", header = T)
dang_nioro.xpehh <- read.table("xpehh/dang_nioro.xpehh.txt", sep = "\t", header = T)

fala_ken.rsb <- read.table("RSB/fala_ken.rsb.txt", sep = "\t", header = T)
fala_ken.xpehh <- read.table("xpehh/fala_ken.xpehh.txt", sep = "\t", header = T)
fala_kol.rsb <- read.table("RSB/fala_kol.rsb.txt", sep = "\t", header = T)
fala_kol.xpehh <- read.table("xpehh/fala_kol.xpehh.txt", sep = "\t", header = T)
fala_nioro.rsb <- read.table("RSB/fala_nioro.rsb.txt", sep = "\t", header = T)
fala_nioro.xpehh <- read.table("xpehh/fala_nioro.xpehh.txt", sep = "\t", header = T)

ken_kol.rsb <- read.table("RSB/ken_kol.rsb.txt", sep = "\t", header = T)
ken_kol.xpehh <- read.table("xpehh/ken_kol.xpehh.txt", sep = "\t", header = T)
ken_nioro.rsb <- read.table("RSB//ken_nioro.rsb.txt", sep = "\t", header = T)
ken_nioro.xpehh <- read.table("xpehh/ken_nioro.xpehh.txt", sep = "\t", header = T)
kol_nioro.rsb <- read.table("RSB/kol_nioro.rsb.txt", sep = "\t", header = T)
kol_nioro.xpehh <- read.table("xpehh/kol_nioro.xpehh.txt", sep = "\t", header = T)

####manhattan plot ihs, xpehh, rsb #############
pdf("./xpehh_mali_pop.pdf",paper = "a4r", width = 10)
par(mfrow=c(4,2),oma=c(3, 3, 2, 2), mar=c(2,1,2,1))
manhattanplot(bko_boug.xpehh,  pval = T, ylim=c(0, 30))
manhattanplot(bko_dang.xpehh, pval = T, ylim=c(0, 30) )
manhattanplot(bko_fala.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(bko_ken.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(bko_kol.xpehh, pval = T, ylim=c(0, 30) )
manhattanplot(bko_nioro.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(boug_dang.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(boug_fala.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(boug_ken.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(boug_kol.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(boug_nioro.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(dang_fala.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(dang_ken.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(dang_kol.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(dang_nioro.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(fala_ken.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(fala_kol.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(fala_nioro.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(ken_kol.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(ken_nioro.xpehh, pval = T,  ylim=c(0, 30))
manhattanplot(kol_nioro.xpehh, pval = T,  ylim=c(0, 30))
dev.off()
####plot of xpehh_vs_rsb
pdf("./rsb_vs_xpehh_ken_pop.pdf",paper = "a4r", width = 10)
par(mfrow=c(4,2),oma=c(3, 3, 2, 2), mar=c(2,1,2,1))
plot(ken_kol.rsb[,3],ken_kol.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
plot(ken_nioro.rsb[,3],ken_nioro.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
dev.off()
plot(fala_nioro.rsb[,3],fala_nioro.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
dev.off()
plot(dang_nioro.rsb[,3],dang_nioro.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
dev.off()
plot(boug_nioro.rsb[,3],boug_nioro.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
dev.off()
plot(bko_nioro.rsb[,3],bko_nioro.xpehh[,3],xlab="Rsb",ylab="XP-EHH",pch=".",
     xlim=c(-10,7.5),ylim=c(-10,7.5))
abline(a = 0, b = 1, lty = 2)
dev.off()
#################ajust p values ##########################
###add an extra column to identify each pair########
bko_boug.rsb$PAIRS <- 1 
bko_dang.rsb$PAIRS <- 2
bko_fala.rsb$PAIRS <- 3
bko_ken.rsb$PAIRS <- 4
bko_kol.rsb$PAIRS <- 5
bko_nioro.rsb$PAIRS <- 6
boug_dang.rsb$PAIRS <- 7
boug_fala.rsb$PAIRS <- 8
boug_ken.rsb$PAIRS <- 9
boug_kol.rsb$PAIRS <- 10
boug_nioro.rsb$PAIRS <- 11
dang_fala.rsb$PAIRS <- 12
dang_ken.rsb$PAIRS <- 13
dang_kol.rsb$PAIRS <- 14
dang_nioro.rsb$PAIRS <- 15
fala_ken.rsb$PAIRS <- 16
fala_kol.rsb$PAIRS <- 17
fala_nioro.rsb$PAIRS <- 18
ken_kol.rsb$PAIRS<- 19
ken_nioro.rsb$PAIRS <- 20 
kol_nioro.rsb$PAIRS<- 21

################concatenate all rsb data into one and extract the pvalues #########
rsb_pairs<- rbind (bko_boug.rsb, bko_dang.rsb, bko_fala.rsb, bko_ken.rsb, bko_kol.rsb, bko_nioro.rsb,
                       boug_dang.rsb, boug_fala.rsb, boug_ken.rsb, boug_kol.rsb, boug_nioro.rsb, 
                        dang_fala.rsb, dang_ken.rsb, dang_kol.rsb, dang_nioro.rsb,
                        fala_ken.rsb, fala_kol.rsb, fala_nioro.rsb,
                        ken_kol.rsb, ken_nioro.rsb, kol_nioro.rsb)
all_rsb.qqman <- data.frame(
  CHR = as.integer(factor(rsb_pairs$CHR, 
                          levels = unique(rsb_pairs$CHR))),
  # chromosomes as integers
  BP = rsb_pairs$POSITION,         # base pairs
  P = rsb_pairs$P,  # transform back to p-values
  SNP = rsb_pairs$SNP     # SNP names
)
library(qqman) man

manhattan(all_rsb.qqman,
          col = c("gray10","gray60"),
          chrlabs = as.character(unique(all_rsb.qqman$CHR)),
          main = "RSB scan between all pairs",
          suggestiveline = 5,
        
           #annotatePval = 0.00001,
         #annotateTop = ,
          genomewideline = F,
          ylim = c(0, 50),
          logp = FALSE, ylab = "-log10(q-value)"
          )

write.table(rsb_pairs, file = "rsb_pairs.txt", row.names = F, quote = F, sep = "\t")


##############p-value ajust#######

##plotting heatmaps
library(heatmap.plus)
library(ComplexHeatmap)
library(circlize)
mat5<-all_pop_xpehh[,4:24]
tmat<-t(mat5)
tmat<-as.matrix(tmat)
myBreaks <- seq(0.0 , 20,  length.out =100)
myCol <- colorRampPalette(c("lightgreen","red"))(100)
pdf('mali_pop_xephh.pdf',width = 20, height = 20)
ha = HeatmapAnnotation(XPEHH=anno_boxplot(tmat, height = unit(1.5, "cm"),
                                        pch = 20,size = unit(1, "mm")))
Heatmap(tmat, name = "XPEHH", column_split=all_pop_xpehh$CHR, heatmap_width = unit(28, "cm"),
        heatmap_height = unit(15, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),show_column_names =F,
        row_names_g = gpar(fontsize = 11,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
        top_annotation = ha,col=colorRamp2(myBreaks, myCol))
dev.off()

###############merge all pop rsb##################

all_pope(merge(merge(merge(merge( merge(merge(bko_boug.xpehh[, c(1,2,5)], bko_dang.xpehh[, c(1,5)], by_xpehh <- merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merg = "SNP", all = TRUE), bko_fala.xpehh[,c(1,5)], by ="SNP", all = TRUE), bko_ken.xpehh[,c(1,5)], by = "SNP", all = TRUE), bko_kol.xpehh[, c(1,5)], by = "SNP", all = TRUE),  bko_nioro.xpehh[, c(1,5)], by = "SNP",all = TRUE), boug_dang.xpehh[, c(1,5)], by = "SNP", all = TRUE), boug_fala.xpehh[, c(1,5)], by= "SNP", all = TRUE), boug_ken.xpehh[, c(1,5)], by = "SNP",all = TRUE),  boug_kol.xpehh[, c(1,5)], by = "SNP", all = TRUE), boug_nioro.xpehh[, c(1,5)], by = "SNP", all = TRUE), dang_fala.xpehh[, c(1,5)], all = TRUE) ,dang_ken.xpehh[, c(1,5)], by ="SNP", all = TRUE), dang_kol.xpehh[, c(1,5)], by="SNP", all = TRUE), dang_nioro.xpehh[, c(1,5)], by = "SNP", all = TRUE), fala_ken.xpehh[, c(1,5)], by = "SNP", all = TRUE ), fala_kol.xpehh[, c(1,5)],by = "SNP", all = TRUE), fala_nioro.xpehh[, c(1,5)], by = "SNP", all = TRUE ), 
                  ken_kol.xpehh[, c(1,5)], by = "SNP", all = TRUE), ken_nioro.xpehh[, c(1,5)], by = "SNP", all = TRUE), kol_nioro.xpehh[, c(1,5)], by = "SNP", all = TRUE)

colnames(merge20)=c("SNP","CHR","bko vs boug", 
                        "bko vs dang","bko vs fala","bko vs ken",
                        "bko vs kol","bko vs nioro","boug vs dang",
                        "boug vs fala","boug vs ken","boug vs kol",
                        "boug vs nioro", "dang vs fala","dang vs ken",
                        "dang vs kol","dang vs nioro","fala vs ken",
                        "fala vs kol","fala vs nioro","ken vs kol",
                        "ken vs nioro","kol vs nioro")
#######################pvalue xephh####################
kol_nioro.xpehh <- na.omit(kol_nioro.xpehh)
kol_nioro.xpehh$P <- 10^(-kol_nioro.xpehh$LOGPVALUE)
fdr = fdrtool(kol_nioro.xpehh$P, statistic = "pvalue")
kol_nioro.xpehh$P <- -log10(fdr$qval)
##plotting heatmaps
library(heatmap.plus)
library(ComplexHeatmap)
library(circlize)
#$sb_all<-na.omit(xpehh_all)
mat5<-rsb_all[,3:23]
tmat<-t(mat5)
tmat<-as.matrix(tmat)
myBreaks <- seq(0.0 , 20,  length.out =100)
myCol <- colorRampPalette(c("lightgreen","red"))(100)
pdf('~/Documents/Aoua_PhD/data/Manuscript1/Figures3/rsb_mali.pdf',width = 20, height = 20)
ha = HeatmapAnnotation("-log10(q-value)"=anno_boxplot(tmat, height = unit(3, "cm"),
                                          pch = 20,size = unit(2, "mm")))
Heatmap(tmat, name = "-log10(q-value)", column_split=rsb_all$CHR, heatmap_width = unit(40, "cm"),
        heatmap_height = unit(25, "cm"), cluster_columns = F,show_column_dend = F,cluster_rows = F, show_row_dend = F,
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),show_column_names =F,
        row_names_g = gpar(fontsize = 16,fontface = "bold"),row_names_rot = 0,column_title_rot = 0,
        top_annotation = ha,col=colorRamp2(myBreaks, myCol))
dev.off()

###add an extra column to identify each pair########
bko_boug.xpehh$PAIRS <- 1 
bko_dang.xpehh$PAIRS <- 2
bko_fala.xpehh$PAIRS <- 3
bko_ken.xpehh$PAIRS <- 4
bko_kol.xpehh$PAIRS <- 5
bko_nioro.xpehh$PAIRS <- 6
boug_dang.xpehh$PAIRS <- 7
boug_fala.xpehh$PAIRS <- 8
boug_ken.xpehh$PAIRS <- 9
boug_kol.xpehh$PAIRS <- 10
boug_nioro.xpehh$PAIRS <- 11
dang_fala.xpehh$PAIRS <- 12
dang_ken.xpehh$PAIRS <- 13
dang_kol.xpehh$PAIRS <- 14
dang_nioro.xpehh$PAIRS <- 15
fala_ken.xpehh$PAIRS <- 16
fala_kol.xpehh$PAIRS <- 17
fala_nioro.xpehh$PAIRS <- 18
ken_kol.xpehh$PAIRS<- 19
ken_nioro.xpehh$PAIRS <- 20 
kol_nioro.xpehh$PAIRS<- 21

################concatenate all rsb data into one and extract the pvalues #########
xpehh_pairs<- rbind (bko_boug.xpehh, bko_dang.xpehh, bko_fala.xpehh, bko_ken.xpehh, bko_kol.xpehh, bko_nioro.xpehh,
                   boug_dang.xpehh, boug_fala.xpehh, boug_ken.xpehh, boug_kol.xpehh, boug_nioro.xpehh, 
                   dang_fala.xpehh, dang_ken.xpehh, dang_kol.xpehh, dang_nioro.xpehh,
                   fala_ken.xpehh, fala_kol.xpehh, fala_nioro.xpehh,
                   ken_kol.xpehh, ken_nioro.xpehh, kol_nioro.xpehh)

xpehh_pairs.qqman <- data.frame(
  CHR = as.integer(factor(xpehh_pairs$CHR, 
                          levels = unique(xpehh_pairs$CHR))),
  # chromosomes as integers
  BP = xpehh_pairs$POSITION,         # base pairs
  P = xpehh_pairs$P,  # transform back to p-values
  SNP = xpehh_pairs$SNP     # SNP names
)
library(qqman)

manhattan(xpehh_pairs.qqman,
          col = c("lightgreen","red"),
          chrlabs = as.character(unique(xpehh_pairs.qqman$CHR)),
          suggestiveline = F,
          #annotatePval = 0.05,
          genomewideline = T,
          ylim = c(0, 50),
          logp = FALSE, ylab = "-log10(q)"
)

#all_pop_xpehh <- merge20
write.table(xpehh_pairs, file = "xpehh_pairs.txt", row.names = F, sep = "\t", quote = F)
#all_pop_xpehh <- read.table("all_pop_xpehh.txt", header = T, sep = "\t")


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
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = "*", cex = 0.5, ...)
    }
  }  
  par(xpd = FALSE)
}

