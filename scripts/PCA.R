setwd("~/Documents/Aoua_PhD/data/PLINK/")
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
ped <- read.table("Maj.ped", sep = "\t", header = FALSE)
bed_maj.fn <- "Pruned_data.bed"
fam_maj.fn <- "Pruned_data.fam"
bim_maj.fn <- "Pruned_data.bim" ####, cvt.chr="char")
#####comvert plink data into gds format############
snpgdsBED2GDS(bed_maj.fn, fam_maj.fn, bim_maj.fn, "Maj.gds", cvt.chr="char")
####summary####
snpgdsSummary("Maj.gds")
genofile_maj <- snpgdsOpen("Maj.gds")

# Get population information
 pop_code <- scan("pop_info.txt", what=character())
sample.id <- read.gdsn(index.gdsn(genofile_maj, "sample.id"))
# Try different LD thresholds for sensitivity analysis
#snpset <- snpgdsLDpruning(genofile_maj, ld.threshold=0.2)
#names(snpset)
#snpset.id <- unlist(snpset)
pca_maj <- snpgdsPCA(genofile_maj, snp.id= NULL,num.thread=2, autosome.only = FALSE )

snpset.id <- read.gdsn(index.gdsn(genofile_maj, "snp.id"))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
# Draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
head(cbind(sample.id, pop_code))
tab <- data.frame(sample.id = pca_maj$sample.id,
                  pop = factor(pop_code)[match(pca_maj$sample.id, sample.id)],
                  EV1 = pca_maj$eigenvect[,1],    # the first eigenvector
                  EV2 = pca_maj$eigenvect[,2],    # the second eigenvector
                  EV3 = pca_maj$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)
#Draw

PC1<- tab$EV1
PC2 <- tab$EV2
PC3 <- tab$EV3
Location <- tab$pop
# variance proportion (%)
pc.percent <- pca_maj$varprop*100
head(round(pc.percent, 2))

ggplot(tab, aes(x=PC1 ,y=PC2))+
  geom_point(aes(fill=Location) ,colour="black", shape=21, size=3)+
  scale_fill_manual(values = c( "Bamako" ="gold4", "Bougoula-Hameau" = "red3", "Dangassa" = "darkmagenta", "Faladje" = "blue2", "Kenieroba" ="darkorange", "Kolle" = "darkgreen", "Nioro" = "black" )) +
 # ggtitle("Principal Componant Analysis")+
 # scale_color_manual(values = Mypalette) +
  xlab("PC1(34%)")+
  ylab("PC2(30%)")+
  theme(legend.position="bottom", axis.text=element_text(size=14),
        axis.title.x = element_text(size =14), axis.title.y = element_text(size =14),legend.text = element_text(size = 12),legend.title = element_text(size =12))

  
ggplot(tab, aes(x=eigenvector_1 ,y=eigenvector_3))+
  geom_point(aes(fill=Location) ,colour="black", shape=21)+
  ggtitle("Principal Componant Analysis")+
  # scale_color_manual(values = Mypalette) +
  theme(legend.position="bottom")
####identity by states####
ibs_maj <- snpgdsIBS(genofile_maj, num.thread=2, remove.monosnp = TRUE)

# individulas in the same population are clustered together
pop.idx <- order(pop_code) 
#, xlab = "ibs", ylab = "ibs"
image(ibs_maj$ibs[pop.idx, pop.idx], col=terrain.colors(21), xlab = "ibs", ylab = "ibs" )
####mds analysis#####
loc_maj <- cmdscale(1 - ibs_maj$ibs)
coordinate_1 <- loc_maj[, 1]; coordinate_2 <- loc_maj[, 2]
race <- as.factor(pop_code)

plot(coordinate_1, coordinate_2, col= tab$pop, xlab = "coordinate1", ylab = "coordinate2",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(tab$pop), text.col=1:nlevels(tab$pop))
ggplot(as.data.frame(loc_maj), aes(x=coordinate_1 ,y=coordinate_2, col = tab$pop ))+
  geom_point(aes(fill=tab$pop) ,colour="black", shape=21)+
  ggtitle( "Multidimensional Scaling Analysis (IBS) Majority allele")+
  scale_color_manual(values = Mypalette) +
  theme(legend.position="bottom") 
