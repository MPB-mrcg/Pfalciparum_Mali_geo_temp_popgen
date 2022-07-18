setwd("~/Documents/Aoua_PhD/data/Popstructure/LD/")
library(sommer)
library(ggplot2)
library(dplyr)
library(stringr)

# import the data
#ld <- read.table("snp_ld.summary", sep="\t",header=)
############split ped file by population and use that as input for plink LD calculation
pedFile <- read.table("~/Documents/Aoua_PhD/data/PLINK/Maj.ped", header = F, sep= "\t")
NEW1_ped <- read.table("~/Documents/Aoua_PhD/data/Selection/IBD/NewPedFile.ped", header = F, sep = "\t" )
write.table(NEW_ped, file = "~/Documents/Aoua_PhD/data/Selection/IBD/NewPedFile.ped", col.names = F, row.names = F, quote = F, sep = "\t")

write.table(bko.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Bamako"), file = "bko.ped", row.names = F, col.names = F, quote = F)
write.table(boug.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Bougoula-Hameau"), file = "boug.ped" , row.names = F , col.names = F, quote =F)
write.table(dang.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Dangassa"), file = "dang.ped", row.names = F, col.names = F, quote = F)
write.table(ken.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Kenieroba"), file = "ken.ped", row.names = F, col.names = F, quote = F)
write.table(kol.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Kolle"), file = "kol.ped", row.names = F, col.names = F, quote = F)
write.table(nioro.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Nioro"), file = "nioro.ped", row.names = F, col.names = F, quote = F)
write.table(fala.ped <-NEW1_ped %>% filter(NEW1_ped$V1 == "Faladje"), file = "fala.ped", row.names = F, col.names = F, quote = F)
###########################
plink --file bko  --allow-extra-chr --ld-window-kb 50000 --r2 --out bko_ld 
#dfr <- read.delim("snp_ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
#colnames(dfr) <- c("dist","rsq")
bko <- read.table("bko_ld.ld", header = T, sep = "\t")
boug <- read.table("boug_ld.ld", header = T, sep = "\t")
dang <- read.table("dang_ld.ld", header = T, sep = "\t")
ken <- read.table("ken_ld.ld", header = T, sep = "\t")
fala <- read.table("fala_ld.ld", header = T, sep = "\t")
kol <- read.table("kol_ld.ld", header = T, sep = "\t")
nioro <- read.table("nioro_ld.ld", header = T, sep = "\t")
###############calculate distance ###########
bko$dist <- bko$BP_B - bko$BP_A
boug$dist <- boug$BP_B - boug$BP_A
dang$dist <- dang$BP_B - dang$BP_A
fala$dist <- fala$BP_B - fala$BP_A
ken$dist <- ken$BP_B - ken$BP_A
kol$dist <- kol$BP_B - kol$BP_A
nioro$dist <- nioro$BP_B - nioro$BP_A

plink --file bko  --allow-extra-chr --ld-window-kb 50000 --r2 --out bko_ld 
##################################################################
#bko1<-bko[order(bko$BP_B-bko$BP_A),]
#boug1<-boug[order(boug$BP_B-boug$BP_A),]
#plot(bko1$BP_B-bko1$BP_A,bko1$R2,type="l"  
    # ,col="red",ylim=c(0,max(bko$R2,bko$R2)),
     #lwd=2,xlab="Distance between SNPs (bp)", 
     #ylab="Correlation")
#points(boug1$BP_B-boug1$BP_A,boug1$R2,type="l",col="blue",lwd=2)
#legend(500,0.15,c("bko","boug"),lwd=c(2,2),col=c("blue","red"))
###############plot ld #######################
bko$Location <- "Bamako"
boug$Location <- "Bougoula-Hameau"
dang$Location <- "Dangassa"
fala$Location <- "Faladje"
ken$Location <- "Kenieroba"
kol$Location <- "Kolle"
nioro$Location <- "Nioro"
################put all ld together######
all_ld_info <- read.table("all_ld_info_new.txt", header = T, sep = "\t")
Location <- all_ld_info$Location
##########correct r2 #######################
intervals <- c(20,30,50,100,300,500,1000,3000,5000,10000,20000)
bins=seq(1,max(all_ld_info$dist), by=1000)
#bins=c(0 + cumsum(c(0, intervals)), max(all_ld_info$dist))
my.means=rep(0, length(bins))
LD.averages=data.frame(bins, my.means, location)
for (i in 1:length(bins)) {
  data.interval=subset(all_ld_info, (all_ld_info$dist >= bins[i] & all_ld_info$dist < (bins[i]+1000))) 
  LD.averages$my.means[i]=mean(data.interval$R2) 
}
head(LD.averages)



Mypalette <-  c("gold4 ", "red3", "darkmagenta ", "blue2","darkorange", "darkgreen" ,"black" )
#all_ld_info <- rbind(bko, boug, dang, fala, ken, kol, nioro)
ggplot(all_ld_info, aes(x=all_ld_info$dist ,y=all_ld_info$R2, group= all_ld_info$Location, color = Location))+
  scale_color_manual(values = c( "Bamako" ="gold4", "Bougoula-Hameau" = "red3", "Dangassa" = "darkmagenta", "Faladje" = "blue2", "Kenieroba" ="darkorange", "Kolle" = "darkgreen", "Nioro" = "black" )) +
   #geom_point()+
 geom_smooth(method = "loess", span = 0.07, size = 1)+
  #geom_smooth(size =1)+
#
  #stat_smooth(method = "loess", formula = y ~ x, size = 1)+
 # xlim(0,5000)+
  #ylim(0, 0.2 )+
  #ggtitle("LD decay")+
  xlab("Distance")+
  ylab("LD(r^2)")+
 # scale_color_manual(values = Mypalette) +
  theme(legend.position="bottom", axis.text=element_text(size=12),
        axis.title.x = element_text(size =12), axis.title.y = element_text(size =12),legend.text = element_text(size = 11),legend.title = element_text(size =11))

  #theme(legend.position="bottom")

devtools::install_github("krlmlr/ulimit")



# first plot
scatter.smooth(x=all_ld_info$dist, y=all_ld_info$R2)
# add plot of second lowess line
lines(loess.smooth(x=all_ld_info$dist, y=all_ld_info$R2, col= all_ld_info$Location))

 Location <- all_ld_info$Location
ggplot(data= all_ld_info, aes(x=all_ld_info$dist, group= all_ld_info$Location, colour = Location)) +
  geom_density(aes(colour = Location))+
  xlab("Distance") +
  ylab("LD(r^2)")
theme_ipsum() 

group.colors <- c("Bamako" == "#333BFF", "Bougoula-Hameau" == "#CC6600", "Dangassa" =="#9633FF", "Faladje" == "#E2FF33", "Kenieroba"== "#E3DB71")
write.table(all_ld_info, file ="all_ld_info.txt", quote = F, row.names = F, sep = "\t")
