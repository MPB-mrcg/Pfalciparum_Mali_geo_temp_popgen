setwd("~/Documents/Aoua_PhD/data/FWS/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("ggpubr")
install.packages("geneplotter")
install.packages("pwr")
library(pwr)
#BiocManager::install("SeqArray", "SeqVarTools", "IRanges", "S4Vectors", "GenomicRanges", "BiocParallel", "BiocStyle"))

metadata <- read.table("sample_pop_info.txt", sep = "\t", header = TRUE) # genotype data
#metadata <- metadata[order(metadata$Samples),]
# convert VCF file to the GDS format 
seqVCF2GDS("Mali_bi_Snps20_indv20_dp5_MQ30_maf2.vcf.recode.vcf", "Mali_genomic.gds") 
#read gds
my_vcf <-seqOpen("Mali_genomic.gds")
seqSummary(my_vcf)
# save sample identifiers
sample.id <- seqGetData(my_vcf, "sample.id")
# get genomic coordinates of all variants
coords <- getCoordinates(my_vcf)
#estimate MOI with fws 
fws <- getFws(my_vcf)
#hist(fws_all)
#plot(fws_all)
#library(ggplot2)
#library(dplyr)
#ggplot(as.data.frame(fws_all))

fws=as.data.frame(fws)
fws$Samples=rownames(fws)
fws$Samples
Final_fws <- merge(fws, metadata, by = "Samples") # add location
Final_fws$
#plot fws
Location<- Final_fws$Location
FWS <-Final_fws$fws
Sample <- Final_fws$Samples
#p<-ggplot(Final_fws, aes(x=Location,  y=FWS, color= Location ))
#p<-p+geom_point()
#p

df(paste("Complexity of infections ", "pdf", sep = ""))
BOX <- boxplot(FWS ~ Location, xlab="Location", ylab="fws") 
stripchart( FWS~ Location , add= TRUE, horizontal = TRUE, method = "jitter", pch = 21)
boxplot.n()
print(BOX)

df(paste("Complexity of infections ", "pdf", sep = ""))

ggplot(Final_fws, aes(x=Location ,y=fws, col = Location ))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2))+
  scale_color_manual(values = c( "Bamako" ="gold4", "Bougoula-Hameau" = "red3", "Dangassa" = "darkmagenta", "Faladje" = "blue2", "Kenieroba" ="darkorange", "Kolle" = "darkgreen", "Nioro" = "black" )) +
  #scale_color_manual(values = Mypalette) +
 
  ylab("Fws")+
  theme(legend.position="bottom", axis.text=element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size =12), axis.title.y = element_text(size =12),legend.text = element_text(size = 11),legend.title = element_text(size =11))

dev.off()

df(paste("Complexity of infections ", "pdf", sep = ""))

ggplot(Final_fws, aes(x=Location, y=fws, col = Location, label = scales::percent(Location)))+
  #geom_boxplot()+
 
  geom_jitter(position=position_jitter(0.2))+
  scale_color_manual(values = c( "Bamako" ="gold4", "Bougoula-Hameau" = "red3", "Dangassa" = "darkmagenta", "Faladje" = "blue2", "Kenieroba" ="darkorange", "Kolle" = "darkgreen", "Nioro" = "black" )) +
  #scale_color_manual(values = Mypalette) +
  
  ylab("Fws")+
  theme(legend.position="bottom", axis.text=element_text(size=12),axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size =12), axis.title.y = element_text(size =12),legend.text = element_text(size = 11),legend.title = element_text(size =11))

dev.off()
myList <-list(c("Bamako", "Nioro"))
ggboxplot(Final_fws, x= "Locations" , y= "fws", color = "Locations",
          palette = c( "gold4", "red3", "darkmagenta", "blue2","darkorange", "darkgreen","black"),
          add = "jitter", size = 1)+
#stat_compare_means(label = "p.signif", method = "wilcox.test")+
 #stat_compare_means(comparisons = myList)+ 
  

  stat_compare_means(method = "kruskal.test",label.y = 1.1 )+
  ylab("Fws")+
theme_pubr( base_size = 16, base_family = "",border = FALSE,
  margin = TRUE ,x.text.angle = 45, legend = "none" )#legend = "top"

 
####if fws >= 0.95 the sample is monoclonal, sample is polyclonal if fws < 0.95 
Final_fws$coil <- ifelse(Final_fws$fws>=0.95, "1", "2")
Final_fws$infection <- ifelse(Final_fws$coil == "1", "single", "multiple")
write.table(Final_fws, file = "fws_file.txt", quote = F, row.names = F, sep = "\t")
Final_fws <- read.table("fws_file.txt", header = TRUE, sep = "\t")
Final_fws
####fws proportion plot ##################
rm(list = ls())
library(tidyverse)

#fws <- read_delim("fws_file.txt", delim = "\t", col_names = c("Samples", "Fws", "Locations", "COI"), skip = 1)
FWS_data <- Final_fws %>% group_by(Locations) %>% arrange(fws) %>% mutate(rw = row_number() ) %>%
  mutate(num = n()) %>% mutate(pop.prop = row_number()/n()) %>% arrange(Locations)
#/n()) %>% arrange(Locations)
####plot fws vs prop ######

#Fws <- fws %>% 
  #group_by(Locations) %>% 
  #arrange(Fws) %>% 
  #mutate(Cum = row_number()/n()) %>% 
  #arrange(Locations)

FWS_data %>% 
  ggplot(aes(x = pop.prop, y = fws, color = Locations)) +
  geom_point() + 
  scale_color_manual(values = c( "Bamako" ="gold4", "Bougoula-Hameau" = "red3", "Dangassa" = "darkmagenta", "Faladje" = "blue2", "Kenieroba" ="darkorange", "Kolle" = "darkgreen", "Nioro" = "black" ))+
  geom_hline(yintercept = 0.95, linetype="dashed", color = "red") +
  labs(x="Proportion", y= "Fws")

## SAVE PLOT
ggsave(filename = "fws_plots.jpeg")
############calculate mean fws in each population ##################
Final_fws %>%
  group_by(Locations) %>% filter(Locations=="Bamako" || Locations=="Bougoula-Hameau") %>% select(Locations, fws) %>% fisher.test()
  dplyr::summarize(Mean = mean(fws, na.rm=TRUE))
###############calculate proportion of fws > 0.95 ##########
Final_fws %>% group_by(coil) %>% summarise(n = n()) %>%
  mutate(freq = n / sum(n))
  ###############paiwaise comparaison using fisher test ########
 data <-  Final_fws %>% filter(Locations=="Bamako" | Locations=="Bougoula-Hameau") %>% select(fws, Locations) %>% 
    pairwise.t.test(data$Locations, data$fws, p.adj = "none")
  
  data <- Final_fws %>% filter(Locations=="Bamako" | Locations=="Faladje") %>% select(Locations, infection) %>%filter(infection =="single") %>% table() %>% fisher.test()
  dplyr::summarize(Mean = mean(fws, na.rm=TRUE))
  
  
  ###calculate p-value###
  monoclonal <- Final_fws %>% group_by(Locations) %>% filter(infection=="multiple") 
  niro_kolle <- monoclonal %>% filter(Locations=="Nioro" || Locations== "Kolle") %>% select(Locations, fws) 
  
  
  