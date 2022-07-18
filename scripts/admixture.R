
library(tidyverse)

validation <- read.delim("data/admixture/Mali.cv.error", header = FALSE, 
                         sep="", col.names = c("K", "Error")) %>% arrange(K)

plot(validation$Error~validation$K , type="b" , lwd=2 , col=rgb(0.1,0.7,0.1,0.8) , 
     ylab="Cross-validation Error" , xlab="K" , bty="l" , pch=20 , cex=2)

# Assign the first argument to prefix
prefix <- "data/admixture/Mali"

# Get individual names in the correct order
labels<-read.table("data/admixture/metadata.txt", header = T)

# Name the columns
names(labels)<-c("ind","pop")

populations <- labels %>% pull(pop) %>% unique()

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n <- factor(labels$pop, levels=sort(populations))
levels(labels$n) <- c(1:length(levels(labels$n)))
labels$n <- as.integer(as.character(labels$n))

# read in the different admixture output files
minK=1
maxK=10
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))

# tbl <- read.table("data/admixture/Mali.3.Q")
# barplot(t(as.matrix(tbl)), col=rainbow(3),
#           xlab="Individual #", ylab="Ancestry", border=NA)

# Prepare spaces to separate the populations/species
rep <- as.vector(table(labels$n))
spaces <- 0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces <- spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1200,res=200)

par(mfrow=c(maxK,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.3,cex.axis=1.2)

# Plot minK
bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=rainbow(n=minK),xaxt="n", 
            border=NA,ylab=paste0("K = ",minK),yaxt="n",space=spaces, font.axis = 4)

# axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)

# Plot higher K values
if(maxK>minK)lapply(2:(maxK), function(x) 
    barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x),xaxt="n", 
            border=NA,ylab=paste0("K = ",x),yaxt="n",space=spaces, cex.axis=2))

axis(1,at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
     labels=sort(populations),font = 4)


dev.off()


