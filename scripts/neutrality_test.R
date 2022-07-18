#!/usr/bin/env Rscript

rm(list = ls())
require(data.table)
library(PopGenome)
require(dplyr)

args = commandArgs(trailingOnly = TRUE)
print(args)
InputFile <- as.character(args[1]) 
gtFile <- as.character(args[2])
pathToOutput <- as.character(args[3])

if(!file.exists(InputFile))
stop(InputFile,": No such file or directory")
if(!file.exists(gtFile))
stop(gtFile,": No such file or directory")
if(!dir.exists(pathToOutput))
stop(pathToOutput,": No such file or directory")

tempDir <- file.path(pathToOutput, gsub(".filtered.vcf.gz", "", basename(InputFile)))
system(paste0("mkdir -p ", tempDir))

## Get VCF file per gene and Remove Invariants
gff <- fread(gtFile, header = TRUE)
inc <- 1
for (i in 1:nrow(gff)) {
	print(paste0("i = ", i))
	region <- paste0(gff$Chromosome[i], ":", gff$Start[i], "-", gff$End [i])
	output <- paste0(tempDir, "/", gff$Gene_ID[i], ".vcf")
	system(paste0("bcftools view -r ", region," -o ", output, " ", InputFile))

	system(paste0("bgzip ", output))
	system(paste0("tabix ", output, ".gz"))

	## Extract header from VCF to get windows size for tajima's D estimation
	header <- system(paste0("zcat ", output, " | grep '##' "), intern = TRUE)
	skipNum <- as.numeric(system(paste("zgrep \"##\"", output," | wc -l"), TRUE))
	vcf <- read.table(output, skip = skipNum, header = TRUE, comment.char = "", 
					  stringsAsFactors = FALSE, check.names=FALSE)

	NbreSNPs <- nrow(vcf)

	if(NbreSNPs >=3){
		GENOME.class <- readVCF(output, NbreSNPs, gff$Chromosome[i], gff$Start[i], gff$End[i], include.unknown = TRUE)
		Stats <- neutrality.stats(GENOME.class, FAST = TRUE)

		if(inc == 1){
			geneName <- gsub(".vcf.gz", "", basename(output))
			ff <- as.data.frame(subset(get.neutrality(Stats)[[1]], select=-c(3,6:9)))

			fuli <- cbind(geneName, gff$Chromosome[i], rownames(ff), ff)
		}
		else{
			geneName <- gsub(".vcf.gz", "", basename(output))
			ff <- as.data.frame(subset(get.neutrality(Stats)[[1]], select=-c(3,6:9)))

			fuli <- rbind(fuli,  cbind(geneName, gff$Chromosome[i], rownames(ff), ff))
		}
		inc <- inc + 1
		file.remove(output, output)
	}
	cat('\n')
	file.remove(output)
}

names(fuli) <- c("Gene", "Chrom", "Pos", "TajimasD", "N.S.S", "Fu.Li.F", "Fu.Li.D")
write.table(fuli, paste0(pathToOutput, gsub(".vcf.gz", "_fuli.xlsx", basename(InputFile))), 
			col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Remove temporaly files
unlink(tempDir, recursive = TRUE)