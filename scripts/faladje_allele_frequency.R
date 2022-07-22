
years <- metadata %>% 
    filter(Location == "Faladje") %>% 
    pull(Year) %>% unique()

# OUTPUT <- "results/tables/Allele_frequency/"
# if(!dir.exists(OUTPUT)) dir.create(OUTPUT)

OUTPUT <- "results/tables/"

frequency_file <- "results/tables/Faladje_drugResistance_allelFreq.xlsx"

VCF <- "data/reference/Allele_frequency/drugResistance_loci.vcf"

for (i in years) {
    # gene_name <- gsub(".vcf", "", basename(VCF))
    
    # Extract sample IDs
    metadata %>% 
        filter(Location == "Faladje" & Year == i) %>%
        pull(Sample) %>% 
        # write.table(paste0(OUTPUT, gene_name, "_", i, ".txt"), col.names = F, row.names = F, quote = F)
        write.table(paste0(OUTPUT, i, ".txt"), col.names = F, row.names = F, quote = F)
    
    # cat("Compute allele frequency for", gene_name, "and year ", i)
    cat("Compute allele frequency for", i)
    
    system(paste0("vcftools --vcf ", VCF,
                  # " --keep ", paste0(OUTPUT, gene_name, "_", i, ".txt"),
                  " --keep ", paste0(OUTPUT, i, ".txt"),
                  # " --freq --out ", paste0(OUTPUT, gene_name, "_", i)))
                  " --freq --out ", paste0(OUTPUT, i)))
    
    ## LOAD ALLELE FREQUENCY FILES
    var_freq <- fread(paste0(OUTPUT, i, ".frq"))
    names(var_freq) = c("chrom", "pos", "nalleles", "nchr", "a1", "a2")
    
    var_freq <- var_freq %>% 
        separate(a1, c("REF"," Freq1"), sep = ":") %>% 
        separate(a2, c("ALT"," Freq2"), sep = ":") %>% 
        dplyr::select(-nalleles, -nchr) %>% 
        mutate(Year = i) %>% 
        relocate(Year, .before = chrom)
    
    ## SAVE MAF FILE
    if(!file.exists(frequency_file)){
        file.create(frequency_file)
        write.table(var_freq, frequency_file, col.names = TRUE,
                    row.names = FALSE, quote = FALSE, sep = '\t')
    }
    else write.table(var_freq, frequency_file, col.names = FALSE, 
                     append = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
    
    file.remove(paste0(OUTPUT, i, ".txt"),
                paste0(OUTPUT, i, ".frq"))
    
    # var_freq <- var_freq %>% 
    #     dplyr::select(a1, a2) %>% 
    #     apply(1, function(z) min(z))
}

data <- read.table(VCF, header = FALSE, stringsAsFactors = FALSE)

gids <- strsplit(as.vector(data$V8), ";")

geneIDs <- NULL
for (i in 1:length(gids))
{
    print(paste0('i= ', i))
    p <- grep("SNPEFF_AMINO_ACID_CHANGE", gids[[i]])
    pp <- grep("SNPEFF_GENE_NAME", gids[[i]])
    
    geneIDs <- rbind(geneIDs, data.frame(Chrom = data$V1[i],
                                         POS = data$V2[i],
                                         AMINO_ACID_CHANGE = gids[[i]][p],
                                         SNPEFF_GENE_NAME = gids[[i]][pp]))
    geneIDs$AMINO_ACID_CHANGE <- gsub("SNPEFF_AMINO_ACID_CHANGE=","", geneIDs$AMINO_ACID_CHANGE)
    geneIDs$SNPEFF_GENE_NAME <- gsub("SNPEFF_GENE_NAME=","", geneIDs$SNPEFF_GENE_NAME)
}

read_tsv(frequency_file) %>% 
    inner_join(geneIDs, by = c("chrom" = "Chrom", "pos" = "POS")) %>% 
    write.table(frequency_file, col.names = TRUE,
                row.names = FALSE, quote = FALSE, sep = '\t')
