
metadata1 <- read_tsv("data/metadata/Mali_metadata.txt") %>% 
    dplyr::select(Sample, Location, Year)

metadata1 <- metadata1[metadata1$Sample %in% sampleIDs,]

data <- read.delim("data/raw_data/Mali_genotype_data.txt", header = T, sep = '\t') %>% as_tibble()

na_count_locus <-sapply(data, function(y) sum(length(which(is.na(y)))))

na_count_locus <- data.frame(na_count_locus)

data_trans <- t(data[,-c(1,2)]) %>% as_tibble()

na_count_sample <-sapply(data_trans, function(y) sum(length(which(is.na(y)))))

na_count_sample <- data.frame(na_count_sample)

cbind(data$Samples, na_count_sample)

na_count_locus %>% filter(na_count > 160)

our_data <- read_tsv("data/out.lmiss") %>% 
    filter(F_MISS <0.2)

aoua_data <- colnames(data)[-c(1,2)] %>% 
    str_split(., pattern = "_(?!.*_)", n = Inf, simplify = FALSE)

aoua_data <- data.frame(matrix(unlist(aoua_data), nrow=length(aoua_data), byrow=TRUE))
aoua_data <- aoua_data %>% mutate_at(vars(X2), as.numeric)

our_data %>% inner_join(aoua_data, by = c("CHR"="X1", "POS"="X2"))
