library(tidyverse)

outliers_detection <- function(data){
    upper_bound <- quantile(data$Fst, 0.99)
    return(upper_bound)
}


fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop1.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))

bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers1 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers1)[2] <- "Pop1"
#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop2.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers2 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers2)[2] <- "Pop2"
#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop3.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers3 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers3)[2] <- "Pop3"
#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop4.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers4 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers4)[2] <- "Pop4"
#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop5.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers5 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers5)[2] <- "Pop5"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop6.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers6 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers6)[2] <- "Pop6"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop7.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers7 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers7)[2] <- "Pop7"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop8.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers8 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers8)[2] <- "Pop8"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop9.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers9 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers9)[2] <- "Pop9"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop10.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers10 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers10)[2] <- "Pop10"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop11.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers11 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers11)[2] <- "Pop11"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop12.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers12 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers12)[2] <- "Pop12"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop13.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers13 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers13)[2] <- "Pop13"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop14.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers14 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers14)[2] <- "Pop14"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop15.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers15 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers15)[2] <- "Pop15"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop16.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers16 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers16)[2] <- "Pop16"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop17.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers17 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers17)[2] <- "Pop17"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop18.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers18 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers18)[2] <- "Pop18"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop19.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers19 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers19)[2] <- "Pop19"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop20.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers20 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers20)[2] <- "Pop20"

#===============================
fst <- read_tsv("/media/Data/Data/Documents_Karim/Fadel/Aoua/pop21.fst.txt") %>% 
    mutate(Fst = if_else(Fst < 0, 0, Fst))
bound <- outliers_detection(fst)

outlier_ind <- which(fst$Fst > bound)

outliers21 <- fst[outlier_ind, ] %>% select(SNP, Fst)
names(outliers21)[2] <- "Pop21"



join_outliers <- full_join(outliers1, outliers2, by = "SNP") %>% 
    full_join(outliers3, by = "SNP") %>% 
    full_join(outliers4, by = "SNP") %>% 
    full_join(outliers5, by = "SNP") %>% 
    full_join(outliers6, by = "SNP") %>% 
    full_join(outliers7, by = "SNP") %>% 
    full_join(outliers8, by = "SNP") %>% 
    full_join(outliers9, by = "SNP") %>% 
    full_join(outliers10, by = "SNP") %>% 
    full_join(outliers11, by = "SNP") %>% 
    full_join(outliers12, by = "SNP") %>% 
    full_join(outliers13, by = "SNP") %>% 
    full_join(outliers14, by = "SNP") %>% 
    full_join(outliers15, by = "SNP") %>% 
    full_join(outliers16, by = "SNP") %>% 
    full_join(outliers17, by = "SNP") %>% 
    full_join(outliers18, by = "SNP") %>% 
    full_join(outliers19, by = "SNP") %>% 
    full_join(outliers20, by = "SNP") %>% 
    full_join(outliers21, by = "SNP")



join_outliers[which(rowMeans(!is.na(join_outliers)) > 0.3), ] %>% 
    arrange(SNP) %>% 
    write.table("/media/Data/Data/Documents_Karim/Fadel/Aoua/outliers_across_pairwise_comparison.xlsx", 
                col.names = T, row.names = F, quote = F, sep = '\t')
