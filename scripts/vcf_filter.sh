

## We will first remove all the missingness: 
# extract the genotypes using bcftools
#for each snp determine the proportion of sample with valide genotype.
Infile=/Users/coulibalya/Documents/Aoua_PhD/data/vcf-6.0/filtering

 #keep only snps remove indels 
 vcftools -gzvcf Mali_Data.vcf.gz --remove-indels  --recode --recode-INFO-all --stdout |  gzip -c >  Mali_snp_data.vcf.gz
####get number of snps
#3844304 snps, 1270 samples
### get only biallelic coding snps
#bcftools view Mali_snp_data.vcf.gz --include 'FILTER="PASS"'  --min-alleles 2 --max-alleles 2 --types snps --output-file Mali_biall_snps.vcf.gz 
#To create a vcf file which contains only PASS bi-allelic coding SNPs with
#VQSLOD > 6: ##&& VQSLOD>3.0 && MQ>20
bcftools view \
--include 'FILTER="PASS" && N_ALT=1 && CDS==1 && TYPE="snp" ' \
--output-type z \
--output-file Mali_biall_snps_qual.vcf.gz \
Mali_snp_data.vcf.gz
### number of snps after filtering 1042185####
#####individuals to keep
bcftools view -S /Users/coulibalya/Documents/Aoua_PhD/data/Metadata/Samples_tokeep_vcf.txt --output-type z --output-file Mali_bi_final_samples.vcf.gz Mali_biall_snps_qual.vcf.gz
####1044 samples. Badiangara sampes were removed because of small sample size. we selected samples from Faladje and Bougoula at time zero

#removing sites with up to 50% missing data
vcftools  --gzvcf Mali_bi_final_samples.vcf.gz --max-missing 0.5   --minDP 5  --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_miss50.vcf.gz
###number of snps  1028391# number of samples 1044

#### individuals with a lot of missing data 50% missing 
vcftools --gzvcf Mali_bi_miss50.vcf.gz  --missing-indv --out lowDP50

#extract individuals with more than 30% missing gdata
#cat miss_indv.imiss| sort -rk5 |   awk 'NR>1 && $5>=0.30 {print $1}' > inviduals_tobe_removed
mawk '!/IN/' miss_indv.imiss | cut -f5 > totalmissing_indv50
mawk '$5 > 0.5' lowDP50.imiss | cut -f1 > lowDP.indv50 ####  119 individuals to be removed 
# removing these individuals
#vcftools  --gzvcf Mali_bi_qual_final_samples.vcf.gz  --remove inviduals_tobe_removed   --recode --recode-INFO-all --out Mali_bi_qual_indv_miss30
vcftools  --gzvcf  Mali_bi_miss50.vcf.gz --remove lowDP.indv50  --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_Snps_indv50.vcf.gz 

### samples size 926####

####removing site with 40%  missingness 
vcftools --gvcf  Mali_bi_Snps_indv50.vcf.gz  --max-missing 0.60  --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_Snps40_indv50.vcf.gz 
###number of snps 1019257  and 926 samples################

####removing individuals with 40% missingness 
vcftools --gzvcf  Mali_bi_Snps40_indv50.vcf.gz --missing-indv --out lowDP40
mawk '$5 > 0.4' lowDP40.imiss | cut -f1 > lowDP.indv40 

# removing these individuals with 40% missingness
vcftools  --gzvcf  Mali_bi_Snps40_indv50.vcf.gz --remove lowDP.indv40  --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_Snps40_indv40.vcf.gz 
#### samples size 897 #### snps 1019257 
####remove snps with 30% missingness

vcftools --gvcf  Mali_bi_Snps40_indv40.vcf.gz --max-missing 0.70     --recode --recode-INFO-all --stdout |  gzip -c Mali_bi_Snps30_indv40.vcf.gz 
####number od snps  990741
####remocving individuals with 30% missingness

vcftools --gzvcf  Mali_bi_Snps30_indv40.vcf.gz --missing-indv --out lowDP30
mawk '$5 > 0.3' lowDP30.imiss | cut -f1 > lowDP.indv30 

# removing these individuals with 30% missingness
vcftools  --gzvcf  Mali_bi_Snps30_indv40.vcf.gz --remove lowDP.indv30  --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_Snps30_indv30.vcf.gz
###number of samples 865 number of snps  990741
####removing site with 25% missingness 
vcftools --gzvcf  Mali_bi_Snps30_indv30.vcf.gz --max-missing 0.75    --recode --recode-INFO-all --stdout |  gzip -c > Mali_bi_Snps25_indv30.vcf.gz
#number of snps 975874

### number of snps 
####remocving individuals with 25% missingness

vcftools --gzvcf  Mali_bi_Snps25_indv30.vcf.gz --missing-indv --out lowDP25
mawk '$5 > 0.25' lowDP25.imiss | cut -f1 > lowDP.indv25 

# removing these individuals with 25% missingness
vcftools  --gzvcf  Mali_bi_Snps25_indv30.vcf.gz --remove lowDP.indv25  --recode --recode-INFO-all --out  Mali_bi_Snps25_indv25_dp5
#number of samples 844 number of snps 975874 

####20% missingness
vcftools --vcf  Mali_bi_Snps25_indv25_dp5.recode.vcf --max-missing 0.8   --recode --recode-INFO-all --out Mali_bi_Snps20_indv25_dp5
# of snps 953289 
####remocving individuals with 20% missingness

vcftools --vcf  Mali_bi_Snps20_indv25_dp5.recode.vcf --missing-indv --out lowDP20
mawk '$5 > 0.20' lowDP20.imiss | cut -f1 > lowDP.indv20

# removing these individuals with 20% missingness
vcftools  --vcf  Mali_bi_Snps20_indv25_dp5.recode.vcf --remove lowDP.indv20  --recode --recode-INFO-all  --stdout |  gzip -c >  Mali_bi_Snps20_indv20_dp5

### quality filtering 

bcftools view \
--include 'MQ>=30' \
--output-file  Mali_bi_Snps20_indv20_dp5_MQ30.vcf  \
 Mali_bi_Snps20_indv20_dp5.vcf
### number of samples 830 
#### number of Snps : 151541
#number of variant before and after filtering
#### minor allele freq
vcftools  --vcf Mali_bi_Snps20_indv20_dp5_MQ30.vcf  --maf 0.02 --recode --recode-INFO-all --out Mali_bi_Snps20_indv20_dp5_MQ30_maf2.vcf
### number of snps 12177 
bcftools view -H snp_only.vcf | wc -l
  607339
Aouas-MacBook-Pro:data coulibalya$ bcftools view -H snp_missing_50_sites.vcf | wc -l
    471767

  #### individuals with a lot of missing data
vcftools --vcf snp_missing_50_sites.vcf --missing-inv --out missing_geno
  #extract individuals with 50% missing gdata
cat missing_geno.imiss | sort -rk5 |   awk 'NR>1 && $5>=0.50 {print $1}' > inviduals_50_missing_data


#removing individuals with more than 50% missing data.

vcftools --vcf snp_missing_50_sites.vcf --remove  inviduals_50_missing_data  --recode --recode-INFO-all --stdout > Genotype_file.imiss50.vcf 


#calculate the number of individual:
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' Genotype_file.imiss50.vcf 

130 out of 184 samples
 #