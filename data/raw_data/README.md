# Data collection

### Study sites
Human malaria parasite infected blood samples were collected between 2007 and 2017 from seven malaria endemic sites in southern Mali, where the high malaria transmission season overlaps with the annual rainy season that occurs from June to December. 

The sites included:
Bamako: 162 samples analysed were collected mostly in 2013. Three samples were from 2012.
Bougoula-Hameau (Bougoula): samples from this site were collected in 2015 (22) and 2017 (19).
Dangassa: Most of the samples from here were from 2016 (158), with 8 from 2017
Faladje: Samples (216) spanned the decade from 2007 to 2017, with the largest number (125) from 2013.
Kenieroba: 89 samples from this site were all collected in 2016  
Kolle: This site presented the earliest sample set (47) collected in 2007. 
Finally, Nioro-du-Sahel (Nioro) 109 samples were collected in 2015 and 2016.

### Data pre-processing 
From each participant, a dried blood spot (DBS) sample was collected by spotting whole blood obtained by finger prick onto filter paper cards. 
The collection was done according to the SPOTMalaria protocol (www.malariagen.net/projects/SPOTMalaria), a project led by the Malaria Genetic Epidemiology Network (MalariaGen) to harness genomic technologies in monitoring the evolution of malaria parasites. 
The filter papers were subjected to DNA extraction using the Qiagen DNA extraction protocol and the genome of P. *falciparum* in each sample was sequenced at the Wellcome Sanger Institute as previously described. 
The short sequence read (fastq) files were mapped to the P. *falciparum* 3D7 reference genome with BWA-MEM and the resultant BAM files were used in variant calling following GATK best practice custom pipelines by the *MalariaGEN P. falciparum community* and **[Pf3K](https://www.malariagen.net/parasite/pf3k)** projects. 
For downstream analysis, vcf files of variants were filtered with vcftools (vcftools.github.io) to retain only biallelic coding SNP variants supported by at least 5 reads and with a mapping quality greater than 30.
Using R scripts, we iterated through filtration steps to remove SNPs or samples with excessive missingness and retaining only SNPs genotyped for at least 80% of individuals and individuals with less than 20% of missing SNP calls. 
We then retained only high-quality PASS bi-allelic SNPs with a minor allele frequency of at least 0.02.
The final dataset included __12,177__ SNPs and __830__ samples, from an original set of 3,844,304 variants and 1270 samples.