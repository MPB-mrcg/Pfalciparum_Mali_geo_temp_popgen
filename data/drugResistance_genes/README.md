# Drug resistance genes

SNPs in the following drug resistance genes were extracted; chloroquine resistance transporter (Pfcrt), dihydrofolate reductase (Pfdhfr), dihydropteroate synthase (Pfdhps), kelch 13 (Pfk13), 
multidrug-resistance protein 1 (Pfmdr1), multidrug-resistance 2 (Pfmdr2) and Plasmepsin II (Pfpm2) using this following command:

```(sh )
bcftools view \
-r Pf3D7_04_v3:747897-750065,Pf3D7_05_v3:955955-963095,Pf3D7_07_v3:402385-406341,Pf3D7_08_v3:547896-551057,Pf3D7_13_v3:1724600-1727877,Pf3D7_14_v3:1953474-1959147 \
-o results/Allele_frequency/drugResistance_loci.vcf \
data/raw_data/Mali_bi_Snps20_indv20_dp5_MQ30_maf2.recode.vcf.gz
```
