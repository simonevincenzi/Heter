January 26, 2017
# Data and code for the manuscript "Stronger effects of heterozygosity on survival in harsher environments"

## 1. Goals of the study
In this work, we test the hypothesis of effects of heterozygosity on survival varying with environmental conditions using as a model system 6 marble trout (Salmo mamoratus) populations living in Western Slovenia. Average survival probabilities of marble trout largely differ among populations and 3 of them have been affected by flash floods causing massive mortalities. 

Our goal was to estimate for each marble trout population the effects of heterozygosity on variation in survival while accounting for year of birth, and season and occurrence of flash floods where applicable. The hypothesis are: (a) stronger (positive) effects of heterozygosity on survival in populations with lower survival probabilities (Zadla) or with marble trout living in sympatry with rainbow trout (LIdri), (b) no effect of heterozygosity on surviving flash floods and (c) greater effects of heterozygosity for fish born after the occurrence of flash floods (Zak, Lipo after 2009) than born before the flash floods, due to higher production of young post-flood and consequent higher competition for resources (food, shelter). For (c), we also predicted stronger effect in Zak than in Lipo due to higher density/stronger competition of fish born after the extreme event in Zak.  For Zadla, the sample size was too small for testing hypothesis (c).

## 2. Data

Tag-recapture and genotyping data are used. Both are in the folder `raw_data`. Tag data are as `NameOfPopulation_pieced_df.csv` (e.g., for Lipovscek `lipo_pieced_df.csv`) and genotyping data are as `NameOfPopulation_genotype_repo.csv` (e.g., `lipo_genotype_repo.csv`)

## 3. Code

R scripts are in the folder `scripts`. The script `data_prep.r` reads the tag-recapture data and uses genotyping data to assign heterozygosity to each individual with at least 50 loci genotyped. The scripts `NameOfPopulation_marked.r` (e.g., `lipo_marked.r`) prepare data and fit models for tag-recapture analysis of survival using the R package `marked`. The script `Plot_survival_heter.r` plots the relationship heterozygosity-survival for the four populations for which models with heterozygosity were strongly supported and the 3-panel plot in Supplementary Material on the relationhsip between heterozygosity and probability of surviving flash floods. The script `maf.r` calculates the minimum allele frequency for each locus for all populations, along with mean and sd across all loci.

## 4. Manuscript

The manuscript (currently under review) can be [downloaded](http://simonevincenzi.com/Publications/Het_paper_Feb2017_complete.pdf)
