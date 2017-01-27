### Estimate minimum allele frequency 

# load libraries
library(readr)

# read the genotype data

lipo_geno_df = read_csv("./raw_data/lipo_genotype_repo.csv")
trebu_geno_df = read_csv("./raw_data/trebu_genotype_repo.csv")
zak_geno_df = read_csv("./raw_data/zak_genotype_repo.csv")
loidri_geno_df = read_csv("./raw_data/loidri_genotype_repo.csv")
uppidri_geno_df = read_csv("./raw_data/uppidri_genotype_repo.csv")
uppidri_geno_df = read_csv("./raw_data/zadla_genotype_repo.csv")

geno.df = trebu_geno_df

col_loci = which(grepl(pattern = "Sma", colnames(geno.df))) # columns with alleles


cont = 1 # count
MAF = rep(0,length(seq(min(col_loci),max(col_loci),2))) # prepare vector of MAF for each locus

for(i in seq(min(col_loci),max(col_loci),2)) {
 # print(i)
  allele = c(as.matrix(geno.df[,i]),as.matrix(geno.df[,i+1]))
  allele = allele[allele!=0]
  MAF[cont] = min(table(allele)/sum(table(allele)))[1]
  cont = cont+1
  
}

mean(MAF[MAF!=1])
sd(MAF[MAF!=1])
