### Data preparation for survival and growth

### Load libraries
library(tidyverse)
library(readr)
library(lubridate)

### function to add heterozygosity

add_heter.f = function(tag_data,geno_data) {
  source("./scripts/heter.r")
  pos.loci.excl = which(grepl("Sma", colnames(geno_data)))
  geno_data$heter = apply(as.matrix(geno_data[,pos.loci.excl]),  1, FUN = heter.f)
  tag_data$heter = geno_data$heter[match(tag_data$Mark_cor,geno_data$SAMPLE_ID)]
  return(tag_data)
}



### Read data

loidri_df = read_csv("./raw_data/loidri_df_pieced.csv", guess_max = 20000)
loidri_df$Date = as.Date(loidri_df$Date,format = "%m/%d/%Y") # Y is year with century
loidri_df = loidri_df %>%
  arrange(.,Mark_cor,Date)
loidri_df = add_heter.f(loidri_df,read_csv("./raw_data/loidri_genotype_repo.csv"))

uppidri_df = read_csv("./raw_data/uppidri_df_pieced.csv", guess_max = 20000)
uppidri_df$Date = as.Date(uppidri_df$Date,format = "%m/%d/%Y") # Y is year with century
uppidri_df = uppidri_df %>%
  arrange(.,Mark_cor,Date)
uppidri_df = add_heter.f(uppidri_df,read_csv("./raw_data/uppidri_genotype_repo.csv"))



zak_df = read_csv("./raw_data/zak_df_pieced.csv", guess_max = 20000)
zak_df$Date = as.Date(zak_df$Date,format = "%m/%d/%Y") # Y is year with century
zak_df = zak_df %>%
  arrange(.,Mark_cor,Date)
zak_df = add_heter.f(zak_df,read_csv("./raw_data/zak_genotype_repo.csv"))

trebu_df = read_csv("./raw_data/trebu_df_pieced.csv", guess_max = 20000)
trebu_df$Date = as.Date(trebu_df$Date,format = "%m/%d/%Y") # Y is year with century
trebu_df = trebu_df %>%
  arrange(.,Mark_cor,Date)
trebu_df = add_heter.f(trebu_df,read_csv("./raw_data/trebu_genotype_repo.csv"))

zadla_df = read_csv("./raw_data/zadla_df_pieced.csv", guess_max = 20000)
zadla_df$Date = as.Date(zadla_df$Date,format = "%m/%d/%Y") # Y is year with century
zadla_df = zadla_df %>%
  arrange(.,Mark_cor,Date)

