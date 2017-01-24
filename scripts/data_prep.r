### Data preparation for survival and growth

### Load libraries
library(tidyverse)
library(readr)
library(lubridate)

### function to add heterozygosity

add_heter.f = function(tag_data,geno_data) {
  source("./scripts/heter.r")
  pos.loci.excl = which(grepl("Sma", colnames(geno_data)))
  geno_data$heter = NA
  geno_data$heter = apply(as.matrix(geno_data[,pos.loci.excl]),  1, FUN = heter.f)
  
 geno_data$n_loci = apply(as.matrix(geno_data[,pos.loci.excl]),  1, FUN = function (x) {length(x[x>0])})
  
 tag_data$heter = NA
  tag_data$heter = geno_data$heter[match(tag_data$Mark_cor,geno_data$SAMPLE_ID)]
  tag_data$n_loci = geno_data$n_loci[match(tag_data$Mark_cor,geno_data$SAMPLE_ID)]
  
  return(tag_data)
}



### Read data

loidri_df = read_csv("./raw_data/loidri_df_pieced.csv", guess_max = 20000)
loidri_df$Date = as.Date(loidri_df$Date,format = "%m/%d/%Y") # Y is year with century
loidri_df = loidri_df %>%
  arrange(.,Mark_cor,Date)
loidri_df = add_heter.f(loidri_df,read_csv("./raw_data/loidri_genotype_repo.csv"))

uppidri_df = read_csv("./raw_data/uppidri_df_pieced.csv", guess_max = 20000)
source("./scripts/correct_uppidri_database.r")
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
zadla_df = add_heter.f(zadla_df,read_csv("./raw_data/zadla_genotype_repo.csv"))

lipo_df = read_csv("./raw_data/lipo_df_pieced.csv", guess_max = 20000)
lipo_df$Date = as.Date(lipo_df$Date,format = "%m/%d/%Y") # Y is year with century
colnames(lipo_df)[which(colnames(lipo_df) == "Mark")] = "Mark_cor"
lipo_df = lipo_df %>%
  arrange(.,Mark_cor,Date)
lipo_df = add_heter.f(lipo_df,read_csv("./raw_data/lipo_genotype_repo.csv"))

lipo_df %>%
  filter(.,!is.na(Mark_cor)) %>%
  group_by(Mark_cor) %>%
  summarise(heter_u = mean(heter)) %>%
  ggplot(.,aes(x = heter_u)) +
  geom_histogram()


 