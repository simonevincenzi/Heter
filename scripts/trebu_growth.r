library(dplyr)
library(grid)
library(ggplot2)


source("./scripts/data_prep.r")



trebu.growth.prep = arrange(trebu_df,Mark_cor,Year,Month)

#trebu.growth.prep = assign.sex.repo.f(
# repo.df = repo.cut, stream.df = trebu.growth.prep,stream = "LIdri")$data.df

### Order cohort levels chronologically  
Cohort.levels = c("U","C00","C01","C02","C03","C04","C05","C06","C07",
                  "C08","C09","C10","C11","C12","C13","C14")



########## September

data_growth = trebu.growth.prep %>% 
  filter(.,!is.na(Mark_cor) & Month == 9 & !Cohort %in% c("C13","C14","U") & !is.na(heter)) %>%
  select(.,-Mark, -Cohort, -Age, -Date)



colnames(data_growth)[which(colnames(data_growth) == "Mark_cor")] = "Mark"
colnames(data_growth)[which(colnames(data_growth) == "Age_cor")] = "Age"
colnames(data_growth)[which(colnames(data_growth) == "Cohort_cor")] = "Cohort"
data_growth$Cohort = factor(data_growth$Cohort)

trebu.AIC_s.df = data.frame(Linf = rep(0,4), k = rep(0,4), AIC = rep(0,4))
caic = 1




linf_var = "Const"
k_var = "Const"
source("./scripts/m_grow3_gen.r")  #model with no predictors
#system("./m_grow3 -est")  ##run model (-est means that standard errors of estimates are 
#not estimated)
trebu_Linf_const_k_const_s.synth = synth.list.3

trebu.AIC_s.df$Linf[caic] = linf_var 
trebu.AIC_s.df$k[caic] = k_var
trebu.AIC_s.df$AIC[caic] = synth.list.3$AIC
caic = caic + 1


linf_var = "Cohort"
k_var = "Cohort"
source("./scripts/m_grow3_gen.r")  #model with no predictors
#system("./m_grow3 -est")  ##run model (-est means that standard errors of estimates are 
#not estimated)
trebu_Linf_coh_k_coh_s.synth = synth.list.3

trebu.AIC_s.df$Linf[caic] = linf_var 
trebu.AIC_s.df$k[caic] = k_var
trebu.AIC_s.df$AIC[caic] = synth.list.3$AIC
caic = caic + 1

linf_var = "Heter + Cohort"
k_var = "Heter + Cohort"
source("./scripts/m_grow3_gen.r")  #model with no predictors
#system("./m_grow3 -est")  ##run model (-est means that standard errors of estimates are 
#not estimated)
trebu_Linf_hetcoh_k_hetcoh_s.synth = synth.list.3

trebu.AIC_s.df$Linf[caic] = linf_var 
trebu.AIC_s.df$k[caic] = k_var
trebu.AIC_s.df$AIC[caic] = synth.list.3$AIC
caic = caic + 1

linf_var = "Heter"
k_var = "Heter"
source("./scripts/m_grow3_gen.r")  #model with no predictors
#system("./m_grow3 -est")  ##run model (-est means that standard errors of estimates are 
#not estimated)
trebu_Linf_het_k_het_s.synth = synth.list.3

trebu.AIC_s.df$Linf[caic] = linf_var 
trebu.AIC_s.df$k[caic] = k_var
trebu.AIC_s.df$AIC[caic] = synth.list.3$AIC

trebu.AIC_s.df = arrange(trebu.AIC_s.df, AIC)