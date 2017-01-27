
library(marked)
library(R2admb)
library(splines)
library(dplyr)

# delete some objects that have been potentially used by other scripts

rm(data_df)
rm(annual)
rm(data.marked)
rm(data.prov)

# annual or monthly time frame

annual = 1

trebu_df = arrange(trebu_df, Mark_cor, Year, Month)

# filter tag-recapture data, excluding NAs and samples of age 0 (not tagged)

data_df = filter(trebu_df, !is.na(heter) & !is.na(Mark_cor) & Age_cor>0 & !is.na(Month) &
                Year >=2006)

# re-factor Cohort

Cohort.levels = c("C98","C99","C00","C01","C02","C03","C04","C05",
                  "C06","C07","C08","C09","C10","C11","C12","C13")

data_df$Cohort <- factor(data_df$Cohort, levels = Cohort.levels)

# first and last year

minyear= 2006
maxyear = 2014

rangeofyears = minyear:maxyear

ncolonne = (maxyear - minyear) + 1 # ncolonne corresponds to the number of sampling occasions. 9 in this case

colonne.names = as.character(minyear:maxyear) # colonne.names are the sampling occasions (format yyyy)

data.marked = as.data.frame(matrix(0,10000,(ncolonne+1))) # prepare data.frame for recapture data

# assign column names

colnames(data.marked) = c("id",colonne.names) #prepare data.frame for recapture data

data.marked$initial.age = rep(0,nrow(data.marked)) # initial age

data.marked$Coh = rep(0,nrow(data.marked)) # Cohort (number)

data.marked$Coh_n = rep(0,nrow(data.marked)) # Cohort (character/factor)

data.marked$heter = 0 # Heterozygosity


unique.mark.data = with(data_df,unique(Mark_cor))  # unique tagged fish

for (i in 1:length(unique.mark.data)) { #loop over unique tagged fish
  
  data.prov = subset(data_df,Mark_cor == unique.mark.data[i]) ## data relative to the tagged fish
  
  year.prov = data.prov$Year ## years in which the fish was sampled
  
  incl.year = rangeofyears %in% year.prov # sampling occasions (true or false) in which the individual was sampled
  
  data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1) # for the tagged fish in the columns
  # of the dataframe I assign 0 if the fish was not sampled, 1 otherwise
  
  if (annual == 0) { # assign age based on months or years
    
    if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
      ### check whether the first capture was in June or September (use <= 7 to be sure)
      
      data.marked[i,"initial.age"] = min(data.prov$Age) *12 ## we are considering months
      
    } else {data.marked[i,"initial.age"] = (min(data.prov$Age) *12) + 3 }} else {
      
      if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
        ### check whether the first capture was in June or September (use <= 7 to be sure)
        
        data.marked[i,"initial.age"] = min(data.prov$Age) *1 ## we are considering months
        
      } else {data.marked[i,"initial.age"] = (min(data.prov$Age) *1) + 0.25 }
      
    }  ### initial.age is a reserved column name
  
  data.marked[i,"Coh"] = data.prov$Cohort[1] ## Cohort of the fish as number
  
  data.marked[i,"Coh_n"] = as.character(data.prov$Cohort[1]) ## Cohort of the fish as character/factor 
  
  data.marked[i,"heter"] = data.prov$heter[1] # heterozygosity
  
  data.marked[i,"id"] = data.prov$Mark_cor[1] # Tag of fish
  
  
}

### I join the columns with the 1/0 for sampled/not sampled

data.marked = filter(data.marked, id!=0) #data.marked[1:sum(data.marked$id>0),]

dataformerge = data.marked[,2:(ncolonne+1)]

dataformerge$ch <- do.call(paste, c(dataformerge, sep=""))

data.marked$ch = dataformerge$ch 

data.marked = data.marked[,-c(2:(ncolonne+1))]

colnames(data.marked)[1] = "Mark"

###

data.marked$Coh_n <- factor(data.marked$Coh_n, levels = Cohort.levels)

# I need to assign the length of sampling occasions in form of a vector (in years or months). Vector has length number of sampling occasions - 1

if (annual == 0) {  
  spacebwtime = rep(12,length(unique(data_df$Year))-1)} else {spacebwtime = rep(1,length(unique(data_df$Year))-1)}


### create the objects for marked


design.Phi=list(static=c("Coh_n","heter"))
design.p=list(static=c("Coh_n","heter"))
design.parameters=list(Phi=design.Phi,p=design.p)
data.proc=process.data(data.marked, model = "CJS", time.interval = spacebwtime)
ddl=make.design.data(data.proc,parameters=design.parameters)

### test best capture model

# fit.models=function()
# {
#   Phi.full=list(formula = ~ Coh_n * heter)
#   p.Coh.=list(formula = ~ Coh_n )
#   p.time = list(formula = ~ time)
#   p.time.Coh = list(formula = ~ time + Coh_n)
#   p.dot = list(formula = ~1)
#   #p.Age.bs=list(formula=~bs(Age))
#   cml=create.model.list(c("Phi","p"))
#   results=crm.wrapper(cml,data=data.proc, ddl=ddl,
#                       external=FALSE,accumulate=FALSE, hessian = T)
#   return(results)
# }
# trebu.mod.heter=fit.models()
# trebu.mod.ddl.heter = ddl

### Fit survival models (bs is splines)

fit.models=function()
{
  Phi.dot=list(formula = ~ 1)
  Phi.int.rel.bs=list(formula = ~bs(heter))
  Phi.int.rel.lm =list(formula = ~ heter)
  Phi.Coh.het.bs.ad=list(formula = ~ Coh_n + bs(heter))
  Phi.Coh.int.rel.mu =list(formula = ~ Coh_n * heter)
  Phi.Coh.=list(formula = ~ Coh_n)
  
  p.time = list(formula = ~ time)
  
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T)
  return(results)
}
trebu.mod.heter=fit.models() # object with models fitted, AIC etc.
trebu.mod.ddl.heter = ddl # design data

test.heter.mod = predict(trebu.mod.heter[[6]], ddl = trebu.mod.ddl.heter, se = T) # best model


### plot the best model

size.title = 20
line.lwd = 1.2
size.label.x = 22
size.text.x = 20
size.point = 6
size.label.y = 22
size.text.y = 20
size.legend.text = 15
size.legend.title = 20
unit.legend.h = 1.8
unit.legend.w = 1.8
size.ann = 10
colour.axis = "gray20"
colour.theme = "black"
colour.axis.line = "gray20"
colour.line = "gray50"
label.T = "Heterozygosity" 


test.gg <- ggplot(test.heter.mod$Phi, aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray40", lty = 2) + 
  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title), 
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.margin = unit(c(1,2,1,1), "cm"),
        axis.title.x = element_text(size=size.label.x,vjust=-1),
        axis.text.x  = element_text(size=size.text.x, vjust = 0.5),
        axis.title.y = element_text(size=size.label.x, vjust = 2),
        axis.text.y  = element_text(size=size.text.x),
        legend.title = element_blank(),
        legend.text = element_text( size = size.legend.text)
  )  +

  labs(y = bquote(phi)) +
  labs(x = label.T)

test.gg