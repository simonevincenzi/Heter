library(marked)
library(R2admb)
library(splines)
library(dplyr)

rm(data_df)
rm(annual)
rm(data.marked)
rm(data.prov)

data_df = arrange(uppidri_df, Mark_cor, Year, Month) %>%
  filter(., !is.na(Mark_cor) & Mark_cor!="NA" & Cohort_cor!="C15" & Age_cor >=1 & !Mark_cor %in% c("BOTTLE","bottle","Dead") & Length >=115 & !grepl("9472A",Mark_cor) & !grepl("dead",Mark_cor) & !is.na(heter))


annual = 1
basic = 1

#####

minyear= 2004
maxyear = 2015

rangeofyears = minyear:maxyear

ncolonne = (length(rangeofyears)*2)  

Cohort.levels = c("C96","C97","C98","C99","C00","C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14")

data_df$Cohort_cor <- factor(data_df$Cohort_cor, levels = Cohort.levels)

colonne.names = as.character(c(62004,92004,62005,92005,62006,92006,62007,92007,
                               62008,92008,62009,92009,62010,92010,62011,92011,62012,92012,62013,92013,62014,92014,62015,92015))

data.marked = as.data.frame(matrix(0,10000,(ncolonne+1)))

colnames(data.marked) = c("id",colonne.names) #prepare data.frame for recapture data

data.marked$initial.age = rep(0,nrow(data.marked)) # initial age

data.marked$Coh = rep(0,nrow(data.marked)) # Cohort (number)

data.marked$Coh_n = rep(0,nrow(data.marked)) # Cohort (character/factor)

data.marked$heter = 0



unique.mark.data = with(data_df,unique(Mark_cor))

for (i in 1:length(unique.mark.data)) {
  
  data.prov = subset(data_df,Mark_cor == unique.mark.data[i])
  
  year.prov = data.prov$Year
  
  monthyear = paste(data.prov$Month,data.prov$Year,sep="")
  
  incl.year = colonne.names %in% monthyear
  
  data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1)
  
  
  if (annual == 0) {
    
    if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
      ### check whether the first capture was in June or September (use <= 7 to be sure)
      
      data.marked[i,"initial.age"] = min(data.prov$Age_cor) *12 ## we are considering months
      
    } else {data.marked[i,"initial.age"] = (min(data.prov$Age_cor) *12) + 3 }} else {
      
      if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
        ### check whether the first capture was in June or September (use <= 7 to be sure)
        
        data.marked[i,"initial.age"] = min(data.prov$Age_cor) *1 ## we are considering months
        
      } else {data.marked[i,"initial.age"] = (min(data.prov$Age_cor) *1) + 0.25 }
      
    }  ### initial.age is a reserved column name
  
  data.marked[i,"Coh"] = data.prov$Cohort_cor[1] ## Cohort of the fish as number
  
  data.marked[i,"Species"] = data.prov$Species[1]
  
  data.marked[i,"Pop"] = data.prov$Pop[1]
  
  data.marked[i,"Stream"] = data.prov$Stream[1]
  
  data.marked[i,"Coh_n"] = as.character(data.prov$Cohort_cor[1]) ## Cohort of the fish as character/factor 
  
  data.marked[i,"heter"] = data.prov$heter[1]
  
  data.marked[i,"id"] = data.prov$Mark_cor[1] # Tag of fish
  
  #data.marked[i,"Sex"] = as.character(data.prov$Sex_gen[1])
}

data.marked = filter(data.marked, id!=0) #data.marked[1:sum(data.marked$id>0),]

dataformerge = data.marked[,2:(ncolonne+1)]

dataformerge$ch <- do.call(paste, c(dataformerge, sep=""))

data.marked$ch = dataformerge$ch 

data.marked = data.marked[,-c(2:(ncolonne+1))]

colnames(data.marked)[1] = "Mark"


data.marked$Coh_n <- factor(data.marked$Coh_n, levels = Cohort.levels)


if (annual == 0) {
  spacebwtime = c(3,rep(c(9,3),11))} else if (annual == 1) 
  {spacebwtime = c(0.25,rep(c(0.75,0.25),11))} else if (annual == 2) {
    spacebwtime = rep(1,23)}


#spacebwtime = rep(1,23)
# Add Season covariate
season.rep = c(1,(rep(c(0,1),11)),0)  # 1 is from June to September, 0 from September to June
# [1] 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0
season.rep = c("J-S",(rep(c("S-J","J-S"),11)),"S-J")

Season=matrix(rep(season.rep,nrow(data.marked)),ncol=length(season.rep),byrow=T)

col.spacebwtime = c(1,cumsum(spacebwtime)+1)

colnames(Season)=paste("Season",col.spacebwtime,sep="")

data.marked=cbind(data.marked,Season)

### add Spawning

design.Phi=list(static=c("Coh_n","heter"),time.varying=c("Season")) 
design.p=list(static=c("Coh_n","heter"),time.varying=c("Season"))


design.parameters=list(Phi=design.Phi,p=design.p)
data.proc=process.data(data.marked, model = "CJS", time.intervals = spacebwtime)
ddl=make.design.data(data.proc,parameters=design.parameters)



fit.models=function()
{
  Phi.dot=list(formula = ~1)
  Phi.het.bs=list(formula = ~bs(heter))
  Phi.het.lm =list(formula = ~ heter)
  Phi.het.Season.mu = list(formula = ~ heter * Season)
  Phi.het.Season.ad = list(formula = ~ heter + Season)
  Phi.het.Season.mu.bs = list(formula = ~ bs(heter) * Season)
  Phi.het.Season.add.bs = list(formula = ~ bs(heter) + Season)
  Phi.Season = list(formula = ~ Season)
  
  
  p.Age= list(formula=~1)
  ##
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T, use.admb = F)
  return(results)
}
test.heter=fit.models()
ddl.heter = ddl


test.heter.mod = predict(test.heter[[2]],ddl = ddl.heter,se = T)
#######
### PLOT

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