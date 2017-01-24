
library(marked)
library(R2admb)
library(splines)
library(dplyr)

rm(data_df)
rm(annual)
rm(data.marked)
rm(data.prov)

annual = 1
#zak_df = readRDS("zak_df.RDS")

zak_df = arrange(zak_df, Mark_cor, Year, Month)

zak_df$Mark_cor[which(zak_df$Mark_cor == "snakelike")] = zak_df$Mark[which(zak_df$Mark_cor == "snakelike")]

data_df = filter(zak_df, !is.na(heter))


#### Cohort pre and post flood

data_df$Coh.p = ifelse(data_df$Cohort == "P", "P", "O")

Coh.p.levels = c("P","O")

data_df$Coh.p = factor(data_df$Coh.p, levels = Coh.p.levels)

data_df$Coh.pflood = rep(0,nrow(data_df))

for (i in 1:nrow(data_df)) {
  
  if (data_df$Cohort[i] == "P") {
    
    data_df$Coh.pflood[i] = "P" } else if (data_df$Cohort[i] %in% c("C98","C99","C00","C01","C02",
                                                              "C03","C04","C05","C06")) {
      
      data_df$Coh.pflood[i] = "Pre" } else { data_df$Coh.pflood[i] = "Post"}}

Coh.pflood.levels = c("P","Pre","Post")

data_df$Coh.pflood = factor(data_df$Coh.pflood,levels = Coh.pflood.levels)

levels(data_df$Coh.pflood) = c("Pre","Pre","Post")



minyear= 1996
maxyear = 2014

rangeofyears = minyear:maxyear

ncolonne = (maxyear - minyear) + 1    ### number of sampling occasions (19)

colonne.names = as.character(minyear:maxyear)

data.marked = as.data.frame(matrix(0,10000,(ncolonne+1)))

colnames(data.marked) = c("id",colonne.names) #prepare data.frame for recapture data

data.marked$initial.age = 0 # initial age

data.marked$Coh = 0 # Cohort (number)

data.marked$heter = 0


unique.mark.data = with(data_df,unique(Mark_cor))

for (i in 1:length(unique.mark.data)) {
  
  data.prov = subset(data_df,Mark_cor == unique.mark.data[i])
  
  year.prov = data.prov$Year
  
  incl.year = rangeofyears %in% year.prov
  
  data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1)
  
  
  data.marked[i,"initial.age"] = min(data.prov$Age) *1 ## we are considering months
  
  data.marked[i,"Coh"] = data.prov$Cohort_cor[1] ## Cohort of the fish as number
  
  data.marked[i,"Coh_n"] = as.character(data.prov$Cohort_cor[1])
  
  data.marked[i,"Coh.pflood"] = as.character(data.prov$Coh.pflood[1])
  
  data.marked[i,"heter"] = data.prov$heter[1]
  
  data.marked[i,"id"] = data.prov$Mark_cor[1] # Tag of fish
  
}


data.marked = filter(data.marked, id!=0) #data.marked[1:sum(data.marked$id>0),]

dataformerge = data.marked[,2:(ncolonne+1)]

dataformerge$ch <- do.call(paste, c(dataformerge, sep=""))

data.marked$ch = dataformerge$ch 

data.marked = data.marked[,-c(2:(ncolonne+1))]

colnames(data.marked)[1] = "Mark"

if (annual == 0) {
  spacebwtime = rep(12,length(unique(data_df$Year))-1)} else {spacebwtime = rep(1,length(unique(data_df$Year))-1)}


col.spacebwtime = c(1,cumsum(spacebwtime)+1)

##add a Flood variable

flood.rep = rep(0,length(unique(data_df$Year)))

flood.rep[12] = 1  # Flood occurred in 2007

Flood=matrix(rep(flood.rep,nrow(data.marked)),ncol=length(flood.rep),byrow=T)

col.spacebwtime = c(1,cumsum(spacebwtime)+1)

colnames(Flood)=paste("Flood",col.spacebwtime,sep="")

data.marked=cbind(data.marked,Flood)

#####

design.Phi=list(static=c("heter","Coh_n","Coh.pflood"), time.varying = ("Flood"))
design.p=list(static=c("heter","Coh_n","Coh.pflood"), time.varying = ("Flood"))
design.parameters=list(Phi=design.Phi,p=design.p)
data.proc=process.data(data.marked, model = "CJS", time.interval = spacebwtime)
ddl=make.design.data(data.proc,parameters=design.parameters)

#data.proc=process.data(data.marked,model="CJS",groups="Coh")
#data.ddl=make.design.data(data.proc)


fit.models=function()
{
  Phi.dot=list(formula = ~1)
  Phi.int.rel.bs=list(formula = ~bs(heter))
  Phi.int.rel.lm =list(formula = ~ heter)
  Phi.het.Flood.add.bs=list(formula = ~bs(heter) + Flood)
  Phi.heter.Flood.add =list(formula = ~ heter + Flood)
  Phi.Coh.p.flood.heter.Flood.ad=list(formula = ~ heter + Coh.pflood + Flood)
  Phi.Coh.p.flood.heter.Flood.mu=list(formula = ~ heter * Coh.pflood + Flood)
  Phi.Coh.p.flood.Flood.heter.ad=list(formula = ~ Coh.pflood + Flood)
  Phi.Coh.p.flood.Flood.mu=list(formula = ~ Coh.pflood * Flood)
  Phi.Flood=list(formula = ~ Flood)
  
  p.Age = list(formula = ~Age)
  
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T,  burnin=1000,iter=5000)
  return(results)
}
zak.mod.heter=fit.models()
zak.mod.ddl.heter = ddl
test.heter.mod = predict(zak.mod.heter[[4]], ddl = zak.mod.ddl.heter , se = T)


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



test.gg <- ggplot(filter(test.heter.mod$Phi, Flood == 0), aes(y = estimate, x = heter, group = Coh.pflood, lty = Coh.pflood)) + 
  #geom_point(size = size.point,lwd = 4) +
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray40", lty = 2) + 
  ggtitle("Zak") +
  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
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
 scale_y_continuous(limits = c(0,1)) +
  labs(y = bquote(phi)) +
  labs(x = label.T)

test.gg








