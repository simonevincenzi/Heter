
library(marked)
library(R2admb)
library(splines)
library(dplyr)



annual = 1
#zak_df = readRDS("zak_df.RDS")

zak_df = arrange(zak_df, Mark_cor, Year, Month)

zak_df$Mark_cor[which(zak_df$Mark_cor == "snakelike")] = zak_df$Mark[which(zak_df$Mark_cor == "snakelike")]

data = filter(zak_df, !is.na(heter))

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


unique.mark.data = with(data,unique(Mark_cor))

for (i in 1:length(unique.mark.data)) {
  
  data.prov = subset(data,Mark_cor == unique.mark.data[i])
  
  year.prov = data.prov$Year
  
  incl.year = rangeofyears %in% year.prov
  
  data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1)
  
  
  data.marked[i,"initial.age"] = min(data.prov$Age) *1 ## we are considering months
  
  data.marked[i,"Coh"] = data.prov$Cohort_cor[1] ## Cohort of the fish as number
  
  data.marked[i,"Coh_n"] = as.character(data.prov$Cohort_cor[1])
  
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
  spacebwtime = rep(12,length(unique(data$Year))-1)} else {spacebwtime = rep(1,length(unique(data$Year))-1)}


col.spacebwtime = c(1,cumsum(spacebwtime)+1)

##add a Flood variable

flood.rep = rep(0,length(unique(data$Year)))

flood.rep[12] = 1  # Flood occurred in 2007

Flood=matrix(rep(flood.rep,nrow(data.marked)),ncol=length(flood.rep),byrow=T)

col.spacebwtime = c(1,cumsum(spacebwtime)+1)

colnames(Flood)=paste("Flood",col.spacebwtime,sep="")

data.marked=cbind(data.marked,Flood)

#####

design.Phi=list(static=c("heter","Coh_n"), time.varying = ("Flood"))
design.p=list(static=c("heter","Coh_n"), time.varying = ("Flood"))
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
  # Phi.int.rel.bs.time=list(formula = ~bs(heter) + time)
  # Phi.int.rel.lm.time =list(formula = ~ heter + time)
  Phi.int.rel.bs.Flood=list(formula = ~bs(heter) + Flood)
  Phi.int.rel.lm.Flood =list(formula = ~ heter + Flood)
  #Phi.Coh.Flood.Temp.mult=list(formula=~Coh_n*Flood*Temperature)
  #Phi.Coh.Flood.Temp.add=list(formula=~Coh_n+Flood+Temperature)
  #Phi.Sex=list(formula=~Sex)
  #Phi.Sex.Age.mult=list(formula=~Sex*Age)
  #Phi.Sex.Age.add=list(formula=~Sex+Age)
  #Phi.linf = list(formula=~linf)
  #Phi.expL3 = list(formula=~expL3)
  #Phi.expL2 = list(formula=~expL2)
  #Phi.expL1 = list(formula=~expL1)
  #Phi.age=list(formula=~age)
  #p.dot=list(formula=~1)
  #p.Density = list(formula = ~ Density)
  #p.Coh=list(formula=~Coh_n)
  #p.Flood=list(formula=~Flood)
  #p.Sex=list(formula=~Sex)
  #p.time=list(formula=~time)
 # p.dot=list(formula=~1)
  p.time = list(formula = ~time)
  #p.Age.bs=list(formula=~bs(Age))
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T,  burnin=1000,iter=5000)
  return(results)
}
test.heter=fit.models()
ddl.heter = ddl
test.heter.mod = predict(test.heter[[5]], ddl = ddl.heter, se = T)


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
  #geom_point(size = size.point,lwd = 4) +
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray40", lty = 2) + 
  #ggtitle("Age") +
  #geom_text(size = size.text) +
  # ggtitle("Mean Length of 0+ in September") +
  #geom_text(vjust = 0.4,hjust = -0.4) + 
  #scale_shape_manual(values=c(16:17), guide = F) +
  #scale_colour_manual(values = c("gray60","black"), guide = F) +
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
  #guides(pch = guide_legend(override.aes = list(size=9)),legend.position = c(0.5,0.5)) +
 # scale_x_continuous() +
#  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  
  #guides(shape = F) +
  labs(y = bquote(phi)) +
  labs(x = label.T)

test.gg

ggsave(file="Plot.zak.surv.heter.pdf", width=11, height=10)





