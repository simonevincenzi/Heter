
library(marked)
library(R2admb)
library(splines)
library(dplyr)

annual = 1

zadla_df = arrange(zadla_df, Mark_cor, Year, Month)


data = filter(zadla_df, !is.na(heter) & !is.na(Mark_cor) & Age_cor>0 & !is.na(Month) &
                Year >=2006)



Cohort.levels = c("C98","C99","C00","C01","C02","C03","C04","C05",
"C06","C07","C08","C09","C10","C11","C12","C13")

data$Cohort <- factor(data$Cohort, levels = Cohort.levels)

minyear= 2006
maxyear = 2014

rangeofyears = minyear:maxyear

ncolonne = (maxyear - minyear) + 1
#ncolonne corresponds to the number of sampling occasions. 9 in this case

colonne.names = as.character(minyear:maxyear)
#colonne.names are the sampling occasions (format yyyy)

data.marked = as.data.frame(matrix(0,10000,(ncolonne+1))) #prepare data.frame for recapture data

colnames(data.marked) = c("id",colonne.names) #prepare data.frame for recapture data

data.marked$initial.age = rep(0,nrow(data.marked)) # initial age

data.marked$Coh = rep(0,nrow(data.marked)) # Cohort (number)

data.marked$Coh_n = rep(0,nrow(data.marked)) # Cohort (character/factor)

data.marked$heter = 0


unique.mark.data = with(data,unique(Mark_cor))  # unique tagged fish

for (i in 1:length(unique.mark.data)) { #loop over unique tagged fish

data.prov = subset(data,Mark_cor == unique.mark.data[i]) ## data relative to the tagged fish

year.prov = data.prov$Year ## years in which the fish was sampled

incl.year = rangeofyears %in% year.prov  

data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1) # for the tagged fish in the columns
# of the dataframe I assign 0 if the fish was not sampled, 1 otherwise

if (annual == 0) {

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

data.marked[i,"heter"] = data.prov$heter[1]

data.marked[i,"id"] = data.prov$Mark_cor[1] # Tag of fish


}

data.marked = filter(data.marked, id!=0) #data.marked[1:sum(data.marked$id>0),]

dataformerge = data.marked[,2:(ncolonne+1)]

dataformerge$ch <- do.call(paste, c(dataformerge, sep=""))

data.marked$ch = dataformerge$ch 

data.marked = data.marked[,-c(2:(ncolonne+1))]

colnames(data.marked)[1] = "Mark"

data.marked$Coh_n <- factor(data.marked$Coh_n, levels = Cohort.levels)



if (annual == 0) {  
spacebwtime = rep(12,length(unique(data$Year))-1)} else {spacebwtime = rep(1,length(unique(data$Year))-1)}

#data.marked = assign.vb.par.f(marked.df = data.marked, vb.pred.df = zadla0_noP.synth$pred_vb)

##add a Flood variable. There was a flood in 2007

flood.rep = rep(0,length(unique(data$Year)))

flood.rep[c(2,7)] = 1 # 2 year is 2007, the year of the flood and 2012 another possible one

Flood=matrix(rep(flood.rep,nrow(data.marked)),ncol=length(flood.rep),byrow=T)
if (annual == 0) {
col.spacebwtime = c(12,cumsum(spacebwtime)+1) } else {col.spacebwtime = c(1,cumsum(spacebwtime)+1)}

colnames(Flood)=paste("Flood",col.spacebwtime,sep="")
# at the end of each column of Flood I need a number starting from 1 indicating the 
# cumulative time

data.marked=cbind(data.marked,Flood)  # bind the flood dataset

#########


design.Phi=list(static=c("Coh_n","heter"),time.varying=c("Flood"))
design.p=list(static=c("Coh_n","heter"),time.varying=c("Flood"))
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
  p.time = list(formula = ~1)
  #p.Age.bs=list(formula=~bs(Age))
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T)
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


test.gg <- ggplot(test.heter.mod$Phi, aes(y = estimate, x = heter, group = Flood, col = Flood)) + 
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





















#best.coh = as.numeric(row.names(data.models.zadla.annual$model.table))[
#  which(grepl("Coh",data.models.zadla.annual$model.table$model))][1]

#predict(data.models.zadla.annual[[best.coh]])$Phi

#with(predict(data.models.zadla.annual[[13]])$Phi, 
#     plot(estimate ~ as.numeric(time), col = "green", type = "p", 
 #         ylim = c(0,0.5), xlab = zadla.time))

#pp = predict(data.models.zadla.annual[[13]])$Phi

#with(pp, 
 #    plot(estimate ~ time, col = "green", type = "p", 
  #        ylim = c(0,0.5)))

#pp$time = zadla.time
#with(predict(data.models.zadla.annual[[14]])$Phi, points(estimate ~ time,col = "red"))


#zadla.time = c("06", "07","08","09","10","11","12","13")

#as.numeric(row.names(data.models.zadla.annual$model.table))


#predPhi = predict(data.models[[3]], se = TRUE)
#pred.sm <- data.frame(spline(predPhi$Age, predPhi$estimate, n=100))
#pred.sm$lcl <- spline(predPhi$Age, predPhi$lcl, n=100)$y
#pred.sm$ucl <- spline(predPhi$Age, predPhi$ucl, n=100)$y

#p <- ggplot(aes(x=x, y=y), data=pred.sm)
#p <- p + geom_line() + geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.25) + 
 # ggtitle("Female dipper survival bs(Time) RMark\n") + xlab("Time") + ylab("Phi(t)")










