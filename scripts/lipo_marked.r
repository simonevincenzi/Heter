library(marked)

library(R2admb)

library(splines)

library(dplyr)

library(ggplot2)

rm(data_df)
rm(annual)
rm(data.marked)
rm(data.prov)


annual = 1



data_df = arrange(lipo_df, Mark_cor, Year, Month) %>%
  filter(.,!is.na(Mark_cor) & !is.na(Cohort) & Age!=0 & !Sector %in% c("F","G","H","I","J","K") & !is.na(heter) & Year < 2015) 



#####

minyear= 2006
maxyear = 2014

rangeofyears = minyear:maxyear

ncolonne = (length(rangeofyears)*2) - 1 # September to September

Cohort.levels = c("C99","C00","C01","C02","C03","C04","C05","C06","C07","C08","C10","C11","C12","C13","C14")

data_df$Cohort <- factor(data_df$Cohort, levels = Cohort.levels)

data$Coh.pflood = rep(0,nrow(data)) # group for cohorts born before and after the flood

# from 2007 on they are after the flood

for (i in 1:nrow(data_df)) {
  
  if (data_df$Cohort[i] %in% c("C99","C00","C01","C02","C03","C04","C05","C06")) {
    
    data_df$Coh.pflood[i] = "Pre" } else { data_df$Coh.pflood[i] = "Post"}}

Coh.pflood.levels = c("Pre","Post")

data_df$Coh.pflood = factor(data_df$Coh.pflood,levels = Coh.pflood.levels)

colonne.names = as.character(c(92006,62007,92007,
                               62008,92008,62009,92009,62010,92010,62011,92011,62012,92012,62013,92013,62014,92014))

data.marked = as.data.frame(matrix(0,10000,(ncolonne+1))) #prepare data.frame for recapture data

colnames(data.marked) = c("id",colonne.names) ## add an id column

data.marked$initial.age = rep(0,nrow(data.marked))  # initial age (1 in June, 1.25 in September)

data.marked$Coh = rep(0,nrow(data.marked)) # Cohort (number)

data.marked$Coh_n = rep(0,nrow(data.marked)) # Cohort (character/factor)

data.marked$Coh.pflood = rep(0,nrow(data.marked)) 

data.marked$heter = 0



unique.mark.data = with(data_df,unique(Mark_cor)) # unique tagged fish


for (i in 1:length(unique.mark.data)) {
  
  data.prov = subset(data_df,Mark_cor == unique.mark.data[i])
  
  year.prov = data.prov$Year
  
  monthyear = paste(data.prov$Month,data.prov$Year,sep="")
  
  incl.year = colonne.names %in% monthyear
  
  data.marked[i,2:(ncolonne+1)] = ifelse(incl.year==F,0,1)
  
  
  if (annual == 0) {
    
    if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
      ### check whether the first capture was in June or September (use <= 7 to be sure)
      
      data.marked[i,"initial.age"] = min(data.prov$Age) *12 ## we are considering months
      
    } else {data.marked[i,"initial.age"] = (min(data.prov$Age) *12) + 3 }} else {
      
      if (min(subset(data.prov, Year == min(data.prov$Year))$Month) <= 7) { 
        ### check whether the first capture was in June or September (use <= 7 to be sure, some were in July)
        
        data.marked[i,"initial.age"] = min(data.prov$Age) *1 ## we are considering months
        
      } else {data.marked[i,"initial.age"] = (min(data.prov$Age) *1) + 0.25 }
      
    }  ### initial.age is a reserved column name
  
  data.marked[i,"Coh"] = data.prov$Cohort[1] ## Cohort of the fish as number
  
  data.marked[i,"Coh_n"] = as.character(data.prov$Cohort[1]) ## Cohort of the fish as character/factor 
  
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

data.marked$Coh_n <- factor(data.marked$Coh_n, levels = Cohort.levels)
data.marked$Coh.pflood <- factor(data.marked$Coh.pflood, levels = Coh.pflood.levels)



if (annual == 0) {  ## if I compute monthly, time between sampling is 9 months, 3 months
  spacebwtime = rep(c(9,3),8)} else {spacebwtime = rep(c(0.75,0.25),8)}

# Add Season covariate
season.rep = c(rep(c(0,1),8),0) # 0 is for Sept to June, 1 from June to Sept

Season=matrix(rep(season.rep,nrow(data.marked)),ncol=length(season.rep),byrow=T)

col.spacebwtime = c(1,cumsum(spacebwtime)+1)

colnames(Season)=paste("Season",col.spacebwtime,sep="")

data.marked=cbind(data.marked,Season)


#data.proc=process.data(data.marked,model="CJS",groups=c("Coh","Season"),
#time.intervals = spacebwtime)
flood.rep = rep(0,17)

flood.rep[c(3,7)] = 1

Flood=matrix(rep(flood.rep,nrow(data.marked)),ncol=length(flood.rep),byrow=T)

col.spacebwtime = c(1,cumsum(spacebwtime)+1)

colnames(Flood)=paste("Flood",col.spacebwtime,sep="")

data.marked=cbind(data.marked,Flood)


#####

design.Phi=list(static=c("Coh_n","Coh.pflood","heter"),time.varying=c("Season","Flood"))
design.p=list(static=c("Coh_n","Coh.pflood","heter"),time.varying=c("Season","Flood"))
design.parameters=list(Phi=design.Phi,p=design.p)
data.proc=process.data(data.marked, model = "CJS", time.interval = spacebwtime)
ddl=make.design.data(data.proc,parameters=design.parameters)


fit.models=function()
{
  Phi.Coh=list(formula=~Coh_n)
  Phi.Coh.pflood=list(formula= ~ Coh.pflood)
  Phi.Coh.pflood.heter.add=list(formula= ~ Coh.pflood + heter + Flood)
  Phi.Coh.pflood.heter.mult=list(formula= ~ Coh.pflood * heter + Flood)
  Phi.Flood=list(formula=~Flood)
  Phi.Flood.heter.add=list(formula=~Flood + heter)
  Phi.Flood.heter.mult=list(formula=~Flood * heter)
 
  Phi.Coh.pflood.Flood = list(formula = ~ Coh.pflood + Flood)
  Phi.dot=list(formula = ~1)
  Phi.Season = list(formula = ~Season)
 
   p.dot=list(formula=~1)
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=data.proc, ddl=ddl,
                      external=FALSE,accumulate=FALSE, hessian = T)
  return(results)
}
test.heter=fit.models()
ddl.heter = ddl
test.heter.mod = predict(test.heter[[8]], ddl = ddl.heter, se = T)

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