library(ggplot2)
library(grid)
library(dplyr)
library(cowplot)


### PLOT

size.title = 20
line.lwd = 1.2
size.label.x = 18
size.text.x = 14
size.point = 6
size.label.y = 18
size.text.y = 14
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

## Theme to be used for all plots

theme.heter =  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
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
                     legend.text = element_text( size = size.legend.text),legend.position = c(0.7, 0.2)
) 
  


## ZAK

zak.heter.best = predict(zak.mod.heter[[4]], ddl = zak.mod.ddl.heter , se = T)

zak.gg <- ggplot(filter(zak.heter.best$Phi, Flood == 0), aes(y = estimate, x = heter, group = Coh.pflood, lty = Coh.pflood)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("Zak") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter
 

# zak.gg

## TREBU

trebu.heter.best = predict(trebu.mod.heter[[6]], ddl = trebu.mod.ddl.heter, se = T)

trebu.gg <- ggplot(filter(trebu.heter.best$Phi), aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2,alpha = 0.5) + 
  ggtitle("Trebu") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter

# trebu.gg

### ZADLA

zadla.heter.best = predict(zadla.heter.mod[[9]], ddl = zadla.ddl.mod.heter, se = T)

zadla.gg <- ggplot(filter(zadla.heter.best$Phi, Flood == 0), aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("Zadla") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter

# zadla.gg

## LIDRI

loidri.heter.best = predict(loidri.mod.heter[[8]],ddl = loidri.mod.ddl.heter,se = T)

loidri.gg <- ggplot(filter(loidri.heter.best$Phi), aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("LIdri") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter

# loidri.gg


### COWPLOT

Plot_heter_all = plot_grid(zak.gg, trebu.gg,zadla.gg, loidri.gg,
                                          labels = c("A", "B","C","D"), 
                                          nrow = 2,ncol = 2, align = "h",hjust = -2.5)

save_plot("Plot_heter_all.pdf", Plot_heter_all,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)


### PLOT of heterozygosity-survival in case of floods (Supplementary Material)


zak.gg.flood <- ggplot(filter(zak.heter.best$Phi, Flood == 1), aes(y = estimate, x = heter)) + 
  #geom_point(size = size.point,lwd = 4) +
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("Zak") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,1,0.1)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter 


zak.gg.flood


zadla.gg.flood <- ggplot(filter(zadla.heter.best$Phi, Flood == 1), aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("Zadla") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,1,0.1)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter

# zadla.gg.flood


lipo.heter.best = predict(lipo.heter.mod[[11]], ddl = lipo.ddl.mod.heter, se = T)


lipo.gg.flood <- ggplot(filter(lipo.heter.best$Phi, Flood == 1), aes(y = estimate, x = heter)) + 
  geom_line(lwd  = line.lwd) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl),width = 0.01, col = "gray80", lty = 2, alpha = 0.5) + 
  ggtitle("Lipo") +
  scale_x_continuous(limits = c(0.1,0.7), breaks = seq(0.1,0.7,0.1)) +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,1,0.1)) +
  labs(y = bquote(phi)) +
  labs(x = label.T) +
  theme.heter

# lipo.gg.flood

### COWPLOT

Plot_heter_flood = plot_grid(zak.gg.flood, zadla.gg.flood, lipo.gg.flood,
                           labels = c("A", "B","C"), 
                           nrow = 3,ncol = 1, align = "h",hjust = -2.5)

save_plot("Plot_heter_flood.pdf", Plot_heter_flood,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)
