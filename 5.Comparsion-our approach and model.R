###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

library(metafor)
library(Matrix)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(ComplexHeatmap)
library("devtools")
library(cowplot)
library(ggrepel)
library(rcartocolor)
library(gridExtra)


rm(list = ls())

temp = seq(0, 6, 0.1)
temp.level = factor(seq(1,6,1))

dat_stock = read.csv("XXXX/Fig data/df_ef_rep9_stock.csv")
dat_stock$layers = factor(dat_stock$layers, levels =c("0-30 cm","30-100 cm", "100-200 cm"),
                          labels = c("0-0.3 m", "0.3-1 m", "1-2 m"))
layer = unique(dat_stock$layers)
tt=seq(0,6,0.1)
rateII0 = c(-0.12,0,0.12)

Q10.fixed = c(1, 1.1, 1.2, 1.4, 1.6, 2.0)
RR = as.data.frame(array(NA, dim = c(length(tt), length(Q10.fixed)+1)))
RR[,1] = tt
for (ii in 1:length(Q10.fixed)) {
  RR[,ii+1]=((0+1)*exp(-0.1*log(Q10.fixed[ii])*tt)-1)*100
}
names(RR) = c("Temp",paste0("Q10 = ", Q10.fixed))
df.RR = melt(RR, id = "Temp")
df.RR$variable = factor(df.RR$variable)

fun_optim = function(x,tt, rateII, dat, SL){
  mu = dat$mu[which(dat$layers == SL)]
  rmse = sqrt(sum(((((1+rateII)*exp(-0.1*log(x)*tt)-1)*100)[which(tt %in% c(1:5))] - mu)**2)/5)
}


Q10II = data.frame(array(NA, dim = c(length(tt),5)))
names(Q10II) = c("Temp","rateII0","Q10.optim", "RR.optim", "Q10.mod")
Q10II[,1] = tt
df.mod.Q10 = list()
table = as.data.frame(array(NA, dim = c(3, 4)))
colnames(table) = c("C_in (%)", paste0("Q10_",layer))
table[,1] = c("-2% °C-1", "0% °C-1", "+2% °C-1")
write.csv(table, "XXXX/new_Fig data/Q10_table.csv")

for (ii in 1:length(layer)) {
  
  Q10.list = list()
  
  for (jj in 1:length(rateII0)) {
    
    if(rateII0[jj] == -0.12){
      rateII = -seq(0,abs(rateII0[jj]), 0.12/(length(tt)-1))
    }else if(rateII0[jj] == 0) {
      rateII = rep(0, length(tt))
    }else{
      rateII = seq(0,abs(rateII0[jj]), 0.12/(length(tt)-1))
    }
    
    Q10II[,2] = rateII0[jj]
    Q10 <- optimize(fun_optim, c(0, 3), tol = 0.0001, tt= tt, rateII = rateII, dat = dat_stock, SL = layer[ii])
    table[jj,ii+1] = round(Q10$minimum, 2)
    print( round(Q10$minimum, 2))
    
    
    Q10II[,3] = paste0("Q10 = ",round(Q10$minimum, 2))
    Q10II[,5] = round(Q10$minimum, 2)
    
    Q10II[,4] = ((rateII+1)*exp(-0.1*log(Q10$minimum)*tt)-1)*100
    
    Q10.list[[jj]] = Q10II
    
  }
  df.Q10II.optim = do.call(rbind.data.frame, Q10.list)
  df.Q10II.optim$Q10.opt = factor(df.Q10II.optim$Q10.optim)
  df.Q10II.optim$layer = layer[ii]
  
  df.mod.Q10[[ii]] = df.Q10II.optim
}
df.mod.Q10 = do.call(rbind.data.frame, df.mod.Q10)
df.mod.Q10$rateII0 = factor(df.mod.Q10$rateII0)


cols = c("#DD1C1A", "#20A4F3", "#248232",
         colorRampPalette((brewer.pal(9, 'Greys'))[-c(1,2,8,9)])(length(Q10.fixed)))
colgrey  = colorRampPalette((brewer.pal(9, 'Greys'))[-c(1,2,8,9)])(length(Q10.fixed))

t1 = ttheme_default(core=list(fg_params = list(col = c("#DD1C1A", "#20A4F3", "#248232")),
                              bg_params = list(fill = "#D9D9D9", alpha = 0.4)), 
                    base_size = 9)

p_stock=
  ggplot()+
  geom_line(data=df.RR, aes(Temp, value, color = variable), size=0.3)+
  geom_point(data=dat_stock, aes(pois.point, mu, color = layers), size=3)+
  geom_errorbar(data=dat_stock,aes(pois.point, mu,ymin = ci.lw, ymax = ci.up,color = layers),
                width = 0, position = position_dodge2(width = 0.2, padding = -0.2), size=0.3, 
                show.legend = F)+
  geom_smooth(data=dat_stock, aes(pois.point, mu, color = layers), method = "lm", se=F, 
              show.legend = F)+
  scale_color_manual(name = "",values = cols)+
  scale_linetype_manual(name = "Carbon input change (%)", values = c(4,2,3),
                        label = c("-2% °C-1", "0% °C-1", "+2% °C-1"))+
  annotate("text",x=5.6, y=-0,  label = expression("Q"[10]*' = 1.0'),  size=3, color = colgrey[1])+
  annotate("text",x=5.6, y=-5,  label = expression("Q"[10]*' = 1.1'),  size=3, color = colgrey[2])+
  annotate("text",x=5.6, y=-10, label = expression("Q"[10]*' = 1.2'),  size=3, color = colgrey[3])+
  annotate("text",x=5.6, y=-18, label = expression("Q"[10]*' = 1.4'),  size=3, color = colgrey[4])+
  annotate("text",x=5.6, y=-23, label = expression("Q"[10]*' = 1.6'),  size=3, color = colgrey[5])+
  annotate("text",x=5.6, y=-33, label = expression("Q"[10]*' = 2.0'),  size=3, color = colgrey[6])+
  # annotate("text",x=5.6, y=-40.5,   label = expression("Q"[10]*' = 2.5'), size=3, color = colgrey[7])+
  scale_y_continuous(expand = c(0,0),limits = c(-35, 5))+
  scale_x_continuous(expand = c(0,0))+
  geom_line(data=df.mod.Q10, aes(Temp, RR.optim, color = layer, linetype = rateII0))+
  ylab("Percentage change of SOC stock (%)")+
  xlab("Warming (°C)")+
  # annotation_custom(grob = tableGrob((table),rows = NULL,theme = t1),
  #                   xmin = 1, xmax = 3,
  #                   ymin = -30, ymax = -25)+
  
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        plot.margin = margin(1,0.15,0.15,0.15, "cm"),
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(vjust = 0.1))

#################################################################################
####SOC content
 
temp = seq(0, 6, 0.1)
temp.level = factor(seq(1,6,1))

dat_content = read.csv("XXXX/Fig data/df_ef_rep10_content.csv")
dat_content$layers = factor(dat_content$layers, levels =c("0-30 cm","30-100 cm", "100-200 cm"),
                          labels = c("0-0.3 m", "0.3-1 m", "1-2 m"))
dat_content.L1 = dat_content[which(dat_content$layers == "0-30 cm"),]

layer = unique(dat_content$layers)
tt=seq(0,6,0.1)
rateII0 = c(-0.12,0,0.12)

Q10II[,1] = tt

RR = as.data.frame(array(NA, dim = c(length(tt), length(Q10.fixed)+1)))
RR[,1] = tt
for (ii in 1:length(Q10.fixed)) {
  RR[,ii+1]=((0+1)*exp(-0.1*log(Q10.fixed[ii])*tt)-1)*100
}
names(RR) = c("Temp",paste0("Q10 = ", Q10.fixed))
df.RR = melt(RR, id = "Temp")
df.RR$variable = factor(df.RR$variable)

fun_optim = function(x,tt, rateII, dat, SL){
  mu = dat$mu[which(dat$layers == SL)]
  rmse = sqrt(sum(((((1+rateII)*exp(-0.1*log(x)*tt)-1)*100)[which(tt %in% c(1:5))] - mu)**2)/5)
}


Q10II = data.frame(array(NA, dim = c(length(tt),5)))
names(Q10II) = c("Temp","rateII0","Q10.optim", "RR.optim", "Q10.mod")
Q10II[,1] = tt
df.mod.Q10 = list()
table = as.data.frame(array(NA, dim = c(3, 4)))
colnames(table) = c("C_in (%)", paste0("Q10_",layer))
table[,1] = c("-2% °C-1", "0% °C-1", "+2% °C-1")

for (ii in 1:length(layer)) {
  
  Q10.list = list()
  
  for (jj in 1:length(rateII0)) {
    
    if(rateII0[jj] == -0.12){
      rateII = -seq(0,abs(rateII0[jj]), 0.12/(length(tt)-1))
    }else if(rateII0[jj] == 0) {
      rateII = rep(0, length(tt))
    }else{
      rateII = seq(0,abs(rateII0[jj]), 0.12/(length(tt)-1))
    }
    
    Q10II[,2] = rateII0[jj]
    Q10 <- optimize(fun_optim, c(0, 3), tol = 0.0001, tt= tt, rateII = rateII, dat = dat_content, SL = layer[ii])
    table[jj,ii+1] = round(Q10$minimum, 2)
    print( round(Q10$minimum, 2))
    
    
    Q10II[,3] = paste0("Q10 = ",round(Q10$minimum, 2))
    Q10II[,5] = round(Q10$minimum, 2)
    
    Q10II[,4] = ((rateII+1)*exp(-0.1*log(Q10$minimum)*tt)-1)*100
    
    Q10.list[[jj]] = Q10II
    
  }
  df.Q10II.optim = do.call(rbind.data.frame, Q10.list)
  df.Q10II.optim$Q10.opt = factor(df.Q10II.optim$Q10.optim)
  df.Q10II.optim$layer = layer[ii]
  
  df.mod.Q10[[ii]] = df.Q10II.optim
}
write.csv(table, "XXXX/Q10_table_content.csv")

df.mod.Q10 = do.call(rbind.data.frame, df.mod.Q10)
df.mod.Q10$rateII0 = factor(df.mod.Q10$rateII0)


cols = c("#DD1C1A", "#20A4F3", "#248232",
         colorRampPalette((brewer.pal(9, 'Greys'))[-c(1,2,8,9)])(length(Q10.fixed)))
colgrey  = colorRampPalette((brewer.pal(9, 'Greys'))[-c(1,2,8,9)])(length(Q10.fixed))

t1 = ttheme_default(core=list(fg_params = list(col = c("#DD1C1A", "#20A4F3", "#248232")),
                              bg_params = list(fill = "#D9D9D9", alpha = 0.4)), 
                    base_size = 9)

#### Fig.3 
p_content<-
ggplot()+
  geom_line(data=df.RR, aes(Temp, value, color = variable), size=0.3)+
  geom_point(data=dat_content, aes(pois.point, mu, color = layers), size=3)+
  geom_errorbar(data=dat_content,aes(pois.point, mu,ymin = ci.lw, ymax = ci.up,color = layers),
                width = 0, position = position_dodge2(width = 0.2, padding = -0.2), size=0.3, 
                show.legend = F)+
  geom_smooth(data=dat_content, aes(pois.point, mu, color = layers), method = "lm", se=F, 
              show.legend = F)+
  scale_color_manual(name = "",values = cols)+
  scale_linetype_manual(name = "Carbon input change (%)", values = c(4,2,3),
                        label = c("-2% °C-1", "0% °C-1", "+2% °C-1"))+
  annotate("text",x=5.6, y=-0,  label = expression("Q"[10]*' = 1.0'),  size=3, color = colgrey[1])+
  annotate("text",x=5.6, y=-5,  label = expression("Q"[10]*' = 1.1'),  size=3, color = colgrey[2])+
  annotate("text",x=5.6, y=-10, label = expression("Q"[10]*' = 1.2'),  size=3, color = colgrey[3])+
  annotate("text",x=5.6, y=-18, label = expression("Q"[10]*' = 1.4'),  size=3, color = colgrey[4])+
  annotate("text",x=5.6, y=-23, label = expression("Q"[10]*' = 1.6'),  size=3, color = colgrey[5])+
  annotate("text",x=5.6, y=-33, label = expression("Q"[10]*' = 2.0'),  size=3, color = colgrey[6])+
  # annotate("text",x=5.6, y=-40.5,   label = expression("Q"[10]*' = 2.5'), size=3, color = colgrey[7])+
  scale_y_continuous(expand = c(0,0),limits = c(-35, 5))+
  scale_x_continuous(expand = c(0,0))+
  geom_line(data=df.mod.Q10, aes(Temp, RR.optim, color = layer, linetype = rateII0))+
  ylab("Percentage change of SOC content (%)")+
  xlab("Warming (°C)")+
  # annotation_custom(grob = tableGrob((table),rows = NULL,theme = t1),
  #                   xmin = 1, xmax = 3,
  #                   ymin = -30, ymax = -25)+
  
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        plot.margin = margin(1,0.15,0.15,0.15, "cm"),
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(vjust = 0.1))

plot_grid(p_stock, p_content, labels = c("a","b"))

