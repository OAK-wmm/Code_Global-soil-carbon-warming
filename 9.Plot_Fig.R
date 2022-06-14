###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

library(cowplot)
library(ggrepel)
library(rcartocolor)
library(gridExtra)
library(ggplot2)
library(metafor)
library(reshape)
library(RColorBrewer)

rm(list=ls())

##Fig.2
#=================================================================================
#=================================================================================
dat_cont = read.csv("XXXX/df_ef_rep10_content.csv")
dat_stock = read.csv("XXXX/Fig data/df_ef_rep9_stock.csv")
# t1 = ttheme_default(core=list(bg_params = list(fill = NULL,col = "grey", lwd = 0.5), 
#                               fg_params = list(col = c("#DD1C1A", "#20A4F3", "#F0C808")),
#                               base_size = 6))

t1 = ttheme_default(core=list(fg_params = list(col = c("#DD1C1A", "#20A4F3", "#248232")),
                              bg_params = list(fill = "#D9D9D9", alpha = 0.4)), 
                    base_size = 9)
dat_cont$layers = factor(dat_cont$layers, 
                         levels =c("0-30 cm","30-100 cm", "100-200 cm"), 
                         labels = c("0-0.3 m", "0.3-1 m", "1-2 m"))

layer = c("0-0.3 m", "0.3-1 m", "1-2 m")
stack_label = as.data.frame(array(NA, dim = c(3,4)))
for (i in 1:3) {
  data = dat_cont[which(dat_cont$layers %in% layer[i]),]
  modlm = lm(mu~temp_ch, data = data)
  suu = summary(modlm)
  stack_label[i,1] = layer[i]
  stack_label[i,2] = round(suu$r.squared,2)
  stack_label[i,3] = round(suu$coefficients[2,1],2)
  stack_label[i,4] = round(suu$coefficients[2,4],2)
}
colnames(stack_label) = c("layers","R2","Slope","p")
# stack_label$'ypois' =  c(-28, -25, -22)
t(stack_label)


dat_cont$pois.text = rep(c(-31, -33, -35), each=5)
dat_cont$"pois.point" = c(rep(c(c(1:5)+0.05), 1), rep(c(1:5), 1), rep(c(1:5)-0.05, 1))

P_RR_Content = ggplot(dat_cont, aes(x=pois.point, y = mu))+
  geom_errorbar(aes(ymin = ci.lw, ymax = ci.up),width = 0, position = position_dodge2(width = 0.2, padding = -0.2))+
  geom_smooth(aes(color = layers), method = "lm", se=F)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", size=1)+
  geom_point(aes(color = layers), size = 3, position = position_dodge2(width = 0.2, padding = -0.2))+
  geom_text(aes(label = (num), x= temp_ch, y=pois.text, color = layers), show.legend = F)+
  # geom_text(data = stack_label, aes(label = paste("R2 = ",round(r2,2)), x=1.3, y=ypois, color =layers))+
  # geom_text(data = stack_label, aes(label = paste("Slope = ",round(slope,2)), x=2.3, y=ypois, color =layers))+
  # geom_text(data = stack_label, aes(label = paste("p_value = ",round(pvalue,2)), x=3.3, y=ypois, color =layers))+
  ylab("Percentage change of SOC content (%)")+
  xlab("Warming (°„C)")+
  scale_color_manual(name = c(""),values = c("#DD1C1A", "#20A4F3", "#248232"))+
  scale_y_continuous(limits = c(-35, 5))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        legend.direction = c("horizontal"), 
        legend.position = c(0.5,0.95), 
        legend.background = element_blank())+
  annotation_custom(grob = tableGrob(stack_label[,-1],rows = NULL,theme = t1),
                    xmin = 1, xmax = 3,
                    ymin = -24, ymax = -20)
P_RR_Content


dat_stock$pois.point = c(rep(c(c(1:5)+0.05), 1), rep(c(1:5), 1), rep(c(1:5)-0.05, 1))
dat_stock$layers = factor(dat_stock$layers, levels =c("0-30 cm","30-100 cm", "100-200 cm"),
                          labels = c("0-0.3 m", "0.3-1 m", "1-2 m"))
dat_stock$pois.text = rep(c(-31, -33, -35), each=5)

stack_label_stack = as.data.frame(array(NA, dim = c(3,4)))
for (i in 1:3) {
  data = dat_stock[which(dat_stock$layers %in% layer[i]),]
  modlm = lm(mu~temp_ch, data = data)
  suu = summary(modlm)
  stack_label_stack[i,1] = layer[i]
  stack_label_stack[i,2] = round(suu$r.squared,2)
  stack_label_stack[i,3] = round(suu$coefficients[2,1],2)
  stack_label_stack[i,4] = round(suu$coefficients[2,4],2)
}
colnames(stack_label_stack) = c("layers","R2","Slope","p")
# stack_label_stack$'ypois' = c(-28, -25, -22)
t(stack_label_stack)

#### Fig.2
P_RR_Stock = ggplot(dat_stock, aes(x=pois.point, y = mu))+
  geom_errorbar(aes(ymin = ci.lw, ymax = ci.up),width = 0, position = position_dodge2(width = 0.5, padding = 0.5))+
  geom_smooth(aes(color = layers), method = "lm", se=F)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", size=1)+
  geom_point(aes(color = layers), size = 3)+
  geom_text(aes(label = (num), x= temp_ch, y=pois.text, color = layers), show.legend = F)+
  # geom_text(data = stack_label_stack, aes(label = paste("R2 = ",round(r2,2)), x=1.3, y=ypois, color =layers))+
  # geom_text(data = stack_label_stack, aes(label = paste("Slope = ",round(slope,2)), x=2.3, y=ypois, color =layers))+
  # geom_text(data = stack_label_stack, aes(label = paste("p_value = ",round(pvalue,2)), x=3.3, y=ypois, color =layers))+
  ylab("Percentage change of SOC stock (%)")+
  xlab("Warming (°„C)")+
  scale_color_manual(name = c(""),values = c("#DD1C1A", "#20A4F3", "#248232"))+
  scale_y_continuous(limits = c(-35, 5))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        legend.direction = c("horizontal"), 
        legend.position = c(0.5,0.95), 
        legend.background = element_blank())+
  annotation_custom(grob = tableGrob(stack_label_stack[,-1],rows = NULL,theme = t1),
                    xmin = 1, xmax = 3,
                    ymin = -24, ymax = -20)

p = cowplot::plot_grid(P_RR_Stock, P_RR_Content, labels = c("a","b"))
p

ggsave(p, filename = "XXXX/Fig1.Percentage_temp.pdf",  width=8.97,height=4.1, unit="in",dpi=300)
ggsave(p, filename = "XXXX/Fig1.Percentage_temp.png",  width=8.97,height=4.1, unit="in",dpi=300)


##==============================================================================
#######################################################
##Fig.4
library(dplyr)
library(cowplot)
Imp_cont = read.csv("E:/Work/World_Soil_Temp/Fig data/new_Fig data/importance.content.csv")
Imp_cont = arrange(Imp_cont, Importance)
Imp_cont$long_Vars = factor((Imp_cont$long_Vars), 
                            levels = (Imp_cont$long_Vars))
Imp_cont$Type = factor(Imp_cont$Type)

P_Importance_cont =
  ggplot(Imp_cont, aes(y=long_Vars,x=Importance))+
  geom_bar(stat = "identity", aes(fill = Type),width = 0.8, show.legend = F)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))+
  geom_text(aes(label = paste("R2 = 0.64"), x=0.6, y=13))+
  geom_text(aes(label = paste("RMSE = 0.59"), x=0.6, y=14))+
  xlab("Relative importance")+
  ylab("")+
  scale_fill_manual(values = c("#5BC0EB","#754F44","#F28F3B","#248232"))

dat.sampsize20.content = read.csv("E:/Work/World_Soil_Temp/Fig data/new_Fig data/dat.sampsize20.content.csv")
id = c(2:15, 18)

dat.sampsize20.content$SoilOrder = factor(dat.sampsize20.content$SoilOrder)
dat.sampsize20.content$Biomes = factor(dat.sampsize20.content$Biomes)
dat.sampsize20.content$Seasonality.pre = factor(dat.sampsize20.content$Seasonality.pre)
dat.sampsize20.content$Soil.layers = factor(dat.sampsize20.content$Soil.layers)

var.names = names(dat.sampsize20.content)
sign = as.data.frame(array(NA, dim = c(length(id),4)))
rr = (exp(dat.sampsize20.content$yi) - 1) *100
for(ii in 1:length(id)){
  corr = cor.test(dat.sampsize20.content[,id[ii]], rr)
  sign[ii,1] = var.names[id[ii]]
  if(corr$p.value < 0.05 & corr$estimate <0){
    sign[ii,2] = "-"
  }else if(corr$p.value < 0.05 & corr$estimate >0){
    sign[ii,2] = "+"
  }else{
    sign[ii,2] = ""
  }
  sign[ii,3] = corr$estimate
  sign[ii,4] = corr$p.value
}

df.pred.val.content = read.csv("XXXX/Fig data/new_Fig data/df.pred.val.content.csv")
df.pred.cal.content = read.csv("XXXX/Fig data/new_Fig data/df.pred.cal.content.csv")
df.pred.val.cal.con = rbind(df.pred.cal.content, df.pred.val.content)
df.pred.val.cal.con$"type" = rep(c("Calibration","Validation"), times = c(nrow(df.pred.cal.content),nrow(df.pred.val.content)))
df.pred.val.cal.con$type = factor(df.pred.val.cal.con$type, 
                                  levels = c("Calibration","Validation"))

p.val.cal=
  ggplot(df.pred.val.cal.con, aes(obs, pred))+
  geom_point(aes(col = type),shape = 1, size=3, alpha = 0.6)+
  geom_smooth(aes(col = type), se = F, method = "lm")+
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey")+
  scale_color_manual(name = "",values = c("#FCCA46", "#8093F1"))+
  scale_y_continuous(limits = c(-5, 5))+
  scale_x_continuous(limits = c(-5, 5))+
  theme_bw(base_size = 12)+
  ylab(c("Predicted effect size"))+
  xlab(c("Observed effect size"))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"), 
        legend.position = c(0.35, 0.9),
        legend.background = element_blank())

p.imp.val.cal.con=
  ggdraw()+
  draw_plot(P_Importance_cont, 0,0,1,1)+
  draw_plot(p.val.cal, 0.52,0.11, 0.45,0.6)



Imp_stock = read.csv("XXXX/new_Fig data/importance.stock.csv")
Imp_stock = arrange(Imp_stock, Importance)

Imp_stock$long_Vars = factor(Imp_stock$long_Vars, 
                             levels = (Imp_stock$long_Vars))
P_Importance_stock = 
  ggplot(Imp_stock, aes(y=long_Vars,x=Importance))+
  geom_bar(stat = "identity", aes(fill = Type), width = 0.8, show.legend = F)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))+
  geom_text(aes(label = paste("R2 = 0.69"), x=0.6, y=13))+
  geom_text(aes(label = paste("RMSE = 0.51"), x=0.6, y=14))+
  scale_fill_manual(values = c("#5BC0EB","#754F44","#F28F3B","#248232"))+
  xlab("Relative importance")+
  ylab("")

df.pred.val.stock = read.csv("XXXX/Fig data/new_Fig data/df.pred.val.stock.csv")
df.pred.cal.stock = read.csv("XXXX/Fig data/new_Fig data/df.pred.cal.stock.csv")
df.pred.val.cal.sto = rbind(df.pred.cal.stock, df.pred.val.stock)
df.pred.val.cal.sto$"type" = rep(c("Calibration","Validation"), times = c( nrow(df.pred.cal.stock), nrow(df.pred.val.stock)))
df.pred.val.cal.sto$type = factor(df.pred.val.cal.sto$type,levels = c("Calibration","Validation"))


p.val.cal.sto=
  ggplot(df.pred.val.cal.sto, aes(obs, pred))+
  geom_point(aes(col = type),shape = 1, size=3, alpha=0.6)+
  geom_smooth(aes(col = type), se = F, method = "lm")+
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey")+
  scale_color_manual(name = "",values = c("#FCCA46", "#8093F1"))+
  scale_y_continuous(limits = c(-4.5, 4))+
  scale_x_continuous(limits = c(-4.5, 4))+
  theme_bw(base_size = 12)+
  ylab(c("Predicted effct size"))+
  xlab(c("Observed effect size"))+
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"), 
        legend.position = c(0.35, 0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank())

p.imp.val.cal.sto=
  ggdraw()+
  draw_plot(P_Importance_stock, 0,0,1,1)+
  draw_plot(p.val.cal.sto, 0.52,0.11, 0.45,0.6)

P_Importance = cowplot::plot_grid(p.imp.val.cal.sto, p.imp.val.cal.con, labels = c("a","b"))
P_Importance
ggsave(P_Importance, filename = "XXXX/Fig/Fig4.P_Importance.pdf",  width=8.26,height=4.4, unit="in",dpi=300) 
ggsave(P_Importance, filename = "XXXX/Fig/Fig4.P_Importance.png",  width=8.26,height=4.4, unit="in",dpi=300) 

#===============================================================================
##Figure: Extended Data Fig. 7. 

Exp = read.csv("XXXX/Fig data/Meat data artical/Meta_toal_final.csv")
# Exp = Exp[which(Exp$T.only == T),]
Exp = Exp[,c(7:13, 15)]
names(Exp)
colnames(Exp) = c("Duration","Temp_ch","soc_con","soc_tre","nrep",
                      "sd_con","sd_tre", "Ecotype")
Exp$sd_con[which(is.na(Exp$sd_con) == T)] = Exp$soc_con[which(is.na(Exp$sd_con) == T)] * 0.2
Exp$sd_tre[which(is.na(Exp$sd_tre) == T)] = Exp$soc_tre[which(is.na(Exp$sd_tre) == T)] * 0.2

# mean(Exp$sd_con[which(!is.na(Exp$sd_con) == T)]/Exp$soc_con[which(!is.na(Exp$sd_con) == T)])
# mean(Exp$sd_tre[which(!is.na(Exp$sd_tre) == T)]/Exp$soc_tre[which(!is.na(Exp$sd_tre) == T)])


exp.yivi = escalc(data = Exp, measure = "ROM",
                  m1i = soc_tre, sd1i = sd_tre, n1i =nrep,
                  m2i = soc_con, sd2i = sd_con, n2i = nrep)

exp.yivi$"Temp_level" = NA
exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch <= 1)] = 1
exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 1 & exp.yivi$Temp_ch <= 2)] = 2
exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 2 & exp.yivi$Temp_ch <= 3)] = 3
exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 3 & exp.yivi$Temp_ch <= 4)] = 4
# exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 4)] = 5

exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 4 & exp.yivi$Temp_ch <= 5)] = 5
exp.yivi$"Temp_level"[which(exp.yivi$Temp_ch > 5)] = 6
exp.yivi$Temp_level = factor(exp.yivi$Temp_level)

mod_temp = rma(yi, vi, mods = ~Temp_ch, data = exp.yivi, method = "REML")
mod_temp_level = rma(yi, vi, mods = ~Temp_level-1, data = exp.yivi, method = "REML")
mod_Eco = rma(yi, vi, mods = ~Eco-1, data = exp.yivi, method = "REML")

temp = seq(0, 6, 0.1)
temp.level = factor(seq(1,6,1))
pred = data.frame(temp,predict(mod_temp, temp))
pred[,c(2,4,5)] = (exp(pred[,c(2,4,5)])-1)*100


dat_stock = read.csv("XXXX/Fig data/df_ef_rep9_stock.csv")
dat_stock.L1 = dat_stock[which(dat_stock$layers == "0-30 cm"),]
model_l1 = readRDS("XXXX/Fig data/mod.meta.Temp.L1.4diff_2_4_degree.Rdata")
pred.toal = data.frame(temp_ch=seq(0, 5, 0.1),predict(model_l1, seq(0, 5, 0.1)))
pred.toal[,c(2,4,5)] = (exp(pred.toal[,c(2,4,5)])-1)*100

  
RR = (exp(mod_temp_level$beta)-1)*100
p_value = mod_temp_level$pval
ci.lb = (exp(mod_temp_level$ci.lb)-1)*100
ci.ub = (exp(mod_temp_level$ci.ub)-1)*100
num = table(exp.yivi$Temp_level)
dat_exp = data.frame(RR, p_value, ci.lb, ci.ub, num)
dat_exp$"layers1" = rep("0-30 cm", times = 6)
dat_exp$"temp_ch" = 1:6

stack_label_exp_our = as.data.frame(array(NA, dim = c(2,4)))
stack_label_exp_our[1,1] = "Our approach"
stack_label_exp_our[1,2] = round(model_l1$QM,2)
stack_label_exp_our[1,3] = round((exp(model_l1$beta[2])-1)*100,2)
stack_label_exp_our[1,4] = round(model_l1$pval[2],2)


stack_label_exp_our[2,1] = "Expermrnt"
stack_label_exp_our[2,2] = round(mod_temp$QM[1],2)
stack_label_exp_our[2,3] = round((exp(mod_temp$beta[2])-1)*100,2)
stack_label_exp_our[2,4] = round(mod_temp$pval[2],2)

colnames(stack_label_exp_our) = c("Approach","QM","Slope","p")


####Extended Data Fig. 7.
p_exp_ourapproach = 
  ggplot()+
  geom_hline(yintercept = 0, linetype = 2, size=0.5, color = "#847577")+
   
  ##our study
  geom_errorbar(data  =dat_stock.L1, aes(x=temp_ch, ymin = ci.lw, ymax = ci.up), width=0, color = "#D62828")+
  geom_point(data = dat_stock.L1, aes(x=temp_ch, y=mu, color="a This study"),size=3.5, show.legend = T)+
  
  geom_smooth(data = dat_stock.L1, aes(x=temp_ch, y=mu), method = "lm", se=F, color = "#D62828")+
  geom_text(data = dat_stock.L1, aes(label = paste0(num), x=temp_ch+0.2,  y=mu+3), color = "#D62828")+
  
  ##meta _analysis
  geom_point(data = Metaother, aes(x=Temp, y=lnRR_SOC, shape = Author), size=3)+
  geom_errorbar(data = Metaother, aes(x=Temp, ymin=ci.lb, ymax=ci.ub), width=0, size=0.3)+
  
  geom_errorbar(data=dat_exp, aes(x=temp_ch, ymin = ci.lb, ymax = ci.ub), width=0,color = "#00A8E8")+
  geom_line(data = pred, aes(x=temp, y=pred), color = "#00A8E8", size=1)+
  geom_point(data = dat_exp, aes(x=temp_ch, y=RR), color="#00A8E8",size=3.5, show.legend = T)+
  geom_text(data = dat_exp, aes(label = paste0(Freq), x=temp_ch+0.2,  y=RR-3), color = "#00A8E8")+
  
  xlab("Warming (°„C)")+
  ylab("Percentage change of soil carbon stock (%)")+
  
  scale_y_continuous(expand = c(0.01,0.01))+
  scale_x_continuous(expand = c(0,0))+
  scale_shape_manual(name = "Published meta-analysis results", values = c(0,2,4,5,6))+
  scale_fill_manual(name = "" ,values = "#636363",label = c("Q10 =1.75"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        plot.margin = margin(1,0.15,0.15,0.15, "cm"),
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(vjust = 0.1))+
  annotation_custom(grob = tableGrob(stack_label_exp_our[,-1],rows = NULL,theme = t2),
                    xmin = 3, xmax = 6,
                    ymin = 9, ymax = 12)

p_exp_ourapproach


################################################################
###exp and our approach between rebiomes
dat_cont_exp = read.csv("XXXX/Fig data/Exp_result.csv")
dat_cont_exp$"Re_percet" = (exp(dat_cont_exp$lnRR.1)-1) *100

df.exp_biomes = aggregate(dat_cont_exp$Re_percet, list(dat_cont_exp$Ecosystem,dat_cont_exp$Magnitude), function(x){
  ttest = t.test(x)
  int = as.numeric(ttest$conf.int)
  mean = as.numeric(ttest$estimate)
  num = length(x)
  return(c(mean, int, num))
})


df.exp_biomes = data.frame(df.exp_biomes$Group.1, df.exp_biomes$Group.2, df.exp_biomes$x)
names(df.exp_biomes)  =c("Biomes", "temp_ch", "mean","ci.lw", "ci.up", "Freq")
df.exp_biomes = df.exp_biomes[-which(df.exp_biomes$Biomes == "Wetland"),]
head(df.exp_biomes)

library(stringr)
filenames = list.files("XXXX/Meta_result_0intercpt_biomess_reclass2/")
id = c(1:15)
sol = list()
for (jj in id) {
  dat = readRDS(paste0("XXXX/Meta_result_0intercpt_biomess_reclass2/", filenames[jj]))
  if(dim(dat[[2]]) > dim(dat[[1]])[1]){
    sol0 = cbind(dat[[1]], dat[[2]][-dim(dat[[2]])])
  }else{
    sol0 = cbind(dat[[1]], dat[[2]])
    
  }
  # sol0 = arrange(sol0, effect_value)
  sol[[jj]] = sol0
}
result_meta_rebiomes = do.call(rbind, sol)
# write.csv(result_meta_rebiomes, "Fig data/result_meta_rebiomes.csv")

biomes_RR = read.csv("XXXX/Fig data/result_meta_rebiomes.csv")

biomes_RR_our = biomes_RR[which(biomes_RR$layers == "L1" & !biomes_RR$Var1 == "Other"),][,c(4,5,7,8,11,12)]
biomes_RR_our$ci.lb = (exp(biomes_RR_our$ci.lb)-1)*100
biomes_RR_our$ci.ub = (exp(biomes_RR_our$ci.ub)-1)*100

biomes_RR_our$temp_ch = factor(biomes_RR_our$temp_ch, labels = paste0("+",1:5, "°„C"))
df.exp_biomes$temp_ch = factor(df.exp_biomes$temp_ch, labels = paste0("+",1:5, "°„C"))
biomes_RR_our = data.frame(biomes_RR_our$Var1, biomes_RR_our[,-5])
names(biomes_RR_our) = names(df.exp_biomes)

df_exp_our = rbind(biomes_RR_our, df.exp_biomes)
df_exp_our$"Type" = rep(c("Our Approach", "Experment"), times = c(nrow(biomes_RR_our), nrow(df.exp_biomes)))
df_exp_our$Biomes = factor(df_exp_our$Biomes, labels = c("Forests", "Grasslands", "Shrublands", "Tundras"))

p_exp_our_ecosystem = 
  ggplot()+
  geom_vline(xintercept = 0, linetype = 2, color = "grey")+
  geom_errorbar(data = df_exp_our, aes(y=Biomes , x= mean, xmin = ci.lw, xmax = ci.up, color =Type),
                position = position_dodge(0.4),width = 0, size = 0.3, show.legend = F)+
  geom_point(data= df_exp_our, aes(y=Biomes, x=mean, color =Type),position = position_dodge(0.4), size=3.5, show.legend = F)+
  geom_text(data = df_exp_our, aes(label = (Freq), x=30, y=Biomes, color =Type),
            position = position_dodge(0.4), show.legend = F)+
  facet_wrap(~temp_ch, nrow=1)+
  scale_color_manual(name = "", values = c("#00A8E8", "#D62828"))+
  scale_x_continuous(limits = c(-40, 40))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(color = "black", size = 12), 
        axis.text.y = element_text(angle = 45),
        axis.title = element_text(vjust = 0.))+
  ylab("Ecosystem types")+
  xlab("Percentage change of SOC content (%)")

####Extended Data Fig. 7. 
p.exp_ourapproach_ecosystem=
  ggdraw()+
  draw_plot(p_exp_ourapproach,   0.065,0.5, 0.945, 0.48)+
  draw_plot(p_exp_our_ecosystem, 0,  0, 1, 0.50)+
  draw_label("a", 0.02,0.94)+
  draw_label("b", 0.02,0.39)+
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))

p.exp_ourapproach_ecosystem
ggsave(p.exp_ourapproach_ecosystem, filename = "XXXX/Fig3.Field study and this study_new_final.pdf",  width=8,height=4, unit="in",dpi=300)


#############==============================================================
## Extended Data Fig. 8. 
Exp = read.csv("XXXX/Fig data/Meat data artical/Meta_toal_final.csv")
# Exp = Exp[which(Exp$T.only == T),]
Exp = Exp[,c(7:13, 15)]
names(Exp)
colnames(Exp) = c("Duration","Temp_ch","soc_con","soc_tre","nrep",
                  "sd_con","sd_tre", "Ecotype")
Exp$sd_con[which(is.na(Exp$sd_con) == T)] = Exp$soc_con[which(is.na(Exp$sd_con) == T)] * 0.2
Exp$sd_tre[which(is.na(Exp$sd_tre) == T)] = Exp$soc_tre[which(is.na(Exp$sd_tre) == T)] * 0.2

# mean(Exp$sd_con[which(!is.na(Exp$sd_con) == T)]/Exp$soc_con[which(!is.na(Exp$sd_con) == T)])
# mean(Exp$sd_tre[which(!is.na(Exp$sd_tre) == T)]/Exp$soc_tre[which(!is.na(Exp$sd_tre) == T)])


exp.yivi = escalc(data = Exp, measure = "ROM",
                  m1i = soc_tre, sd1i = sd_tre, n1i =nrep,
                  m2i = soc_con, sd2i = sd_con, n2i = nrep)

mod_Duration = rma(yi, vi, mods = ~Duration, data = exp.yivi, method = "REML")

stack_label_exp_Duration = as.data.frame(array(NA, dim = c(1,4)))
stack_label_exp_Duration[1,1] = "Our approach"
stack_label_exp_Duration[1,2] = round(mod_Duration$QM,2)
stack_label_exp_Duration[1,3] = round((exp(mod_Duration$beta[2])-1)*100,2)
stack_label_exp_Duration[1,4] = round(mod_Duration$pval[2],2)
colnames(stack_label_exp_Duration) = c("Approach","QM","Slope","p")

t3 = ttheme_default(core=list(fg_params = list(col = c("blue")),
                              bg_params = list(fill = "#D8E2DC")), 
                    base_size = 9)

wi<- 1/sqrt(exp.yivi$vi)
size  <-1+ 50.0 * (wi - min(wi,na.rm = T))/(max(wi,na.rm = T) - min(wi,na.rm = T))
exp.yivi$"size" = size
dur = seq(0.4, 25, 0.1)
Pred = data.frame(predict(mod_Duration, dur), dur)
summary(lm(dur~Pred$pred))

###Figure: the change of RR with duration of exp
p_exp_Expermentduration = 
  ggplot(exp.yivi, aes(Duration, (exp(yi)-1)*100))+
  # geom_smooth(method = "lm", se = F)+
  geom_hline(yintercept = 0, linetype = 2, size=1, color = "grey")+
  geom_point(aes(size=size), shape=1, color = "blue", show.legend = F)+
  geom_line(data = Pred, aes(x=dur, (exp(pred)-1)*100),color = "black", size=1)+
  geom_ribbon(data = Pred, aes(x=dur, pred, ymin = (exp(ci.lb)-1)*100, ymax=(exp(ci.ub)-1)*100),
              fill = "blue", alpha=0.4, size=1)+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        plot.margin = margin(1,0.15,0.1,0.1, "cm"))+
  xlab("Experment duration (yrear)")+
  ylab("Percentage change of SOC (%)")+
  annotation_custom(grob = tableGrob(stack_label_exp_Duration[,-1],rows = NULL,theme = t3),
                    xmin = 20, xmax = 22,
                    ymin = 60, ymax = 65)
p_exp_Expermentduration

ggsave(p_exp_our_ecosystem, filename = "XXXX/Fig.S3.Experment and our approach in biomes.pdf",  width=9.5,height=4, unit="in",dpi=300) 
ggsave(p_exp_our_ecosystem, filename = "XXXX/Fig.S3.Experment and our approach in biome.png",  width=9.5,height=4, unit="in",dpi=300) 









