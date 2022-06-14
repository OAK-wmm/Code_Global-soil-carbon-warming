###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

#### Sensitivity analysis of randomly select class, sample size, croplands, 
#### soilorder, precipitation categorization

library(metafor)
library(Matrix)
library(ggplot2)
library(gtools)

#### ===========================================================================
#### Meta-analysis function
fun_meta = function(dat, tau2 = 0.01,change = 1, prec = .000001){
  
  if(length(which(is.na(dat$vi) == T)) == 0){
    dat = dat
  }else{
    dat = dat[-which(is.na(dat$vi) == T), ]
  }
  
  while (change > prec) {
    
    save   <- tau2
    wi     <- 1/(dat$vi + tau2)
    mu     <- sum(wi*dat$yi)/sum(wi)
    adj    <- (sum(wi^2 * (dat$yi - mu)^2) + sum(wi^2)/sum(wi) - sum(wi)) / (sum(wi^2) - 2*sum(wi^3)/sum(wi) + (sum(wi^2)/sum(wi))^2)
    while (tau2 + adj < 0) # use step-halving if necessary
      adj <- adj / 2
    tau2   <- tau2 + adj
    change <- abs(save - tau2)
    
  }
  
  wii = 1/(dat$vi + tau2)
  se = sqrt(1/sum(wii))
  
  CI = c((exp(mu-1.96*se)-1)*100, (exp(mu+1.96*se)-1)*100)
  
  num = dim(dat)[1]
  mu = (exp(mu)-1)*100
  return(c(mu, tau2, CI, num))
}


### effect value of sample size > 20
##==============================================================================
filenames = mixedsort(list.files("dat_Mean_sd_data/"))

id = c(41:45)
effect_value = array(NA, dim = c(3, 5, 5))
for (ii in 1:length(id)) {
  
  data.meta = readRDS(paste0("dat_Mean_sd_data/", filenames[id[ii]]))
  id.na = which(data.meta$sd.soc.100_200cm.con<0 | data.meta$sd.soc.100_200cm.tre<0)
  if(length(id.na)  == 0){}else{
    data.meta = data.meta[-id.na, ]
  }
  
  data.meta = data.meta[which(data.meta$Num.pro.tre >=20 & data.meta$Num.pro.con >=20),]
  
  df.yivi.L1 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.0_30cm.tre, sd1i = sd.soc.0_30cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.0_30cm.con, sd2i = sd.soc.0_30cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L1$yi) > 100 | abs(df.yivi.L1$vi) > 100 |
                   df.yivi.L1$vi ==0 | is.na(df.yivi.L1$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L1 = df.yivi.L1
  }else{
    df.yivi.L1 = df.yivi.L1[-id.del,]
  }
  effect_value[1, ,ii] = fun_meta(df.yivi.L1)  
  
  df.yivi.L2 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.30_100cm.tre, sd1i = sd.soc.30_100cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.30_100cm.con, sd2i = sd.soc.30_100cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L2$yi) > 100 | abs(df.yivi.L2$vi) > 100 |
                   df.yivi.L2$vi ==0 | is.na(df.yivi.L2$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L2 = df.yivi.L2
  }else{
    df.yivi.L2 = df.yivi.L2[-id.del,]
  }
  effect_value[2, ,ii] = fun_meta(df.yivi.L2)  
  
  df.yivi.L3 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.100_200cm.tre, sd1i = sd.soc.100_200cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.100_200cm.con, sd2i = sd.soc.100_200cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L3$yi) > 100 | abs(df.yivi.L3$vi) > 100 |
                   df.yivi.L3$vi ==0 | is.na(df.yivi.L3$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L3 = df.yivi.L3
  }else{
    df.yivi.L3 = df.yivi.L3[-id.del,]
  }
  effect_value[3, ,ii] = fun_meta(df.yivi.L3)
}

df_ef = rbind(effect_value[,,1],effect_value[,,2],effect_value[,,3],effect_value[,,4], effect_value[,,5])
colnames(df_ef) = c("RR", "tau2", "int_low", "int_up","number")
layers = rep(c("0-30 cm", "30-100 cm", "100-200 cm"), times = 5)
layers = factor(layers, levels = c("0-30 cm", "30-100 cm", "100-200 cm"))
temp_ch = rep(c(1:5), each = 3)
df_ef = data.frame(cbind(layers, temp_ch, df_ef))
df_ef$layers = factor(df_ef$layers, labels =  c("0-30 cm", "30-100 cm", "100-200 cm"))
# df_ef$temp_ch = factor(df_ef$temp_ch)
df_ef$"pso.num" = df_ef$temp_ch
df_ef$"pso.num.y" = rep(c(-50, -48, -46),times=5)
# write.csv(df_ef, "Fig data/df_ef_size20_stock.csv")

####============================================================================
### effect value of delete corp land
filenames = mixedsort(list.files("dat_Mean_sd_data/"))

id = c(41:45)
effect_value = array(NA, dim = c(3, 5, 5))
for (ii in 1:length(id)) {
  
  data.meta = readRDS(paste0("dat_Mean_sd_data/", filenames[id[ii]]))
  id.na = which(data.meta$sd.soc.100_200cm.con<0 | data.meta$sd.soc.100_200cm.tre<0)
  if(length(id.na)  == 0){}else{
    data.meta = data.meta[-id.na, ]
  }
  data.meta = data.meta[-which(data.meta$Biomes == "9"),]
  
  # data.meta = data.meta[which(data.meta$Num.pro.tre >=20 & data.meta$Num.pro.con >=20),]
  
  df.yivi.L1 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.0_30cm.tre, sd1i = sd.soc.0_30cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.0_30cm.con, sd2i = sd.soc.0_30cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L1$yi) > 100 | abs(df.yivi.L1$vi) > 100 |
                   df.yivi.L1$vi ==0 | is.na(df.yivi.L1$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L1 = df.yivi.L1
  }else{
    df.yivi.L1 = df.yivi.L1[-id.del,]
  }
  effect_value[1, ,ii] = fun_meta(df.yivi.L1)  
  
  df.yivi.L2 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.30_100cm.tre, sd1i = sd.soc.30_100cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.30_100cm.con, sd2i = sd.soc.30_100cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L2$yi) > 100 | abs(df.yivi.L2$vi) > 100 |
                   df.yivi.L2$vi ==0 | is.na(df.yivi.L2$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L2 = df.yivi.L2
  }else{
    df.yivi.L2 = df.yivi.L2[-id.del,]
  }
  effect_value[2, ,ii] = fun_meta(df.yivi.L2)  
  
  df.yivi.L3 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.100_200cm.tre, sd1i = sd.soc.100_200cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.100_200cm.con, sd2i = sd.soc.100_200cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L3$yi) > 100 | abs(df.yivi.L3$vi) > 100 |
                   df.yivi.L3$vi ==0 | is.na(df.yivi.L3$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L3 = df.yivi.L3
  }else{
    df.yivi.L3 = df.yivi.L3[-id.del,]
  }
  effect_value[3, ,ii] = fun_meta(df.yivi.L3)
}
df_ef = rbind(effect_value[,,1],effect_value[,,2],effect_value[,,3],effect_value[,,4], effect_value[,,5])
colnames(df_ef) = c("RR", "tau2", "int_low", "int_up","number")
layers = rep(c("0-30 cm", "30-100 cm", "100-200 cm"), times = 5)
layers = factor(layers, levels = c("0-30 cm", "30-100 cm", "100-200 cm"))
temp_ch = rep(c(1:5), each = 3)
df_ef = data.frame(cbind(layers, temp_ch, df_ef))
df_ef$layers = factor(df_ef$layers, labels =  c("0-30 cm", "30-100 cm", "100-200 cm"))
# df_ef$temp_ch = factor(df_ef$temp_ch)
df_ef$"pso.num" = df_ef$temp_ch
df_ef$"pso.num.y" = rep(c(-38, -36, -34),times=5)
# write.csv(df_ef, "Fig data/df_ef_del_corp_stock.csv")

#### ===========================================================================
#### effect of 10 randomly selecting class with 50 mm precipitation categorization
filenames = mixedsort(list.files("dat_yivi/"))

eff_val_10rep <- as.data.frame(array(NA, dim = c(150, 8)))
for (jj in 1:150) {
  dat = readRDS(paste0("dat_yivi/",filenames[jj]))
  
  eff_val_10rep[jj, 4:8] = fun_meta(dat)
}
layers = rep(paste0("L", c(1:3)), each  =  50)
rep_times = rep(c(1:10), each = 5, times = 3)
temp_ch = rep(c(1:5), times = 30)
eff_val_10rep[,1] = layers
eff_val_10rep[,2] = rep_times
eff_val_10rep[,3] = temp_ch
eff_val_10rep = as.numeric(data.frame(eff_val_10rep))

colnames(eff_val_10rep) = c("layers","rep_times", "temp_ch", "mu", "tau2", "ci.lw", "ci.up", "num")
eff_val_10rep$'pois.text' = rep(c(-34, -31, -28), each = 50)
eff_val_10rep$"pois.point" = c(rep(c(c(1:5)+0.1), 10), rep(c(1:5), 10), rep(c(1:5)-0.1, 10))
eff_val_10rep$"x.pois" = rep(seq(1.3, 4.8,0.8), times = 10)
head(eff_val_10rep)
eff_val_10rep$layers  =factor(eff_val_10rep$layers, 
                              levels = c("L1", "L2", "L3"),
                              labels = c("0-30 cm", "30-100 cm","100-200 cm"))

P_10rep = ggplot(eff_val_10rep, aes(x=pois.point, y = mu))+
  geom_errorbar(aes(ymin = ci.lw, ymax = ci.up),width = 0, position = position_dodge2(width = 0.5, padding = 0.5))+
  geom_smooth(aes(color = layers), method = "lm", se=F)+
  geom_hline(yintercept = 0, linetype = 2, color = "grey", size=1)+
  geom_point(aes(color = layers), size = 3)+
  facet_wrap(~rep_times)+
  geom_text(aes(label = num, x=x.pois, y=pois.text, color = layers), size=3)+
  ylab("Percentage change (%)")+
  xlab("Warming (¡ãC)")+
  scale_color_manual(name = c(""),values = c("#DD1C1A", "#20A4F3", "#F0C808"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.8, 0.2))
P_10rep
ggsave(P_10rep, filename = "Fig/Fig_new/P_10rep_Stock.pdf",  width=8.27,height=7.2, unit="in",dpi=300) 

######Figure:  cropland, sample size>20
dat_cont = read.csv("XXXX/df_ef_rep10_content.csv")
dat_stock = read.csv("XXXX/df_ef_rep9_stock.csv")
dat_cont_size20 = read.csv("XXXX/df_ef_size20_content.csv")
dat_stock_size20 = read.csv("XXXX/df_ef_size20_stock.csv")
dat_cont_delcorp = read.csv("XXXX/df_ef_del_corp_content.csv")
dat_stock_delcorp = read.csv("XXXX/df_ef_del_corp_stock.csv")

dat_cont = dat_cont[,!names(dat_cont)%in%c("rep_times", "x.pois")]
dat_stock = dat_stock[,!names(dat_stock)%in%c("rep_times", "x.pois")]
names(dat_cont)=names(dat_stock)=names(dat_cont_suborder) = names(dat_stock_suborder) = names(dat_stock_size20)

dat_sen = rbind(dat_stock,dat_stock_delcorp, dat_stock_size20,   
                dat_cont, dat_cont_delcorp, dat_cont_size20)

dat_sen$"Ctype" = factor(rep(c("Stock","Content"), each = 45), levels = c("Stock","Content"))
dat_sen$"sentype" = factor(rep(c("All data","non-corplans","Sample size > 20"), each = 15, times = 2),
                           levels = c("All data","non-corplans","Sample size > 20"))
layer = c("0-0.3 m","0.3-1 m", "1-2 m")
dat_sen$layers = factor(dat_sen$layers, 
                        levels =c("0-30 cm","30-100 cm", "100-200 cm"),
                        labels = c("0-0.3 m","0.3-1 m", "1-2 m"))
dat_sen$"pois.point" = c(rep(c(1:5)-0.3, times = 3), rep(c(1:5)-0.15, each = 3), rep(c(1:5),each = 3))
dat_sen$pois.text = rep(c(-47,-50, -53), each=15, times = 2)


#### Extended Data Fig. 2
P_cropland_samplesize<-
  ggplot(dat_sen, aes(pois.point, RR))+
  geom_point(aes(color = sentype), size = 3)+
  geom_smooth(aes(color = sentype), method = "lm", se = F)+
  geom_errorbar(aes(ymin = int_low, ymax = int_up,color = sentype),width = 0, position = position_dodge2(width = 0.5, padding = 0.5))+
  geom_text(aes(label = (number), x= temp_ch, y=pois.text, color = sentype), show.legend = F, size=4)+
  # scale_y_reverse()+
  scale_color_manual(name = "", values = c("#E71D36", "#568259", "#3E92CC", "#F7CB15"))+
  scale_y_continuous(limits = c(-55, 6), expand = c(0,0))+
  ylab("Percentage change of SOC (%)")+
  xlab("Warming (¡ãC)")+
  facet_grid(layers~Ctype)+
  theme_bw(base_size = 12)+
  theme(legend.position = c(0.88,0.95), 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(color = "black", size = 12))

ggsave("XXXX/Fig.S3.Percentage_temp_samplesize20_suborder.pdf",P_cropland_samplesize,  width = 8.26, height=9, 
       units = c("in"))

####============================================================================
#### effect of soil orders on results


#### Extended Data Fig. 5
dat_contUSDA = read.csv("XXXX/df_ef_rep10_content.csv")
dat_stockUSDA = read.csv("XXXX/Fig data/df_ef_rep9_stock.csv")

dat_eff_soiloFAO = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_FAO_stock.csv")
dat_eff_soiloWRB  = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_WRB_stock.csv")

dat_eff_soiloFAO_conttent = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_FAO_content.csv")
dat_eff_soiloWRB_content  = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_WRB_content.csv")

dat_cont_suborder = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_content.csv")
dat_stock_suborder = read.csv("XXXX/Fig data/new_Fig data/eff_suborder_stock.csv")

dat = rbind(dat_stockUSDA,dat_eff_soiloFAO, dat_eff_soiloWRB, dat_stock_suborder,
            dat_contUSDA, dat_eff_soiloFAO_conttent, dat_eff_soiloWRB_content, dat_cont_suborder)
dat$"Ctype" = factor(rep(c("Stock","Content"), each = 60), levels = c("Stock","Content"))
dat$"Soilclass" = factor(rep(c("USDA","FAO","WRB", "Soil suborder (USDA)"), each = 15, times = 2), 
                         levels = c("USDA","FAO","WRB", "Soil suborder (USDA)"))
layer = c("0-0.3 m","0.3-1 m", "1-2 m")
dat$layers = factor(dat$layers, 
                    levels =c("0-30 cm","30-100 cm", "100-200 cm"),
                    labels = c("0-0.3 m","0.3-1 m", "1-2 m"))
dat$pois.point = c(rep(c(1:5)+0.195, 3), rep(c(1:5)+0.065,3), rep(c(1:5)-0.065, 3), rep(c(1:5)-0.195, 3))
dat$pois.text = rep(c(-31,-34,-37, -40), each=15, times = 2)

p_soilorder<- ggplot(dat, aes(pois.point, mu))+
  geom_point(aes(color = Soilclass), size = 3)+
  geom_smooth(aes(color = Soilclass), method = "lm", se = F)+
  geom_errorbar(aes(ymin = ci.lw, ymax = ci.up,color = Soilclass),width = 0, position = position_dodge2(width = 0.5, padding = 0.5))+
  geom_text(aes(label = (num), x= temp_ch, y=pois.text, color = Soilclass), show.legend = F)+
  # scale_y_reverse()+
  scale_color_manual(name = "Soil classfications",values = c("#D80032", "#4A7C59", "#35A7FF", "#523F38"))+
  scale_y_continuous(limits = c(-42, 6), expand = c(0,0))+
  ylab("Percentage change of SOC (%)")+
  xlab("Warming (¡ãC)")+
  facet_grid(layers~Ctype)+
  theme_bw(base_size = 12)+
  theme(legend.position = c(0.38,0.93), 
        legend.background = element_blank(),
        axis.text = element_text(color = "black", size = 12))

ggsave("XXXX/Fig.S7.new-FAO_WRB_soilorder.pdf",p_soilorder,  width = 8.26, height=9, 
       units = c("in"))


####============================================================================
#### effect of 10 randomly selecting class with 25 mm precipitation categorization

filenames = list.files("Mean_sd_25P/")

for (ii in 1:length(filenames)) {
  
  data.meta = readRDS(paste0("Mean_sd_25P/", filenames[ii]))
  
  df.yivi.L1 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.0_30cm.tre, sd1i = sd.soc.0_30cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.0_30cm.con, sd2i = sd.soc.0_30cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L1$yi) > 100 | abs(df.yivi.L1$vi) > 100 |
                   df.yivi.L1$vi ==0 | is.na(df.yivi.L1$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L1 = df.yivi.L1
  }else{
    df.yivi.L1 = df.yivi.L1[-id.del,]
  }
  
  saveRDS(df.yivi.L1, paste0("dat_yivi_25p/df.yivi.L1_", gsub("Mean_sd_Stock_diff_", "", filenames[ii])))
  
  df.yivi.L2 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.30_100cm.tre, sd1i = sd.soc.30_100cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.30_100cm.con, sd2i = sd.soc.30_100cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L2$yi) > 100 | abs(df.yivi.L2$vi) > 100 |
                   abs(df.yivi.L2$vi)==0 | is.na(df.yivi.L2$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L2 = df.yivi.L2
  }else{
    df.yivi.L2 = df.yivi.L2[-id.del,]
  }
  saveRDS(df.yivi.L2, paste0("dat_yivi_25p/df.yivi.L2_",gsub("Mean_sd_Stock_diff_", "", filenames[ii])))
  
  df.yivi.L3 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.100_200cm.tre, sd1i = sd.soc.100_200cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.100_200cm.con, sd2i = sd.soc.100_200cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L3$yi) > 100 | abs(df.yivi.L3$vi) > 100 |
                   abs(df.yivi.L3$vi)==0 | is.na(df.yivi.L3$vi) == T)
  
  if(length(id.del) == 0){
    df.yivi.L3 = df.yivi.L3
  }else{
    df.yivi.L3 = df.yivi.L3[-id.del,]
  }
  saveRDS(df.yivi.L3, paste0("dat_yivi_25p/df.yivi.L3_", gsub("Mean_sd_Stock_diff_", "", filenames[ii])))
  
}

#####make sensitivity test of precitation gap
eff_val_10rep$"MAPgap" = "¡À50 mm"
eff_val_10rep_25P$"MAPgap" = "¡À25 mm"
eff_val = rbind(eff_val_10rep, eff_val_10rep_25P)

eff_Ptest_Content = eff_val[which(eff_val$rep_times == 1),]
eff_Ptest_Content$"SOC.type" = "Content"
# write.csv(eff_Ptest_Content, "Fig data/eff_Ptest_Content.csv")

#####plot 
eff_Ptest_Content = read.csv("Fig data/eff_Ptest_Content.csv")
eff_Ptest_Stock = read.csv("XXXX/eff_Ptest_stock.csv")
eff_Ptest = rbind(eff_Ptest_Content, eff_Ptest_Stock)
eff_Ptest$x.pois[which(eff_Ptest$MAPgap == "¡À25 mm")] = eff_Ptest$temp_ch[which(eff_Ptest$MAPgap == "¡À25 mm")]-0.1
eff_Ptest$x.pois[which(eff_Ptest$MAPgap == "¡À50 mm")] = eff_Ptest$temp_ch[which(eff_Ptest$MAPgap == "¡À50 mm")]+0.1

eff_Ptest$pois.text[which(eff_Ptest$MAPgap == "¡À25 mm")] = -32
eff_Ptest$pois.text[which(eff_Ptest$MAPgap == "¡À50 mm")] = -36
eff_Ptest$SOC.type = factor(eff_Ptest$SOC.type, levels = c("Stock", "Content"))

eff_Ptest$layers  =factor(eff_Ptest$layers, 
                          levels = c("0-30 cm", "30-100 cm", "100-200 cm"),
                          labels = c("0-0.3 m", "0.3-1 m","1-2 m"))

#### Figure: Extended Data Fig. 11
P_eff_Ptest<-
  ggplot(eff_Ptest, aes(x=x.pois, y=mu))+
  geom_point(aes(color = MAPgap), size=3)+
  geom_smooth(aes(color = MAPgap), method = "lm", se = F)+
  geom_errorbar(aes(ymin = ci.lw, ymax = ci.up, color = MAPgap),width = 0)+
  geom_text(aes(label = (num), x= temp_ch, y=pois.text, color = MAPgap), show.legend = F)+
  facet_grid(layers~SOC.type)+
  ylab("Percentage change of SOC stock (%)")+
  xlab("Warming (¡ãC)")+
  scale_color_manual(name = c(""),values = c("#DD1C1A", "#20A4F3"), label = c("25 mm categorization","50 mm categorization"))+
  scale_y_continuous(limits = c(-40, 8), expand = c(0,0))+
  theme_bw(base_size = 12)+
  theme(#panel.grid = element_blank(), 
    legend.direction = c("vertical"), 
    legend.position = c(0.85,0.95), 
    legend.background = element_blank(),
    axis.text = element_text(colour = "black"))
P_eff_Ptest

ggsave("XXXX/P_eff_Ptest.pdf",P_eff_Ptest,  width = 6.75, height=7.72, 
       units = c("in"))

