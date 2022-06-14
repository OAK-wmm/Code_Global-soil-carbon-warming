###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

library(foreach)
library(doParallel)
library(emdbook)
library(doSNOW)
library(maps)
library(stringr)
library(gtools)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

#### Load soil profile SOC stock data
setwd("XXXX")
Stock_SOC = readRDS("Stock_SOC.Rdata")


### ============================================================================
####Selecting “ambient�? and “warm�? class
Sel.profilr = function(del_temp){
  
  Stock_SOC.con = Stock_SOC_class
  Stock_SOC.tre = Stock_SOC_class
  
  ii = 1
  # SOC_tre = as.data.frame(array(NA, dim = c(dim(Stock_SOC)[1], 20)))
  # SOC_con = as.data.frame(array(NA, dim = c(dim(Stock_SOC)[1], 20)))
  
  SOC_tre = list()
  SOC_con = list()
  
  while (dim(Stock_SOC.con)[1] > 0) {
    
    nrow.sample = sample(1:dim(Stock_SOC.con)[1], 1)
    Data0 = Stock_SOC.con[nrow.sample, ]
    
    index.con = which(Stock_SOC.con$iTemp >= (Data0$iTemp - sei) & Stock_SOC.con$iTemp <= (Data0$iTemp + sei) & 
                        Stock_SOC.con$Soilorder_wrb ==  Data0$Soilorder_wrb & 
                        Stock_SOC.con$iPre >= (Data0$iPre - del_pre) & Stock_SOC.con$iPre <= (Data0$iPre + del_pre) &
                        Stock_SOC.con$Seasonality.pre == Data0$Seasonality.pre &
                        Stock_SOC.con$Landform == Data0$Landform)
    
    index.tre = which(Stock_SOC.tre$iTemp >= (Data0$iTemp + del_temp - sei) &
                        Stock_SOC.tre$iTemp <= (Data0$iTemp + del_temp + sei) & 
                        Stock_SOC.tre$Soilorder_wrb ==  Data0$Soilorder_wrb & 
                        Stock_SOC.tre$iPre >= (Data0$iPre - del_pre) & Stock_SOC.tre$iPre <= (Data0$iPre + del_pre) &
                        Stock_SOC.tre$Seasonality.pre == Data0$Seasonality.pre &
                        Stock_SOC.tre$Landform == Data0$Landform)
    
    if(length(index.con) == 0){
      Data.con = NA
      Stock_SOC.con = Stock_SOC.con[-nrow.sample,]
    }else{
      
      Data.con = Stock_SOC.con[index.con,]
      
      Data.con$"Num.pro" = length(index.con)
      Data.con$"iProfile.base" = Data0$iProfile
      Data.con$"Lon.base" = Data0$longitude
      Data.con$"Lat.base" = Data0$latitude
      colnames(Data.con) <- paste0(colnames(Data.con), ".con")
      
      Stock_SOC.con = Stock_SOC.con[-which(Stock_SOC.con$iProfile %in% Data.con$iProfile.con), ]
      # dim_con = dim(Data.con)
      # id.na= which.min(is.na(SOC_con[,1])==T)
      # SOC_con[id.na:(id.na+dim_con[1]-1),] = Data.con
    }
    SOC_con[[ii]] = Data.con
    
    if(length(index.tre) == 0){
      Data.tre = NA
      Stock_SOC.tre = Stock_SOC.tre
    }else{
      
      Data.tre = Stock_SOC.tre[index.tre,]
      
      Data.tre$"Num.pro" = length(index.tre)
      Data.tre$"iProfile.base" = Data0$iProfile
      Data.tre$"Lon.base" = Data0$longitude
      Data.tre$"Lat.base" = Data0$latitude
      colnames(Data.tre) <- paste0(colnames(Data.tre), ".tre")
      
      Stock_SOC.tre = Stock_SOC.tre[-which(Stock_SOC.tre$iProfile %in% Data.tre$iProfile.tre), ]
      
      # dim_tre = dim(Data.tre)
      # id.na= which.min(is.na(SOC_tre[,1])==T)
      # SOC_tre[id.na:(id.na+dim_tre[1]-1),] = Data.tre
    }
    SOC_tre[[ii]] = Data.tre
    
    ii=ii+1
    # print(ii)
  }
  # colnames(SOC_tre) = paste0(c(names(Stock_SOC), "Num.pro","iProfile.base", "Lon.base","Lat.base"),"_tre")
  # colnames(SOC_con) = paste0(c(names(Stock_SOC), "Num.pro","iProfile.base", "Lon.base","Lat.base"),"_con")
  # 
  return(list(SOC_con,SOC_tre))
}

seq_temp = c(1,2,3,4,5)
del_pre = 20
sei = 0.5

cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)
foreach(xx = 1:10) %dopar% {
  for (jj in seq_temp) {
    saveRDS(assign(paste0("SOC_tre_con", jj), Sel.profilr(jj)), 
            paste0("Select_pro_25P/Stock_diff_",xx,"_", jj, "_degree.Rdata"))
  }
}
stopCluster(cl)


### ============================================================================
### Calculating the mean and sd of each “ambient “warm�? class
fun_mean.sd = function(jj){
  
  fun_amount.max = function(x){
    frent.table = data.frame(table(x))
    sds = frent.table[,1][which.max(frent.table$Freq)]
    return(as.character(sds))
  }
  
  SOC.2T <- readRDS(paste0("Select_pro_25P/", filenames[jj]))
  
  SOC.con = SOC.2T[[1]]
  SOC.tre = SOC.2T[[2]]
  
  df.SOC.con = do.call(rbind.data.frame, SOC.con)
  names(df.SOC.con)[4:6] = c("soc.0_30cm.con", "soc.30_100cm.con", "soc.100_200cm.con")
  df.SOC.con$"soc_whole.con" = rowSums(df.SOC.con[,4:6])
  id.sel = match(df.SOC.con$iProfile.con, Stock_SOC$iProfile )
  # data.frame(Stock_SOC$iProfile[id.sel],  df.SOC.con$iProfile.con)
  
  
  vars_del = c("Landuse","Litho")
  vars_factor = c("Landuse","Biomes", "SoilOrder","Landform","Seasonality.pre","suborder","sublandform", "xx", "Soilorder_fao", "Soilorder_wrb")
  vars_local = c("iProfile", "latitude", "longitude")
  vars_soc = c("soc.0-30 cm", "soc.30-100 cm", "soc.100-200 cm")
  
  
  env.factor = Stock_SOC[, names(Stock_SOC) %in% vars_factor, ][id.sel, ]
  env.factor.aggerate = aggregate(env.factor, list(df.SOC.con$iProfile.base.con), fun_amount.max)
  
  env.num = Stock_SOC[, !names(Stock_SOC) %in% c(vars_factor,vars_local, vars_soc), ][id.sel, ]
  env.num.aggerate= aggregate(env.num,list(df.SOC.con$iProfile.base.con), function(x)mean(x, na.rm = T))
  
  SOC.sel.con = df.SOC.con[,c(4:6,7,8, 14, 16)]
  df.SOC.mean.con = aggregate(SOC.sel.con,list(df.SOC.con$iProfile.base.con), function(x)mean(x, na.rm=T))
  df.SOC.sd.con = aggregate(SOC.sel.con,  list(df.SOC.con$iProfile.base.con), function(x)sd(x, na.rm=T))
  colnames(df.SOC.sd.con) = paste0("sd.",colnames(df.SOC.sd.con))
  
  # CV = colMeans(df.SOC.sd.con[,c(2:4,8)]/df.SOC.mean.con[,c(2:4,8)], na.rm = T)
  ID.1pro = which(df.SOC.mean.con$Num.pro.con == 1)
  df.SOC.sd.con[ID.1pro,c(2:4,8)] = t(t(df.SOC.mean.con[ID.1pro,c(2:4,8)]) * CV)
  
  
  #####tre
  df.SOC.tre = do.call(rbind.data.frame, SOC.tre)
  df.SOC.tre$"soc_whole.tre" = rowSums(df.SOC.tre[,4:6])
  
  names(df.SOC.tre)[4:6] = c("soc.0_30cm.tre", "soc.30_100cm.tre", "soc.100_200cm.tre")
  SOC.sel.tre =  df.SOC.tre[,c(4:6,7,8, 14,16)]
  df.SOC.mean.tre = aggregate(SOC.sel.tre,list(df.SOC.tre$iProfile.base.tre), function(x)mean(x, na.rm=T))
  df.SOC.sd.tre = aggregate(SOC.sel.tre,  list(df.SOC.tre$iProfile.base.tre), function(x)sd(x, na.rm=T))
  colnames(df.SOC.sd.tre) = paste0("sd.",colnames(df.SOC.sd.tre))
  
  # CV = colMeans(df.SOC.sd.tre[,c(2:4,8)]/df.SOC.mean.tre[,c(2:4,8)], na.rm = T)
  ID.1pro = which(df.SOC.mean.tre$Num.pro.tre == 1)
  df.SOC.sd.tre[ID.1pro,c(2:4,8)] = t(t(df.SOC.mean.tre[ID.1pro,c(2:4,8)]) * CV)
  
  
  id.base = match(df.SOC.mean.tre$Group.1, df.SOC.mean.con$Group.1 )
  iProfile.base = df.SOC.mean.con$Group.1[id.base]
  
  df.mean.sd = cbind(iProfile.base, 
                     df.SOC.mean.con[id.base,c(2:4,8)], 
                     df.SOC.sd.con[id.base,c(2:4,8)], 
                     Num.pro.con = df.SOC.mean.con$Num.pro.con[id.base],
                     df.SOC.mean.tre[,c(2:4,8)],
                     df.SOC.sd.tre[,c(2:4,8)],
                     Num.pro.tre = df.SOC.mean.tre$Num.pro.tre,
                     env.num.aggerate[id.base,], 
                     env.factor.aggerate[id.base,])
  
  df.mean.sd$'Pre.diff' = df.SOC.mean.tre$iPre.tre - df.SOC.mean.con[id.base, ]$iPre.con
  df.mean.sd$'Temp.diff' = df.SOC.mean.tre$iTemp.tre - df.SOC.mean.con[id.base, ]$iTemp.con
  
  saveRDS(df.mean.sd, paste0("Mean_sd_25P/Mean_sd_", gsub(".Rdata", "", filenames[jj]), ".Rdata"))
}

filenames = list.files("Select_pro_25P/")
cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl)

foreach(jj = 1:length(filenames)) %dopar% {
  fun_mean.sd(jj)
}
stopCluster(cl)


###=============================================================================
##### comparison between 10 repeats 
filenames = mixedsort(list.files("Mean_sd_10rep/"))

for (ii in 1:length(filenames)) {
  
  data.meta = readRDS(paste0("Mean_sd_10rep/", filenames[ii]))
  
  id.na = which(data.meta$sd.soc.100_200cm.con<0 | data.meta$sd.soc.100_200cm.tre<0)
  if(length(id.na)  == 0){}else{
    data.meta = data.meta[-id.na, ]
  }
  
  
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
  
  saveRDS(df.yivi.L1, paste0("yivi data/df.yivi.L1_", gsub("Mean_sd_SOC_", "", filenames[ii])))
  
  
  
  df.yivi.L2 = escalc(data =data.meta, measure = "ROM", 
                      m1i = soc.3_100cm.tre, sd1i = sd.soc.3_100cm.tre, n1i = Num.pro.tre, 
                      m2i = soc.3_100cm.con, sd2i = sd.soc.3_100cm.con,  n2i = Num.pro.con)
  id.del = which(abs(df.yivi.L2$yi) > 100 | abs(df.yivi.L2$vi) > 100 |
                   abs(df.yivi.L2$vi)==0 | is.na(df.yivi.L2$vi) == T)
  if(length(id.del) == 0){
    df.yivi.L2 = df.yivi.L2
  }else{
    df.yivi.L2 = df.yivi.L2[-id.del,]
  }
  saveRDS(df.yivi.L2, paste0("yivi data/df.yivi.L2_",gsub("Mean_sd_SOC_", "", filenames[ii])))
  
  
  
  
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
  saveRDS(df.yivi.L3, paste0("yivi data/df.yivi.L3_", gsub("Mean_sd_SOC_", "", filenames[ii])))
  
}

library(gtools)
filenames = mixedsort(list.files("yivi data/"))
eff_val_10rep <- as.data.frame(array(NA, dim = c(length(filenames), 8)))
for (jj in 1:length(filenames)) {
  dat = readRDS(paste0("yivi data/",filenames[jj]))
  
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
eff_val_10rep$'pois.text' = rep(c(-40, -37, -34), each = 50)
eff_val_10rep$"pois.point" = c(rep(c(c(1:5)+0.05), 10), rep(c(1:5), 10), rep(c(1:5)-0.05, 10))
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
  xlab("Warming (°C)")+
  scale_color_manual(name = c(""),values = c("#DD1C1A", "#20A4F3", "#523F38"))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        legend.position = c(0.8, 0.2))
ggsave(P_10rep, filename = "Fig/P_10rep_Content.pdf",  width=8.27,height=7.2, unit="in",dpi=300) 
ggsave(P_10rep, filename = "Fig/P_10rep_Content.png", width=8.27,height=7.2, unit="in",dpi=300) 
##==============================================================================

###=============================================================================
### effect value between biomes, soilorder, landform
######
filenames = list.files("dat_meat_0intrcp/")
id = 1:60


sol = list()
for (jj in id) {
  dat = readRDS(paste0("dat_meat_0intrcp/", filenames[jj]))
  if(dim(dat[[2]]) > dim(dat[[1]])[1]){
    sol0 = cbind(dat[[1]], dat[[2]][-dim(dat[[2]])])
  }else{
    sol0 = cbind(dat[[1]], dat[[2]])
    
  }
  # sol0 = arrange(sol0, effect_value)
  sol[[jj]] = sol0
}

filenames_qm = list.files("dat_meat_0intrcp/")
id = 1:60
filenames_qm[id]
sol_qm = list()
for (jj in id ) {
  dat = readRDS(paste0("dat_meat_0intrcp/", filenames_qm[jj]))
  if(dim(dat[[2]]) > dim(dat[[1]])[1]){
    sol0 = cbind(dat[[1]], dat[[2]][-dim(dat[[2]])])
  }else{
    sol0 = cbind(dat[[1]], dat[[2]])
    
  }
  # sol0 = arrange(sol0, effect_value)
  sol_qm[[jj]] = sol0
}
sol_qm = do.call(rbind.data.frame,sol_qm)
result_meta = do.call(rbind.data.frame, sol)
result_meta$temp_ch = as.numeric(substr(as.character(result_meta$temp_ch), str_length(as.character(result_meta$temp_ch)),str_length(as.character(result_meta$temp_ch))))
result_meta$QM = sol_qm$QM
result_meta$p_QM = sol_qm$p_QM

# result_meta = result_meta[-which(result_meta$Var1 %in% 	"BiomesNA"), -1]
result_meta$ci.lb = (exp(result_meta$ci.lb)-1)*100
result_meta$ci.ub = (exp(result_meta$ci.ub)-1)*100
result_meta$"sing" = NA
result_meta$sing[which(result_meta$ci.lb< 0 & result_meta$ci.ub > 0 )]="UNSIG"
result_meta$sing[which(result_meta$ci.lb< 0 & result_meta$ci.ub < 0 & result_meta$effect_value<0)]="SIG1"
result_meta$sing[which(result_meta$ci.lb> 0 & result_meta$ci.ub > 0 & result_meta$effect_value<0)]="SIG2"

result_meta$sing = factor(result_meta$sing)
result_meta$names_long = paste0(result_meta$Var1, result_meta$layers)
levelss = c(result_meta$names_long[which(result_meta$temp_ch == 1 & result_meta$layers == "L1")], 
            result_meta$names_long[which(result_meta$temp_ch == 1 & result_meta$layers == "L2")],
            result_meta$names_long[which(result_meta$temp_ch == 1 & result_meta$layers == "L3")])

result_meta$"Labels" = NA
result_meta$Labels[which(result_meta$Var1 == "Soilorder111")] = "Alfisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder112")] = "Andisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder113")] = "Aridisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder114")] = "Entisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder115")] = "Gelisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder116")] = "Histosols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder117")] = "Inceptisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder118")] = "Mollisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder119")] = "Oxisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder121")] = "Spodosols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder122")] = "Ultisols"
result_meta$Labels[which(result_meta$Var1 == "Soilorder123")] = "Vertisols"
result_meta$Labels[which(result_meta$Var1 == "Landform21")] = "Plains"
result_meta$Labels[which(result_meta$Var1 == "Landform22")] = "Plateaus"
result_meta$Labels[which(result_meta$Var1 == "Landform23")] = "Mountains"
result_meta$Labels[which(result_meta$Var1 == "Biomes1")] = "TS forests"
result_meta$Labels[which(result_meta$Var1 == "Biomes2")] = "TS grasslands/savannas"
result_meta$Labels[which(result_meta$Var1 == "Biomes3")] = "Temperate forests"
result_meta$Labels[which(result_meta$Var1 == "Biomes4")] = "Temperate grasslands"
result_meta$Labels[which(result_meta$Var1 == "Biomes5")] = "Med/Mon shrublands"
result_meta$Labels[which(result_meta$Var1 == "Biomes6")] = "Boreal forest"
result_meta$Labels[which(result_meta$Var1 == "Biomes7")] = "Tundra"
result_meta$Labels[which(result_meta$Var1 == "Biomes8")] = "Deserts"
result_meta$Labels[which(result_meta$Var1 == "Biomes9")] = "Croplands"
result_meta$Labels[which(result_meta$Var1 == "R.winter")] = "Winter type"
result_meta$Labels[which(result_meta$Var1 == "R.Summer")] = "Summer type"
result_meta$Labels[which(result_meta$Var1 == "uniformity")] = "Uniform"

albelss = c(result_meta$Labels[which(result_meta$temp_ch == 1 & result_meta$layers == "L1")], 
            result_meta$Labels[which(result_meta$temp_ch == 1 & result_meta$layers == "L2")],
            result_meta$Labels[which(result_meta$temp_ch == 1 & result_meta$layers == "L3")])

result_meta$names_long = factor(result_meta$names_long, 
                                levels = levelss, labels = albelss)
# write.csv(result_meta, "XXXX/result_meta.csv")

###
result_meta = read.csv("Fig data/result_meta.csv")
result_meta$layers  =factor(result_meta$layers, 
                            levels = c("L1", "L2", "L3"),
                            labels = c("0-0.3 m", "0.3-1 m","1-2 m"))

result_meta$sing[which(result_meta$ci.lb< 0 & result_meta$ci.ub > 0 )]="UNSIG"
result_meta$sing[which(result_meta$ci.lb< 0 & result_meta$ci.ub < 0 & result_meta$effect_value<0)]="SIG1"
result_meta$sing[which(result_meta$ci.lb> 0 & result_meta$ci.ub > 0 & result_meta$effect_value>0)]="SIG2"
result_meta$sing = factor(result_meta$sing)

result_meta$Freq[which(result_meta$Freq > 100)]=NA
result_meta$Type = factor(result_meta$Type)
result_meta$Labels = factor(result_meta$Labels, 
                            levels = rev(c("Boreal forest",          "Crops",                 
                                           "Deserts",                "Med shrublands",        
                                           "Temperate forests",      "Temperatre grasslands", 
                                           "TS forests",             "TS grasslands/savannas",
                                           "Tundra",                 "Summer type"  , 
                                           "Winter type" ,           "Uniform",           
                                           "Alfisols",    "Andisols",    "Aridisols",   "Entisols",  
                                           "Gelisols",    "Histosols",   "Inceptisols", "Mollisols",  
                                           "Oxisols",     "Spodosols",   "Ultisols",    "Vertisols",             
                                           "Mountains", "Plateaus", "Plains" , "Overall effect")))



result_meta$temp_ch = factor(result_meta$temp_ch, labels = paste0("+",1:5, "°C"))

# write.csv(result_meta, "result_meta.csv")

LABS = rev(c("Temperatre grasslands","Boreal forest","TS forests",
             "TS grasslands/savannas","Crops","Deserts",
             "Temperate forests","Med shrublands","Tundra",
             "Summer type","Winter type","Uniform",
             "Alfisols","Andisols","Aridisols",  
             "Entisols","Gelisols","Histosols","Inceptisols",
             "Mollisols","Oxisols","Spodosols","Ultisols",
             "Plateaus","Plains","Mountains","Overall effect"))

mixedsort(c("TS forests"  ,           "Temperate forests",        
            "Boreal forest",          "Temperatre grasslands",
            "TS grasslands/savannas", "Med shrublands", 
            "Deserts",                "Tundra"  ,             
            "Crops"))




id.big = which(abs(result_meta$ci.ub) > 50 | abs(result_meta$ci.lb) > 50)
out.result_meta = result_meta[id.big,]
xstart = c()
xend = c()
direct = c()
for (ii in 1:length(out.result_meta$vars)) {
  if(out.result_meta$ci.lb[ii]< -50 & out.result_meta$ci.ub[ii] > 50){
    xstart[ii] = -44
    xend[ii] = 42
    direct[ii] = "both"
  }else if(out.result_meta$ci.lb[ii]< -50 & out.result_meta$ci.ub[ii] < 50){
    xstart[ii] = -42
    xend[ii] = out.result_meta$ci.ub[ii]
    direct[ii] = "first"
  }else if(out.result_meta$ci.lb[ii]> -50 & out.result_meta$ci.ub[ii] > 50){
    xstart[ii] = out.result_meta$ci.lb[ii]
    xend[ii] = 42
    direct[ii] = "last"
  }
}
out.result_meta$"xstart" = xstart
out.result_meta$"xend" = xend
out.result_meta$"direct" = direct

col.text = rev(rep(c("#248232","#20A4F3","#F2BB05", "#685044","#011627"), c(9,3,12,3,1)))
face.text = rev(rep(c("plain","bold"), c(27,1)))


#### Extended Data Fig.3 
P_RR_all=
  ggplot(result_meta, aes(y=Labels, x= effect_value))+
  geom_errorbar(data = result_meta, aes(y=Labels, x= effect_value, xmin = ci.lb, xmax = ci.ub),
                position = position_dodge(0.4),width = 0, size = 0.3)+
  geom_segment(data=out.result_meta[which(out.result_meta$direct == "both"),], 
               mapping = aes(x=xstart, y=Labels, xend=xend, yend=Labels),
               arrow = arrow(ends = "both",length=unit(0.2, "cm")), size = 0.3)+
  
  geom_segment(data=out.result_meta[which(out.result_meta$direct == "first"),], 
               mapping = aes(x=xstart, y=Labels, xend=xend, yend=Labels),
               arrow = arrow(ends = "first",length=unit(0.2, "cm")), size = 0.3)+
  
  geom_segment(data=out.result_meta[which(out.result_meta$direct == "last"),], 
               mapping = aes(x=xstart, y=Labels, xend=xend, yend=Labels),
               arrow = arrow(ends = "last",length=unit(0.2, "cm")), size = 0.3)+
  
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey")+
  geom_point(aes(color = sing, shape = Type), show.legend = F, size = 2.5)+
  geom_text(aes(label = Freq, x=48), size=3)+
  xlab("Warming effect on soil carbon stock (%)")+
  ylab("")+
  facet_grid(layers~temp_ch, scale="free_y")+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(angle = 0, color = col.text, face = face.text), 
        axis.text.x = element_text(color = "black"))+
  scale_x_continuous(limits = c(-50, 50))+
  scale_color_manual(values = c( "#D62828","blue", "grey"))+
  scale_shape_manual(values = c( 16, 18))
P_RR_all

ggsave(P_RR_all, filename = "Fig/Fig_new/Fig.S.4.RR_all_stock_new.pdf",  width=8.27,height=11.69, unit="in",dpi=300) 
ggsave(P_RR_all, filename = "Fig/Fig_new/Fig1b.RR_all_new.png", width =8.27,height=11.69, unit="in",dpi=300) 

################################################################################
#### The code of analyzing SOC content is same with SOC stock







