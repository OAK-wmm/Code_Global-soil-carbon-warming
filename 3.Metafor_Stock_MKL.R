###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

rm(list = ls())
setwd("XXXX/")

#==========================meta analysis===================================
library(metafor)
library(Matrix)

filenames = list.files("dat_Mean_sd_data/")

setMKLthreads(40)

id = c(46:50)

for (ii in 1:length(id)) { 
  
  order.rep = gsub("dat_Mean_sd_Stock_","",filenames[id[ii]])
  data.meta = readRDS(paste0("dat_Mean_sd_data/", filenames[id[ii]]))
  dat.sampsize_soilorder112 = data.meta[which(data.meta$SoilOrder == "Soilorder112" &
                                                data.meta$Num.pro.con >=10 & data.meta$Num.pro.tre >= 10), ]
  

  data.meta$Biomes = paste0("Biomes",data.meta$Biomes)
  data.meta$SoilOrder = paste0("Soilorder", data.meta$SoilOrder)
  data.meta$Landform = paste0("Landform", data.meta$Landform)

  df.yivi.L1 = escalc(data =data.meta, measure = "ROM",
                      m1i = soc.0_30cm.tre, sd1i = sd.soc.0_30cm.tre, n1i = Num.pro.tre,
                      m2i = soc.0_30cm.con, sd2i = sd.soc.0_30cm.con,  n2i = Num.pro.con)
  df.yivi.L1$id <- 1:nrow(df.yivi.L1)
  id.del = which(abs(df.yivi.L1$yi) > 100 | abs(df.yivi.L1$vi) > 100 |
                   df.yivi.L1$vi ==0 | is.na(df.yivi.L1$vi) == T |
                   is.na(df.yivi.L1$yi) == T)
  if(length(id.del) == 0){
    df.yivi.L1 = df.yivi.L1
  }else{
    df.yivi.L1 = df.yivi.L1[-id.del,]
  }
  saveRDS(df.yivi.L1, paste0("dat_yivi_select9/df.yivi.L1_", gsub("Mean_sd_Stock_diff_", "", filenames[id[ii]])))


  df.yivi.L2 = escalc(data =data.meta, measure = "ROM",
                      m1i = soc.30_100cm.tre, sd1i = sd.soc.30_100cm.tre, n1i = Num.pro.tre,
                      m2i = soc.30_100cm.con, sd2i = sd.soc.30_100cm.con,  n2i = Num.pro.con)
  df.yivi.L2$id <- 1:nrow(df.yivi.L2)
  id.del = which(abs(df.yivi.L2$yi) > 100 | abs(df.yivi.L2$vi) > 100 |
                   abs(df.yivi.L2$vi)==0 | is.na(df.yivi.L2$vi) == T |
                   is.na(df.yivi.L2$yi) == T)
  if(length(id.del) == 0){
    df.yivi.L2 = df.yivi.L2
  }else{
    df.yivi.L2 = df.yivi.L2[-id.del,]
  }
  saveRDS(df.yivi.L2, paste0("dat_yivi_select9/df.yivi.L2_",gsub("Mean_sd_Stock_diff_", "", filenames[id[ii]])))


  df.yivi.L3 = escalc(data =data.meta, measure = "ROM",
                      m1i = soc.100_200cm.tre, sd1i = sd.soc.100_200cm.tre, n1i = Num.pro.tre,
                      m2i = soc.100_200cm.con, sd2i = sd.soc.100_200cm.con,  n2i = Num.pro.con)
  df.yivi.L3$id <- 1:nrow(df.yivi.L3)
  id.del = which(abs(df.yivi.L3$yi) > 100 | abs(df.yivi.L3$vi) > 100 |
                   abs(df.yivi.L3$vi)==0 | is.na(df.yivi.L3$vi) == T |
                   is.na(df.yivi.L3$yi) == T)

  if(length(id.del) == 0){
    df.yivi.L3 = df.yivi.L3
  }else{
    df.yivi.L3 = df.yivi.L3[-id.del,]
  }
  saveRDS(df.yivi.L3, paste0("dat_yivi_select9/df.yivi.L3_",gsub("Mean_sd_Stock_diff_", "", filenames[id[ii]])))
  
  
  
  mod.meta.Biomes.L1 <-rma.mv(yi, vi, mods = ~biomes_re-1, random = ~ 1 | id, data=df.yivi.L1, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L1, paste0("dat_meat_0intrcp/mod.meta.Biomes.L1.",order.rep))

  effect_value = (exp(mod.meta.Biomes.L1$beta)-1)*100
  p_value = mod.meta.Biomes.L1$pval
  ci.lb = mod.meta.Biomes.L1$ci.lb
  ci.ub = mod.meta.Biomes.L1$ci.ub
  QM = rep(mod.meta.Biomes.L1$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L1$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L1$biomes_re)
  layers = rep("L1", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("test_dat_meat_0intrcp_20samplesize/Meta.sol.Biomes.L1.",order.rep))


  mod.meta.Biomes.L2 <-rma.mv(yi, vi,  mods = ~biomes_re-1, random = ~ 1 | id, data=df.yivi.L2, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L2, paste0("dat_meat_0intrcp/mod.meta.Biomes.L2.",order.rep))

  effect_value = (exp(mod.meta.Biomes.L2$beta)-1)*100
  p_value = mod.meta.Biomes.L2$pval
  ci.lb = mod.meta.Biomes.L2$ci.lb
  ci.ub = mod.meta.Biomes.L2$ci.ub
  QM = rep(mod.meta.Biomes.L2$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L2$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L2$biomes_re)
  layers = rep("L2", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("test_dat_meat_0intrcp_20samplesize/Meta.sol.Biomes.L2.",order.rep))

  
  mod.meta.Biomes.L3 <-rma.mv(yi, vi,  mods = ~biomes_re-1, random = ~ 1 | id, data=df.yivi.L3, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L3, paste0("dat_meat_0intrcp/mod.meta.Biomes.L3.",order.rep))

  effect_value = (exp(mod.meta.Biomes.L3$beta)-1)*100
  p_value = mod.meta.Biomes.L3$pval
  ci.lb = mod.meta.Biomes.L3$ci.lb
  ci.ub = mod.meta.Biomes.L3$ci.ub
  QM = rep(mod.meta.Biomes.L3$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L3$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L3$biomes_re)
  layers = rep("L3", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("test_dat_meat_0intrcp_20samplesize/Meta.sol.Biomes.L3.",order.rep))

  ###########################
  ### the relationship between RR and Biomes
  mod.meta.Biomes.L1 <-rma.mv(yi, vi, mods = ~Biomes-1, random = ~ 1 | id, data=df.yivi.L1, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L1, paste0("dat_meat_0intrcp/mod.meta.Biomes.L1.",ii, order.rep))

  effect_value = (exp(mod.meta.Biomes.L1$beta)-1)*100
  p_value = mod.meta.Biomes.L1$pval
  ci.lb = mod.meta.Biomes.L1$ci.lb
  ci.ub = mod.meta.Biomes.L1$ci.ub
  QM = rep(mod.meta.Biomes.L1$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L1$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L1$Biomes)
  layers = rep("L1", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Biomes.L1.",order.rep))
  #
  mod.meta.Biomes.L2 <-rma.mv(yi, vi,  mods = ~Biomes-1, random = ~ 1 | id, data=df.yivi.L2, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L2, paste0("dat_meat_0intrcp/mod.meta.Biomes.L2.",ii, order.rep))
  effect_value = (exp(mod.meta.Biomes.L2$beta)-1)*100
  p_value = mod.meta.Biomes.L2$pval
  ci.lb = mod.meta.Biomes.L2$ci.lb
  ci.ub = mod.meta.Biomes.L2$ci.ub
  QM = rep(mod.meta.Biomes.L2$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L2$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L2$Biomes)
  layers = rep("L2", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Biomes.L2.",order.rep))
  #
  mod.meta.Biomes.L3 <-rma.mv(yi, vi,  mods = ~Biomes-1, random = ~ 1 | id, data=df.yivi.L3, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Biomes.L3, paste0("dat_meat_0intrcp/mod.meta.Biomes.L3.",ii, order.rep))
  effect_value = (exp(mod.meta.Biomes.L3$beta)-1)*100
  p_value = mod.meta.Biomes.L3$pval
  ci.lb = mod.meta.Biomes.L3$ci.lb
  ci.ub = mod.meta.Biomes.L3$ci.ub
  QM = rep(mod.meta.Biomes.L3$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Biomes.L3$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L3$Biomes)
  layers = rep("L3", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Biomes.L3.",order.rep))
 
  
  # ###########################
  # ### the relationship between RR and SoilOrder
  mod.meta.SoilOrder.L1 <-rma.mv(yi, vi, mods = ~SoilOrder-1, random = ~ 1 | id, data=df.yivi.L1, method="REML", sparse=TRUE)
  saveRDS(mod.meta.SoilOrder.L1, paste0("dat_meat_0intrcp/mod.meta.SoilOrder.L1.",ii, order.rep))
  effect_value = (exp(mod.meta.SoilOrder.L1$beta)-1)*100
  p_value = mod.meta.SoilOrder.L1$pval
  ci.lb = mod.meta.SoilOrder.L1$ci.lb
  ci.ub = mod.meta.SoilOrder.L1$ci.ub
  QM = rep(mod.meta.SoilOrder.L1$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.SoilOrder.L1$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L1$SoilOrder)
  layers = rep("L1", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.soilorder.L1.",order.rep))
  #
  mod.meta.SoilOrder.L2 <-rma.mv(yi, vi,  mods = ~SoilOrder-1, random = ~ 1 | id, data=df.yivi.L2, method="REML", sparse=TRUE)
  saveRDS(mod.meta.SoilOrder.L2, paste0("dat_meat_0intrcp/mod.meta.SoilOrder.L2.",ii, order.rep))
  effect_value = (exp(mod.meta.SoilOrder.L2$beta)-1)*100
  p_value = mod.meta.SoilOrder.L2$pval
  ci.lb = mod.meta.SoilOrder.L2$ci.lb
  ci.ub = mod.meta.SoilOrder.L2$ci.ub
  QM = rep(mod.meta.SoilOrder.L2$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.SoilOrder.L2$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L2$SoilOrder)
  layers = rep("L2", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.soilorder.L2.",order.rep))
  #
  mod.meta.SoilOrder.L3 <-rma.mv(yi, vi,  mods = ~SoilOrder-1, random = ~ 1 | id, data=df.yivi.L3, method="REML", sparse=TRUE)
  saveRDS(mod.meta.SoilOrder.L3, paste0("dat_meat_0intrcp/mod.meta.SoilOrder.L3.",ii, order.rep))
  effect_value = (exp(mod.meta.SoilOrder.L3$beta)-1)*100
  p_value = mod.meta.SoilOrder.L3$pval
  ci.lb = mod.meta.SoilOrder.L3$ci.lb
  ci.ub = mod.meta.SoilOrder.L3$ci.ub
  QM = rep(mod.meta.SoilOrder.L3$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.SoilOrder.L3$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L3$SoilOrder)
  layers = rep("L3", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.soilorder.L3.",order.rep))
 
  # ###########################
  # ### the relationship between RR and Landform
  mod.meta.Landform.L1 <-rma.mv(yi, vi, mods = ~Landform-1, random = ~ 1 | id, data=df.yivi.L1, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Landform.L1, paste0("dat_meat_0intrcp/mod.meta.Landform.L1.",ii, order.rep))
  effect_value = (exp(mod.meta.Landform.L1$beta)-1)*100
  p_value = mod.meta.Landform.L1$pval
  ci.lb = mod.meta.Landform.L1$ci.lb
  ci.ub = mod.meta.Landform.L1$ci.ub
  QM = rep(mod.meta.Landform.L1$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Landform.L1$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L1$Landform)
  layers = rep("L1", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Landform.L1.",order.rep))
  #
  mod.meta.Landform.L2 <-rma.mv(yi, vi,  mods = ~Landform-1, random = ~ 1 | id, data=df.yivi.L2, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Landform.L2, paste0("dat_meat_0intrcp/mod.meta.Landform.L2.",ii, order.rep))
  effect_value = (exp(mod.meta.Landform.L2$beta)-1)*100
  p_value = mod.meta.Landform.L2$pval
  ci.lb = mod.meta.Landform.L2$ci.lb
  ci.ub = mod.meta.Landform.L2$ci.ub
  QM = rep(mod.meta.Landform.L2$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Landform.L2$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L2$Landform)
  layers = rep("L2", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Landform.L2.",order.rep))

  mod.meta.Landform.L3 <-rma.mv(yi, vi,  mods = ~Landform-1, random = ~ 1 | id, data=df.yivi.L3, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Landform.L3, paste0("dat_meat_0intrcp/mod.meta.Landform.L3.",ii, order.rep))

  effect_value = (exp(mod.meta.Landform.L3$beta)-1)*100
  p_value = mod.meta.Landform.L3$pval
  ci.lb = mod.meta.Landform.L3$ci.lb
  ci.ub = mod.meta.Landform.L3$ci.ub
  QM = rep(mod.meta.Landform.L3$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Landform.L3$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L3$Landform)
  layers = rep("L3", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Landform.L3.",order.rep))
  
  
  # ###########################
  # ### the relationship between RR and Seasonality.pre
  mod.meta.Seasonality.pre.L1 <-rma.mv(yi, vi, mods = ~Seasonality.pre-1, random = ~ 1 | id, data=df.yivi.L1, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Seasonality.pre.L1, paste0("dat_meat_0intrcp/mod.meta.Seasonality.pre.L1",ii, order.rep))
  effect_value = (exp(mod.meta.Seasonality.pre.L1$beta)-1)*100
  p_value = mod.meta.Seasonality.pre.L1$pval
  ci.lb = mod.meta.Seasonality.pre.L1$ci.lb
  ci.ub = mod.meta.Seasonality.pre.L1$ci.ub
  QM = rep(mod.meta.Seasonality.pre.L1$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Seasonality.pre.L1$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L1$Seasonality.pre)
  layers = rep("L1", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Seasonality.pre.L1.",order.rep))
  
  
  mod.meta.Seasonality.pre.L2 <-rma.mv(yi, vi,  mods = ~Seasonality.pre-1, random = ~ 1 | id, data=df.yivi.L2, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Seasonality.pre.L2, paste0("dat_meat_0intrcp/mod.meta.Seasonality.pre.L2",ii, "degres.Rdata"))
  effect_value = (exp(mod.meta.Seasonality.pre.L2$beta)-1)*100
  p_value = mod.meta.Seasonality.pre.L2$pval
  ci.lb = mod.meta.Seasonality.pre.L2$ci.lb
  ci.ub = mod.meta.Seasonality.pre.L2$ci.ub
  QM = rep(mod.meta.Seasonality.pre.L2$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Seasonality.pre.L2$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L2$Seasonality.pre)
  layers = rep("L2", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Seasonality.pre.L2.",order.rep))

  
  mod.meta.Seasonality.pre.L3 <-rma.mv(yi, vi,  mods = ~Seasonality.pre-1, random = ~ 1 | id, data=df.yivi.L3, method="REML", sparse=TRUE)
  saveRDS(mod.meta.Seasonality.pre.L3, paste0("dat_meat_0intrcp/mod.meta.Seasonality.pre.L3",ii, order.rep))
  effect_value = (exp(mod.meta.Seasonality.pre.L3$beta)-1)*100
  p_value = mod.meta.Seasonality.pre.L3$pval
  ci.lb = mod.meta.Seasonality.pre.L3$ci.lb
  ci.ub = mod.meta.Seasonality.pre.L3$ci.ub
  QM = rep(mod.meta.Seasonality.pre.L3$QM, times = dim(effect_value)[1])
  vars = rownames(effect_value)
  p_QM = rep(mod.meta.Seasonality.pre.L3$QMp, times = dim(effect_value)[1])
  num = table(df.yivi.L3$Seasonality.pre)
  layers = rep("L3", length(effect_value))
  temp_ch = rep(ii, length(effect_value))
  sol = list(data.frame(vars,layers, temp_ch, effect_value, p_value, ci.lb, ci.ub, QM, p_QM), num)
  saveRDS(sol, paste0("dat_meat_0intrcp/Meta.sol.Seasonality.pre.L3.",order.rep))

}














