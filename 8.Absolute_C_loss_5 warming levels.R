###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

###############################################################################
#library(sf)
library(foreach)
library(doSNOW)
library(parallel)
library(terra)
library(stringi)

rm(list = ls())

Input.dat = readRDS("XXXX/Sel_face_data/inputdata_datatable_face_WISE_sg_HWSD1.Rdata")
XY = readRDS("XXXX/Sel_face_data/XY_face_WISE_sg_HWSD1.Rdata")

path = "XXXX"
Temp <- rast(paste0(path, "wc2.1_30s_bio_1.tif"))
NPPi = rast(paste0(path,"MOD17A3_Science_NPP_mean_00_15.tif"))
Temp = crop(Temp, ext(NPPi))

###Soilgrids
sgocs = rast("XXXX/Source_data/ocs_gl1.tif")##unit:t/ha
sgocs = crop(sgocs, ext(NPPi))
sgocs = tapp(sgocs, index = c(1,1,1,2,2,3), sum)
sgocs = resample(sgocs, Temp, method = "near")

###HWSD
hwsd030 = rast("XXXX/finalgrid/t_oc/w001001.adf") ###unit:%
hwsd3010 = rast("XXXX/finalgrid/s_oc/w001000.adf")###unit:%
hwsd030.bd = rast("XXXX/finalgrid/t_bd/w001000.adf")###unit:kg dm-3
hwsd3010.bd = rast("XXXX/finalgrid/s_bd/w001000.adf")###unit:kg dm-3
hwsd030.gravel = rast("XXXX/finalgrid/t_gravel/w001000.adf")###unit£º%
hwsd3010.gravel = rast("XXXX/finalgrid/s_gravel/w001000.adf")###unit: %

hwsd.t.soc = hwsd030*hwsd030.bd*(1-hwsd030.gravel/100)*3*100*10000/1000/100 #unit: t/ha
hwsd.s.soc = hwsd3010*hwsd3010.bd*(1-hwsd3010.gravel/100)*7*100*10000/1000/100

hwsd = c(hwsd.t.soc, hwsd.s.soc)
hwsd = resample(hwsd, Temp, method = "near")
hwsd = crop(hwsd, ext(NPPi))

ocs = c(sgocs, hwsd)
names(ocs) = c("sg.L1", "sg.L2","sg.L3","hwsd.L1","hwsd.L2")
rm(sgocs, hwsd, hwsd.t.soc, hwsd.s.soc, hwsd030, hwsd3010, hwsd030.bd, hwsd3010.bd, hwsd030.gravel, hwsd3010.gravel)


###=============================================================================

#### Table.1 and Extended Data Table 3.


mod = readRDS("XXXX/mod_Stock_20sample_unscale.Rdata")
mod.final = mod$finalModel

Input.dat = readRDS("XXXX/inputdata_datatable_face_WISE_sg_HWSD1.Rdata")
df.Input.dat = as.data.frame(do.call(rbind, Input.dat))

socnames = c("wise.soc.L1","wise.soc.L2","wise.soc.L3", "sg.L1","sg.L2","sg.L3","hwsd.L1", "hwsd.L2")
tt = c(1:5)
Varspace = expand.grid( SOClyr = socnames, deltaT = tt)
Varspace$SOClyr = as.character(Varspace$SOClyr)

input.base = df.Input.dat[,names(df.Input.dat)%in%c("iPre","iTemp","ai_et0","iPC","P.summer","elevation","NPP","Biomes","SoilOrder",
                                                    "Landform", "Seasonality.pre","Pre.diff","Temp.diff","Soil.layers")]
fun_pred_T5 = function(hh){

  input.soc = df.Input.dat[,names(df.Input.dat)%in%Varspace[hh, 1]] 
  pred.df = data.frame(input.soc, input.base)
  names(pred.df) = c( "soc","iTemp","iPre","iPC","P.summer","ai_et0","elevation","NPP","SoilOrder","Landform",
                      "Seasonality.pre","Biomes","Temp.diff","Soil.layers")

  pred.df$Soil.layers = as.numeric(substr(Varspace[hh, 1], nchar(Varspace[hh, 1]), nchar(Varspace[hh, 1])))

  pred.df$Biomes = factor(pred.df$Biomes)
  pred.df$SoilOrder = factor(pred.df$SoilOrder)
  pred.df$Landform = factor(pred.df$Landform)
  pred.df$Seasonality.pre = factor(pred.df$Seasonality.pre)
  pred.df$Soil.layers = factor(pred.df$Soil.layers)
  pred.df$Temp.diff = Varspace[hh, 2]
  pred.df$"Pre.diff" = 0
  
  pred.df.rNA = pred.df[complete.cases(pred.df), ]
  pred.df.rNA$soc = pred.df.rNA$soc
  
  pred = predict(mod.final, pred.df.rNA, type = "se")
  pred.rr = array(NA, dim = c(nrow(pred.df), 2))
  
  pred.rr[complete.cases(pred.df), 1] = (exp(pred$predictions) - 1) * 100
  pred.rr[complete.cases(pred.df), 2] = (exp(pred$se) -1) *100
  colnames(pred.rr) = c("rr", "se")
  # summary(pred.rr$rr)
  
  sol = cbind(pred.rr, soc = pred.df$soc, biomes = as.numeric(as.character(pred.df$Biomes)))
  saveRDS(sol, paste0("I:/Select profile_source stock data/Pred_df_XY_SOCs/Pred_3dataset/new_Pred_", 
                      Varspace[hh, 1],Varspace[hh, 2],".Rdata"))
}

cl = makeCluster(15)
registerDoSNOW(cl)

foreach(hh = 1:nrow(Varspace), .combine = list) %dopar% {
  fun_pred_T5(hh)
}
stopCluster(cl)










