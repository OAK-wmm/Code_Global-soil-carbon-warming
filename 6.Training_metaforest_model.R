###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

rm(list = ls())
setwd("XXXX")


library(metaforest)
library(caret)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(rcartocolor)
options(bitmapType="cairo")
library(parallel)
library(doParallel)
library(e1071)
library(doSNOW)
library(RColorBrewer)

dat_all = readRDS("df.yivi.env.stock.Rdata")
dat_all = na.omit(dat_all)
yi = dat_all$yi
id.out = which(yi %in% boxplot.stats(yi)$out)
dat_all = dat_all[-id.out,]

vars_del = c("Num.pro.con", "Num.pro.tre", "Temp_ch")

dat_all$Soil.layers[which(dat_all$Soil.layers == "0-30 cm")] = 1
dat_all$Soil.layers[which(dat_all$Soil.layers == "30-100 cm")] = 2
dat_all$Soil.layers[which(dat_all$Soil.layers == "100-200 cm")] = 3
dat_all$Soil.layers = factor(dat_all$Soil.layers)

dat_all$Seasonality.pre[which(dat_all$Seasonality.pre =="R.Summer")] = 11
dat_all$Seasonality.pre[which(dat_all$Seasonality.pre =="R.winter")] = 12
dat_all$Seasonality.pre[which(dat_all$Seasonality.pre =="uniformity")] = 13
dat_all$Seasonality.pre = factor(dat_all$Seasonality.pre)


dat_all$Biomes = factor(as.numeric(gsub("Biomes","",dat_all$Biomes)))
dat_all$Landform = factor(as.numeric(gsub("Landform", "",dat_all$Landform)))
dat_all$SoilOrder = factor(as.numeric(gsub("Soilorder","",dat_all$SoilOrder)))

set.seed(7564)
id.test = sample(1:nrow(dat_all), ceiling(0.2*nrow(dat_all)))
dat_all.test = dat_all[id.test,]
dat_all.train = dat_all[-id.test,]


cpu.num = 50
cl = makeCluster(cpu.num)
registerDoParallel(cl) # Run in parallel -> 31 cores
#
set.seed(6839)
cv_folds <- trainControl(method = "cv", number=10, allowParallel = T)
#
# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid <- expand.grid(whichweights = c("random", "fixed"),
                           mtry = 2:6,
                           min.node.size = 2:6)
#
metaforest_10cv_Stock <- train(y=dat_all.train$yi,
                               x = dat_all.train[,!names(dat_all.train) %in% "yi"],
                               method = ModelInfo_mf(),
                               trControl = cv_folds,
                               tuneGrid = tuning_grid,
                               keep.inbag = TRUE,
                               num.threads = cpu.num,
                               verbose=TRUE,
                               num.trees = 2000)
stopCluster(cl)
saveRDS(metaforest_10cv_Stock, "mod_metaforest/mod_Stock_20sample_unscale.Rdata")
metaforest_10cv_Stock$finalModel$forest$r.squared
VarImpPlot(metaforest_10cv_Stock$finalModel)

mod = readRDS("mod_metaforest/mod_Stock_20sample_unscale.Rdata")
mod.final = mod$finalModel
importance.stock = mod.final$forest$variable.importance

pred.validation = predict(mod.final, dat_all.test, type = "se")
pred.calibation = predict(mod.final, dat_all.train, type = "se")

df.pred.val = data.frame(pred = pred.validation$predictions, obs = dat_all.test$yi)
df.pred.cal = data.frame(pred = pred.calibation$predictions, obs = dat_all.train$yi)

write.csv(importance.stock, "XXXX/importance.stock.csv")
write.csv(df.pred.val, "XXXX/df.pred.val.stock.csv")
write.csv(df.pred.cal, "XXXX/df.pred.cal.stock.csv")
write.csv(dat_all, "XXXX/dat_all.csv")

plot(metaforest_10cv_Stock$finalModel)
varImp(metaforest_10cv_Stock)
mod_final = metaforest_10cv_Stock$finalModel
# saveRDS(mod_final, "mod_final_Stock.Rdata")

importance_stock = sort(mod_final$forest$variable.importance)
write.csv(importance_stock, "Fig data/importance_stock_20sample_unscale_new.csv")

p_PartialDependence_stock = PartialDependence(mod_final, pi=0.95)

ggsave("Fig/p_PartialDependence_stock.pdf",p_PartialDependence_stock, width =  12, height = 8.3, units = "in")

plot_conv_stock <- plot(mod_final)
ggsave("Fig/plot_conv_stock.pdf",plot_conv_stock, width =  6.875, units = "in")

