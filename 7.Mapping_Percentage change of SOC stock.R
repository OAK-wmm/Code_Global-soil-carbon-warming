###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

library(sf)

rm(list = ls())

Input.dat = readRDS("XXXX/inputdata_datatable_face_1.Rdata")
XY = readRDS("XXXX/XY_face_1.Rdata")
mod = readRDS("XXXX/mod_Stock_20sample_unscale.Rdata")
mod.final = mod$finalModel
mod$finalModel$forest$variable.importance

df.Input.dat = as.data.frame(do.call(rbind, Input.dat))

socnames = c("soc.L1","soc.L2","soc.L3")
lay = c("L1","L2","L3")
input.base = df.Input.dat[,names(df.Input.dat)%in%c("iPre","iTemp","ai_et0","iPC","P.summer","elevation","NPP","Biomes","SoilOrder",
                                                    "Landform", "Seasonality.pre","Pre.diff","Temp.diff","Soil.layers")]

for(jj in 1:3) {
  input.soc = df.Input.dat[,names(df.Input.dat)%in%socnames[jj]] 
  pred.df = data.frame(input.soc, input.base)
  names(pred.df) = c( "soc","iTemp","iPre","iPC","P.summer","ai_et0","elevation","NPP","SoilOrder","Landform",
                      "Seasonality.pre","Biomes","Temp.diff","Soil.layers")
  if(jj==1){
    pred.df$Soil.layers = 1
  }else if(jj==2){
    pred.df$Soil.layers = 2
  }else{
    pred.df$Soil.layers = 3
  }
  
  pred.df$Biomes = factor(pred.df$Biomes)
  pred.df$SoilOrder = factor(pred.df$SoilOrder)
  pred.df$Landform = factor(pred.df$Landform)
  pred.df$Seasonality.pre = factor(pred.df$Seasonality.pre)
  pred.df$Soil.layers = factor(pred.df$Soil.layers)
  pred.df$Temp.diff = 2
  pred.df$"Pre.diff" = 0
  
  pred.df.rNA = pred.df[complete.cases(pred.df), ]
  pred.df.rNA$soc = pred.df.rNA$soc
  
  pred = predict(mod.final, pred.df.rNA, type = "se")
  pred.rr = as.data.frame(array(NA, dim = c(nrow(pred.df), 2)))
  names(pred.rr) = c("rr", "se")
  
  pred.rr[complete.cases(pred.df), 1] = (exp(pred$predictions) - 1) * 100
  pred.rr[complete.cases(pred.df), 2] = (exp(pred$se) -1) *100
  summary(pred.rr$rr)
  
  rr = paste0("rr", lay[jj])
  se = paste0("se", lay[jj])
  
  if(jj==1){
    for (ii in 1:length(XY)){
      XY[[ii]]$"rrL1" = pred.rr$rr[ii]
      XY[[ii]]$"seL1" = pred.rr$se[ii]
      XY[[ii]]$"biomess" = as.numeric(as.character(pred.df$Biomes[ii]))
    }
  }else if (jj==2){
    for (ii in 1:length(XY)){
      XY[[ii]]$"rrL2" = pred.rr$rr[ii]
      XY[[ii]]$"seL2" = pred.rr$se[ii]
    }
  }else{
    for (ii in 1:length(XY)){
      XY[[ii]]$"rrL3" = pred.rr$rr[ii]
      XY[[ii]]$"seL3" = pred.rr$se[ii]
    }
  }
}
XY.pred = do.call(rbind, XY)
saveRDS(XY.pred,"XXXX/XY.pred.stock.Rdata")

################################################################################
library(terra)
library(RColorBrewer)
library(rgdal)
library(dplyr)
library(maps)

rm(list = ls())

bbox <- readOGR(dsn="XXXX/ne_50m_wgs84_bounding_box.shp", layer="ne_50m_wgs84_bounding_box")
bbox <- spTransform(bbox, CRS("+proj=robin"))

ocean <- readOGR(dsn="XXXX/ne_50m_ocean.shp", layer="ne_50m_ocean")
ocean <- spTransform(ocean, CRS("+proj=robin"))

wmap <- readOGR(dsn="XXXX/ne_50m_land.shp", layer="ne_50m_land")
wmap_wgs <- spTransform(wmap, CRS("+proj=robin"))


XY.pred = readRDS("XXXX/XY.pred.stock.Rdata")
crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

maps.l1 = rast(XY.pred[,c("x","y", "rrL1")], crs=crs, type = "xyz")
maps.l2 = rast(XY.pred[,c("x","y", "rrL2")], crs=crs, type = "xyz")
maps.l3 = rast(XY.pred[,c("x","y", "rrL3")], crs=crs, type = "xyz")
maps.rr = c(maps.l1, maps.l2, maps.l3)
rm(maps.l1, maps.l2, maps.l3)

maps.robin.rr = terra::project(maps.rr, bbox)
rm(maps.rr)
# writeRaster(maps.robin.rr,
#             filename = file.path("H:/Select profile_source stock data/Maps_ribon/maps.robin.stock.rr.tif"))


# cols = c(colorRampPalette(rev(brewer.pal(9, 'YlOrRd')))(6), "black")
cols = c((brewer.pal(11, 'RdYlBu'))[c(1:6, 8:10)] )
cols.eco = c(rev(brewer.pal(9, "BuGn")[c(1:7)]), "#F0C808", "#EE964B")
breaks = c(-100, seq(-50,20, 10), 1000)
leg.txt = c("<-50","-50~-40","-40~-30","-30~-20","-20~-10","-10~0","0~10", "10~20",">20")
eco.type = c("TS forests",
             "Temperate forests",        
             "Temperate grasslands",
             "TS grasslands/savannas",
             "Med/Mon shrublands", 
             "Boreal forests", 
             "Tundra", 
             "Deserts",            
             "Croplands")
maxpixels=5000000

fun_cols.eco = function(XX, breaks, cols, ga){
  cols.eco = c()
  ccs = max(XX)
  for (ii in 1:length(XX)) {
    id = which(breaks >= XX[ii] & breaks < XX[ii]+ga)
    while(length(id) == 0){
      id = which(breaks >= XX[ii] & breaks < ccs+ga)
      ccs = ccs + ga
    }
    cols.eco[ii] = cols[id-1]
  }
  return(cols.eco)
}

set.seed(2395)
XY.pred.sample = XY.pred[sample(1:dim(XY.pred)[1], maxpixels),]
XY.pred.sample$biomess = factor(XY.pred.sample$biomess)
rrL1.mean = tapply(XY.pred.sample$rrL1, list(XY.pred.sample$biomess), function(x)mean(x, na.rm = T))
rrL2.mean = tapply(XY.pred.sample$rrL2, list(XY.pred.sample$biomess), function(x)mean(x, na.rm = T))
rrL3.mean = tapply(XY.pred.sample$rrL3, list(XY.pred.sample$biomess), function(x)mean(x, na.rm = T))
###BLUE : #A6DEF6  grey:#DCDCDD

###Mapping the change of SOC stock to 2°C warming in 3 soil layers and among biomes
nf <- layout(matrix(c(1,2,3,4,5,6,7,8),4), widths = c(1,1), heights = c(2,2,2,1))
layout.show(nf)

par(pin=c(4, 2))
plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$rrL1,  useRaster=TRUE, maxcell=maxpixels, bty = "n",
      axes = F, col = cols, breaks = breaks, add=T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("a"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("0 - 0.3 m"), cex = 1.3, font=1)


plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$rrL2,  useRaster=TRUE, maxcell=maxpixels,bty = "n", axef = F,
      ylim = c(-90, 90), col = cols, legend = F, breaks = breaks, add=T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("c"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("0.3 - 1 m"), cex = 1.3, font=1)


plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$rrL3,  useRaster=TRUE, maxcell=maxpixels,bty = "n", axef = F,
      ylim = c(-90, 90), col = cols, legend = F, breaks = breaks, add = T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("e"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("1 - 2 m"), cex = 1.3, font=1)


par(pin=c(4, 0.55))
plot(c(1,0), c(0,1), ylim = c(0,10), xlim = c(0.5,10.5),
     type = "n", ylab = "",
     xaxt="n",yaxt="n",bty="n",
     xlab= ("Percentage changes of soil carbon stock under 2C°warming (%)"))

leg.txt1 = c("<-50","-50","-40","-30","-20","-10","0","10", "20",">20")
axis(1,at=seq(1, 10,1), labels = leg.txt1,  pos=8)
dw=1
sag <-seq(1,9, by=1)
rect(0+sag, 8,dw+sag, 13,  col = cols, border = NA)
text(x=5.6, y=1, labels = c("Percentage changes of soil carbon stock under 2°C warming (%)"))


par(pin=c(3, 1.8))
boxplot(rrL1~biomess, data = XY.pred.sample, ylab=c("Percentage change (%)"),
        xaxt = "n", col = fun_cols.eco(rrL1.mean, breaks, cols, ga=10), ylim = c(-80, 210))
points(x=c(1:9), y=rrL1.mean, col = "red", pch = 16, cex = 1.7)
abline(a=0, b=0, lty=2, col = "grey")
# text(x=0.5, y=210, labels = c("b"), cex = 1.5, font=2)
text(x=1.2, y=208, labels = c("0 - 0.3 m"), cex = 1.3, font=1)


cols.eco = data.frame(cols, ord= order(rrL2.mean)) %>% arrange( ord)
boxplot(rrL2~biomess, data = XY.pred.sample, ylab=c("Percentage change (%)"),
        xaxt="n", col = fun_cols.eco(rrL2.mean, breaks, cols, ga=10), ylim = c(-80, 210))
points(x=c(1:9), y=rrL2.mean, col = "red", pch = 16, cex = 1.7)
abline(a=0, b=0, lty=2, col = "grey")
# text(x=0.5, y=210, labels = c("d"), cex = 1.5, font=2)
text(x=1.2, y=208, labels = c("0.3 - 1 m"), cex = 1.3, font=1)


cols.eco = data.frame(cols, ord= order(rrL2.mean)) %>% arrange( ord)
boxplot(rrL3~biomess, data = XY.pred.sample, ylab=c("Percentage change (%)"),
        xaxt="n", col = fun_cols.eco(rrL3.mean, breaks, cols, ga=10), ylim = c(-80, 210))
abline(a=0, b=0, lty=2, col = "grey")
points(x=c(1:9), y=rrL3.mean, col = "red", pch = 16, cex = 1.7)
# text(x=0.5, y=210, labels = c("f"), cex = 1.5, font=2)
text(x=1, y=208, labels = c("1 - 2 m"), cex = 1.3, font=1)

par(pin=c(3, 0.6))
plot(c(1,0), c(0,1), ylim = c(0,3), xlim = c(0.5,9.5),
     type = "n", ylab = "",
     xaxt="n",yaxt="n",bty="n",
     xlab= ("Biomes"))

axis(1,at=seq(1, 9,1), labels = c("","","","","","","","",""),  las = 2, pos=4.5)
text(x=1:9, y = 4, srt = 45, adj=1, labels = eco.type, xpd = T)
##size:8.02*7.01
# dev.off()

######################################################
##Mapping uncetinty
maps.se.l1 = rast(XY.pred[,c("x","y", "seL1")], crs=crs, type = "xyz")
maps.se.l2 = rast(XY.pred[,c("x","y", "seL2")], crs=crs, type = "xyz")
maps.se.l3 = rast(XY.pred[,c("x","y", "seL3")], crs=crs, type = "xyz")
rm(XY.pred)
maps.se = c(maps.se.l1, maps.se.l2, maps.se.l3)
rm(maps.se.l1, maps.se.l2, maps.se.l3)
maps.robin.rr = terra::project(maps.se, bbox)
rm(maps.se)
# writeRaster(maps.robin.rr,filename = file.path("D:/MM/Mpas_robin/maps.robin.rr.tif"))

cols.unc = colorRampPalette(rev((brewer.pal(11, 'RdYlBu'))[c(1:6)]))(11)
cols.eco = c(rev(brewer.pal(9, "BuGn")[c(1:7)]), "#F0C808", "#EE964B")
breaks = c(seq(0,20, 2),70)
leg.txt = c("<-50","-50~-40","-40~-30","-30~-20","-20~-10","-10~0","0~10", "10~20",">20")
eco.type = c("TS forests",
             "Temperate forests",        
             "Temperate grasslands",
             "TS grasslands/savannas",
             "Med/Mon shrublands", 
             "Boreal forests", 
             "Tundra", 
             "Deserts",            
             "Croplands")
maxpixels=5000000

XY.pred.unc.sample = XY.pred[sample(1:dim(XY.pred)[1], maxpixels),]
XY.pred.unc.sample$biomess = factor(XY.pred.unc.sample$biomess)
seL1.mean = tapply(XY.pred.unc.sample$seL1, list(XY.pred.unc.sample$biomess), function(x)mean(x, na.rm = T))
seL2.mean = tapply(XY.pred.unc.sample$seL2, list(XY.pred.unc.sample$biomess), function(x)mean(x, na.rm = T))
seL3.mean = tapply(XY.pred.unc.sample$seL3, list(XY.pred.unc.sample$biomess), function(x)mean(x, na.rm = T))
#

nf <- layout(matrix(c(1,2,3,4,5,6,7,8),4), widths = c(1,1), heights = c(2,2,2,1))
layout.show(nf)

###BLUE : #A6DEF6  grey:#DCDCDD

par(pin=c(4, 2))
plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$seL1,  useRaster=TRUE, maxcell=maxpixels, bty = "n",
      axes = F, col = cols.unc, breaks = breaks, add=T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("a"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("0 - 0.3 m"), cex = 1.3, font=1)


plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$seL2,  useRaster=TRUE, maxcell=maxpixels,bty = "n", axef = F,
      ylim = c(-90, 90), col = cols.unc, legend = F, breaks = breaks, add=T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("c"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("0.3 - 1 m"), cex = 1.3, font=1)


plot(ocean,col = "#DCDCDD", border = NA)
plot(bbox, add = T)
image(maps.robin.rr$seL3,  useRaster=TRUE, maxcell=maxpixels,bty = "n", axef = F,
      ylim = c(-90, 90), col = cols.unc, legend = F, breaks = breaks, add = T)
plot(wmap_wgs, add = T, lwd = 0.5)
text(x=-16000000, y=8000000, labels = c("e"), cex = 1.5, font=2)
text(x=-11000000, y=-5000000, labels = c("1 - 2 m"), cex = 1.3, font=1)


par(pin=c(4, 0.55))
plot(c(1,0), c(0,1), ylim = c(0,10), xlim = c(0,13),
     type = "n", ylab = "",
     xaxt="n",yaxt="n",bty="n",
     xlab= ("Uncertainity in 2C°warming effect estimate (%)"))

leg.txt1 = c(seq(0,20,2),">20")
axis(1,at=seq(1, 12,1), labels = leg.txt1,  pos=8)
dw=1
sag <-seq(1,11, by=1)
rect(0+sag, 8,dw+sag, 13,  col = cols.unc, border = NA)
text(x=9, y = 1, adj=1, labels = c("Prediction Uncertainity (%)"), xpd = T)


par(pin=c(3, 1.8))
boxplot(seL1~biomess, data = XY.pred.unc.sample, ylab=c("Uncertainity (%)"),
        xaxt = "n", col = fun_cols.eco(XX=seL1.mean, breaks, cols=cols.unc, ga=2),
        ylim = c(3, 55))
points(x=c(1:9), y=seL1.mean, col = "red", pch = 16, cex = 1.7)
abline(a=0, b=0, lty=2, col = "grey")
# text(x=0.5, y=55, labels = c("b"), cex = 1.5, font=2)
text(x=1.2, y=55, labels = c("0 - 0.3 m"), cex = 1.3, font=1)

boxplot(seL2~biomess, data = XY.pred.unc.sample, ylab=c("Uncertainity (%)"),
        xaxt="n", col = fun_cols.eco(XX=seL2.mean, breaks, cols=cols.unc, ga=2),
        ylim = c(3, 55))
points(x=c(1:9), y=seL2.mean, col = "red", pch = 16, cex = 1.7)
abline(a=0, b=0, lty=2, col = "grey")
# text(x=0.5, y=55, labels = c("d"), cex = 1.5, font=2)
text(x=1.2, y=55, labels = c("0.3 - 1 m"), cex = 1.3, font=1)

boxplot(seL3~biomess, data = XY.pred.unc.sample, ylab=c("Uncertainity (%)"),
        xaxt="n", col = fun_cols.eco(XX=seL3.mean, breaks, cols=cols.unc, ga=2),
        ylim = c(3, 55))
points(x=c(1:9), y=seL2.mean, col = "red", pch = 16, cex = 1.7)
abline(a=0, b=0, lty=2, col = "grey")
# text(x=0.5, y=55, labels = c("f"), cex = 1.5, font = 2)
text(x=1.2, y=55, labels = c("1 - 2 m"), cex = 1.3, font=1)

par(pin=c(3, 0.6))
plot(c(1,0), c(0,1), ylim = c(0,3), xlim = c(0.5,9.5),
     type = "n", ylab = "",
     xaxt="n",yaxt="n",bty="n",
     xlab= ("Ecosystem type"))

axis(1,at=seq(1, 9,1), labels = c("","","","","","","","",""),  las = 2, pos=4.5)
text(x=1:9, y = 4, srt = 45, adj=1, labels = eco.type, xpd = T)
##size:8.02*7.01
# dev.off()











