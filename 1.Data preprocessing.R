###############################################################################
##### Depth-dependent soil carbon losses under warmer climate revealed ########
###################### by assessing global soil profiles ######################
###############################################################################

#### Enconding with UTF-8

## Imputing the missing bulk density and gravl content in WoSIS 

library(dplyr)
library(matrixStats)
library(Xgboost)

# Load WoSIS data
attributes = read.csv("XXXX\\wosis_201909_attributes.tsv",header = T, sep = "\t")
profiles = read.csv("XXXX\\wosis_201909_profiles.tsv",header = T, sep = "\t")
profiles <- profiles[,c("profile_id","latitude","longitude","dsds")]
chemical = read.csv("XXXX\\wosis_201909_layers_chemical.tsv",header = T, sep = "\t")
physcial = read.csv("XXX\\wosis_201909_layers_physical.tsv",header = T, sep = "\t")

variables = c("profile_id","upper_depth","lower_depth","dsds",paste(tolower(attributes$code),"_value_avg",sep=""))
mer <- merge(x = chemical, y = physcial, by = c("profile_id","upper_depth","lower_depth"), all = TRUE)
mer <- mer[,na.omit(match(variables,colnames(mer)))]
mer_profiles <- merge(x = profiles, y = mer, by = c("profile_id"), all = TRUE)

variables = c("profile_id","latitude","longitude","dsds","upper_depth","lower_depth" ,tolower(attributes$code))
variables = variables[-c(17,21:23)]    #鍒犲�?"cfao"   "cstx"   "cwrb"   "dsds"
colnames(mer_profiles) = variables
mer_profiles = mer_profiles[,c(1:6, 26, 10, 18, 7:9,11:17,19:25,27:54)]
summary(mer_profiles[,1:9])
mer_sdata <- mer_profiles


### Training XG-Boost model of BD and G
df <- mer_sdata[,c(5,6,8:54)]
df$depth <- (df$upper_depth + df$lower_depth)/2
df <- df[,c(-1,-2)]
df <- df[-which(df$depth<=0),]

df_bdfiod <- df[-which(apply(df[,-1],1,function(x) all(is.na(x)))),]
df_cfvo <- df[-which(apply(df[,-2],1,function(x) all(is.na(x)))),]


fold10 <- function(dataframe, y, n, nrounds, params){
  r2_total <- list()
  for (i in 1:n){
    set.seed(i+227)
    samp <- sample(1:length(dataframe[,1]),length(dataframe[,1])*0.8)
    train <- dataframe[samp,]
    test <- dataframe[-samp,]
    matrix_train = DP_xgb(train, y)
    model <- xgboost(data = matrix_train, nrounds = nrounds, silent = 1, params = params)
    t=test[, y]
    p=pred_xgboost(test[,-y], 1:length(test[,-y]), model)
    r2 = summary(lm(t~p))$adj.r.squared
    # print(r2)
    r2_total = c(r2_total, r2)
  }
  r2_total = unlist(r2_total)
  # print(mean(r2_total))
  # print(CV(r2_total))
  return(r2_total)
}

uncertainty <- function(train, train_y, test, model_nrounds, model_params, n_bootstrap, n_repeat){
  prediction = data.frame(num = 1:length(test[,1]))
  for (i in 1:n_repeat){
    num = sample(1:length(train[,1]), n_bootstrap, T)
    boot_train = train[num,]
    matrix_train = DP_xgb(boot_train, train_y)
    model <- xgboost(data = matrix_train, nrounds = model_nrounds, params = model_params)
    p = pred_xgboost(test, 1:length(test[1,]), model)
    prediction = cbind(prediction, p)
  }
  library(dplyr)
  library(matrixStats)
  prediction$cv = rowSds(as.matrix(prediction[2:n_repeat+1]))/rowMeans(prediction[2:n_repeat+1])
  prediction$mean = rowMeans(prediction[2:n_repeat+1])
  prediction = prediction[,-1]
  return(prediction)
}

# Missing BD prediction
nrounds.bd = 1400
params = list(eta = 0.18, max_depth = 8, min_child_weight = 1, gamma = 0)
num <- which(is.na(df_bdfiod$bdfiod) == T)
train.bd <- df_bdfiod[-num,]
unknown <- df_bdfiod[num,]
bd_10fold = fold10(train, 1, 10, 1400, params)
bd_uncer = uncertainty(train=train.bd, train_y, test, model_nrounds=nrounds.bd, model_params=params, n_bootstrap=50, n_repeat=50)
df.uncertainty_bdfi = data.frame(bd_10fold, CV.bd=bd_uncer)

# Missing G prediction
nrounds.g = 320
parmas = list(eta = 0.2, max_depth = 8, min_child_weight = 3, gamma = 0)
num <- which(is.na(df_cfvo$cfvo) == T)
train.cf <- df_cfvo[-num,]
unknown <- df_cfvo[num,]
cf_10fold = fold10(train, 2, 10, 320, params)
cf_uncer = uncertainty(train=train.cf, train_y, test, model_nrounds=nrounds.bd, model_params=params, n_bootstrap=50, n_repeat=50)
df.uncertainty_cfvo = data.frame(cf_10fold, CV.cf=cf_uncer)

# Selecting CV<1
df.uncertainty_bdfi[which(abs(df.uncertainty_bdfi$CV.bd) > 1), 52] = NA
df.uncertainty_bdfi[which(df.uncertainty_bdfi$mean < 0), 52] = NA
summary(df.uncertainty_bdfi$mean)

df.uncertainty_cfvo[which(abs(df.uncertainty_cfvo$CV.cf) > 1), 52] = NA
df.uncertainty_cfvo[which(df.uncertainty_cfvo$mean < 0), 52] = NA
df.uncertainty_cfvo[which(df.uncertainty_cfvo$mean > 100), 52] = NA












