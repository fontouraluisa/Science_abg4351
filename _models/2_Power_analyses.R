# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: September 2020
# updated: June 2021
# outputs: Power of the coefficients


# packages
library(dotwhisker)
library(broom)
library(dplyr)
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)
library(brms)

# load original data
  # use the 0.Load_DATA_biomass.R script to load data

# load global Model 1 (ACTIVE 1) and ACTIVE 1 for functional groups
Global.active1.biom <- readRDS("/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/GLOBAL/ACTIVE 1/BIOMASS/MOD_BIOM_1_run_ACTIVE.1_V2_May.Rds")
TRANSIENT.active1.biom <- readRDS("/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/ACTIVE_1/Transient/Biomass/MOD_BIOM_1_run_TRANSIENT_V2_May.Rds")
RESID.active1.biom <- readRDS("/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/ACTIVE_1/Resident/Biomass/MOD_BIOM_1_run_RESID_V2_May.Rds")
PARENTAL.active1.biom <- readRDS("/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/ACTIVE_1/Parental/Biomass/MOD_BIOM_1_run_PARENTAL_V2_May.Rds")
CRYPTO.active1.biom <- readRDS("/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/ACTIVE_1/Cryptic/Biomass/MOD_BIOM_1_run_CRYPTIC_V2_May.Rds")

# newdata 
post.pred.GLOBAL.biomass.mod1.ACT1 <- posterior_predict(Global.active1.biom) # 8000*263 (biomass)
post.pred.TRANSIENT <- posterior_predict(TRANSIENT.active1.biom) # 1800*256 (biomass)
post.pred.RESID <- posterior_predict(RESID.active1.biom ) # 8000*263 (biomass)
post.pred.PARENTAL <- posterior_predict(PARENTAL.active1.biom) # 8000*262 (biomass)
post.pred.CRYPTIC <- posterior_predict(CRYPTO.active1.biom) # 8000*264 (biomass)

# Number of posterior to draw from each matrix = 20% = 1600
#n <- 8000*20/100 ; n # ideally 1600 but could be very time consuming
# try first with 200
n <- 200
post.draws <- sample(1:8000,n) # sample randomly 1600 draws from the 8000 from each model - 1st test with 200 and check computation time

# post.pred for biomass
Biom.GLOBAL.ACT1 <- post.pred.GLOBAL.biomass.mod1.ACT1
Biom.TRANSIENT <- post.pred.TRANSIENT
Biom.RESID <- post.pred.RESID
Biom.PARENTAL <- post.pred.PARENTAL
Biom.CRYPTIC <- post.pred.CRYPTIC 

# newdata
rm(newdata.GLOBAL.Biom.ACT1)
newdata.GLOBAL.Biom.ACT1 <-Global.active1.biom$data # rm NA
rm(newdata.TRANSIENT)
newdata.TRANSIENT <-TRANSIENT.active1.biom$data # rm NA
rm(newdata.RESID)
newdata.RESID <- RESID.active1.biom$data # rm NA
rm(newdata.PARENTAL)
newdata.PARENTAL <- PARENTAL.active1.biom$data # rm NA
rm(newdata.CRYPTIC)
newdata.CRYPTIC <- CRYPTIC.active1.biom$data # rm NA

# updates the model with each raw of post.pred and extract the coefficient values
  # create list to store the coefficient for each of the 30 model
rm(coef.biom.GLOBAL)
coef.biom.GLOBAL <- list()
rm(coef.biom.TRANSIENT)
coef.biom.TRANSIENT <- list()
rm(coef.biom.RESID)
coef.biom.RESID <- list()
rm(coef.biom.PARENTAL)
coef.biom.PARENTAL <- list()
rm(coef.biom.CRYPTIC)
coef.biom.CRYPTIC <- list()
  
for (i in 1:length(post.draws)){
  paste("i=",i,sep="_")
  
  Todraw <- post.draws[i]
  
  # GLOBAL
  newdata.GLOBAL.Biom.ACT1$log_biomassarea <- Biom.GLOBAL.ACT1[Todraw,]
  temp.fit.GLOBAL.Biom.ACT1 <- update(Global.active1.biom, newdata=newdata.GLOBAL.Biom.ACT1)
  a.GLOBAL <- mcmc_intervals(temp.fit.GLOBAL.Biom.ACT1)
  biom.var.GLOBAL <- as.data.frame(a.GLOBAL$data)[grep("b_",as.data.frame(a.GLOBAL$data)[,"parameter"]),"parameter"]
  coef.biom.GLOBAL[[i]]  <- mcmc_intervals_data(temp.fit.GLOBAL.Biom.ACT1,pars = as.character(biom.var.GLOBAL))
  
  # TRANSIENT
  newdata.TRANSIENT$log_biomassarea <- Biom.TRANSIENT[Todraw,]
  temp.fit.TRANSIENT <- update(TRANSIENT.active1.biom, newdata=newdata.TRANSIENT)
    #extract coefficients for each TRANSIENT model
  a.TRANSIENT <- mcmc_intervals(temp.fit.TRANSIENT)
  biom.var.TRANSIENT <- as.data.frame(a.TRANSIENT$data)[grep("b_",as.data.frame(a.TRANSIENT$data)[,"parameter"]),"parameter"]
  coef.biom.TRANSIENT[[i]]  <- mcmc_intervals_data(temp.fit.TRANSIENT,pars = as.character(biom.var.TRANSIENT))
  
  # RESIDENT
  newdata.RESID$log_biomassarea <- Biom.RESID[Todraw,]
  temp.fit.RESID <- update(RESID.active1.biom, newdata=newdata.RESID)
  # extract coefficients for each TRANSIENT model
  a.RESID <- mcmc_intervals(temp.fit.RESID)
  biom.var.RESID <- as.data.frame(a.RESID$data)[grep("b_",as.data.frame(a.RESID$data)[,"parameter"]),"parameter"]
  coef.biom.RESID[[i]]  <- mcmc_intervals_data(temp.fit.RESID,pars = as.character(biom.var.RESID))
  
  # PARENTAL
  newdata.PARENTAL$log_biomassarea <- Biom.PARENTAL[Todraw,]
  temp.fit.PARENTAL <- update(all_fit_brms.tot.PARENTAL.intr.extr, newdata=newdata.PARENTAL)
  # extract coefficients for each TRANSIENT model
  a.PARENTAL <- mcmc_intervals(temp.fit.PARENTAL)
  biom.var.PARENTAL <- as.data.frame(a.PARENTAL$data)[grep("b_",as.data.frame(a.PARENTAL$data)[,"parameter"]),"parameter"]
  coef.biom.PARENTAL[[i]]  <- mcmc_intervals_data(temp.fit.PARENTAL,pars = as.character(biom.var.PARENTAL))
  
  # CRYPTIC
  newdata.CRYPTIC$log_biomassarea <- Biom.CRYPTIC[Todraw,]
  temp.fit.CRYPTIC <- update(CRYPTO.active1.biom, newdata=newdata.CRYPTIC)
  # extract coefficients for each TRANSIENT model
  a.CRYPTIC <- mcmc_intervals(temp.fit.CRYPTIC)
  biom.var.CRYPTIC <- as.data.frame(a.CRYPTIC$data)[grep("b_",as.data.frame(a.CRYPTIC$data)[,"parameter"]),"parameter"]
  
  coef.biom.CRYPTIC[[i]]  <- mcmc_intervals_data(temp.fit.CRYPTIC,pars = as.character(biom.var.CRYPTIC))
  
  # clean temporary object at the end of the loop
  rm(temp.fit.TRANSIENT); rm(temp.fit.RESID);rm(temp.fit.PARENTAL);rm(temp.fit.CRYPTIC)
  rm(a.TRANSIENT); rm(a.RESID); rm(a.PARENTAL); rm(a.CRYPTIC)
  rm(rich.var.TRANSIENT);  rm(rich.var.RESID);  rm(rich.var.PARENTAL);  rm(rich.var.CRYPTIC)
  rm(biom.var.TRANSIENT);  rm(biom.var.RESID);  rm(biom.var.PARENTAL);  rm(biom.var.CRYPTIC)

}

coef.biom.GLOBAL <- coef.biom.GLOBAL[-which(sapply(coef.biom.GLOBAL, is.null))]
coef.biom.TRANSIENT <- coef.biom.TRANSIENT[-which(sapply(coef.biom.TRANSIENT, is.null))]
coef.biom.RESID <- coef.biom.RESID[-which(sapply(coef.biom.RESID, is.null))]
coef.biom.PARENTAL <- coef.biom.PARENTAL[-which(sapply(coef.biom.PARENTAL, is.null))]
coef.biom.CRYPTIC <- coef.biom.CRYPTIC[-which(sapply(coef.biom.CRYPTIC, is.null))]

rm(temp.biom.GLOBAL)
temp.biom.GLOBAL <-  coef.biom.GLOBAL[[1]]$parameter
rm(temp.biom.TRANSIENT)
temp.biom.TRANSIENT <-  coef.biom.TRANSIENT[[1]]$parameter
rm(temp.biom.RESID)
temp.biom.RESID <-  coef.biom.RESID[[1]]$parameter
rm(temp.biom.PARENTAL)
temp.biom.PARENTAL <-  coef.biom.PARENTAL[[1]]$parameter
rm(temp.biom.CRYPTIC)
temp.biom.CRYPTIC <-  coef.biom.CRYPTIC[[1]]$parameter

for (j in 1:length(coef.biom.GLOBAL)){
  temp.biom.GLOBAL <- cbind(temp.biom.GLOBAL,as.vector(as.data.frame(coef.biom.GLOBAL[[j]][,"m"])))
  temp.biom.TRANSIENT <- cbind(temp.biom.TRANSIENT,as.vector(as.data.frame(coef.biom.TRANSIENT[[j]][,"m"])))
  temp.biom.PARENTAL <- cbind(temp.biom.PARENTAL,as.vector(as.data.frame(coef.biom.PARENTAL[[j]][,"m"])))
  temp.biom.RESID <- cbind(temp.biom.RESID,as.vector(as.data.frame(coef.biom.RESID[[j]][,"m"])))
  temp.biom.CRYPTIC <- cbind(temp.biom.CRYPTIC,as.vector(as.data.frame(coef.biom.CRYPTIC[[j]][,"m"])))
  
}

# GLOBAL
  # biomass
  rm(df.biom.GLOBAL)
  df.biom.GLOBAL <- temp.biom.GLOBAL
  rownames(df.biom.GLOBAL) <- temp.biom.GLOBAL[,1]; df.biom.GLOBAL$temp.biom <- NULL
  saveRDS(df.biom.GLOBAL,"/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/Power analysis/Biomass.power.pred.GLOBAL.rds")

# TRANSIENT
 # biomass
 rm(df.biom.TRANSIENT)
 df.biom.TRANSIENT <- temp.biom.TRANSIENT 
 rownames(df.biom.TRANSIENT) <- temp.biom.TRANSIENT[,1]; df.biom.TRANSIENT$temp.biom <- NULL
 saveRDS(df.biom.TRANSIENT,"/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/Power analysis/Biomass.power.pred.TRANSIENT.rds")

# RESIDENT
 # biomass
 rm(df.biom.RESIDENT)
 df.biom.RESIDENT <- temp.biom.RESID
 rownames(df.biom.RESIDENT) <- temp.biom.RESID[,1]; df.biom.RESIDENT$temp.biom <- NULL
 saveRDS(df.biom.RESIDENT,"/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/Power analysis/Biomass.power.pred.RESIDENT.rds")
 
# PARENTAL
 # biomass
 rm(df.biom.PARENTAL)
 df.biom.PARENTAL <- temp.biom.PARENTAL 
 rownames(df.biom.PARENTAL) <- temp.biom.PARENTAL[,1]; df.biom.PARENTAL$temp.biom <- NULL
 saveRDS(df.biom.PARENTAL,"/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/Power analysis/Biomass.power.pred.PARENTAL.rds")
 
# CRYPTIC
 # biomass
 rm(df.biom.CRYPTIC)
 df.biom.CRYPTIC <- temp.biom.CRYPTIC 
 rownames(df.biom.CRYPTIC) <- temp.biom.CRYPTIC[,1]; df.biom.CRYPTIC$temp.biom <- NULL
 saveRDS(df.biom.CRYPTIC,"/Users/stephdagata/Dropbox/Science_MS_Connectivity_Review/Git_and_SI/GITHUB/Models output/Power analysis/Biomass.power.pred.CRYPTIC.rds")
 
 
 #### END #####
 
 

  
 
 
 
 
 

