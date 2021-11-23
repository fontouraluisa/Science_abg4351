require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)
### marginal plots
library(ggeffects)
library(here)
library(viridis)


#Load data ----
#Active Global
ACTIVE.1.sub.V2 <- readRDS(here("_data","Datasets","GLOBAL","ACTIVE.1.sub.V2_May.Rds"))
ACTIVE.1.sub.V2.std <- readRDS(here("_data","Datasets","GLOBAL","ACTIVE.1.sub.V2.std_May.Rds"))

MOD_BIOM_1_run_ACTIVE.1 <- readRDS(here("Models output","MOD_BIOM_1_run_ACTIVE.1_V2_May.Rds"))
MOD_S_1_run_ACTIVE.1 <- readRDS(here("Models output","MOD_S_1_run_ACTIVE.1_V2_May.Rds"))

#cryptic 5
MOD_BIOM_1_run_CRYPTIC <- readRDS(here("Models output","ACTIVE_1","CRYPTIC","Biomass","MOD_BIOM_1_run_CRYPTIC_V2_May.Rds"))
MOD_S_1_run_CRYPTIC <- readRDS(here("Models output","ACTIVE_1","CRYPTIC","Richness","MOD_S_1_run_CRYPTIC_V2_May.Rds"))

#Parental 5
MOD_BIOM_1_run_PARENTAL <- readRDS(here("Models output","ACTIVE_1","PARENTAL","Biomass","MOD_BIOM_1_run_PARENTAL_V2_May.Rds"))
MOD_S_1_run_PARENTAL <- readRDS(here("Models output","ACTIVE_1","PARENTAL","Richness","MOD_S_1_run_PARENTAL_V2_May.Rds"))

#Resident 15
MOD_BIOM_1_run_RESID <- readRDS(here("Models output","ACTIVE_1","Resident","Biomass","MOD_BIOM_1_run_RESID_V2_May.Rds"))
MOD_S_1_run_RESID <- readRDS(here("Models output","ACTIVE_1","Resident","Richness","MOD_S_1_run_RESID_V2_May.Rds"))

#Transient 15
MOD_BIOM_1_run_TRANSIENT <- readRDS(here("Models output","ACTIVE_1","Transient","Biomass","MOD_BIOM_1_run_TRANSIENT_V2_May.Rds"))
MOD_S_1_run_TRANSIENT <- readRDS(here("Models output","ACTIVE_1","Transient","Richness","MOD_S_1_run_TRANSIENT_V2_May.Rds"))

#Figure 1B and Extended (FIG S1) Fig 1 Biomass and Species richness----
# TRANSIENT TOT
rm(a)
a <- mcmc_intervals(MOD_BIOM_1_run_ACTIVE.1)
# no classification
rm(biom.var)
biom.var <- a$data[1:10,]$parameter

Biomass_inflow_GLOBAL_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_ACTIVE.1,pars = as.character(biom.var))
Biomass_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_TRANSIENT,pars = as.character(biom.var))
Biomass_inflow_CRYPTO_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_CRYPTIC,pars = as.character(biom.var))
Biomass_inflow_RESID_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_RESID,pars = as.character(biom.var))
Biomass_inflow_PARENTAL_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_PARENTAL,pars = as.character(biom.var))


# Classification
INT <- "Intercept"
CON <- c("log_SelfR","Netflow","log_InflowNei")
HUMAN <- c("log_grav_total","_ClassClosed","_ClassRestricted")
ENV <- c("Richness","temp","prod.annual")

PARAM <- c("Intercept","Richness","Temperature","Productivity",
           "Tot.Gravity","No-Take","Restricted gears",
           "Netflow","Local Recruit.","Exogenous Inflow")

CAT <- c("Intercept","Human/Env.","Human/Env.","Human/Env.","Human/Env.","Human/Env.",
         "Human/Env.",
         "Connectivity",
         "Connectivity",
         "Connectivity")


DesiredOrder <- c("Intercept",
                  "Richness",
                  "Temperature",
                  "Productivity",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Netflow",
                  "Local Recruit.",
                  "Exogenous Inflow")

#GLOBAL
rm(m2_df_global)
m2_df_global <- Biomass_inflow_GLOBAL_simp_data %>%  mutate(FE = "GLOBAL") %>%  mutate(model = "Biomass")
m2_df_global <- cbind(PARAM,m2_df_global,CAT)
m2_df_global$PARAM <- factor(m2_df_global$PARAM, levels = rev(DesiredOrder))
m2_df_global$CAT <- factor(m2_df_global$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#TRANSIENT
rm(m2_df_transient)
m2_df_transient <- Biomass_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Biomass")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#CRYPTIC
rm(m2_df_crypto)
m2_df_crypto <- Biomass_inflow_CRYPTO_simp_data %>%  mutate(FE = "CRYPTO") %>%  mutate(model = "Biomass")
m2_df_crypto <- cbind(PARAM,m2_df_crypto,CAT)
m2_df_crypto$PARAM <- factor(m2_df_crypto$PARAM, levels = rev(DesiredOrder))
m2_df_crypto$CAT <- factor(m2_df_crypto$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#PARENTAL
rm(m2_df_PARENTAL)
m2_df_PARENTAL <- Biomass_inflow_PARENTAL_simp_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Biomass")
m2_df_PARENTAL <- cbind(PARAM,m2_df_PARENTAL,CAT)
m2_df_PARENTAL$PARAM <- factor(m2_df_PARENTAL$PARAM, levels = rev(DesiredOrder))
m2_df_PARENTAL$CAT <- factor(m2_df_PARENTAL$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#RESID
rm(m2_df_RESID)
m2_df_RESID <- Biomass_inflow_RESID_simp_data %>%  mutate(FE = "RESID") %>%  mutate(model = "Biomass")
m2_df_RESID <- cbind(PARAM,m2_df_RESID,CAT)
m2_df_RESID$PARAM <- factor(m2_df_RESID$PARAM, levels = rev(DesiredOrder))
m2_df_RESID$CAT <- factor(m2_df_RESID$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)


# biomass data.frame
rm(four_FE_Biomass)
four_FE_Biomass <- rbind(m2_df_global,m2_df_transient,m2_df_crypto,m2_df_RESID,m2_df_PARENTAL)
four_FE_Biomass$FE <- factor(four_FE_Biomass$FE, levels=c("CRYPTO","PARENTAL","RESID","TRANSIENT","GLOBAL"))
four_FE_Biomass$CAT <- as.factor(four_FE_Biomass$CAT)
four_FE_Biomass$model <- as.factor(four_FE_Biomass$model)
summary(four_FE_Biomass)

# plot biomass
rm(Biomass.FE)
Biomass.FE <- ggplot(four_FE_Biomass %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE,size=6), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21,fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +  scale_color_viridis(discrete=TRUE) +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Biomass") 
Biomass.FE  # The trick to these is position_dodge()

#Richness
rm(d)
d <- mcmc_intervals(MOD_S_1_run_ACTIVE.1)
# no classification
rm(rich.var)
rich.var <- d$data[1:10,]$parameter

Richness_inflow_GLOBAL_simp_data <- mcmc_intervals_data(MOD_S_1_run_ACTIVE.1,pars = as.character(rich.var))
Richness_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(MOD_S_1_run_TRANSIENT,pars = as.character(rich.var))
Richness_inflow_CRYPTO_simp_data <- mcmc_intervals_data(MOD_S_1_run_CRYPTIC,pars = as.character(rich.var))
Richness_inflow_RESID_simp_data <- mcmc_intervals_data(MOD_S_1_run_RESID,pars = as.character(rich.var))
Richness_inflow_PARENTAL_simp_data <- mcmc_intervals_data(MOD_S_1_run_PARENTAL,pars = as.character(rich.var))


# Classification
INT <- "Intercept"
CON <- c("log_btwdegree","Netflow","log_SelfR","log_Indegree_Neigh")
HUMAN <- c("log_grav_total","_ClassClosed","_ClassRestricted")
ENV <- c("temp","prod.annual")

PARAM <- c("Intercept","Temperature","Productivity",
           "Betweenness","Netflow","Local Recruit.","Exogenous Indegree",
           "Tot.Gravity","No-Take","Restricted gears")

CAT <- c("Intercept","Human/Env.","Human/Env.",
         "Connectivity",
         "Connectivity",
         "Connectivity",
         "Connectivity","Human/Env.","Human/Env.",
         "Human/Env.")

DesiredOrder <- c("Intercept",
                  "Temperature",
                  "Productivity",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Betweenness",
                  "Netflow",
                  "Local Recruit.",
                  "Exogenous Indegree")
#GLOBAL
rm(m2_df_global)
m2_df_global <- Richness_inflow_GLOBAL_simp_data %>%  mutate(FE = "GLOBAL") %>%  mutate(model = "Richness")
m2_df_global <- cbind(PARAM,m2_df_global,CAT)
m2_df_global$PARAM <- factor(m2_df_global$PARAM, levels = rev(DesiredOrder))
m2_df_global$CAT <- factor(m2_df_global$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#TRANSIENT
rm(m2_df_transient)
m2_df_transient <- Richness_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Richness")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#CRYPTIC
rm(m2_df_crypto)
m2_df_crypto <- Richness_inflow_CRYPTO_simp_data %>%  mutate(FE = "CRYPTO") %>%  mutate(model = "Richness")
m2_df_crypto <- cbind(PARAM,m2_df_crypto,CAT)
m2_df_crypto$PARAM <- factor(m2_df_crypto$PARAM, levels = rev(DesiredOrder))
m2_df_crypto$CAT <- factor(m2_df_crypto$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#PARENTAL
rm(m2_df_PARENTAL)
m2_df_PARENTAL <- Richness_inflow_PARENTAL_simp_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Richness")
m2_df_PARENTAL <- cbind(PARAM,m2_df_PARENTAL,CAT)
m2_df_PARENTAL$PARAM <- factor(m2_df_PARENTAL$PARAM, levels = rev(DesiredOrder))
m2_df_PARENTAL$CAT <- factor(m2_df_PARENTAL$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#RESID
rm(m2_df_RESID)
m2_df_RESID <- Richness_inflow_RESID_simp_data %>%  mutate(FE = "RESID") %>%  mutate(model = "Richness")
m2_df_RESID <- cbind(PARAM,m2_df_RESID,CAT)
m2_df_RESID$PARAM <- factor(m2_df_RESID$PARAM, levels = rev(DesiredOrder))
m2_df_RESID$CAT <- factor(m2_df_RESID$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)


# biomass data.frame
rm(four_FE_Richness)
four_FE_Richness <- rbind(m2_df_global,m2_df_transient,m2_df_crypto,m2_df_RESID,m2_df_PARENTAL)
four_FE_Richness$FE <- factor(four_FE_Richness$FE, levels=c("CRYPTO","PARENTAL","RESID","TRANSIENT","GLOBAL"))
four_FE_Richness$CAT <- as.factor(four_FE_Richness$CAT)
four_FE_Richness$model <- as.factor(four_FE_Richness$model)
summary(four_FE_Richness)

# plot richness
rm(Richness.FE)
Richness.FE <- ggplot(four_FE_Richness %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE,size=6), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21,fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +  scale_color_viridis(discrete=TRUE) +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Richness") 
Richness.FE  # The trick to these is position_dodge()

BAYE_coef_tot_rich_biom_FULL <- ggarrange(Biomass.FE,Richness.FE,
                                     ncol=2,nrow=1,labels = c("A","B"),align="hv",common.legend = T,legend="bottom")


#Fig.1D&E Marginal plots ----

#Netflow and indegree of neighbors #back transform netflow, species richness and indegree of neighbors
conditional_effects(MOD_S_1_run_ACTIVE.1, "log_Indegree_Neigh")
##back transform richness
meanRich<-mean(ACTIVE.1.sub.V2$Richness,na.rm=T)
sdRich<- 1*sd(ACTIVE.1.sub.V2$Richness,na.rm=T)
FUNINVLU = function(x){(x * sdRich) + meanRich}
##back transform log indegree
meanINN<-mean(ACTIVE.1.sub.V2$log_Indegree_Neigh,na.rm=T)
sdINN<- 1*sd(ACTIVE.1.sub.V2$log_Indegree_Neigh,na.rm=T)
FUNINVIN = function(x){(x * sdINN) + meanINN}
#marginal plot (richness vs extrinsic indegree )
margCry<-plot(conditional_effects(MOD_S_1_run_ACTIVE.1, "log_Indegree_Neigh"))
marRich_Cor<-margCry$log_Indegree_Neigh
try2<-margCry$log_Indegree_Neigh$plot_env$plots$log_Indegree_Neigh$plot_env$x$log_Indegree_Neigh
rm(margPlot)
margPlot<-try2[,c("log_Indegree_Neigh","estimate__","lower__","upper__")]
margPlot$rwRich<-FUNINVLU(margPlot$estimate__) #unscale richness
margPlot$rwRichlow<-FUNINVLU(margPlot$lower__)
margPlot$rwRichupp<-FUNINVLU(margPlot$upper__)
margPlot$log_Indegree_NeighBT <- FUNINVIN(margPlot$log_Indegree_Neigh)
###

PartC<-ggplot(margPlot, aes(y = rwRich, x = log_Indegree_NeighBT)) +
  geom_point(size=0.01) + geom_point(data=ACTIVE.1.sub.V2, mapping=aes(y=Richness,x=log_Indegree_Neigh),shape=21 ,size=2,alpha=0.35, color="gray",fill="gray25") +
  geom_ribbon( aes(ymin = rwRichlow, ymax = rwRichupp), alpha = .3) +
  geom_line( aes(y = rwRich), size = 1.5, color="black") + theme_classic() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))
PartC


##Netflow
##back transform Netflow
meanINfN<-mean(ACTIVE.1.sub.V2$Netflow,na.rm=T)
sdINfN<- 1*sd(ACTIVE.1.sub.V2$Netflow,na.rm=T)
FUNINVLUIF = function(x){(x * sdINfN) + meanINfN}
#marginal plot (log biomass vs Netflow - back transform Netflow )
rm(margCry)
rm(try2)
margCry<-plot(conditional_effects(MOD_BIOM_1_run_ACTIVE.1, "Netflow"))
try2<-margCry$Netflow$plot_env$plots$Netflow$plot_env$x$Netflow
rm(margPlot)
margPlot<-try2[,c("Netflow","estimate__","lower__","upper__")]
margPlot$NetflowBT <- FUNINVLUIF(margPlot$Netflow)
summary(margPlot)

rm(PartD)
PartD<-ggplot(margPlot, aes(y = estimate__, x = NetflowBT)) +
  geom_point(size=0.01) +  geom_point(data=ACTIVE.1.sub.V2, mapping=aes(y=log_biomassarea, x=Netflow), shape=21 ,size=2,alpha=0.35, color="gray", fill="gray25") +
  geom_ribbon( aes(ymin = lower__, ymax = upper__), alpha = .3) +
  geom_line( aes(y = lestimate__), size = 1.3, color="black") + theme_classic() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
PartD

Netflow_Indegree <- ggarrange(PartD,PartC,
                              ncol=2,nrow=1,labels = c("A","B"),align="hv")



#Fig. S1 ----
# TRANSIENT TOT
rm(a)
a <- mcmc_intervals(MOD_BIOM_1_run_ACTIVE.1)
# no classification
rm(biom.var)
biom.var <- a$data[1:17,]$parameter

Biomass_inflow_GLOBAL_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_ACTIVE.1,pars = as.character(biom.var))
Biomass_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_TRANSIENT,pars = as.character(biom.var))
Biomass_inflow_CRYPTO_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_CRYPTIC,pars = as.character(biom.var))
Biomass_inflow_RESID_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_RESID,pars = as.character(biom.var))
Biomass_inflow_PARENTAL_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_PARENTAL,pars = as.character(biom.var))


# Classification
INT <- "Intercept"
CON <- c("log_SelfR","Netflow","log_InflowNei")
HUMAN <- c("log_grav_total","_ClassClosed","_ClassRestricted")
ENV <- c("Richness","temp","prod.annual")
INT <- c("log_grav_total:ClassClosed","log_grav_total:ClassRestricted","log_grav_total:Netflow",
         "ClassClosed:Netflow","ClassRestricted:Netflow","log_grav_total:ClassClosed:Netflow",
         "log_grav_total:ClassRestricted:Netflow")

PARAM <- c("Intercept","Richness","Temperature","Productivity",
           "Tot.Gravity","No-Take","Restricted gears",
           "Netflow","Local Recruit.","Exogenous Inflow",
           "Tot.Gravity x No-Take","Tot.Gravity x Restricted gears",
           "Tot.Gravity x Netflow","No-Take x Netflow","Restricted gears x Netflow",
           "Tot.Gravity x No-Take x Netflow", "Tot.Gravity x Restricted gears x Netflow")

CAT <- c("Intercept","Human/Env.","Human/Env.","Human/Env.","Human/Env.","Human/Env.",
         "Human/Env.",
         "Connectivity",
         "Connectivity",
         "Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity",
         "Management x Human pressure x Connectivity")


DesiredOrder <- c("Intercept",
                  "Richness",
                  "Temperature",
                  "Productivity",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Netflow",
                  "Local Recruit.",
                  "Exogenous Inflow",
                  "Tot.Gravity x No-Take",
                  "Tot.Gravity x Restricted gears",
                  "Tot.Gravity x Netflow",
                  "No-Take x Netflow",
                  "Restricted gears x Netflow",
                  "Tot.Gravity x No-Take x Netflow", 
                  "Tot.Gravity x Restricted gears x Netflow")

#GLOBAL
rm(m2_df_global)
m2_df_global <- Biomass_inflow_GLOBAL_simp_data %>%  mutate(FE = "GLOBAL") %>%  mutate(model = "Biomass")
m2_df_global <- cbind(PARAM,m2_df_global,CAT)
m2_df_global$PARAM <- factor(m2_df_global$PARAM, levels = rev(DesiredOrder))
m2_df_global$CAT <- factor(m2_df_global$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#TRANSIENT
rm(m2_df_transient)
m2_df_transient <- Biomass_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Biomass")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#CRYPTIC
rm(m2_df_crypto)
m2_df_crypto <- Biomass_inflow_CRYPTO_simp_data %>%  mutate(FE = "CRYPTO") %>%  mutate(model = "Biomass")
m2_df_crypto <- cbind(PARAM,m2_df_crypto,CAT)
m2_df_crypto$PARAM <- factor(m2_df_crypto$PARAM, levels = rev(DesiredOrder))
m2_df_crypto$CAT <- factor(m2_df_crypto$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#PARENTAL
rm(m2_df_PARENTAL)
m2_df_PARENTAL <- Biomass_inflow_PARENTAL_simp_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Biomass")
m2_df_PARENTAL <- cbind(PARAM,m2_df_PARENTAL,CAT)
m2_df_PARENTAL$PARAM <- factor(m2_df_PARENTAL$PARAM, levels = rev(DesiredOrder))
m2_df_PARENTAL$CAT <- factor(m2_df_PARENTAL$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)

#RESID
rm(m2_df_RESID)
m2_df_RESID <- Biomass_inflow_RESID_simp_data %>%  mutate(FE = "RESID") %>%  mutate(model = "Biomass")
m2_df_RESID <- cbind(PARAM,m2_df_RESID,CAT)
m2_df_RESID$PARAM <- factor(m2_df_RESID$PARAM, levels = rev(DesiredOrder))
m2_df_RESID$CAT <- factor(m2_df_RESID$CAT , levels = c("Intercept","Human/Env.","Connectivity","Management x Human pressure x Connectivity"),ordered = TRUE)


# biomass data.frame
rm(four_FE_Biomass)
four_FE_Biomass <- rbind(m2_df_global,m2_df_transient,m2_df_crypto,m2_df_RESID,m2_df_PARENTAL)
four_FE_Biomass$FE <- factor(four_FE_Biomass$FE, levels=c("CRYPTO","TRANSIENT","RESID","PARENTAL","GLOBAL"))
four_FE_Biomass$CAT <- as.factor(four_FE_Biomass$CAT)
four_FE_Biomass$model <- as.factor(four_FE_Biomass$model)
summary(four_FE_Biomass)

# plot biomass
rm(Biomass.FE)
Biomass.FE <- ggplot(four_FE_Biomass %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE,size=6), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21,fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +  scale_color_viridis(discrete=TRUE) +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Biomass") 
Biomass.FE  # The trick to these is position_dodge()

#Fig 2.  3-way interaction netflow x gravity x management----
# class  
rm(glob.model)
#Global Active - total connectivity shaped by the main four fish groups

glob.model<-MOD_BIOM_1_run_ACTIVE.1
conditions <- make_conditions(glob.model, "Class")

rm(global_pred_class)
global_pred_class <- conditional_effects(glob.model, "log_grav_total:Netflow", conditions = conditions,rug=T,points=T, re_formula = NULL)
rm(df)
df <- as.data.frame(global_pred_class$`log_grav_total:Netflow`)

rm(pred_dat)
Netflow <- c(min(ACTIVE.1.sub.V2.std$Netflow),mean(ACTIVE.1.sub.V2.std$Netflow),max(ACTIVE.1.sub.V2.std$Netflow))
Class <- df$Class %>% unique() %>% as.character()
log_grav_total <- seq(min(df$log_grav_total), max(df$log_grav_total), length.out =  100)
Richness <- df$Richness %>% mean()
temp <- df$temp %>% mean()
prod.annual <- df$prod.annual %>% mean()
log_SelfR <- df$log_SelfR %>% mean()
log_InflowNei <- df$log_InflowNei %>% mean()
region <- df$region %>% unique()

pred_dat <- expand.grid(Netflow, Class, log_grav_total, Richness, temp, prod.annual, log_SelfR, log_InflowNei, region)
names(pred_dat) <- c("Netflow", "Class", "log_grav_total", "Richness", "temp", "prod.annual", "log_SelfR", "log_InflowNei", "region")
y_hat <- predict(object = glob.model, newdata = pred_dat) %>% as_tibble()

pred_dat[['log_biomassarea']] <- y_hat$Estimate
pred_dat[['log_biomassarea_2.5']] <- y_hat$Q2.5
pred_dat[['log_biomassarea_97.5']] <- y_hat$Q97.5

#unscale Netflow
rm(meanINfN)
meanINfN<-mean(ACTIVE.1.sub.V2$Netflow,na.rm=T)
rm(sdINfN)
sdINfN<- 1*sd(ACTIVE.1.sub.V2$Netflow,na.rm=T)
FUNINVLUIF = function(x){(x * sdINfN) + meanINfN}
pred_dat$Netflowun<- round(FUNINVLUIF(pred_dat$Netflow))
head(pred_dat)
#unscale Graviy
rm(meanINfN)
meanINfN<-mean(ACTIVE.1.sub.V2$log_grav_total,na.rm=T)
rm(sdINfN)
sdINfN<- 1*sd(ACTIVE.1.sub.V2$log_grav_total,na.rm=T)
FUNINVLUIF = function(x){(x * sdINfN) + meanINfN}
pred_dat$log_grav_totalUn<- FUNINVLUIF(pred_dat$log_grav_total)
head(pred_dat)

pred_dat$Class <- factor(pred_dat$Class, levels=c("Closed","Restricted","Fished"))

threewayGlobalActive <- ggplot(NULL,aes(x,y)) +
  geom_smooth(data=pred_dat,aes(x=log_grav_totalUn,y=log_biomassarea,colour=as.factor(Netflowun)),method="lm") +
  geom_ribbon(data=pred_dat,aes(x=log_grav_totalUn,y=log_biomassarea,ymin=log_biomassarea_2.5,ymax=log_biomassarea_97.5,fill=as.factor(Netflowun),colour=NULL), alpha = .15)+
  scale_colour_manual(values = c("#39568CFF","#1F968BFF","#95D840FF")) +
  scale_fill_manual(values = c("#39568CFF","#1F968BFF","#95D840FF")) +
  geom_rug(data=ACTIVE.1.sub.V2,aes(y=log_biomassarea,x=log_grav_total)) +
  theme_bw() +  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  ggtitle("") + labs(colour="Netflow",fill=F) +
  ylab("Log Standing Biomass (g/m2)") +theme(text = element_text(size = 15, family="Helvetica")) +
  xlab("Log Human Gravity") + facet_grid(~Class)

threewayGlobalActive
