### script to load data for the connectivity project when necessary


# brms

#if (!requireNamespace("remotes")) {
#  install.packages("remotes")
#
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan", build_opts = "")

# stan
#remove.packages("rstan")
#if (file.exists(".RData")) file.remove(".RData")
#pkgbuild::has_build_tools(debug = TRUE)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# packagesa
#install.packages(c("dplyr","here","forcats","corrgram","ggpubr"))

require(dplyr)
require(here)
require(forcats)
library(corrgram) # for corrgram
library(ggpubr)


# load data
rm(all.data)
all.data<-read.csv(here("_data","Connectivity_Biomass_data.csv"),h=T, stringsAsFactors = F,dec=".")

# clean first column
all.data$X <- NULL

# check 
head(all.data)
summary(all.data)
dim(all.data)
summary(all.data)
str(all.data)
names(all.data)
apply(all.data,2,class)

# log all the data
all.data$log_grav_total <- log(all.data$grav_total+1)
all.data$log_biomassarea <-log(all.data$biomassare+1)

# chage to factor
all.data$region <- as.factor(all.data$region)
all.data$locality <- as.factor(all.data$locality)
all.data$sites <- as.factor(all.data$sites)
all.data$Class <- as.factor(all.data$Class)
all.data$ModelMode <- as.factor(all.data$ModelMode)
all.data$Larval_behaviour <- as.factor(all.data$Larval_beh)
all.data$FE <- as.factor(all.data$FE)


  # summarize productivity to annual productivity
all.data <- all.data%>%   
  rowwise() %>%
  mutate(prod.annual = mean(c(prod.Jan:prod.Dec),na.rm=T)) %>%
  as.data.frame()

# rm NAs for netflow
all.data <- all.data %>% filter(!is.na(Netflow))
head(all.data)
# sub dataframe for the global models
  # ACTIVE.1
ACTIVE.1 <- all.data %>% filter(Larval_behaviour == "active" & ModelMode %in% c("pare5","crypto5","resid15","transi15")) %>% droplevels()
rm(ACTIVE.1.sub)
ACTIVE.1.sub <- ACTIVE.1 %>% dplyr::group_by(ID) %>%
  summarise(temp = mean(temp),Richness = mean(Richness),grav_total=mean(grav_total),
            biomassare = mean(biomassare), OutFlow = sum(OutFlow),Outdegree = sum(Outdegree),btwdegree = sum(btwdegree),
            SelfR = sum(SelfR),Inflow = sum(Inflow),Indegree = sum(Indegree),
            IndegreeNe = sum(IndegreeNe), InflowNei = sum(InflowNei), Netflow = mean(Netflow),prod.annual = mean(prod.annual))
ACTIVE.1.sub.V2 <- left_join(ACTIVE.1 %>% dplyr::select(ID,region,locality,sites,Kulbicki,        
PROVINCE,Class,Lon,Lat),ACTIVE.1.sub,by="ID") %>% distinct()
head(ACTIVE.1.sub.V2)           
dim(ACTIVE.1.sub.V2) 

# log transformed skewed data
ACTIVE.1.sub.V2$log_biomassarea <- log(ACTIVE.1.sub.V2$biomassare+1)
ACTIVE.1.sub.V2$log_btwdegree <- log(ACTIVE.1.sub.V2$btwdegree+1)
ACTIVE.1.sub.V2$log_SelfR <- log(ACTIVE.1.sub.V2$SelfR+1)
ACTIVE.1.sub.V2$log_InflowNei <- log(ACTIVE.1.sub.V2$InflowNei+1)
ACTIVE.1.sub.V2$log_Inflow <- log(ACTIVE.1.sub.V2$Inflow+1)
ACTIVE.1.sub.V2$log_annual_prod <- log(ACTIVE.1.sub.V2$prod.annual+1)
ACTIVE.1.sub.V2$log_Indegree <- log(ACTIVE.1.sub.V2$Indegree+1)
ACTIVE.1.sub.V2$log_Indegree_Neigh <- log(ACTIVE.1.sub.V2$IndegreeNe+1)
ACTIVE.1.sub.V2$log_Outdegree <- log(ACTIVE.1.sub.V2$Outdegree+1)
ACTIVE.1.sub.V2$log_Outflow <- log(ACTIVE.1.sub.V2$OutFlow+1)
ACTIVE.1.sub.V2$log_grav_total <- log(ACTIVE.1.sub.V2$grav_total+1)

rm(ACTIVE.1.sub.V2.std)
ACTIVE.1.sub.V2.std<-data.frame(apply(X = ACTIVE.1.sub.V2[,c("Richness","temp","prod.annual",
                                                             "Netflow","log_grav_total",
                                                             "log_btwdegree","log_SelfR","log_Inflow",  
                                                             "log_InflowNei","log_Indegree","log_Indegree_Neigh","log_Outdegree","log_Outflow")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
ACTIVE.1.sub.V2.std <- cbind(ACTIVE.1.sub.V2$region,ACTIVE.1.sub.V2.std)
ACTIVE.1.sub.V2.std <- cbind(ACTIVE.1.sub.V2.std,ACTIVE.1.sub.V2$Class)
ACTIVE.1.sub.V2.std$logbiomassarea <- log(ACTIVE.1.sub.V2$biomassare+1)
colnames(ACTIVE.1.sub.V2.std)[1] <- "region"
colnames(ACTIVE.1.sub.V2.std)[15] <- "Class"
head(ACTIVE.1.sub.V2.std)

# Fished as the reference
ACTIVE.1.sub.V2.std$Class <- relevel(ACTIVE.1.sub.V2.std$Class, ref="Fished")
summary(ACTIVE.1.sub.V2.std)

# save all.data file
#saveRDS(ACTIVE.1.sub.V2.std,here::here("_data","Datasets","Global","ACTIVE.1.sub.V2.std_May.rds"))
#saveRDS(ACTIVE.1.sub.V2,here::here("_data","Datasets","Global","ACTIVE.1.sub.V2_May.rds"))


# ACTIVE.2
rm(ACTIVE.2)
ACTIVE.2 <- all.data %>% filter(Larval_behaviour == "active" & ModelMode %in% c("pare15","crypto15","resid35","transi35")) %>% droplevels()
rm(ACTIVE.2.sub)
ACTIVE.2.sub <- ACTIVE.2 %>% dplyr::group_by(ID) %>%
  summarise(temp = mean(temp),Richness = mean(Richness),grav_total=mean(grav_total),
            biomassare = mean(biomassare), OutFlow = sum(OutFlow),Outdegree = sum(Outdegree),btwdegree = sum(btwdegree),
            SelfR = sum(SelfR),Inflow = sum(Inflow),Indegree = sum(Indegree),
            IndegreeNe = sum(IndegreeNe), InflowNei = sum(InflowNei), Netflow = mean(Netflow),prod.annual = mean(prod.annual))
ACTIVE.2.sub.V2 <- left_join(ACTIVE.2 %>% dplyr::select(ID,region,locality,sites,Kulbicki,        
                                                        PROVINCE,Class,Lon,Lat),ACTIVE.2.sub,by="ID") %>% distinct()
#head(ACTIVE.2.sub.V2)           

# log transformed skewed data
ACTIVE.2.sub.V2$log_biomassarea <- log(ACTIVE.2.sub.V2$biomassare+1)
ACTIVE.2.sub.V2$log_btwdegree <- log(ACTIVE.2.sub.V2$btwdegree+1)
ACTIVE.2.sub.V2$log_SelfR <- log(ACTIVE.2.sub.V2$SelfR+1)
ACTIVE.2.sub.V2$log_InflowNei <- log(ACTIVE.2.sub.V2$InflowNei+1)
ACTIVE.2.sub.V2$log_Inflow <- log(ACTIVE.2.sub.V2$Inflow+1)
ACTIVE.2.sub.V2$log_annual_prod <- log(ACTIVE.2.sub.V2$prod.annual+1)
ACTIVE.2.sub.V2$log_Indegree <- log(ACTIVE.2.sub.V2$Indegree+1)
ACTIVE.2.sub.V2$log_Indegree_Neigh <- log(ACTIVE.2.sub.V2$IndegreeNe+1)
ACTIVE.2.sub.V2$log_Outdegree <- log(ACTIVE.2.sub.V2$Outdegree+1)
ACTIVE.2.sub.V2$log_Outflow <- log(ACTIVE.2.sub.V2$OutFlow+1)
ACTIVE.2.sub.V2$log_grav_total <- log(ACTIVE.2.sub.V2$grav_total+1)

rm(ACTIVE.2.sub.V2.std)
ACTIVE.2.sub.V2.std<-data.frame(apply(X = ACTIVE.2.sub.V2[,c("Richness","temp","prod.annual",
                                                             "Netflow","log_grav_total",
                                                             "log_btwdegree","log_SelfR","log_Inflow",  
                                                             "log_InflowNei","log_Indegree","log_Indegree_Neigh","log_Outdegree","log_Outflow")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
ACTIVE.2.sub.V2.std <- cbind(ACTIVE.2.sub.V2$region,ACTIVE.2.sub.V2.std)
ACTIVE.2.sub.V2.std <- cbind(ACTIVE.2.sub.V2.std,ACTIVE.2.sub.V2$Class)
ACTIVE.2.sub.V2.std$logbiomassarea <- log(ACTIVE.2.sub.V2$biomassare+1)
colnames(ACTIVE.2.sub.V2.std)[1] <- "region"
colnames(ACTIVE.2.sub.V2.std)[15] <- "Class"
head(ACTIVE.2.sub.V2.std)

ACTIVE.2.sub.V2.std$Class <- relevel(ACTIVE.2.sub.V2.std$Class, ref="Fished")
summary(ACTIVE.2.sub.V2.std)

# save all.data file
#saveRDS(ACTIVE.2.sub.V2.std,here::here("_data","Connectivity_Biomass_SEMGLMMDATA_April2021_GLOBAL_ACTIVE2_std.rds"))


# PASSIVE
rm(PASSIVE)
PASSIVE <- all.data %>% filter(Larval_behaviour == "passive") %>% droplevels()
rm(PASSIVE.sub)
PASSIVE.sub <- PASSIVE %>% dplyr::group_by(ID) %>%
  summarise(temp = mean(temp),Richness = mean(Richness),grav_total=mean(grav_total),
            biomassare = mean(biomassare), OutFlow = sum(OutFlow),Outdegree = sum(Outdegree),btwdegree = sum(btwdegree),
            SelfR = sum(SelfR),Inflow = sum(Inflow),Indegree = sum(Indegree),
            IndegreeNe = sum(IndegreeNe), InflowNei = sum(InflowNei), Netflow = mean(Netflow),prod.annual = mean(prod.annual))
PASSIVE.sub.V2 <- left_join(PASSIVE %>% dplyr::select(ID,region,locality,sites,Kulbicki,        
                                                        PROVINCE,Class,Lon,Lat),PASSIVE.sub,by="ID") %>% distinct()
head(PASSIVE.sub.V2)           

# log transformed skewed data
PASSIVE.sub.V2$log_biomassarea <- log(PASSIVE.sub.V2$biomassare+1)
PASSIVE.sub.V2$log_btwdegree <- log(PASSIVE.sub.V2$btwdegree+1)
PASSIVE.sub.V2$log_SelfR <- log(PASSIVE.sub.V2$SelfR+1)
PASSIVE.sub.V2$log_InflowNei <- log(PASSIVE.sub.V2$InflowNei+1)
PASSIVE.sub.V2$log_Inflow <- log(PASSIVE.sub.V2$Inflow+1)
PASSIVE.sub.V2$log_annual_prod <- log(PASSIVE.sub.V2$prod.annual+1)
PASSIVE.sub.V2$log_Indegree <- log(PASSIVE.sub.V2$Indegree+1)
PASSIVE.sub.V2$log_Indegree_Neigh <- log(PASSIVE.sub.V2$IndegreeNe+1)
PASSIVE.sub.V2$log_Outdegree <- log(PASSIVE.sub.V2$Outdegree+1)
PASSIVE.sub.V2$log_Outflow <- log(PASSIVE.sub.V2$OutFlow+1)
PASSIVE.sub.V2$log_grav_total <- log(PASSIVE.sub.V2$grav_total+1)

rm(PASSIVE.sub.V2.std)
PASSIVE.sub.V2.std<-data.frame(apply(X = PASSIVE.sub.V2[,c("Richness","temp","prod.annual",
                                                             "Netflow","log_grav_total",
                                                             "log_btwdegree","log_SelfR","log_Inflow",  
                                                             "log_InflowNei","log_Indegree","log_Indegree_Neigh","log_Outdegree","log_Outflow")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
PASSIVE.sub.V2.std <- cbind(PASSIVE.sub.V2$region,PASSIVE.sub.V2.std)
PASSIVE.sub.V2.std <- cbind(PASSIVE.sub.V2.std,PASSIVE.sub.V2$Class)
PASSIVE.sub.V2.std$logbiomassarea <- log(PASSIVE.sub.V2$biomassare+1)
colnames(PASSIVE.sub.V2.std)[1] <- "region"
colnames(PASSIVE.sub.V2.std)[15] <- "Class"
head(PASSIVE.sub.V2.std)

# reference level
PASSIVE.sub.V2.std$Class <- relevel(PASSIVE.sub.V2.std$Class, ref="Fished")
summary(PASSIVE.sub.V2.std)


# save all.data file
#saveRDS(PASSIVE.sub.V2.std,here::here("_data","Connectivity_Biomass_SEMGLMMDATA_April2021_GLOBAL_PASSIVE_std.rds"))


