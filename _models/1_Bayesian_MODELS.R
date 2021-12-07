# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# updated: May 2021
# outputs: SEM coefficients


#rm(list=ls())

# brms

#if (!requireNamespace("remotes")) {
#  install.packages("remotes")
#}
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan", build_opts = "")
 
# stan
#remove.packages("rstan")
#if (file.exists(".RData")) file.remove(".RData")

#Sys.setenv(MAKEFLAGS = "-j10") # 10 cores used
#install.packages(c("Rcpp", "RcppEigen", "RcppParallel", "StanHeaders"), type = "source")
#install.packages("rstan", type = "source")

# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)
require(Rserve)
require(rJava)
library(performance) # for r2_bayes
library(here)


#########################
#####BUILDING MODELS#####
#########################----
#### SPECIES RICHNESS = Environment + Connectivity
rm(MOD_1_S_null) # NULL
MOD_1_S_null <- bf(Richness ~ 1 + (1 |region))

#### SPECIES RICHNESS = Environment
S_mod_nocon <- bf(Richness ~ temp +  prod.annual + Class + log_grav_total + (1 +temp +  prod.annual + Class + log_grav_total |region))

# models with connectivity
rm(MOD_1_S) # ENV + CON (btw,netflow,indegreeNei,SelfR)
MOD_1_S <- bf(Richness ~ temp + prod.annual + 
                log_btwdegree + Netflow +
                log_SelfR + log_Indegree_Neigh + 
                log_grav_total + Class + (1 + temp + prod.annual + 
                                            log_btwdegree + Netflow +
                                            log_SelfR + log_Indegree_Neigh + 
                                            log_grav_total |region))


#ALTERNATIVE MODELS WITH INFLOW/OUTDEGREE/INFLOWNEI/OUTFLOW/INDEGREE

rm(MOD_2_S) # ENV + CON (selfR, btw, indegree, outflow)
MOD_2_S <- bf(Richness ~ temp + prod.annual + log_SelfR + 
                log_btwdegree + log_Indegree + log_Outflow + 
                log_grav_total + Class + 
                (1 + temp + prod.annual + log_SelfR + 
                   log_btwdegree + log_Indegree + Netflow + 
                   log_grav_total |region))

rm(MOD_3_S) # ENV + CON (btw, inflow, outdegree)
MOD_3_S <- bf(Richness ~ temp + prod.annual + log_Outdegree + 
                log_btwdegree + log_Inflow +
                log_grav_total + Class + (1 + temp + prod.annual + 
                                            log_Outdegree + log_btwdegree +
                                            log_Inflow +
                                            log_grav_total |region))

rm(MOD_4_S) # ENV + CON (selfR, btw, inflowNEI, netflow)
MOD_4_S <- bf(Richness ~ temp + prod.annual + log_SelfR + log_btwdegree +
                log_InflowNei + Netflow +
                log_grav_total + Class + (1 + temp + prod.annual + 
                                            log_SelfR + log_btwdegree +
                                            log_InflowNei + Netflow + 
                                            log_grav_total |region))


#### BIOMASS = Environment + Connectivity
# null model
rm(MOD_B.null) # ENV + CON
MOD_B.null <- bf(log_biomassarea ~ 1 + (1 |region))

#### Biomass = Environment
B_mod_nocon <- bf(log_biomassarea ~ temp +  prod.annual + Class*log_grav_total + 
                    (1 + temp +  prod.annual + log_grav_total  |region))

# full model #ENV + CON
rm(MOD_1_B)
MOD_1_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class*Netflow + 
                log_SelfR +  log_InflowNei +
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total+Netflow + log_InflowNei +
                   log_SelfR |region))

#Alternative models with inflow/indegree/outflow/outdegree)
rm(MOD_2_B) # ENV + CON
MOD_2_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class + 
                log_SelfR + 
                log_Indegree + log_Outflow +
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total  + 
                   log_SelfR + log_Outflow +
                   log_Indegree |region))

rm(MOD_3_B) # ENV + CON
MOD_3_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class + 
                log_SelfR + 
                log_Outdegree  + log_Inflow +
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total + 
                   log_SelfR + 
                   log_Outdegree  + log_Inflow|region))


####Load .std datasets for to rerun models here(_data, "Datasets","Active 1") for each functional group and 
#here(_data, "Datasets","Global") for total connectivity ----

##########################
###      TRANSIENT    ####
##########################----

### run species richness models
  # null model
MOD_S_1_run_TRANSIENT_null %>% rm()  #  
MOD_S_1_run_TRANSIENT_null <-brm(MOD_1_S_null , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_1_run_TRANSIENT_null,"Models/ACTIVE V1/MOD_S_1_run_TRANSIENT_null_V2_May.Rds")

  # model 1
MOD_S_1_run_TRANSIENT %>% rm()  #  
MOD_S_1_run_TRANSIENT <-brm(MOD_1_S , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_1_run_TRANSIENT,"Models/ACTIVE V1/MOD_S_1_run_TRANSIENT_V2_May.Rds")

  # model 2
MOD_S_2_run_TRANSIENT %>% rm()  #  
MOD_S_2_run_TRANSIENT <-brm(MOD_2_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_2_run_TRANSIENT,"Models/ACTIVE V1/MOD_S_2_run_TRANSIENT_V2_May.Rds")

# model 3
MOD_S_3_run_TRANSIENT %>% rm()  #  
MOD_S_3_run_TRANSIENT <-brm(MOD_3_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_3_run_TRANSIENT,"Models/ACTIVE V1/MOD_S_3_run_TRANSIENT_V2_May.Rds")

# model 4
MOD_S_4_run_TRANSIENT %>% rm()  #  
MOD_S_4_run_TRANSIENT <-brm(MOD_4_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_4_run_TRANSIENT,"Models/ACTIVE V1/MOD_S_4_run_TRANSIENT_V2_May.Rds")

# model S env
MOD_S_env.FE %>% rm()  #  
MOD_S_env.FE <-brm(S_mod_nocon , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_env.FE,"Models/ACTIVE V1/MOD_S_env_V2_May.Rds")

#### test all S models transient
r2_bayes.transient.S <- rbind(c(r2_bayes(MOD_S_1_run_TRANSIENT_null)[1],
r2_bayes(MOD_S_1_run_TRANSIENT)[1],
r2_bayes(MOD_S_2_run_TRANSIENT)[1],
r2_bayes(MOD_S_3_run_TRANSIENT)[1],
r2_bayes(MOD_S_4_run_TRANSIENT)[1],
r2_bayes(MOD_S_env.FE)[1]),
c(r2_bayes(MOD_S_1_run_TRANSIENT_null)[2],
  r2_bayes(MOD_S_1_run_TRANSIENT)[2],
  r2_bayes(MOD_S_2_run_TRANSIENT)[2],
  r2_bayes(MOD_S_3_run_TRANSIENT)[2],
  r2_bayes(MOD_S_4_run_TRANSIENT)[2],
  r2_bayes(MOD_S_env.FE)[2]))
colnames(r2_bayes.transient.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.transient.S) <- rep("TRANSIENT",2)
#weights
model.W.transient.S <- model_weights(MOD_S_1_run_TRANSIENT_null,MOD_S_1_run_TRANSIENT,MOD_S_env.FE,weights="loo")
names(model.W.transient.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_TRANSIENT, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=TRANSIENT.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_TRANSIENT_null,MOD_S_1_run_TRANSIENT,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_TRANSIENT)
#pp_check(MOD_S_1_run_TRANSIENT, resp="log_biomassarea")

#brms::pp_check(MOD_S_1_run_TRANSIENT, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_TRANSIENT, resp="Richness", nsamples = 100)

#r2_bayes(MOD_S_1_run_TRANSIENT_null)

### run BIOMASS models
#NULL
MOD_BIOM_1_run_TRANSIENT.null %>% rm()  #  
MOD_BIOM_1_run_TRANSIENT.null <-brm(MOD_B.null , data=TRANSIENT.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
#saveRDS(MOD_BIOM_1_run_TRANSIENT.null,"Models/ACTIVE V1/MOD_BIOM_1_run_TRANSIENT_null_V2.May.Rds")

# model B env
MOD_B_env %>% rm()  #  
MOD_B_env.FE <-brm(B_mod_nocon , data=TRANSIENT.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_B_env,"Models/ACTIVE V1/MOD_B_env_V2_May.Rds")

#MODEL 1
MOD_BIOM_1_run_TRANSIENT %>% rm()  #  
MOD_BIOM_1_run_TRANSIENT <-brm(MOD_1_B , data=TRANSIENT.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_1_run_TRANSIENT,"Models/ACTIVE V1/MOD_BIOM_1_run_TRANSIENT_V2_May.Rds")

#MODEL 2
MOD_BIOM_2_run_TRANSIENT %>% rm()  #  
MOD_BIOM_2_run_TRANSIENT <-brm(MOD_2_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_2_run_TRANSIENT,"Models/ACTIVE V1/MOD_BIOM_2_run_TRANSIENT_V2_May.Rds")


#MODEL 3
MOD_BIOM_3_run_TRANSIENT %>% rm()  #  
MOD_BIOM_3_run_TRANSIENT <-brm(MOD_3_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_3_run_TRANSIENT,"Models/ACTIVE V1/MOD_BIOM_3_run_TRANSIENT_V2_May.Rds")


#### test all B models
r2_bayes.transient.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_TRANSIENT.null)[1],
                                r2_bayes(MOD_BIOM_1_run_TRANSIENT)[1],
                                r2_bayes(MOD_BIOM_2_run_TRANSIENT)[1],
                                r2_bayes(MOD_BIOM_3_run_TRANSIENT)[1],
                                r2_bayes(MOD_B_env.FE)[1]),
                              c(r2_bayes(MOD_BIOM_1_run_TRANSIENT.null)[2],
                                r2_bayes(MOD_BIOM_1_run_TRANSIENT)[2],
                                r2_bayes(MOD_BIOM_2_run_TRANSIENT)[2],
                                r2_bayes(MOD_BIOM_3_run_TRANSIENT)[2],
                                r2_bayes(MOD_B_env.FE)[2]))
colnames(r2_bayes.transient.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.transient.B) <- rep("TRANSIENT",2)

model.W.transient.B <- model_weights(MOD_BIOM_1_run_TRANSIENT.null,
                                     MOD_BIOM_1_run_TRANSIENT,
                                     MOD_B_env.FE,weights="loo")
names(model.W.transient.B) <- c("B_model_null","B_model_1","B_model_ENV")

#rownames(model.W.transient.B) <- rep("TRANSIENT",1)

##########################
###      CRYPTIC      ####
##########################----

### run species models
# null model
MOD_S_1_run_CRYPTIC_null %>% rm()  #  
MOD_S_1_run_CRYPTIC_null <-brm(MOD_1_S_null , data=CRYPTIC.std,cores=4,chains = 4,
                                 iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_1_run_CRYPTIC_null,"Models/ACTIVE V1/MOD_S_1_run_CRYPTIC_null_V2_May.Rds")

# model S env
MOD_S_env.FE %>% rm()  #  
MOD_S_env.FE <-brm(S_mod_nocon , data=CRYPTIC.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_env.FE,here("Models output","ACTIVE 1","Cryptic","Richness","MOD_S_env.FE_June.Rds"))

# model 1
MOD_S_1_run_CRYPTIC %>% rm()  #  
MOD_S_1_run_CRYPTIC <-brm(MOD_1_S , data=CRYPTIC.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_1_run_CRYPTIC,"Models/ACTIVE V1/MOD_S_1_run_CRYPTIC_V2_May.Rds")

# model 2
MOD_S_2_run_CRYPTIC %>% rm()  #  
MOD_S_2_run_CRYPTIC <-brm(MOD_2_S , data=CRYPTIC.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_2_run_CRYPTIC,"Models/ACTIVE V1/MOD_S_2_run_CRYPTIC_V2_May.Rds")

# model 3
MOD_S_3_run_CRYPTIC %>% rm()  #  
MOD_S_3_run_CRYPTIC <-brm(MOD_3_S , data=CRYPTIC.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_3_run_CRYPTIC,"Models/ACTIVE V1/MOD_S_3_run_CRYPTIC_V2_May.Rds")

# model 4
MOD_S_4_run_CRYPTIC %>% rm()  #  
MOD_S_4_run_CRYPTIC <-brm(MOD_4_S , data=CRYPTIC.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_4_run_CRYPTIC,"Models/ACTIVE V1/MOD_S_4_run_CRYPTIC_V2_May.Rds")


#### test all S models CRYPTIC
r2_bayes.CRYPTIC.S <- rbind(c(r2_bayes(MOD_S_1_run_CRYPTIC_null)[1],
                                r2_bayes(MOD_S_1_run_CRYPTIC)[1],
                                r2_bayes(MOD_S_2_run_CRYPTIC)[1],
                                r2_bayes(MOD_S_3_run_CRYPTIC)[1],
                                r2_bayes(MOD_S_4_run_CRYPTIC)[1],
                                r2_bayes(MOD_S_env.FE)[1]),
                              c(r2_bayes(MOD_S_1_run_CRYPTIC_null)[2],
                                r2_bayes(MOD_S_1_run_CRYPTIC)[2],
                                r2_bayes(MOD_S_2_run_CRYPTIC)[2],
                                r2_bayes(MOD_S_3_run_CRYPTIC)[2],
                                r2_bayes(MOD_S_4_run_CRYPTIC)[2],
                                r2_bayes(MOD_S_env.FE)[2]))
colnames(r2_bayes.CRYPTIC.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.CRYPTIC.S) <- rep("CRYPTIC",2)
model.W.CRYPTIC.S <- model_weights(MOD_S_1_run_CRYPTIC_null,
                                   MOD_S_1_run_CRYPTIC,
                                   MOD_S_env.FE,weights="loo")
names(model.W.CRYPTIC.S) <- c("S_model_null","S_model_1","S_model_ENV")
model.W.CRYPTIC.S


#fit3 <- update(MOD_S_1_run_CRYPTIC, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=CRYPTIC.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_CRYPTIC_null,MOD_S_1_run_CRYPTIC,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_CRYPTIC)
#pp_check(MOD_S_1_run_CRYPTIC, resp="log_biomassarea")
#brms::pp_check(MOD_S_1_run_CRYPTIC, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_CRYPTIC, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_CRYPTIC_null)

### run biomass models
## NULL
MOD_BIOM_1_run_CRYPTIC.null %>% rm()  # 
MOD_BIOM_1_run_CRYPTIC.null <-brm(MOD_B.null , data=CRYPTIC.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
#saveRDS(MOD_BIOM_1_run_CRYPTIC.null,"Models/ACTIVE V1/MOD_BIOM_1_run_CRYPTIC_null_V2.May.Rds")

# model B env
MOD_B_env %>% rm()  #  
MOD_B_env.FE <-brm(B_mod_nocon , data=CRYPTIC.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_B_env.FE,here("Models output","ACTIVE 1","Cryptic","Biomass","MOD_B_env.FE_June.Rds"))

#MODEL 1
MOD_BIOM_1_run_CRYPTIC %>% rm()  #
MOD_BIOM_1_run_CRYPTIC <-brm(MOD_1_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_1_run_CRYPTIC,"Models/ACTIVE V1/MOD_BIOM_1_run_CRYPTIC_V2_May.Rds")

#MODEL 2
MOD_BIOM_2_run_CRYPTIC %>% rm()  # 
MOD_BIOM_2_run_CRYPTIC <-brm(MOD_2_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_2_run_CRYPTIC,"Models/ACTIVE V1/MOD_BIOM_2_run_CRYPTIC_V2_May.Rds")

#MODEL 3
MOD_BIOM_3_run_CRYPTIC %>% rm()  #  
MOD_BIOM_3_run_CRYPTIC <-brm(MOD_3_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_3_run_CRYPTIC,"Models/ACTIVE V1/MOD_BIOM_3_run_CRYPTIC_V2_May.Rds")

#### test all B models
r2_bayes.CRYPTIC.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_CRYPTIC.null)[1],
                                r2_bayes(MOD_BIOM_1_run_CRYPTIC)[1],
                                r2_bayes(MOD_BIOM_2_run_CRYPTIC)[1],
                                r2_bayes(MOD_BIOM_3_run_CRYPTIC)[1],
                                #r2_bayes(MOD_BIOM_4_run_CRYPTIC)[1],
                                r2_bayes(MOD_B_env.FE)[1]),
                              c(r2_bayes(MOD_BIOM_1_run_CRYPTIC.null)[2],
                                r2_bayes(MOD_BIOM_1_run_CRYPTIC)[2],
                                r2_bayes(MOD_BIOM_2_run_CRYPTIC)[2],
                                r2_bayes(MOD_BIOM_3_run_CRYPTIC)[2],
                                #r2_bayes(MOD_BIOM_4_run_CRYPTIC)[2],
                                r2_bayes(MOD_B_env.FE)[2]))
colnames(r2_bayes.CRYPTIC.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.CRYPTIC.B) <- rep("CRYPTIC",2)
model.W.CRYPTIC.B <- model_weights(MOD_BIOM_1_run_CRYPTIC.null,
                                   MOD_BIOM_1_run_CRYPTIC,
                                     MOD_B_env.FE,weights="loo")
names(model.W.CRYPTIC.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.CRYPTIC.B) <- rep("CRYPTIC",1)

##########################
###     RESIDENT      ####
##########################----

### run species models
# null model
MOD_S_1_run_RESID_null %>% rm()  #  
MOD_S_1_run_RESID_null <-brm(MOD_1_S_null , data=RESID.std,cores=4,chains = 4,
                                 iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_1_run_RESID_null,"Models/ACTIVE V1/MOD_S_1_run_RESID_null_V2_May.Rds")

# model S env
MOD_S_env.FE %>% rm()  #  
MOD_S_env.FE <-brm(S_mod_nocon , data=RESID.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_env.FE,here("Models output","ACTIVE 1","Resident","Richness","MOD_B_env.FE_June.Rds"))

# model 1
MOD_S_1_run_RESID %>% rm()  #  
MOD_S_1_run_RESID <-brm(MOD_1_S , data=RESID.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_1_run_RESID,"Models/ACTIVE V1/MOD_S_1_run_RESID_V2_May.Rds")

# model 2
MOD_S_2_run_RESID %>% rm()  #  
MOD_S_2_run_RESID <-brm(MOD_2_S , data=RESID.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_2_run_RESID,"Models/ACTIVE V1/MOD_S_2_run_RESID_V2_May.Rds")

# model 3
MOD_S_3_run_RESID %>% rm()  #  
MOD_S_3_run_RESID <-brm(MOD_3_S , data=RESID.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_3_run_RESID,"Models/ACTIVE V1/MOD_S_3_run_RESID_V2_May.Rds")

# model 4
MOD_S_4_run_RESID %>% rm()  #  
MOD_S_4_run_RESID <-brm(MOD_4_S , data=RESID.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_4_run_RESID,"Models/ACTIVE V1/MOD_S_4_run_RESID_V2_May.Rds")

#### test all S models RESID
r2_bayes.RESID.S <- rbind(c(r2_bayes(MOD_S_1_run_RESID_null)[1],
                                r2_bayes(MOD_S_1_run_RESID)[1],
                                r2_bayes(MOD_S_2_run_RESID)[1],
                                r2_bayes(MOD_S_3_run_RESID)[1],
                                r2_bayes(MOD_S_4_run_RESID)[1],
                                r2_bayes(MOD_S_env.FE)[1]),
                              c(r2_bayes(MOD_S_1_run_RESID_null)[2],
                                r2_bayes(MOD_S_1_run_RESID)[2],
                                r2_bayes(MOD_S_2_run_RESID)[2],
                                r2_bayes(MOD_S_3_run_RESID)[2],
                                r2_bayes(MOD_S_4_run_RESID)[2],
                                r2_bayes(MOD_S_env.FE)[2]))
colnames(r2_bayes.RESID.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.RESID.S) <- rep("RESID",2)
model.W.RESID.S <- model_weights(MOD_S_1_run_RESID_null,
                                 MOD_S_1_run_RESID,
                                  MOD_S_env.FE,weights="loo")
names(model.W.RESID.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_RESID, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=RESID.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_RESID_null,MOD_S_1_run_RESID,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_RESID)
#pp_check(MOD_S_1_run_RESID, resp="log_biomassarea")
#brms::pp_check(MOD_S_1_run_RESID, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_RESID, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_RESID_null)


### run biomass models
## NULL
MOD_BIOM_1_run_RESID.null %>% rm()  #  
MOD_BIOM_1_run_RESID.null <-brm(MOD_B.null , data=RESID.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
#saveRDS(MOD_BIOM_1_run_RESID.null,"Models/ACTIVE V1/MOD_BIOM_1_run_RESID_null_V2.May.Rds")

# model B env
MOD_B_env %>% rm()  #  
MOD_B_env.FE <-brm(B_mod_nocon , data=RESID.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_B_env.FE,here("Models output","ACTIVE 1","Resident","Biomass","MOD_B_env.FE_June.Rds"))

#MODEL 1
MOD_BIOM_1_run_RESID %>% rm()  #  
MOD_BIOM_1_run_RESID <-brm(MOD_1_B , data=RESID.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_1_run_RESID,"Models/ACTIVE V1/MOD_BIOM_1_run_RESID_V2_May.Rds")

#MODEL 2
MOD_BIOM_2_run_RESID %>% rm()  #  
MOD_BIOM_2_run_RESID <-brm(MOD_2_B , data=RESID.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_2_run_RESID,"Models/ACTIVE V1/MOD_BIOM_2_run_RESID_V2_May.Rds")

#MODEL 3
MOD_BIOM_3_run_RESID %>% rm()  #  
MOD_BIOM_3_run_RESID <-brm(MOD_3_B , data=RESID.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_3_run_RESID,"Models/ACTIVE V1/MOD_BIOM_3_run_RESID_V2_May.Rds")

#### test all B models
r2_bayes.RESID.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_RESID.null)[1],
                                r2_bayes(MOD_BIOM_1_run_RESID)[1],
                                r2_bayes(MOD_BIOM_2_run_RESID)[1],
                                r2_bayes(MOD_BIOM_3_run_RESID)[1],
                                #r2_bayes(MOD_BIOM_4_run_RESID)[1],
                                r2_bayes(MOD_B_env.FE)[1]),
                              c(r2_bayes(MOD_BIOM_1_run_RESID.null)[2],
                                r2_bayes(MOD_BIOM_1_run_RESID)[2],
                                r2_bayes(MOD_BIOM_2_run_RESID)[2],
                                r2_bayes(MOD_BIOM_3_run_RESID)[2],
                                #r2_bayes(MOD_BIOM_4_run_RESID)[2],
                                r2_bayes(MOD_B_env.FE)[2]))
colnames(r2_bayes.RESID.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.RESID.B) <- rep("RESID",2)
model.W.RESID.B <- model_weights(MOD_BIOM_1_run_RESID.null,
                                 MOD_BIOM_1_run_RESID,
                                MOD_B_env.FE,weights="loo")
names(model.W.RESID.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.RESID.B) <- rep("RESID",1)


##########################
###     PARENTAL      ####
##########################----


### run species models
# null model
MOD_S_1_run_PARENTAL_null %>% rm()  #  
MOD_S_1_run_PARENTAL_null <-brm(MOD_1_S_null , data=PARENTAL.std,cores=4,chains = 4,
                             iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_1_run_PARENTAL_null,"Models/ACTIVE V1/MOD_S_1_run_PARENTAL_null_V2_May.Rds")

# model S env
MOD_S_env.FE %>% rm()  #  
MOD_S_env.FE <-brm(S_mod_nocon , data=PARENTAL.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_env.FE,here("Models output","ACTIVE 1","Parental","Richness","MOD_B_env.FE_June.Rds"))


# model 1
MOD_S_1_run_PARENTAL %>% rm()  #  
MOD_S_1_run_PARENTAL <-brm(MOD_1_S , data=PARENTAL.std,cores=4,chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                        prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_1_run_PARENTAL,"Models/ACTIVE V1/MOD_S_1_run_PARENTAL_V2_May.Rds")

# model 2
MOD_S_2_run_PARENTAL %>% rm()  #  
MOD_S_2_run_PARENTAL <-brm(MOD_2_S , data=PARENTAL.std,cores=4,chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                        prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_2_run_PARENTAL,"Models/ACTIVE V1/MOD_S_2_run_PARENTAL_V2_May.Rds")

# model 3
MOD_S_3_run_PARENTAL %>% rm()  #  
MOD_S_3_run_PARENTAL <-brm(MOD_3_S , data=PARENTAL.std,cores=4,chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                        prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_3_run_PARENTAL,"Models/ACTIVE V1/MOD_S_3_run_PARENTAL_V2_May.Rds")

# model 4
MOD_S_4_run_PARENTAL %>% rm()  #  
MOD_S_4_run_PARENTAL <-brm(MOD_4_S , data=PARENTAL.std,cores=4,chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                        prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_4_run_PARENTAL,"Models/ACTIVE V1/MOD_S_4_run_PARENTAL_V2_May.Rds")

#### test all S models PARENTAL
r2_bayes.PARENTAL.S <- rbind(c(r2_bayes(MOD_S_1_run_PARENTAL_null)[1],
                                r2_bayes(MOD_S_1_run_PARENTAL)[1],
                                r2_bayes(MOD_S_2_run_PARENTAL)[1],
                                r2_bayes(MOD_S_3_run_PARENTAL)[1],
                                r2_bayes(MOD_S_4_run_PARENTAL)[1],
                                r2_bayes(MOD_S_env.FE)[1]),
                              c(r2_bayes(MOD_S_1_run_PARENTAL_null)[2],
                                r2_bayes(MOD_S_1_run_PARENTAL)[2],
                                r2_bayes(MOD_S_2_run_PARENTAL)[2],
                                r2_bayes(MOD_S_3_run_PARENTAL)[2],
                                r2_bayes(MOD_S_4_run_PARENTAL)[2],
                                r2_bayes(MOD_S_env.FE)[2]))
colnames(r2_bayes.PARENTAL.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.PARENTAL.S) <- rep("PARENTAL",2)
model.W.PARENTAL.S <- model_weights(MOD_S_1_run_PARENTAL_null,
                                    MOD_S_1_run_PARENTAL,
                            MOD_S_env.FE,weights="loo")
names(model.W.PARENTAL.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_PARENTAL, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=PARENTAL.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_PARENTAL_null,MOD_S_1_run_PARENTAL,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_PARENTAL)
#pp_check(MOD_S_1_run_PARENTAL, resp="log_biomassarea")
#brms::pp_check(MOD_S_1_run_PARENTAL, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_PARENTAL, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_PARENTAL_null)

### run biomass models
## NULL
MOD_BIOM_1_run_PARENTAL.null %>% rm()  #  
MOD_BIOM_1_run_PARENTAL.null <-brm(MOD_B.null , data=PARENTAL.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
#saveRDS(MOD_BIOM_1_run_PARENTAL.null,"Models/ACTIVE V1/MOD_BIOM_1_run_PARENTAL_null_V2.May.Rds")

# model B env
MOD_B_env %>% rm()  #  
MOD_B_env.FE <-brm(B_mod_nocon , data=PARENTAL.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_B_env.FE,here("Models output","ACTIVE 1","Parental","Biomass","MOD_B_env.FE_June.Rds"))

#MODEL 1
MOD_BIOM_1_run_PARENTAL %>% rm()  #  
MOD_BIOM_1_run_PARENTAL <-brm(MOD_1_B , data=PARENTAL.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_1_run_PARENTAL,"Models/ACTIVE V1/MOD_BIOM_1_run_PARENTAL_V2_May.Rds")

#MODEL 2
MOD_BIOM_2_run_PARENTAL %>% rm()  #  
MOD_BIOM_2_run_PARENTAL <-brm(MOD_2_B , data=PARENTAL.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_2_run_PARENTAL,"Models/ACTIVE V1/MOD_BIOM_2_run_PARENTAL_V2_May.Rds")

#MODEL 3
MOD_BIOM_3_run_PARENTAL %>% rm()  #  
MOD_BIOM_3_run_PARENTAL <-brm(MOD_3_B , data=PARENTAL.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_3_run_PARENTAL,"Models/ACTIVE V1/MOD_BIOM_3_run_PARENTAL_V2_May.Rds")

# model B env
MOD_B_env %>% rm()  #  
MOD_B_env <-brm(B_mod_nocon , data=PARENTAL.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_B_env,"Models/ACTIVE V1/MOD_B_env_V2_May.Rds")

#### test all B models
r2_bayes.PARENTAL.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_PARENTAL.null)[1],
                                r2_bayes(MOD_BIOM_1_run_PARENTAL)[1],
                                r2_bayes(MOD_BIOM_2_run_PARENTAL)[1],
                                r2_bayes(MOD_BIOM_3_run_PARENTAL)[1],
                                #r2_bayes(MOD_BIOM_4_run_PARENTAL)[1],
                                r2_bayes(MOD_B_env.FE)[1]),
                              c(r2_bayes(MOD_BIOM_1_run_PARENTAL.null)[2],
                                r2_bayes(MOD_BIOM_1_run_PARENTAL)[2],
                                r2_bayes(MOD_BIOM_2_run_PARENTAL)[2],
                                r2_bayes(MOD_BIOM_3_run_PARENTAL)[2],
                                #r2_bayes(MOD_BIOM_4_run_PARENTAL)[2],
                                r2_bayes(MOD_B_env.FE)[2]))
colnames(r2_bayes.PARENTAL.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.PARENTAL.B) <- rep("PARENTAL",2)
model.W.PARENTAL.B <- model_weights(MOD_BIOM_1_run_PARENTAL.null,MOD_BIOM_1_run_PARENTAL,
                                     MOD_B_env,weights="loo")
names(model.W.PARENTAL.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.PARENTAL.B) <- rep("PARENTAL",1)



##########################
### GLOBAL ACTIVE V1  ####
##########################----

### run species models
# null model
MOD_S_1_run_ACTIVE.1_null %>% rm()  #  
MOD_S_1_run_ACTIVE.1_null <-brm(MOD_1_S_null , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_1_run_ACTIVE_null,"Models/GLOBAL/MOD_S_1_run_ACTIVE.1_null_V2_May.Rds")

# model 1
MOD_S_1_run_ACTIVE.1 %>% rm()  #  
MOD_S_1_run_ACTIVE.1 <-brm(MOD_1_S , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_1_run_ACTIVE.1,"Models/GLOBAL/MOD_S_1_run_ACTIVE.1_V2_May.Rds")

# model 2
MOD_S_2_run_ACTIVE.1 %>% rm()  #  
MOD_S_2_run_ACTIVE.1 <-brm(MOD_2_S , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_2_run_ACTIVE.1,"Models/GLOBAL/MOD_S_2_run_ACTIVE.1_V2_May.Rds")

# model 3
MOD_S_3_run_ACTIVE.1 %>% rm()  #  
MOD_S_3_run_ACTIVE.1 <-brm(MOD_3_S , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_3_run_ACTIVE.1,"Models/GLOBAL/MOD_S_3_run_ACTIVE.1_V2_May.Rds")

# model 4
MOD_S_4_run_ACTIVE.1 %>% rm()  #  
MOD_S_4_run_ACTIVE.1 <-brm(MOD_4_S , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_4_run_ACTIVE.1,"Models/GLOBAL/MOD_S_4_run_ACTIVE.1_V2_May.Rds")

# model S env
MOD_S_env %>% rm()  #  
MOD_S_env <-brm(S_mod_nocon , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_S_env,"Models/GLOBAL/MOD_S_env_V2_May.Rds")

#### test all S models ACTIVE.1.sub.V2
#fix MOD_S_1_run_ACTIVE.1_null 
r2_bayes.ACTIVE.1.sub.V2.S <- rbind(c(r2_bayes(MOD_S_1_run_ACTIVE.1_null)[1],
                              r2_bayes(MOD_S_1_run_ACTIVE.1)[1],
                              r2_bayes(MOD_S_2_run_ACTIVE.1)[1],
                              r2_bayes(MOD_S_3_run_ACTIVE.1)[1],
                              r2_bayes(MOD_S_4_run_ACTIVE.1)[1],
                              r2_bayes(MOD_S_env)[1]),
                            c(r2_bayes(MOD_S_1_run_ACTIVE.1_null)[2],
                              r2_bayes(MOD_S_1_run_ACTIVE.1)[2],
                              r2_bayes(MOD_S_2_run_ACTIVE.1)[2],
                              r2_bayes(MOD_S_3_run_ACTIVE.1)[2],
                              r2_bayes(MOD_S_4_run_ACTIVE.1)[2],
                              r2_bayes(MOD_S_env)[2]))
colnames(r2_bayes.ACTIVE.1.sub.V2.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.ACTIVE.1.sub.V2.S) <- rep("ACTIVE.1.sub.V2",2)
model.W.ACTIVE.1.sub.V2.S <- model_weights(MOD_S_1_run_ACTIVE.1_null,MOD_S_1_run_ACTIVE.1,
                                 MOD_S_env,weights="loo")
names(model.W.ACTIVE.1.sub.V2.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_ACTIVE.1, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=ACTIVE.1.sub.V2.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_ACTIVE.1_null,MOD_S_1_run_ACTIVE.1,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_ACTIVE.1)
#pp_check(MOD_S_1_run_ACTIVE.1, resp="log_biomassarea")
brms::pp_check(MOD_S_1_run_ACTIVE.1, resp="Richness", nsamples = 100)
brms::pp_check(MOD_BIOM_1_run_ACTIVE.1, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_ACTIVE.1_null)

### run biomass models
ACTIVE.1.sub.V2.std$log_biomassarea <- ACTIVE.1.sub.V2.std$logbiomassarea
## ACTIVE.1.sub.V2
#null
MOD_BIOM_1_run_ACTIVE.1.null %>% rm()  #  
MOD_BIOM_1_run_ACTIVE.1.null <-brm(MOD_B.null , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                                  iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
#saveRDS(MOD_BIOM_1_run_ACTIVE.1.null,"Models/GLOBAL/MOD_BIOM_1_run_ACTIVE.1_null_V2.May.Rds")

# model S env
MOD_B_env %>% rm()  #  
MOD_B_env <-brm(B_mod_nocon , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#MODEL 1
MOD_BIOM_1_run_ACTIVE.1 %>% rm()  #  
MOD_BIOM_1_run_ACTIVE.1 <-brm(MOD_1_B , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                             iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                             prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_1_run_ACTIVE.1,"Models/GLOBAL/MOD_BIOM_1_run_ACTIVE.1_V2_May.Rds")

#MODEL 2
MOD_BIOM_2_run_ACTIVE.1 %>% rm()  #  
MOD_BIOM_2_run_ACTIVE.1 <-brm(MOD_2_B , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                             iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                             prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_2_run_ACTIVE.1,"Models/GLOBAL/MOD_BIOM_2_run_ACTIVE.1_V2_May.Rds")

#MODEL 3
MOD_BIOM_3_run_ACTIVE.1 %>% rm()  #  
MOD_BIOM_3_run_ACTIVE.1 <-brm(MOD_3_B , data=ACTIVE.1.sub.V2.std,cores=4,chains = 4,
                             iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                             prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#saveRDS(MOD_BIOM_3_run_ACTIVE.1,"Models/GLOBAL/MOD_BIOM_3_run_ACTIVE.1_V2_May.Rds")

#### test all B models
r2_bayes.ACTIVE.1.sub.V2.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_ACTIVE.1.null)[1],
                              r2_bayes(MOD_BIOM_1_run_ACTIVE.1)[1],
                              r2_bayes(MOD_BIOM_2_run_ACTIVE.1)[1],
                              r2_bayes(MOD_BIOM_3_run_ACTIVE.1)[1],
                              #r2_bayes(MOD_BIOM_4_run_ACTIVE.1)[1],
                              r2_bayes(MOD_B_env)[1]),
                            c(r2_bayes(MOD_BIOM_1_run_ACTIVE.1.null)[2],
                              r2_bayes(MOD_BIOM_1_run_ACTIVE.1)[2],
                              r2_bayes(MOD_BIOM_2_run_ACTIVE.1)[2],
                              r2_bayes(MOD_BIOM_3_run_ACTIVE.1)[2],
                              #r2_bayes(MOD_BIOM_4_run_ACTIVE.1)[2],
                              r2_bayes(MOD_B_env)[2]))
colnames(r2_bayes.ACTIVE.1.sub.V2.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.ACTIVE.1.sub.V2.B) <- rep("ACTIVE.1.sub.V2",2)
model.W.ACTIVE.1.sub.V2.B <- model_weights(MOD_BIOM_1_run_ACTIVE.1.null,MOD_BIOM_1_run_ACTIVE.1,
                                   MOD_B_env,weights="loo")
names(model.W.ACTIVE.1.sub.V2.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.ACTIVE.1.sub.V2.B) <- rep("ACTIVE.1.sub.V2",1)


##########################
### GLOBAL ACTIVE V2  ####
##########################----

### run species models
# null model
MOD_S_1_run_ACTIVE.2_null %>% rm()  #  
MOD_S_1_run_ACTIVE.2_null <- brm(MOD_1_S_null , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2)

# model 1
MOD_S_1_run_ACTIVE.2 %>% rm()  #  
MOD_S_1_run_ACTIVE.2 <-brm(MOD_1_S , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# model 2
MOD_S_2_run_ACTIVE.2 %>% rm()  #  
MOD_S_2_run_ACTIVE.2 <-brm(MOD_2_S , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# model 3
MOD_S_3_run_ACTIVE.2 %>% rm()  #  
MOD_S_3_run_ACTIVE.2 <-brm(MOD_3_S , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# model 4
MOD_S_4_run_ACTIVE.2 %>% rm()  #  
MOD_S_4_run_ACTIVE.2 <-brm(MOD_4_S , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# model S env
MOD_S_env %>% rm()  #  
MOD_S_env <-brm(S_mod_nocon , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#### test all S models ACTIVE.2.sub.V2
r2_bayes.ACTIVE.2.sub.V2.S <- rbind(c(r2_bayes(MOD_S_1_run_ACTIVE.2_null)[1],
                                      r2_bayes(MOD_S_1_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_S_2_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_S_3_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_S_4_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_S_env)[1]),
                                    c(r2_bayes(MOD_S_1_run_ACTIVE.2_null)[2],
                                      r2_bayes(MOD_S_1_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_S_2_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_S_3_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_S_4_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_S_env)[2]))
colnames(r2_bayes.ACTIVE.2.sub.V2.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.ACTIVE.2.sub.V2.S) <- rep("ACTIVE.2.sub.V2",2)
model.W.ACTIVE.2.sub.V2.S <- model_weights(MOD_S_1_run_ACTIVE.2_null,MOD_S_1_run_ACTIVE.2,
                                        MOD_S_env,weights="loo")
names(model.W.ACTIVE.2.sub.V2.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_ACTIVE.2, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=ACTIVE.2.sub.V2.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_ACTIVE.2_null,MOD_S_1_run_ACTIVE.2,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_ACTIVE.2)
#pp_check(MOD_S_1_run_ACTIVE.2, resp="log_biomassarea")
#brms::pp_check(MOD_S_1_run_ACTIVE.2, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_ACTIVE.2, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_ACTIVE.2_null)

### run biomass models
ACTIVE.2.sub.V2.std$log_biomassarea <- ACTIVE.2.sub.V2.std$logbiomassarea
## ACTIVE.2.sub.V2
#NULL
MOD_BIOM_1_run_ACTIVE.2.null %>% rm()  #  
MOD_BIOM_1_run_ACTIVE.2.null <-brm(MOD_B.null , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))

#MODEL 1
MOD_BIOM_1_run_ACTIVE.2 %>% rm()  #  
MOD_BIOM_1_run_ACTIVE.2 <-brm(MOD_1_B , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#MODEL 2
MOD_BIOM_2_run_ACTIVE.2 %>% rm()  #  
MOD_BIOM_2_run_ACTIVE.2 <-brm(MOD_2_B , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#MODEL 3
MOD_BIOM_3_run_ACTIVE.2 %>% rm()  #  
MOD_BIOM_3_run_ACTIVE.2 <-brm(MOD_3_B , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#Env. variables
MOD_B_env %>% rm()  #  
MOD_B_env <-brm(B_mod_nocon , data=ACTIVE.2.sub.V2.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


#### test all B models
r2_bayes.ACTIVE.2.sub.V2.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_ACTIVE.2.null)[1],
                                      r2_bayes(MOD_BIOM_1_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_BIOM_2_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_BIOM_3_run_ACTIVE.2)[1],
                                      #r2_bayes(MOD_BIOM_4_run_ACTIVE.2)[1],
                                      r2_bayes(MOD_B_env)[1]),
                                    c(r2_bayes(MOD_BIOM_1_run_ACTIVE.2.null)[2],
                                      r2_bayes(MOD_BIOM_1_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_BIOM_2_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_BIOM_3_run_ACTIVE.2)[2],
                                      #r2_bayes(MOD_BIOM_4_run_ACTIVE.2)[2],
                                      r2_bayes(MOD_B_env)[2]))
colnames(r2_bayes.ACTIVE.2.sub.V2.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.ACTIVE.2.sub.V2.B) <- rep("ACTIVE.2.sub.V2",2)
model.W.ACTIVE.2.sub.V2.B <- model_weights(MOD_BIOM_1_run_ACTIVE.2.null,MOD_BIOM_1_run_ACTIVE.2,
                                           MOD_B_env,weights="loo")
names(model.W.ACTIVE.2.sub.V2.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.ACTIVE.2.sub.V2.B) <- rep("ACTIVE.2.sub.V2",1)


##########################
### GLOBAL PASSIVE  ####
##########################----

### run species models
# null model
MOD_S_1_run_PASSIVE_null %>% rm()  #  
MOD_S_1_run_PASSIVE_null <-brm(MOD_1_S_null , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2)
# model 1
MOD_S_1_run_PASSIVE %>% rm()  #  
MOD_S_1_run_PASSIVE <-brm(MOD_1_S , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
# model 2
MOD_S_2_run_PASSIVE %>% rm()  #  
MOD_S_2_run_PASSIVE <-brm(MOD_2_S , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
# model 3
MOD_S_3_run_PASSIVE %>% rm()  #  
MOD_S_3_run_PASSIVE <-brm(MOD_3_S , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
# model 4
MOD_S_4_run_PASSIVE %>% rm()  #  
MOD_S_4_run_PASSIVE <-brm(MOD_4_S , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# model S env
MOD_S_env %>% rm()  #  
MOD_S_env <-brm(S_mod_nocon , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#### test all S models PASSIVE.sub.V2
r2_bayes.PASSIVE.sub.V2.S <- rbind(c(r2_bayes(MOD_S_1_run_PASSIVE_null)[1],
                                      r2_bayes(MOD_S_1_run_PASSIVE)[1],
                                      r2_bayes(MOD_S_2_run_PASSIVE)[1],
                                      r2_bayes(MOD_S_3_run_PASSIVE)[1],
                                      r2_bayes(MOD_S_4_run_PASSIVE)[1],
                                      r2_bayes(MOD_S_env)[1]),
                                    c(r2_bayes(MOD_S_1_run_PASSIVE_null)[2],
                                      r2_bayes(MOD_S_1_run_PASSIVE)[2],
                                      r2_bayes(MOD_S_2_run_PASSIVE)[2],
                                      r2_bayes(MOD_S_3_run_PASSIVE)[2],
                                      r2_bayes(MOD_S_4_run_PASSIVE)[2],
                                      r2_bayes(MOD_S_env)[2]))
colnames(r2_bayes.PASSIVE.sub.V2.S) <- c("S_model_null","S_model_1","S_model_2","S_model_3","S_model_4","S_model_ENV")
rownames(r2_bayes.PASSIVE.sub.V2.S) <- rep("PASSIVE.sub.V2",2)
model.W.PASSIVE.sub.V2.S <- model_weights(MOD_S_1_run_PASSIVE_null,MOD_S_1_run_PASSIVE,MOD_S_env,weights="loo")
names(model.W.PASSIVE.sub.V2.S) <- c("S_model_null","S_model_1","S_model_ENV")

#fit3 <- update(MOD_S_1_run_PASSIVE, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=PASSIVE.sub.V2.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_PASSIVE_null,MOD_S_1_run_PASSIVE,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_PASSIVE)
#pp_check(MOD_S_1_run_PASSIVE, resp="log_biomassarea")
#brms::pp_check(MOD_S_1_run_PASSIVE, resp="Richness", nsamples = 100)
#brms::pp_check(MOD_BIOM_1_run_PASSIVE, resp="Richness", nsamples = 100)
#r2_bayes(MOD_S_1_run_PASSIVE_null)

### run biomass models
PASSIVE.sub.V2.std$log_biomassarea <- PASSIVE.sub.V2.std$logbiomassarea
## PASSIVE.sub.V2
#NULL
MOD_BIOM_1_run_PASSIVE.null %>% rm()  #  
MOD_BIOM_1_run_PASSIVE.null <-brm(MOD_B.null , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))

#MODEL 1
MOD_BIOM_1_run_PASSIVE %>% rm()  #  
MOD_BIOM_1_run_PASSIVE <-brm(MOD_1_B , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#MODEL 2
MOD_BIOM_2_run_PASSIVE %>% rm()  #  
MOD_BIOM_2_run_PASSIVE <-brm(MOD_2_B , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#MODEL 3
MOD_BIOM_3_run_PASSIVE %>% rm()  #  
MOD_BIOM_3_run_PASSIVE <-brm(MOD_3_B , data=PASSIVE.sub.V2.std,cores=4,chains = 4,
                              iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                              prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

#### test all B models
r2_bayes.PASSIVE.sub.V2.B <- rbind(c(r2_bayes(MOD_BIOM_1_run_PASSIVE.null)[1],
                                      r2_bayes(MOD_BIOM_1_run_PASSIVE)[1],
                                      r2_bayes(MOD_BIOM_2_run_PASSIVE)[1],
                                      r2_bayes(MOD_BIOM_3_run_PASSIVE)[1],
                                      #r2_bayes(MOD_BIOM_4_run_PASSIVE)[1],
                                      r2_bayes(MOD_B_env)[1]),
                                    c(r2_bayes(MOD_BIOM_1_run_PASSIVE.null)[2],
                                      r2_bayes(MOD_BIOM_1_run_PASSIVE)[2],
                                      r2_bayes(MOD_BIOM_2_run_PASSIVE)[2],
                                      r2_bayes(MOD_BIOM_3_run_PASSIVE)[2],
                                      #r2_bayes(MOD_BIOM_4_run_PASSIVE)[2],
                                      r2_bayes(MOD_B_env)[2]))
colnames(r2_bayes.PASSIVE.sub.V2.B) <- c("B_model_null","B_model_1","B_model_2","B_model_3","B_model_ENV")
rownames(r2_bayes.PASSIVE.sub.V2.B) <- rep("PASSIVE.sub.V2",2)

model.W.PASSIVE.sub.V2.B <- model_weights(MOD_BIOM_1_run_PASSIVE.null,MOD_BIOM_1_run_PASSIVE,
                                           MOD_B_env,weights="loo")
names(model.W.PASSIVE.sub.V2.B) <- c("B_model_null","B_model_1","B_model_ENV")
#rownames(model.W.PASSIVE.sub.V2.B) <- rep("PASSIVE.sub.V2",1)

####################3

###Merging database##----

#biomass r2
r2_bayes.FE.B<-rbind(r2_bayes.ACTIVE.1.sub.V2.B,r2_bayes.ACTIVE.2.sub.V2.B,r2_bayes.PASSIVE.sub.V2.B,
                     r2_bayes.transient.B,r2_bayes.CRYPTIC.B, r2_bayes.PARENTAL.B,r2_bayes.RESID.B)
rm(r2_Global)
r2_Global<-cbind(r2_bayes.FE.B[c(1:14),1],
r2_bayes.FE.B[c(1:14),2],
r2_bayes.FE.B[c(1:14),3],
r2_bayes.FE.B[c(1:14),4],
r2_bayes.FE.B[c(1:14),5])
colnames(r2_Global) <- colnames(r2_bayes.FE.B)
#write.csv(r2_Global,"Models output/GLOBAL_B_R2_Jun21.csv")

#richness
r2_bayes.FE.S<-rbind(r2_bayes.ACTIVE.1.sub.V2.S,r2_bayes.ACTIVE.2.sub.V2.S,r2_bayes.PASSIVE.sub.V2.S,
                     r2_bayes.transient.S,r2_bayes.CRYPTIC.S, r2_bayes.PARENTAL.S,r2_bayes.RESID.S)
rm(r2_Global)
r2_Global<-cbind(r2_bayes.FE.S[c(1:14),1],
                 r2_bayes.FE.S[c(1:14),2],
                 r2_bayes.FE.S[c(1:14),3],
                 r2_bayes.FE.S[c(1:14),4],
                 r2_bayes.FE.S[c(1:14),5],
                 r2_bayes.FE.S[c(1:14),6])
colnames(r2_Global) <- colnames(r2_bayes.FE.S)
#write.csv(r2_Global,"Models output/GLOBAL_S_R2_Jun21.csv")


#model weights

model.W.FE.B<-rbind(model.W.ACTIVE.1.sub.V2.B,model.W.ACTIVE.2.sub.V2.B,model.W.PASSIVE.sub.V2.B,
                     model.W.transient.B,model.W.CRYPTIC.B, model.W.PARENTAL.B,model.W.RESID.B)

rm(model.W.GLobal)
model.W.GLobal<- cbind(model.W.FE.B[c(1:7),1],
                  model.W.FE.B[c(1:7),2],
                  model.W.FE.B[c(1:7),3])
colnames(model.W.GLobal) <- colnames(model.W.FE.B)
write.csv(model.W.GLobal,"Models output/GLOBAL_B_modelweight_Jun21.csv")


#richness
model.W.FE.S<-rbind(model.W.ACTIVE.1.sub.V2.S,model.W.ACTIVE.2.sub.V2.S,model.W.PASSIVE.sub.V2.S,
                    model.W.transient.S,model.W.CRYPTIC.S, model.W.PARENTAL.S,model.W.RESID.S)

rm(model.W.GLobal)
model.W.GLobal<- cbind(model.W.FE.S[c(1:7),1],
                       model.W.FE.S[c(1:7),2],
                       model.W.FE.S[c(1:7),3])
colnames(model.W.GLobal) <- colnames(model.W.FE.S)
write.csv(model.W.GLobal,"Models output/GLOBAL_S_modelweight_Jun21.csv")

###END - save outputs
