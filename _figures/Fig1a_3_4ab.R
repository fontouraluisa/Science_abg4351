# Figures - Global connectivity 
# author: Luisa Fontoura and Majambo Gamoyo
# date: April 2020
# outputs: Maps for Fig 1A and Fig4A.B; Fig 3B (biplot percentiles ), % in MPA for Table S5 and Fig 4C-E 
# obs # global connections were spatially reproduced using Gephi (input probability matrices) and overlaid to the map(s) a posteriori

#load library ----
library(raster)
library(rgdal)
library(sp)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
#devtools::install_github("ropenscilabs/rnaturalearthdata")
#install.packages("rnaturalearthhires",
#                 repos = "http://packages.ropensci.org",
#                 type = "source")
library(hrbrthemes)
library(ggplot2)
library(maps)
library(colorBlindness)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(extrafont)
#font_import()
#loadfonts(device = "pdf")
library(here)
library(stats)



# Load global connectivity data ----
ID<-read.csv(here("_data","Data_Global_Connectivity","Data_global_connectivity.csv"),h=T, stringsAsFactors = F,dec=".")
#load data required for plotting figures
ID$LonTeste<-ifelse(ID$Lon <0, ID$Lon+360,ID$Lon)
#load all.data.Rds - Plot biomass/richness site points
all.data <- readRDS(here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.Rds"))
all.data$LonTeste<-ifelse(all.data$Lon <0, all.data$Lon+360,all.data$Lon)

#Calculate netflow for all reefs
dataNet <- ID
dataNet$crypto5Netflow <- (dataNet$crypto5BROF - dataNet$crypto5BRIF)/(dataNet$crypto5BROF + dataNet$crypto5BRIF)
dataNet$pare5Netflow <- (dataNet$pare5BROF - dataNet$pare5BRIF)/(dataNet$pare5BROF + dataNet$pare5BRIF)
dataNet$resid15Netflow <- (dataNet$resid15BROF - dataNet$resid15BRIF)/(dataNet$resid15BROF + dataNet$resid15BRIF)
dataNet$transi15Netflow <- (dataNet$transi15BROF - dataNet$transi15BRIF)/(dataNet$transi15BROF + dataNet$transi15BRIF)

# rm NAs for netflow for plotting purposes
dataNet <- dataNet %>% filter(!is.na(crypto5Netflow))

#Set color for top 10% sinks and sources RETHINK COLOURS
dataNet$colorSinkSource<- ifelse(dataNet$crypto5Netflow > quantile(dataNet$crypto5Netflow,0.9),"darkorange",
                                 ifelse(dataNet$crypto5Netflow < quantile(dataNet$crypto5Netflow,0.1),"navyblue","gray25"))
SinkSourceData<-dataNet %>% filter(colorSinkSource %in% c("darkorange","navyblue"))
colourSinks<-as.character(SinkSourceData$colorSinkSource)

#Set World Map
world <- ne_coastline(scale = "large", returnclass = "sf")
world2 = maps::map(wrap=c(0,360), proj="mollweide",plot=FALSE, fill=TRUE,interior = FALSE)
class(world)


# Fig 1A Map####----
theme_set(theme_classic(base_size = 22))
t <- ggplot() +
  geom_sf() +
  borders("world2",fill="gray90", col="gray90", bg="lightblue") +
  #geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group)) +
  geom_point(data = ID,
             aes(x = LonTeste, y = Lat), color="gray15",fill="gray55",
             alpha = 0.5, size=0.75) +
  theme_bw() +
  #scale_colour_gradientn(colours = rev(c("DarkRed","orange","white","cyan3","Blue2"))) +
  #scale_fill_discrete(values=Blue2DarkRed12Steps) +
  theme(legend.title = element_blank(),legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm')) +
  labs(title = "") +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  theme(text = element_text(size = 14, family="Helvetica")) +
  coord_sf(xlim = c(0, 360), ylim = c(-40, 40), expand = FALSE) +
  scale_y_continuous(breaks = seq(-40, 40, by=20))
#Fig 1
fin<-t +
  geom_point(data = all.data, 
             aes(x = LonTeste, y = Lat), color="white",fill="#336699",
             alpha = 0.35, size=2, shape=21)

fin 

##Create transparent
fin+theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))

#ggsave("Fig1A_basemap.png",device = "png",bg = "transparent",dpi = 1000,width=10,height=5)


# Fig 4 Corridors and Source/sink maps----
t4 <- ggplot() +
  geom_sf() +
  borders("world2",fill="gray95", col="gray95") +
  #geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group)) +
  geom_point(data = ID,
             aes(x = LonTeste, y = Lat), color="#999999",fill="#CCCCCC",
             alpha = 0.8, size=0.75) +
  theme_bw() +
  #scale_colour_gradientn(colours = rev(c("DarkRed","orange","white","cyan3","Blue2"))) +
  #scale_fill_discrete(values=Blue2DarkRed12Steps) +
  theme(legend.title = element_blank(),legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm')) +
  labs(title = "") +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  theme(text = element_text(size = 14, family="Helvetica")) +
  coord_sf(xlim = c(0, 360), ylim = c(-40, 40), expand = FALSE) +
  scale_y_continuous(breaks = seq(-40, 40, by=20))
t4


#Netflow # Top 10% sinks and 10% sources dataset (SinkSourceData)
NetTop<-t +
  geom_point(data = SinkSourceData, 
             aes(x = LonTeste, y = Lat, group = colorSinkSource),
             color=SinkSourceData$colorSinkSource,
             alpha = 0.8, size=1, shape=19)

##Create transparent
NetTop+theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))

##Indegree of upstream reefs (outdegree >0) - top dispersal corridors for each biogeographical regions
quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$pare5BRin,0.9) #62
quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$pare5BRin,0.9) #60
quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$pare5BRin,0.9) #118
quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$pare5BRin,0.9) #96

InDTop<-t4 +
  geom_point(data = dataNet[dataNet$Kulbicki == "Western Indian" &
                              dataNet$pare5BRin>61,], 
             aes(x = LonTeste, y = Lat),
             color="#FFCC00",
            alpha=0.45,size=1, shape=19) + 
  geom_point(data = dataNet[dataNet$Kulbicki == "Western Atlantic" &
                              dataNet$pare5BRin>59,], 
             aes(x = LonTeste, y = Lat),
             color="#FFCC00",
             alpha=0.45,size=1, shape=19) +
  geom_point(data = dataNet[dataNet$Kulbicki == "Central Indo_Pacific" &
                              dataNet$pare5BRin>117,], 
             aes(x = LonTeste, y = Lat),
             color="#FFCC00",
             alpha=0.45,size=1, shape=19) +
  geom_point(data = dataNet[dataNet$Kulbicki == "Central Pacific" &
                             dataNet$pare5BRin>95,], 
            aes(x = LonTeste, y = Lat),
            color="#FFCC00",
            alpha=0.45,size=1, shape=19)
  

##Create transparent
InDTop+theme(
  panel.background = element_rect(fill ="transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))

#ggsave("IndegreeSources_pare5_basemap.png",device = "png",bg = "transparent",dpi = 1000,width=10,height=5)


#colour lightblue bg
bgworld <- ggplot() +
  geom_sf() +
  borders("world2",fill="gray95", col="gray95") +
  theme_bw() +
  #scale_colour_gradientn(colours = rev(c("DarkRed","orange","white","cyan3","Blue2"))) +
  #scale_fill_discrete(values=Blue2DarkRed12Steps) +
  theme(legend.title = element_blank(),legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm')) +
  labs(title = "") +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  theme(text = element_text(size = 14, family="Helvetica")) +
  coord_sf(xlim = c(0, 360), ylim = c(-40, 40), expand = FALSE) +
  scale_y_continuous(breaks = seq(-40, 40, by=20)) +
theme(
  panel.background = element_rect(fill ="#D9ECF6",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) 

#ggsave("blue_bg_basemap.png",device = "png",bg = "transparent",dpi = 1000,width=10,height=5)


#ggsave("Netflow_crypto5_basemap.png",type = "png",bg = "transparent",dpi = 1000,width=10,height=5)
#ggsave("netflowmap.pdf",dpi = 600,width=12,height=7)
#ggsave("nflowmap.png",type = "png",dpi = 1000,width=10,height=5)





# Fig 3A,B - conceptual framework and netflow/indegree biplot ----
##obs - we applied connectivity values for crypto5 and pare5 models - 
#models' choice was made based on the strongest predictors of biomass and richness, respectively (refer to Figure 1B,C).
head(dataNet)
#filter 

#ggplot() +geom_point(data = dataNet,
#                     aes(y=globalNetflow, x=log(globalBRin), color=General), alpha=0.5) +
#  scale_colour_manual(values = c("seashell","#1F968BFF")) + theme_minimal()


##The stat for adding vetical lines at pecentile locations 
StatPercentileX <- ggproto("StatPercentileX", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$x, probs=probs)
                             data.frame(xintercept=percentiles)
                           },
                           required_aes = c("x")
)

stat_percentile_x <- function(mapping = NULL, data = NULL, geom = "vline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileX, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileXLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$x, probs=probs)
                                   data.frame(x=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           round(percentiles, digits=3)))
                                 },
                                 required_aes = c("x")
)

stat_percentile_xlab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
##The stat for adding vetical lines at pecentile locations 
StatPercentileY <- ggproto("StatPercentileY", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$y, probs=probs)
                             data.frame(yintercept=percentiles)
                           },
                           required_aes = c("y")
)

stat_percentile_y <- function(mapping = NULL, data = NULL, geom = "hline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileY, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileYLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$y, probs=probs)
                                   data.frame(y=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           round(percentiles, digits=3)))
                                 },
                                 required_aes = c("y")
)

stat_percentile_ylab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

percentile <- ecdf(dataNet$crypto5Netflow)
dataNet$crypto5netflow_percentile<-percentile(dataNet$crypto5Netflow)
percentile <- ecdf(dataNet$pare5BRin)
dataNet$pare5BRin_percentile<-percentile(dataNet$pare5BRin)


#All coral reefs 
allRe<- ggplot(dataNet, aes(y=crypto5netflow_percentile, x=pare5BRin_percentile)) +
  geom_point(data=dataNet, 
             aes(y=crypto5netflow_percentile, x=pare5BRin_percentile),alpha=0.45, shape=21, color="gray50") + 
  geom_point(data=dataNet %>% filter(crypto5netflow_percentile > 0.9  & General == "MPA"), 
             aes(y=crypto5netflow_percentile, x=pare5BRin_percentile),alpha=0.45, shape=21, fill="#339966") + 
  geom_point(data=dataNet %>% filter(crypto5netflow_percentile < 0.1 & General == "MPA"), 
             aes(y=crypto5netflow_percentile, x=pare5BRin_percentile),alpha=0.45, shape=21,fill="#339966") + 
  geom_point(data=dataNet %>% filter(pare5BRin_percentile > 0.9 & pare5BROF > 1 & General == "MPA"), 
             aes(y=crypto5netflow_percentile, x=pare5BRin_percentile),alpha=0.45, shape=21,fill="#339966") + 
  theme_bw() +
  stat_percentile_x(probs=c(0.9), linetype=2) +
  stat_percentile_xlab(probs=c(0.9), hjust=1, vjust=1.5, angle=180) +
  stat_percentile_y(probs=c(0.1,0.90), linetype=2) +
  stat_percentile_ylab(probs=c(0.1,0.90), hjust=1, vjust=1.5, angle=180) +
  scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent) + 
  labs(y="Net larval flow",x="Number of incoming connections")

allRe

allRe_CF<- ggplot(dataNet, aes(y=crypto5netflow_percentile, x=pare5BRin_percentile)) +
  #geom_point(data=dataNet, 
  #           aes(y=crypto5netflow_percentile, x=pare5BRin_percentile),alpha=0, color="gray50") + 
  theme_bw()  +
  stat_percentile_x(probs=c(0.9), linetype=2) +
  stat_percentile_xlab(probs=c(0.9), hjust=1, vjust=1.5, angle=180) +
  stat_percentile_y(probs=c(0.1,0.90), linetype=2) +
  stat_percentile_ylab(probs=c(0.1,0.90), hjust=1, vjust=1.5, angle=180) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) + 
  labs(y="Net larval flow",x="Number of incoming connections")

relatinnetfin <- ggarrange(allRe_CF, allRe, nrow=1,common.legend = T,legend="right",labels = c("A","B"))
relatinnetfin




#Proportions sources, sinks and corridors in MPA for all 4 groups (applied on Fig 4 and Table S5) ----
#Table S5, Fig.3 & 4 were built based on crypto5 and pare5 
#(these models were chosen based on their higher predictive power)

rm(percPro)
percPro<- as.data.frame(c("Sources","Sinks","Dispersal Corridors"))
colnames(percPro) <- c("ReefType")

propData <- dataNet
propData <- propData%>% filter(!is.na(crypto5Netflow))

#Crypto 5
percentile <- ecdf(propData$crypto5Netflow)
propData$crypto5netflow_percentile<-percentile(propData$crypto5Netflow)
percentile <- ecdf(propData$crypto5BRin)
propData$crypto5BRin_percentile<-percentile(propData$crypto5BRin)

#thresholds within regions given large disparities in indegree values 
corWIO<- quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$crypto5BRin,0.9,na.rm=T) #126
corWA<-quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$crypto5BRin,0.9,na.rm=T) #101
corCIP<- quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$crypto5BRin,0.9,na.rm=T) #214
corCP<-quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$crypto5BRin,0.9,na.rm=T) #163

percPro$globalcrypto5 <- c(length(propData[propData$crypto5netflow_percentile  >  0.9 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile > 0.9,]$ID)*100,
                           length(propData[propData$crypto5netflow_percentile < 0.1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile < 0.1,]$ID)*100,
                           length(propData[propData$crypto5BRin_percentile > 0.9 & propData$crypto5BROF > 1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$crypto5BRin_percentile > 0.9 & propData$crypto5BROF > 10,]$ID)*100)

percPro$WAcrypto5 <- c(length(propData[propData$crypto5netflow_percentile  >  0.9 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile > 0.9 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #sinks
                       length(propData[propData$crypto5netflow_percentile < 0.1 &propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile < 0.1 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #corridors
                       length(propData[propData$crypto5BRin > corWA & propData$crypto5BROF > 1 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$crypto5BRin > corWA & propData$crypto5BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)*100)

percPro$WIOcrypto5 <- c(length(propData[propData$crypto5netflow_percentile  >  0.9 & propData$Kulbicki == "Western Indian" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile > 0.9 & propData$Kulbicki == "Western Indian",]$ID)*100,
                       #sinks
                       length(propData[propData$crypto5netflow_percentile < 0.1 &propData$Kulbicki == "Western Indian" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile < 0.1 & propData$Kulbicki == "Western Indian",]$ID)*100,
                       #corridors
                       length(propData[propData$crypto5BRin > corWIO & propData$crypto5BROF > 1 & propData$Kulbicki == "Western Indian" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$crypto5BRin > corWIO & propData$crypto5BROF > 1 & propData$Kulbicki == "Western Indian",]$ID)*100)

percPro$CPcrypto5 <- c(length(propData[propData$crypto5netflow_percentile  >  0.9 & propData$Kulbicki == "Central Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile > 0.9 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                        #sinks
                        length(propData[propData$crypto5netflow_percentile < 0.1 &propData$Kulbicki == "Central Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile < 0.1 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                        #corridors
                        length(propData[propData$crypto5BRin > corCP & propData$crypto5BROF > 1 & propData$Kulbicki == "Central Pacific" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$crypto5BRin >corCP & propData$crypto5BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID)*100)

percPro$CIPcrypto5 <- c(length(propData[propData$crypto5netflow_percentile  >  0.9 & propData$Kulbicki == "Central Indo_Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile > 0.9 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                       #sinks
                       length(propData[propData$crypto5netflow_percentile < 0.1 &propData$Kulbicki == "Central Indo_Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$crypto5netflow_percentile < 0.1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                       #corridors
                       length(propData[propData$crypto5BRin > corCIP & propData$crypto5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$crypto5BRin > corCIP & propData$crypto5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100)

#Parental 5
percentile <- ecdf(propData$pare5Netflow)
propData$pare5netflow_percentile<-percentile(propData$pare5Netflow)
percentile <- ecdf(propData$pare5BRin)
propData$pare5BRin_percentile<-percentile(propData$pare5BRin)

#thresholds within regions given large disparities in indegree values 
corWIO<- quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$pare5BRin,0.9,na.rm=T)
corWA<-quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$pare5BRin,0.9,na.rm=T)
corCIP<- quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$pare5BRin,0.9,na.rm=T)
corCP<-quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$pare5BRin,0.9,na.rm=T)

percPro$globalpare5 <- c(length(propData[propData$pare5netflow_percentile  >  0.9 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile > 0.9,]$ID)*100,
                           length(propData[propData$pare5netflow_percentile < 0.1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile < 0.1,]$ID)*100,
                           length(propData[propData$pare5BRin_percentile > 0.9 & propData$pare5BROF > 1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$pare5BRin_percentile > 0.9 & propData$pare5BROF > 10,]$ID)*100)

percPro$WApare5 <- c(length(propData[propData$pare5netflow_percentile  >  0.9 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile > 0.9 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #sinks
                       length(propData[propData$pare5netflow_percentile < 0.1 &propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile < 0.1 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #corridors
                       length(propData[propData$pare5BRin > corWA & propData$pare5BROF > 1 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$pare5BRin > corWA & propData$pare5BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)*100)

percPro$WIOpare5 <- c(length(propData[propData$pare5netflow_percentile  >  0.9 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile > 0.9 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #sinks
                        length(propData[propData$pare5netflow_percentile < 0.1 &propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile < 0.1 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #corridors
                        length(propData[propData$pare5BRin > corWIO & propData$pare5BROF > 1 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$pare5BRin > corWIO & propData$pare5BROF > 1 & propData$Kulbicki == "Western Indian",]$ID)*100)

percPro$CPpare5 <- c(length(propData[propData$pare5netflow_percentile  >  0.9 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile > 0.9 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #sinks
                       length(propData[propData$pare5netflow_percentile < 0.1 &propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile < 0.1 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #corridors
                       length(propData[propData$pare5BRin > corCP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$pare5BRin >corCP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID)*100)

percPro$CIPpare5 <- c(length(propData[propData$pare5netflow_percentile  >  0.9 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile > 0.9 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #sinks
                        length(propData[propData$pare5netflow_percentile < 0.1 &propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$pare5netflow_percentile < 0.1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #corridors
                        length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100)

#Resident 15
percentile <- ecdf(propData$resid15Netflow)
propData$resid15netflow_percentile<-percentile(propData$resid15Netflow)
percentile <- ecdf(propData$resid15BRin)
propData$resid15BRin_percentile<-percentile(propData$resid15BRin)

#thresholds within regions given large disparities in indegree values 
corWIO<- quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$resid15BRin,0.9,na.rm=T) #126
corWA<-quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$resid15BRin,0.9,na.rm=T) #101
corCIP<- quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$resid15BRin,0.9,na.rm=T) #214
corCP<-quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$resid15BRin,0.9,na.rm=T) #163

percPro$globalresid15 <- c(length(propData[propData$resid15netflow_percentile  >  0.9 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile > 0.9,]$ID)*100,
                           length(propData[propData$resid15netflow_percentile < 0.1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile < 0.1,]$ID)*100,
                           length(propData[propData$resid15BRin_percentile > 0.9 & propData$resid15BROF > 1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$resid15BRin_percentile > 0.9 & propData$resid15BROF > 10,]$ID)*100)

percPro$WAresid15 <- c(length(propData[propData$resid15netflow_percentile  >  0.9 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile > 0.9 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #sinks
                       length(propData[propData$resid15netflow_percentile < 0.1 &propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile < 0.1 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #corridors
                       length(propData[propData$resid15BRin > corWA & propData$resid15BROF > 1 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$resid15BRin > corWA & propData$resid15BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)*100)

percPro$WIOresid15 <- c(length(propData[propData$resid15netflow_percentile  >  0.9 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile > 0.9 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #sinks
                        length(propData[propData$resid15netflow_percentile < 0.1 &propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile < 0.1 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #corridors
                        length(propData[propData$resid15BRin > corWIO & propData$resid15BROF > 1 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$resid15BRin > corWIO & propData$resid15BROF > 1 & propData$Kulbicki == "Western Indian",]$ID)*100)

percPro$CPresid15 <- c(length(propData[propData$resid15netflow_percentile  >  0.9 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile > 0.9 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #sinks
                       length(propData[propData$resid15netflow_percentile < 0.1 &propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile < 0.1 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #corridors
                       length(propData[propData$resid15BRin > corCP & propData$resid15BROF > 1 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$resid15BRin >corCP & propData$resid15BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID)*100)

percPro$CIPresid15 <- c(length(propData[propData$resid15netflow_percentile  >  0.9 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile > 0.9 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #sinks
                        length(propData[propData$resid15netflow_percentile < 0.1 &propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$resid15netflow_percentile < 0.1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #corridors
                        length(propData[propData$resid15BRin > corCIP & propData$resid15BROF > 1 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$resid15BRin > corCIP & propData$resid15BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100)
#Transient 15
percentile <- ecdf(propData$transi15Netflow)
propData$transi15netflow_percentile<-percentile(propData$transi15Netflow)
percentile <- ecdf(propData$transi15BRin)
propData$transi15BRin_percentile<-percentile(propData$transi15BRin)

#thresholds within regions given large disparities in indegree values 
corWIO<- quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$transi15BRin,0.9,na.rm=T) #126
corWA<-quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$transi15BRin,0.9,na.rm=T) #101
corCIP<- quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$transi15BRin,0.9,na.rm=T) #214
corCP<-quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$transi15BRin,0.9,na.rm=T) #163

percPro$globaltransi15 <- c(length(propData[propData$transi15netflow_percentile  >  0.9 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile > 0.9,]$ID)*100,
                           length(propData[propData$transi15netflow_percentile < 0.1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile < 0.1,]$ID)*100,
                           length(propData[propData$transi15BRin_percentile > 0.9 & propData$transi15BROF > 1 &
                                             propData$General == "MPA",]$ID)/ length(propData[propData$transi15BRin_percentile > 0.9 & propData$transi15BROF > 10,]$ID)*100)

percPro$WAtransi15 <- c(length(propData[propData$transi15netflow_percentile  >  0.9 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile > 0.9 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #sinks
                       length(propData[propData$transi15netflow_percentile < 0.1 &propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile < 0.1 & propData$Kulbicki == "Western Atlantic",]$ID)*100,
                       #corridors
                       length(propData[propData$transi15BRin > corWA & propData$transi15BROF > 1 & propData$Kulbicki == "Western Atlantic" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$transi15BRin > corWA & propData$transi15BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)*100)

percPro$WIOtransi15 <- c(length(propData[propData$transi15netflow_percentile  >  0.9 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile > 0.9 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #sinks
                        length(propData[propData$transi15netflow_percentile < 0.1 &propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile < 0.1 & propData$Kulbicki == "Western Indian",]$ID)*100,
                        #corridors
                        length(propData[propData$transi15BRin > corWIO & propData$transi15BROF > 1 & propData$Kulbicki == "Western Indian" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$transi15BRin > corWIO & propData$transi15BROF > 1 & propData$Kulbicki == "Western Indian",]$ID)*100)

percPro$CPtransi15 <- c(length(propData[propData$transi15netflow_percentile  >  0.9 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile > 0.9 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #sinks
                       length(propData[propData$transi15netflow_percentile < 0.1 &propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile < 0.1 & propData$Kulbicki == "Central Pacific",]$ID)*100,
                       #corridors
                       length(propData[propData$transi15BRin > corCP & propData$transi15BROF > 1 & propData$Kulbicki == "Central Pacific" &
                                         propData$General == "MPA",]$ID)/length(propData[propData$transi15BRin >corCP & propData$transi15BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID)*100)

percPro$CIPtransi15 <- c(length(propData[propData$transi15netflow_percentile  >  0.9 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile > 0.9 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #sinks
                        length(propData[propData$transi15netflow_percentile < 0.1 &propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/ length(propData[propData$transi15netflow_percentile < 0.1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100,
                        #corridors
                        length(propData[propData$transi15BRin > corCIP & propData$transi15BROF > 1 & propData$Kulbicki == "Central Indo_Pacific" &
                                          propData$General == "MPA",]$ID)/length(propData[propData$transi15BRin > corCIP & propData$transi15BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100)

percPro

percPro_TableSM <- gather(percPro,key="Fish group",value= "percentage",2:21) 
#write.csv(percPro_TableSM,"percPro.csv")





#Proportions of functionally important reefs across bio regions (applied on Fig.4C and Table S5)----
#sources
length(propData[propData$crypto5netflow_percentile > 0.9  & propData$Kulbicki == "Central Indo_Pacific",]$ID)/
    length(propData[propData$crypto5netflow_percentile > 0.9,]$ID) #43.6%
    
length(propData[propData$crypto5netflow_percentile > 0.9  & propData$Kulbicki == "Central Pacific",]$ID)/
  length(propData[propData$crypto5netflow_percentile > 0.9,]$ID) #29.9

length(propData[propData$crypto5netflow_percentile > 0.9  & propData$Kulbicki == "Western Atlantic",]$ID)/
  length(propData[propData$crypto5netflow_percentile > 0.9,]$ID) #15.9

length(propData[propData$crypto5netflow_percentile > 0.9  & propData$Kulbicki == "Western Indian",]$ID)/
  length(propData[propData$crypto5netflow_percentile > 0.9,]$ID) #9.8%

#sinks
length(propData[propData$crypto5netflow_percentile < 0.1  & propData$Kulbicki == "Central Indo_Pacific",]$ID)/
  length(propData[propData$crypto5netflow_percentile < 0.1,]$ID) #39.6

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$Kulbicki == "Central Pacific",]$ID)/
  length(propData[propData$crypto5netflow_percentile < 0.1,]$ID) #35.4

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$Kulbicki == "Western Indian",]$ID)/
  length(propData[propData$crypto5netflow_percentile < 0.1,]$ID) #11.8

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$Kulbicki == "Western Atlantic",]$ID)/
  length(propData[propData$crypto5netflow_percentile < 0.1,]$ID) #12.2


#corridors relative to each biogeo threshold (given the disparities on number of reefs)

allDP<- sum(c(length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID),
length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID),
length(propData[propData$pare5BRin > corWIO & propData$pare5BROF > 1 & propData$Kulbicki == "Western Indian",]$ID),
length(propData[propData$pare5BRin > corWA & propData$pare5BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)))

length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific",]$ID)/allDP #54.7
length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Pacific",]$ID)/allDP #20.8
length(propData[propData$pare5BRin > corWIO & propData$pare5BROF > 1 & propData$Kulbicki == "Western Indian",]$ID)/allDP #13.8%
length(propData[propData$pare5BRin > corWA & propData$pare5BROF > 1 & propData$Kulbicki == "Western Atlantic",]$ID)/allDP #10.6%


#Representativeness within MPAs
#sources
length(propData[propData$crypto5netflow_percentile > 0.9  & propData$General %in% "MPA" & propData$Kulbicki == "Central Indo_Pacific",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100

length(propData[propData$crypto5netflow_percentile > 0.9  & propData$General %in% "MPA" & propData$Kulbicki == "Central Pacific",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Pacific",]$ID)*100

length(propData[propData$crypto5netflow_percentile > 0.9  & propData$General %in% "MPA" & propData$Kulbicki == "Western Indian",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Indian",]$ID)*100

length(propData[propData$crypto5netflow_percentile > 0.9  & propData$General %in% "MPA" & propData$Kulbicki == "Western Atlantic",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Atlantic",]$ID)*100 

#sinks
length(propData[propData$crypto5netflow_percentile < 0.1  & propData$General %in% "MPA" & propData$Kulbicki == "Central Indo_Pacific",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$General %in% "MPA" & propData$Kulbicki == "Central Pacific",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Pacific",]$ID)*100

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$General %in% "MPA" & propData$Kulbicki == "Western Indian",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Indian",]$ID)*100

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$General %in% "MPA" & propData$Kulbicki == "Western Atlantic",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Atlantic",]$ID)*100 


#corridors 
length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Indo_Pacific" & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Indo_Pacific",]$ID)*100

length(propData[propData$pare5BRin > corCIP & propData$pare5BROF > 1 & propData$Kulbicki == "Central Pacific" & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Central Pacific",]$ID)*100

length(propData[propData$pare5BRin > corWIO & propData$pare5BROF > 1 & propData$Kulbicki == "Western Indian" & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Indian",]$ID)*100

length(propData[propData$pare5BRin > corWA & propData$pare5BROF > 1 & propData$Kulbicki == "Western Atlantic" & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA" & propData$Kulbicki == "Western Atlantic",]$ID)*100 


#Global
length(propData[propData$crypto5netflow_percentile > 0.9 & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA",]$ID)*100

length(propData[propData$crypto5netflow_percentile < 0.1  & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA",]$ID)*100

length(propData[propData$pare5BRin_percentile > 0.9 & propData$pare5BROF > 1 & propData$General %in% "MPA",]$ID)/
  length(propData[propData$General %in% "MPA",]$ID)*100

