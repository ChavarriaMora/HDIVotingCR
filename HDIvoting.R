#Final project for time series/spatial statistics class

library(ggplot2)
library(sp) #this is better than sf, more up to date and uses an R style shp file
library(tigris) #to join the spatial data with the dataset
library(spdep) #can be used to create spatial weight matrices
library(gridExtra)
library(spatialreg)
library(dplyr)

setwd("C:/Elias/1 Serious/Academia/University of Pittsburgh/3rd Year/Time Series/FinalPaper")

CRmap<-readRDS("gadm36_CRI_2_sf.rds")  #rds is a type of file for R
#It is better to use this type of file than the shapefile. 
head(CRmap) #Variable para el nombre de canton NAME_2


CRdata<-read.csv("CantonesCr.csv", encoding = "UTF-8") #encoding importante por tildes
CRdata<-as.data.frame(CRdata)
head(CRdata)

#nest is joining the dataset and the shapefile
CRFull<-geo_join(
  spatial_data=CRmap,
  data_frame=CRdata,
  by_sp="NAME_2",
  by_df="X.U.FEFF.Canton",
  how="inner"
)
head(CRFull)
#OJO: hay unos huecos en los cantones, el mapa tiene 76, no pegó 5


#make a map, for example, human development index for 2018
HDIc<-ggplot(CRFull, aes(fill=idh_2018))+
               geom_sf()+
               ggtitle("IDH de CR, por cantones, 2018")+
               theme(plot.title = element_text(hjust=.5))+
               theme(plot.title=element_text(size=24, face="bold"))
HDIc

PACMAP<-ggplot(CRFull, aes(fill=porc_pac_2))+
  geom_sf()+
  theme(plot.title = element_text(hjust=.5))+
  theme(plot.title=element_text(size=24, face="bold"))+
  scale_fill_gradient(low = "yellow", high = "gold", name = "PAC support")
PACMAP

PRNMAP<-ggplot(CRFull, aes(fill=porc_prn_2))+
  geom_sf()+
  theme(plot.title = element_text(hjust=.5))+
  theme(plot.title=element_text(size=24, face="bold"))+
  scale_fill_gradient(low = "lightskyblue1", high = "navy", name = "PRN support")
PRNMAP

grid.arrange(PACMAP, PRNMAP, ncol = 2)

#Create a weight matrix
wneig<- poly2nb(CRmap)            #create a list of neighbors based on shapefiles
w<-nb2listw(wneig, style="W") #Zero policy calculates the matrix even if there are islands
class(wneig)

#check for spatial clustering in Y, use Moran's I
#Remember: Moran's I H0 is that the distribution is random, Ha is that it is clustered. 
#you have to interpret based on the z-score and p-value: 
#positive z-score: clustering
#negative z-score: competitive clustering, high values repel high values

#Voting for PAC
moran.test(CRFull$porc_pac_2, w) #alternative hypotesis is greater: there is clustering. Moran I= 0.59

olspac<-lm(porc_pac_2~idh_2018, data=CRFull)
moran.test((resid(olspac)), w) #alternative hypothesis is greater, adding covariates does not eliminate clustering

lmtestspac<-lm.LMtests(olspac, w, test=c("LMerr","LMlag", "RLMerr","RLMlag")) #lagrange multiplier test
summary(lmtestspac) #NONE is significant, use SLX model

#Voting for PRN
moran.test(CRFull$porc_prn_2, w) #alternative hypothesis is greater, Moran's I=0.41
olsprn<-lm(porc_prn_2~idh_2018, data=CRFull) 
moran.test((resid(olsprn)), w) #alternative hypothesis is greater
lmtestsprn<-lm.LMtests(olsprn, w, test=c("LMerr","LMlag", "RLMerr","RLMlag")) #lagrange multiplier test
summary(lmtestsprn) #RLMlag is significant, use SAR model


class(CRMat)
class(wmat)

CRdata2<- select(CRdata,-1)
  
CRMat<-as.matrix(CRdata2) #Al fin funciona, pero tiene una fila que deno borrar
CRMat<-CRMat[-82,] 
wmat<-nb2mat(wneig, style = "W")
SLX<-wmat%*%CRMat


#SAR model for PAC
lagpac<-lagsarlm(porc_pac_2~idh_2018, data=CRFull,  listw=w, method = "eigen", zero.policy=T, tol.solve = 1e-11)
summary(lagpac)

impacts.sarlm(lagpac, data=CRFull, listw=w, method="eigen", zero.policy=T, tol.sove=1e-11)

#SLX model for PAC
lagxprn<-lm(porc_prn_2~idh_2018+SLX[,1], data=CRFull)
summary(lagxprn)