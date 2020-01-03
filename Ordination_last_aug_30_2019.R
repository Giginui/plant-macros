rm(list=ls())
# setwd("/Users/Gastorga/Dropbox/chapter2_25_may_2025/data_second")
setwd("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros")
# setwd("/Users/giselleastorga/Documents/Thesis_GA/second_chapter/data_second")
getwd()

library(vegan)
library(ggplot2)
library(MASS)
library(paleoMAS)
library(analogue)
library(rioja)
library(moments)
library(stats)
library(nlme)
# library(MVA)
library(ggplot2)
library(rgl)

#----------------------------DCA and PCA on species data-------------------------------------------#
macros5 <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/combined_macros_shorted(5).csv")
# macros5 <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/combined_macros_shorted(5).csv")
# macros_shorted_5 only includes species higher than 5 %

macros.pa <- decostand(macros5, method = "pa")
macros.per <- (decostand(macros5, method = "total")*100)

rowtotals <- rowSums(macros.per)
columntotals <- colSums(macros.per)

str(macros5)
summary(macros5)
macros.log <- log1p(macros5)# computes log(1+x)
macros.sqrt <- sqrt(macros5)

dca.raw <- decorana(macros5, iweigh = 1, ira = 0)
dca.raw
dca.log <- decorana(macros.log, iweigh = 1, ira = 0)
dca.log
dca.sqrt <- decorana(macros.sqrt, iweigh = 1, ira = 0)
dca.sqrt
# DCA performed on log(raw count data). 80 samples 52 species including individual conifer abundances
# argument ira = 0 for detrending and iweight = 1 for downweighting of rare species
# DCA2 2.9 DCA2 2.2
# If raw abundances are used DCA2 = 3.7; DCA2 = 2.2

dca.dob1 <- decorana(macros.log, iweigh = 1, ira = 0)
ca.dob <- cca(macros.log)
ca.dob$tot.chi

(dca.sqrt$evals / 2.69518)*100

dca1.samplescores <- scores(dca.dob, display = c("sites"), choices = 1)
is.data.frame(dca1.samplescores)
dca1.samplescores <- as.data.frame(dca1.samplescores)

#-------------------------------------Interpolate---------------------------------------------------#
# a try with antartic data with AICC2022
dobson <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/dobson.csv")
# dobson <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/dobson.csv")

head(dobson)
dob.depth <- dobson$depth
dob.age <- dobson$age

# All the records are presented in the AICC2022 chronology
# Composite CO2 curve first row correspond to EDML core from Monin et al 2004; 
# Siegenthaier et al 2005 (in antarctica2025co2.xls) 
# The following ages are from EDC CO2 composite curve of Lüthi et al., 2008 transfered on 
# AICC2022 (icecores_data_on_AICC2022.xlsx). I could use the whole chronology proposed for EDC in 
# 2025 version of the AICC2022 chronology, but we will see the fit with the composite 
# chronology is better.
# antartic <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Antartic_data.csv")
antartic <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Antartic_data.csv")
head(antartic)
age.dome <- antartic$age_dome
temp.dome1 <- antartic$temp_dome
temp.dome2 <- antartic$temp_dome
dome.temp <- cbind(temp.dome1, temp.dome2)
dome.interp <-
  interp.dataset(
    y = dome.temp,
    x = age.dome,
    xout = dob.age,
    rep.negt = F
  )
head(dome.interp)
dome.interp <- subset(dome.interp, select = -c(2))
head(dome.interp)
dome.interp <- as.data.frame(dome.interp)
dome.stand <- decostand(dome.interp, "standardize")
dome.interp <-dome.interp$temp.dome1
dome.stand <- dome.stand$temp.dome1

age.vos <- antartic$age_vos
temp.vos1 <- antartic$deltaTS
temp.vos2 <- antartic$deltaTS
vos.temp <- cbind(temp.vos1, temp.vos2)
vos.interp <-
  interp.dataset(
    y = vos.temp,
    x = age.vos,
    xout = dob.age,
    rep.negt = F
  )
head(vos.interp)
vos.interp <- subset(vos.interp, select = -c(2))
head(vos.interp)
vos.interp <- as.data.frame(vos.interp)
vos.stand <- decostand(vos.interp, "standardize")
vos.interp <- vos.interp$temp.vos1
vos.stand <- vos.stand$temp.vos1

age.co2 <- antartic$age_composite
co21 <- antartic$CO2_composite
co22 <- antartic$CO2_composite
co2.dome <- cbind(co21, co22)
co2.interp <-
    interp.dataset(
        y = co2.dome,
        x = age.co2,
        xout = dob.age,
        rep.negt = F
    )
head(co2.interp)
co2.interp <- subset(co2.interp, select = -c(2))
co2.interp <- as.data.frame(co2.interp)
head(co2.interp)
co2.stand <- decostand(co2.interp, method = "standardize")
co2.interp <- co2.interp$co21
co2.stand <- co2.stand$co21

#---------------------------------Eagle Tarn interpolation-----------------------------------------#
# interpolated to Dobson ages
# Eagle <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Eagle_TWARM.csv", header=TRUE, sep=";")
Eagle <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Eagle.csv")
head(Eagle)
age.et <- Eagle$age
temp.et1 <- Eagle$TWARM
temp.et2 <- Eagle$TWARM
temp.et <- cbind(temp.et1, temp.et2)
et.interp <- interp.dataset(y = temp.et, x = age.et, xout = dob.age)
head(et.interp)
et.interp <- subset(et.interp, select = -c(2))
head(et.interp)
et.interp <- as.data.frame(et.interp)
et.stand <- decostand(et.interp, method = "standardize")
et.interp <- et.interp$temp.et1
et.stand <- et.stand$temp.et1

#------------------------------Lake Dobson Chironomids interpolation-------------------------------#
# I should not interpolated to the LD age, because chironomids are from the same core as the macros
# In the original data from A. Rees depth 945 had an age of 24879, but I calculate my own 
# age-depth model. I could interpolated to the LD depths of plant macros
# chiro_dob <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Chiro_Dob.csv", sep = ",")
chiro_dob <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Chiro_Dob.csv", sep = ";")
head(chiro_dob)
age.chiro <- chiro_dob$age
depth.chiro <- chiro_dob$depth
dca1.dob <- chiro_dob$dca1.chiro
dca2.dob <- chiro_dob$dca1.chiro
dca.dob <- cbind(dca1.dob, dca2.dob)
chiro.interp <- interp.dataset(y = dca.dob , x = age.chiro, xout = dob.age)
head(chiro.interp)
chiro.interp <- subset(chiro.interp, select = -c(2))
head(chiro.interp)
chiro.interp <- as.data.frame(chiro.interp)
chiro.stand <- decostand(chiro.interp, method = "standardize")
chiro.interp <- chiro.interp$dca1.dob
chiro.stand <- chiro.stand$dca1.dob

#---------------------------------------NZ isotopes interpolation----------------------------------#
# Speleothem data from NW South Island NZ Williams 2005
# Updated 20/03/2008
# NZW <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/NZW_both_isotopes.csv")
NZW <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/NZW_both_isotopes.csv")
head(NZW)
age.NZW = NZW$age
O18.nzw = NZW$nzw_18O
C13.nzw = NZW$nzw_13C
iso.NZW <- cbind(O18.nzw, C13.nzw)

NZW.interp <-
    interp.dataset(
        y = iso.NZW,
        x = age.NZW,
        xout = dob.age,
        rep.negt = F
    )
head(NZW.interp)
NZW.interp <-as.data.frame(NZW.interp)
O18.nzw.interp <- NZW.interp$O18.nzw
C13.nzw.interp <- NZW.interp$C13.nzw
NZW.stand <- decostand(NZW.interp, method = "standardize")
O18.nzw.stand <- NZW.stand$O18.nzw
C13.nzw.stand <- NZW.stand$C13.nzw

# Oxygen isotope variation in Mt. Arthur speleothems primarily represents changes in 
# meteoric waters falling above the caves, possibly responding to latitudinal changes
# in the position of the Subtropical Front in the Tasman Sea. How is the STF in Tasmania compared
# to the NZ one? 
# Carbon isotope variations in the speleothems record, represent changes in forest productivity,
# closely matching existing paleovegetation records in NZ according to Hellstrom 2998.
# NZA <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/NZ.arthur.csv")
NZA <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/NZ.arthur.csv")
head(NZA)
age.nza = NZA$age
O18.nza = NZA$d18O
C13.nza = NZA$d13C
iso.nza <- cbind(O18.nza, C13.nza)
NZA.interp <-
    interp.dataset(
        y = iso.nza,
        x = age.nza,
        xout = dob.age,
        rep.negt = F
    )
head(NZA.interp)
NZA.interp <- as.data.frame(NZA.interp)
O18.nza.interp <- NZA.interp$O18.nza
C13.nza.interp <- NZA.interp$C13.nza
NZA.stand <- decostand(NZA.interp, method = "standardize")
O18.nza.stand <- NZA.stand$O18.nza
C13.nza.stand <- NZA.stand$C13.nza

# ----------------------------------------charcoal interpolation------------------------------------#
# Charcoal data is from the same core as tha plant macros so the 2 records should have same age.
# However, counts of charcoal particles were made every 2 cm starting at 0.5 cm and plant macros 
# starting at 2 cm and counted every 20 cm. I interpoletated to LD depth.
# char <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/charcoal.csv", sep =  ";")
char <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/charcoal.csv", sep = ";")
head(char)
age.char <- char$age
age.char
depth.char <- char$depth
depth.char
char1 = char$influx
char2 = char$influx
chars <- cbind(char1, char2)
chars <- subset(chars, select = -c(2))
char.interp <-
    interp.dataset(
        y = chars,
        x = age.char,
        xout = dob.age,
        rep.negt = F
    )
head(char.interp)
char.interp <- as.data.frame(char.interp)
char.stand <- decostand(char.interp, method = "standardize")
char.interp <- char.interp$char1
char.stand <-char.stand$char1

#-------------------------------Loss on Ignision interpolation-------------------------------------#
# loi data interpolated to plant macros depth
# loi <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/loi.csv", sep = ";")
loi <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/loi.csv",sep = ";")
head(loi)
depth.loi <- loi$depth
age.loi <- loi$age
loi1 = loi$LOI
loi2 = loi$LOI
lois <- cbind(loi1, loi2)
head(lois)
lois <- subset(lois, select = -c(2))
loi.interp <-
    interp.dataset(
        y = lois,
        x = age.loi,
        xout = dob.age,
        rep.negt = F
    )
# loi_interpolated <- cbind(loi.interp, dob.age)
head(loi.interp)
loi.interp <- as.data.frame(loi.interp)
loi.stand <- decostand(loi.interp, method = "standardize")
loi.stand <- loi.stand$loi1
loi.interp <- loi.interp$loi1

#-------------------------------------Data exploration---------------------------------------------#
# assembling environmental variables based on Lake Dobson age scale
# clim.complete with antartic and NZ records (NZW and NZA) not joined in PCA

clim.complete <- cbind(dome.interp, vos.interp, 
                       char.interp, et.interp, chiro.interp, loi.interp, co2.interp, 
                       O18.nza.interp, C13.nza.interp, O18.nzw.interp, C13.nzw.interp)
head(clim.complete)
clim.complete.stand <- cbind(dome.stand, vos.stand, char.stand,
                             et.stand, chiro.stand, loi.stand, co2.stand,
                             O18.nza.stand, C13.nza.stand, O18.nzw.stand, C13.nzw.stand)
head(clim.complete.stand)

# Patterns of correlations among env.var
(cor(clim.complete, use = "pairwise.complete.obs"))

(cor(clim.complete.stand, use = "pairwise.complete.obs"))

library(corrplot)
cor.env <- cor(clim.complete)
dev.new()
corrplot(cor.env, method = "circle")

antar.interp <- cbind(vos.interp, dome.interp)
head(antar.interp)
antar.stand <- cbind(vos.stand, dome.stand)
head(antar.stand)

# PCA based on a covariance matrix because variables are in the same scale 
# not sure if I should center
# Arguments scale = TRUE and centre = TRUE calls for a standardization of the variables
# antar.temp include the Vostok and DomeC variacion of temperature (delta TS)
head(antar.temp)
cor(antar.temp)
antar.pca <- prcomp(antar.interp, scale. = TRUE, center = TRUE) # temperature data of Antarctica 
summary(antar.pca) 
# values based on antar.interp PC1 = 90.54%; PC2 = 9.457%
# values based on antar.standPC1 = 89.97%; PC2 = 10.003%
antar.rda <- rda(antar.stand, scale = TRUE, center = TRUE)
antar.rda
summary(antar.rda)
names(antar.rda)

# calculate axis-importance and draw the barplots:
ev.antar <- antar.rda$CA$eig
source("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/evplot.R")
# source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/evplot.R")
dev.new()
# pdf("ev.pdf")
evplot(ev.antar)
# dev.off()
# only PC1 is important 

sites.antar <- scores(antar.rda, display = "sites", choices = 1)
head(sites.antar)
temp.antar<- as.data.frame(sites.antar)
head(temp.antar)   #only PC1
PC1.antar <-temp.antar$PC1
# PCA based on a covariance matrix because variables are in the same scale 
# not sure if I should center
# Arguments scale = TRUE and centre = TRUE calls for a standardization of the variables

source("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/cleanplot.pca.R")
# source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/cleanplot.pca.R")
dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(antar.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(antar.rda, scaling = 2, mar.percent = 0.04)

# PCA based on a covariance matrix because variables are in the same scale 
# not sure if I should center
# Arguments scale = TRUE and centre = TRUE calls for a standardization of the variables
# antar.temp include the Vostok and DomeC variacion de deuterium data

#---------------------Antartic & NZW.O18 together like in the Thesis-----------------------------#
# non-including 13C
southern.temp <- cbind(antar.interp, NZW.interp)
head(southern.temp)
southern.temp <- subset(southern.temp, select = -c(4))
head(southern.temp)
south.temp.pca <- prcomp(southern.temp, scale. = TRUE, center = TRUE) 
south.temp.pca
summary(south.temp.pca) 
# PC1 = 74.9%; PC2 = 18.4% based on temperature data
# Using the correlation matrix is equivalent to standardizing each of the variables 
# (to mean 0 and standard deviation 1).

south.temp.rda <- rda(southern.temp, scale = TRUE, center = TRUE)
south.temp.rda
summary(south.temp.rda)
# percentage of variance explained by each  component
(temp.eig <- south.temp.rda$CA$eig/south.temp.rda$tot.chi)

dev.new()
screeplot(south.temp.rda, bstick = TRUE)

ev.south.temp <- south.temp.pca$sdev^2
# pdf("ev.pdf")
dev.new()
evplot(ev.south.temp)

# only PC1 is significant
sites.south.temp.rda <- scores(south.temp.rda, choices = c(1), display = "sites")
sh.temp <-sites.south.temp.rda
head(sh.temp)
sh.temp <- as.data.frame(sh.temp)
PC1.sh.temp <- sh.temp$PC1 # PC1 of antarctic records + O18.nzw

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp.rda, scaling = 2, mar.percent = 0.04)

#----------------------------sh.temp1: Antartic + NZW complete------------------------------#

southern.temp1 <- cbind(antar.interp, NZW.interp)
head(southern.temp1)
head(southern.temp1)
south.temp.pca1 <- prcomp(southern.temp1, scale. = TRUE, center = TRUE) 
south.temp.pca1
summary(south.temp.pca1) 
# PC1 = 59.66%; PC2 = 23.86%; PC3 = 22.60%

south.temp.rda1 <- rda(southern.temp1, scale = TRUE, center = TRUE)
south.temp.rda1
summary(south.temp.rda1)
# percentage of variance explained by each  component
(temp.eig <- south.temp.rda1$CA$eig/south.temp.rda1$tot.chi)

dev.new()
screeplot(south.temp.rda1, bstick = TRUE)

ev.south.temp1 <- south.temp.pca1$sdev^2
# pdf("ev.pdf")
dev.new()
evplot(ev.south.temp1)
# still only PC1 is significant

sites.south.temp.rda1 <- scores(south.temp.rda1, choices = c(1), display = "sites")
sh.temp1 <-sites.south.temp.rda1
head(sh.temp1)
sh.temp1 <- as.data.frame(sh.temp1)
PC1.sh.temp1 <- sh.temp1$PC1 # PC1 of antarctic records + nzw complete

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp.rda1, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp.rda1, scaling = 2, mar.percent = 0.04)

#----------------------------Antartic records + O18 from Port Arthur-------------------------------#
# I dont think this is correct, becasue O18 from Port Arthur
# does not represent temerature, but instead meteoric water 
# 13C represent changes in forest productivitym matching paleovegetation records
southern.temp2 <- cbind(antar.interp, NZA.interp)
head(southern.temp2)
southern.temp2 <- subset(southern.temp2, select = -c(4))
head(southern.temp2)
south.temp2.pca <- prcomp(southern.temp2, scale. = TRUE, center = TRUE) 
south.temp2.pca
summary(south.temp2.pca)
south.temp2.rda <- rda(southern.temp2, scale = TRUE, center = TRUE)
south.temp2.rda
summary(south.temp2.rda)

# percentage of variance explained by each  component
(temp.eig <- south.temp2.rda$CA$eig/south.temp2.rda$tot.chi)
summary(south.temp2.pca) 
# PC1 = 78.13%; PC2 = 15.37% 

dev.new()
screeplot(south.temp2.rda, bstick = TRUE)

# select the data frame with eigenvalues of particular axes:
ev.south.temp2 <- south.temp2.pca$sdev^2
# pdf("ev.pdf")
dev.new()
evplot(ev.south.temp2)
# only PC1 is significant

sites.south.temp2.rda <- scores(south.temp2.rda, display = "sites", choices = 1)
sh.temp2 <- sites.south.temp2.rda
head(sh.temp2)
# temp.comp  contains only PC1
sh.temp2 <- as.data.frame(sh.temp2)
PC1.sh.temp2 <- sh.temp2$PC1 # PC1 of antarctic + O18 nza

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp2.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp2.rda, scaling = 2, mar.percent = 0.04)

#  it is weird how NZW and NZA points towards the oldest sites from Lake Dobson
#  could it be because the 2 records represent precipitation values instead of temperature?

#--------------------------------Antartic records + O18 from nza and nzw---------------------------# 

NZ.isotopes <- cbind(NZW.interp, NZA.interp)
head(NZ.isotopes)
O18.nz <- subset(NZ.isotopes, select = -c(2,4))
head(O18.nz)

O18.pca <- prcomp(O18.nz, scale. = FALSE, center = TRUE) 
O18.pca 
summary(O18.pca) # PC1 = 73.86%; PC2 = 26.24%
# percentage of variance explained by each  component

O18.rda <- rda(O18.nz, scale = FALSE, center = TRUE)
O18.rda
summary(O18.rda)
(O18.eig <- O18.rda$CA$eig/O18.rda$tot.chi)
# PC1       PC2 
# 0.7386376 0.2613624 

dev.new()
screeplot(O18.rda, bstick = TRUE)
# source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/Scree.Plot.R")
source("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/Scree.Plot.R")

R <- cor(O18.nz[,1:2])
Scree.Plot(R, main ="Scree Plot (Environmental variables)") # only PC1 is important

# pdf("ev.pdf")
dev.new()
ev.O18 <- O18.pca$sdev^2
evplot(ev.O18)
# only PC1 is signidicant

sites.O18.rda <- scores(O18.rda, display = "sites", choices = 1)
O18.sh <- sites.O18.rda
head(O18.sh)
# temp.comp  contains only PC2
O18.sh<- as.data.frame(O18.sh)
PC1.sh.O18 <- O18.sh$PC1

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(O18.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(O18.rda, scaling = 2, mar.percent = 0.04)

head(NZ.isotopes)
C13.nz <-subset(NZ.isotopes, select = -c(1,3))
head(C13.nz)
C13.pca <- prcomp(C13.nz, scale. = FALSE, center = TRUE) 
C13.pca 
summary(C13.pca) 
# PC1 = 53.99%; PC2 = 46.01%

C13.rda <- rda(C13.nz, scale = FALSE, center = TRUE)
C13.rda
summary(C13.rda)
(C13.eig <- C13.rda$CA$eig/C13.rda$tot.chi)

dev.new()
screeplot(C13.rda, bstick = TRUE)

# pdf("ev.pdf")
dev.new()
ev.C13 <- C13.pca$sdev^2
evplot(ev.C13)
# Only PC1 
sites.C13.rda <- scores(C13.rda, display = "sites", choices = 1)
C13.sh <- sites.C13.rda
head(C13.sh)
C13.sh<- as.data.frame(C13.sh)
PC1.sh.C13 <- C13.sh$PC1

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(C13.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(C13.rda, scaling = 2, mar.percent = 0.04)

# ------------------------Antartic data and O18 isotopes from nzw and nza--------------------------#

southern.climate1 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate1)
southern.climate1 <- subset(southern.climate1, select = -c(4,6))# removing C13 from both records
head(southern.climate1)
south.climate1.pca <- prcomp(southern.climate1, scale. = TRUE, center = TRUE) 
south.climate1.pca
summary(south.climate1.pca)
# PC1 = 69.33%; PC2 = 14.36% with O18 isotopes from both records

south.climate1.rda <- rda(southern.climate1, scale = TRUE, center = TRUE)
south.climate1.rda
summary(south.climate1.rda)
(temp.eig <- south.climate1.rda$CA$eig/south.climate1.rda$tot.chi)

dev.new()
screeplot(south.climate1.rda, bstick = TRUE)

R <- cor(southern.climate1[,1:4])
Scree.Plot(R, main ="Scree Plot (Environmental variables)") # only PC1 is important,

# select the data frame with eigenvalues of particular axes:
ev <- south.climate1.rda$CA$eig
# pdf("ev.pdf")
dev.new()
evplot(ev)
# only PC1 is significant when running the analysis with O18 data from nzw and nza

sites.south.climate1.rda <- scores(south.climate1.rda, choices = c(1), display = "sites")
sh.climate1 <- sites.south.climate1.rda
head(sh.climate1)
sh.climate1 <- as.data.frame(sh.climate1)
PC1.sh.clim1 <- sh.climate1$PC1 # antactic records + O18 nzw & O18 nza

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate1.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate1.rda, scaling = 2, mar.percent = 0.04)

#-------------------------------Antartic + nzw & nza complete-------------------------------------#

southern.climate2 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate2)
south.climate2.pca <- prcomp(southern.climate2, scale. = TRUE, center = TRUE) 
south.climate2.pca
summary(south.climate2.pca)
# PC1 = 48.55%; PC2 = 19.78%; PC3 = 14.66% with O18 and C13 isotopes from both records

south.climate2.rda <- rda(southern.climate2, scale = TRUE, center = TRUE)
south.climate2.rda
summary(south.climate2.rda)
(temp.eig <- south.climate2.rda$CA$eig/south.climate2.rda$tot.chi)

dev.new()
screeplot(south.climate2.rda, bstick = TRUE)

R <- cor(southern.climate2[,1:6])
Scree.Plot(R,main ="Scree Plot (Environmental variables)")
# PC1 & PC2 are significant

# select the data frame with eigenvalues of particular axes:
ev <- south.climate2.rda$CA$eig
# pdf("ev.pdf")
dev.new()
evplot(ev)
# first two PC are significant when running the analysis with O18 and C13 
# isotopes data from nzw and nza

sites.south.climate2.rda <- scores(south.climate2.rda, choices = c(1:2), display = "sites")
sh.climate2 <- sites.south.climate2.rda
head(sh.climate2)
sh.climate2 <- as.data.frame(sh.climate2)
PC1.sh.clim2 <- sh.climate2$PC1
PC2.sh.clim2 <- sh.climate2$PC2

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate2.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate2.rda, scaling = 2, mar.percent = 0.04)
# scaling 1 = distancia entre los sitios o muestras (scaling sites)
# scaling 1, incluye circulo de contribucion de equilibrio # However we can clearly see that 
# scaling 2 = correlacion entre variables representada por el angulo de los vectores (scaling environmental variables)
# I have the impression that I have multicollinearity between co2 and loi
# also perhaps collinearity (correlation) between antartic records  of temperature

# ------------------------antarctic + O18 from nzw & nza + C13 nzw---------------------------------#

southern.climate3 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate3)
southern.climate3 <- subset(southern.climate3, select = -c(6))# removing C13.nza from both records
head(southern.climate3)
south.climate3.pca <- prcomp(southern.climate3, scale. = TRUE, center = TRUE) 
south.climate3.pca
summary(south.climate3.pca)
# PC1 = 56.97%; PC2 = 21.93% with O18 from nzw and nza + C13 isotopes from nzw

south.climate3.rda <- rda(southern.climate3, scale = TRUE, center = TRUE)
south.climate3.rda
summary(south.climate3.rda)
(temp.eig <- south.climate3.rda$CA$eig/south.climate3.rda$tot.chi)

dev.new()
screeplot(south.climate3.rda, bstick = TRUE)

R <- cor(southern.climate3[,1:5])
Scree.Plot(R,main ="Scree Plot (Environmental variables)")
# Two first PCA are significant, althoug PC2 is just above the red line

# select the data frame with eigenvalues of particular axes:
ev <- south.climate3.rda$CA$eig

# pdf("ev.pdf")
dev.new()
evplot(ev)
# first two PC are significant when running the analysis with O18 and C13 
# isotopes data from nzw and nza

sites.south.climate3.rda <- scores(south.climate3.rda, choices = c(1:2), display = "sites")
sh.climate3 <- sites.south.climate3.rda
head(sh.climate3)
sh.climate3 <- as.data.frame(sh.climate3)
PC1.sh.clim3 <- sh.climate3$PC1
PC2.sh.clim3 <- sh.climate3$PC2

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate3.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate3.rda, scaling = 2, mar.percent = 0.03)

#-------------------------------------O18 from NZW & NZA + C13 nzw---------------------------------#

nz.climate <- cbind(NZW.interp, NZA.interp)
head(nz.climate)
nz.clim <- subset(nz.climate, select = -c(4))# removing C13.nza 
head(nz.clim)
nz.clim.pca <- prcomp(nz.clim, scale. = TRUE, center = TRUE) 
nz.clim.pca
summary(nz.clim.pca)
# PC1 = 52.89%; PC2 = 34.05% with O18 from nzw and nza + C13 isotopes from nzw

nz.clim.rda <- rda(nz.clim, scale = TRUE, center = TRUE)
nz.clim.rda
summary(nz.clim.rda)
(temp.eig <- nz.clim.rda$CA$eig/nz.clim.rda$tot.chi)

dev.new()
screeplot(nz.clim.rda, bstick = TRUE)

R <- cor(nz.clim[,1:3])
Scree.Plot(R,main ="Scree Plot (Environmental variables)")
# Two first PCA are significant, althoug PC2 is just above the red line

# select the data frame with eigenvalues of particular axes:
ev <- nz.clim.rda$CA$eig

# pdf("ev.pdf")
dev.new()
evplot(ev)
# first two PC are significant when running the analysis with O18 and C13 
# isotopes data from nzw and nza

sites.nz.clim.rda <- scores(nz.clim.rda, choices = c(1:2), display = "sites")

head(sites.nz.clim.rda)
NZ.clim <- as.data.frame(sites.nz.clim.rda)
PC1.nz.clim <- NZ.clim$PC1
PC2.nz.clim <- NZ.clim$PC2

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(nz.clim.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(nz.clim.rda, scaling = 2, mar.percent = 0.03)

# --------------------------------------local climate LD-------------------------------------------#
local.clim <- cbind(chiro.interp, et.interp)
head(local.clim)
local.clim.pca <- prcomp(local.clim, scale. = TRUE, center = TRUE) 
local.clim.pca
summary(local.clim.pca)
# PC1 = 79.6%; PC2 = 20.41%

local.clim.rda <- rda(local.clim, scale = TRUE, center = TRUE)
local.clim.rda
summary(local.clim.rda)
(temp.eig <- local.clim.rda$CA$eig/local.clim.rda$tot.chi)

dev.new()
screeplot(local.clim.rda, bstick = TRUE)

R <- cor(local.clim[,1:2])
Scree.Plot(R, main ="Scree Plot (Environmental variables)")
# PC1 is significant

# select the data frame with eigenvalues of particular axes:
ev <- local.clim.rda$CA$eig
# pdf("ev.pdf")
dev.new()
evplot(ev)
# PC1 is significant 
sites.local.clim.rda <- scores(local.clim.rda, choices = c(1), display = "sites")
local.climate <- sites.local.clim.rda
head(local.climate)
local.climate <- as.data.frame(local.climate)
PC1.local.clim <- local.climate$PC1

head(clim.complete)
clim.complete1 <- cbind(PC1.ant.temp, PC1.sh.temp, PC1.sh.temp1, PC1.sh.temp2, 
                        PC1.sh.clim1, PC1.sh.clim2, PC2.sh.clim2, PC1.sh.clim3, PC2.sh.clim3, 
                        co2.interp, loi.interp, char.interp, et.interp, chiro.interp, PC1.sh.O18, 
                        PC1.sh.C13, PC1.nz.clim, PC2.nz.clim, PC1.local.clim)
head(clim.complete1)
is.data.frame(clim.complete1)
clim.complete1 <- as.data.frame(clim.complete1)

#----------------------------------------individual cca--------------------------------------------#
# Individual CCA with each variable log1p(macros); 80 columns 52 variables
# macros5 <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/combined_macros_shorted(5).csv")
macros5 <- read.csv("/Users/giselleastorga/Google Drive/Lake Dobson/plant-macros/DataLD/combined_macros_shorted(5).csv")
str(macros5)
macros.log <- log1p(macros5)# computes log(1+x)
macros.sqrt <- sqrt(macros5)

# dataset containing 47 species and 80 samples, alpine conifers are clumped together
macros.short <- read.csv("combined_macros_shorted.csv")
str(macros.short)
macros.short.log <- log2p(macros.short)

head(clim.complete1)

antar.mod <- cca (macros.log ~ PC1.ant.temp, data = clim.complete1)
antar.mod
summary(antar.mod)
# PC1.antar.temp explains about 12% (11.740%) of the total variability
control <- how(within = Within(type = "series", mirror = TRUE))
(check2 <- check(macros.log, control))
summary(check2)
(anova(antar.mod,  permutations = control))

#           Df ChiSquare      F Pr(>F)    
# Model     1   0.31643 10.376 0.00625 **
# Residual 78   2.37875                        

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.5 ‘ ’ 1

head(clim.complete1)
sh.temp.mod <- cca (macros.log ~ PC1.sh.temp, data = clim.complete1)
sh.temp.mod
summary(sh.temp.mod)
# sh.temp (antartic + O18.nzw) explains about 12% (12.35%) of the total variability
(anova(sh.temp.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.33298 10.995 0.00625 **
# Residual 78   2.36220                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)
sh.temp1.mod <- cca (macros.log ~ PC1.sh.temp1, data = clim.complete1)
sh.temp1.mod
summary(sh.temp1.mod)
# sh.temp1 (antartic + nzw complete) explains about 12% (11.54%) of the total variability
(anova(sh.temp1.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.31115 10.18 0.00625 **
# Residual 78   2.38403                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

sh.temp2.mod <- cca(macros.log ~ PC1.sh.temp2, data = clim.complete1)  
sh.temp2.mod 
summary(sh.temp2.mod )
# sh.temp2 PC1 explains about ~13% (12.52%) of the constrained variability
(anova(sh.temp2.mod,  permutations = control))
#           Df ChiSquare      F Pr(>F)    
# Model     2   0.39971 6.704 0.00625 **
# Residual 77   2.29547                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
sh.clim1.mod <- cca(macros.log ~ PC1.sh.clim1, data = clim.complete1)  
sh.clim1.mod
summary(sh.clim1.mod)
# sh.clim1 CCA1 explains about ~13% (12.91%) of the constrained variability
(anova(sh.clim1.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.34783 11.558 0.00625 **
# Residual 78   2.34736                       

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
sh.clim2.mod <- cca(macros.log ~ PC1.sh.clim2 + PC2.sh.clim2, data = clim.complete1)  
sh.clim2.mod
summary(sh.clim2.mod)
# sh.clim2 using PCA1 & PCA2 explains about ~15% (14.83%) of the constrained variability
# PC1.sh.clim2 explains 12.51% of the constrained variability

(anova(sh.clim2.mod,  permutations = control))
# Df ChiSquare     F  Pr(>F)   
# Model     2   0.39971 6.704 0.00625 **
# Residual 77   2.29547                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
sh.clim3.mod <- cca(macros.log ~ PC1.sh.clim3 + PC2.sh.clim3, data = clim.complete1)  
sh.clim3.mod
summary(sh.clim3.mod)

# sh.clim3 using PCA1 & PCA2 explains about ~15% (14.57%) of the constrained variability
# PC1.sh.clim3 explains 12.63% of the constrained variability

(anova(sh.clim3.mod,  permutations = control))
# Df ChiSquare      F  Pr(>F)   
# Model     1   0.34042 11.276 0.00625 **
# Residual 78   2.35476                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


head(clim.complete1)  
local.clim.mod <- cca(macros.log ~ PC1.local.clim, data = clim.complete1)  
local.clim.mod
#                Inertia Proportion Rank
# Total         2.69518    1.00000     
# Constrained   0.26491    0.09829    1
# Unconstrained 2.43027    0.90171   51
# Inertia is scaled Chi-square 

(anova(local.clim.mod,  permutations = control))
# Df ChiSquare      F  Pr(>F)   
# Model     1   0.26491 8.5025 0.00625 **
# Residual 78   2.43027            

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
nz.clim.mod <- cca(macros.log ~ PC1.nz.clim + PC2.nz.clim, data = clim.complete1)  
nz.clim.mod
#                 Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.3419     0.1268    2
# Unconstrained  2.3533     0.8732   51
# Inertia is scaled Chi-square 
summary(nz.clim.mod)

# Accumulated constrained eigenvalues
# Importance of components:
#                 CCA1      CCA2
# Total          2.6952     1.0000     
# Constrained    0.3419     0.1268    2
# Unconstrained  2.3533     0.8732   51

# nz.clim.mod using PCA1 & PCA2 explains about ~13% (12.68%) of the constrained variability
summary(nz.clim.mod)

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1    CCA2
# Eigenvalue            0.2891 0.05274
# Proportion Explained  0.8457 0.15429
# Cumulative Proportion 0.8457 1.00000

(anova(nz.clim.mod,  permutations = control))
# Df ChiSquare      F  Pr(>F)   
# Model     2   0.34186 5.5927 0.00625 **
# Residual 77   2.35332               

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)
co2.mod <- cca(macros.log ~ co2.stand, data = clim.complete1)  
co2.mod

# co2 explains about ~13% (12.51%) of the total variability
(anova(co2.mod,  permutation = control))
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.33711 11.151  0.001 ***
# Residual 78   2.35807                

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

chiro.mod <- cca(macros.log ~ chiro.stand, data = clim.complete1)  
chiro.mod
# Chironomids DCA2analysis explains about ~13% (12.71) of the total variability
(anova(chiro.mod,  permutations = control))
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.34247 11.354 0.00625 **
# Residual 78   2.35271  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

loi.mod <- cca(macros.log ~ loi.stand, data = clim.complete1)  
loi.mod
summary(loi.mod)
# loi explains about 9% (9.266) of the total variability
(anova(loi.mod,  permutations = control))
#          Df ChiSquare      F Pr(>F)    
# Model     1   0.24972 7.9652  0.025 *
# Residual 78   2.44546                      

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

char.mod <- cca(macros.log ~ char.stand, data = clim.complete1)  
char.mod
summary(char.mod)
# char explains just about 3 % (2.543%) of the total variability
(anova(char.mod,  permutations = control))
# model is marginally significant at the 0.2 level (p-value = 0.08125) 

et.mod <- cca(macros.log ~ et.stand, data = clim.complete1)  
et.mod
summary(et.mod)
# Eagle Tarn temperature explains about 4% (4.447%) of total variability
(anova(et.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.11985 3.6298  0.025 *
# Residual 78   2.57533                   

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

O18.mod <- cca(macros.log ~ PC1.sh.O18, data = clim.complete1)  
O18.mod

# O18 isotopes from nza and nzw explain together 10.66% (CCA1)
(anova(O18.mod,  permutations = control))
#          Df ChiSquare      F  Pr(>F)   
# Model     1    0.2874 9.3102 0.00625 **
# Residual 78    2.4078                  

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

C13.mod <- cca(macros.log ~ PC1.sh.C13, data = clim.complete1)  
C13.mod
# C13 isotopes from nza and nzw explain together 2.207% (CCA1)
(anova(C13.mod,  permutations = control))
# non-significant    

local.mod <- cca(macros.log ~ PC1.local.clim, data = clim.complete1)  
local.mod
summary(local.mod)
# local explains just about 10% (9.829%) of the total variability
(anova(local.mod,  permutations = control))
#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.26491 8.5025 0.00625 **
# Residual 78   2.43027                

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

set.seed(60)
head(clim.complete1)
clim.complete2 <-cbind(clim.complete1, dob.age)
head(clim.complete2)
cond.model <- cca(macros.log ~ dob.age, data = clim.complete2)
cond.model
#                Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.3311     0.1228    1
# Unconstrained  2.3641     0.8772   51
# Inertia is scaled Chi-square 
summary(cond.model)
(anova(cond.model,  permutations = control))
#          Df ChiSquare      F  Pr(>F)   
# Model     1    0.3311 10.924 0.00625 **
# Residual 78    2.3641                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# age explains about 12% (12.28%) of the total variability 
# permutation test showed the sediment age was significantly 
# correlated with the plant macrofossil data at 0.01 level (p-value = 0.00625)
# therefore sediment age is used as covariable in subsequent ordination

# Manual model building
head(clim.complete2)
# clim.complete3 includes PC1.ant.temp, co21, loi1, char1, temp.et1, dca1.dob, PC1.sh.O18
# PC1.sh.O18 does include O18 isotope data from nzw and nza
clim.complete3 <- subset(clim.complete2, select = -c(2:9,17:18))
head(clim.complete3)
clim.complete3 <- as.data.frame(clim.complete3)

head(clim.complete.stand)
head(clim.complete1)
#--------------------------------correlation of environmental variable-----------------------------#

# Function bioenv() finds the best subset of environmental variables, so that the Euclidean 
# distances of scaled environmental variables have the maximum (rank) correlation with 
# community dissimilarities.

sol <- bioenv(wisconsin(macros.log) ~ PC1.ant.temp + co21 +  char1 + loi1 + temp.et1 + dca1.dob +
                PC1.sh.O18, clim.complete3, method = "pearson", index = "bray")
sol
summary(sol)

## Is vegetation related to environment?

veg.dist <- vegdist(macros.log) # Bray-Curtis
env.dist <- vegdist(scale(clim.complete3), "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method = "pearson")


library("FactoMineR")
clim3.pca <- PCA(clim.complete3, scale.unit = TRUE, graph = FALSE)
library("factoextra")
eig.val <- get_eigenvalue(clim3.pca)
eig.val
dev.new()
fviz_eig(clim3.pca, addlabels = TRUE, ylim = c(0, 80))
var.clim3 <- get_pca_var(clim3.pca)
var.clim3
# var$coord: coordinates of variables to create a scatter plot
# var$cos2: represents the quality of representation for variables on the factor map. 
# It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# var$contrib: contains the contributions (in percentage) of the variables to the 
# principal components. The contribution of a variable (var) to a given principal
# component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
head(var.clim3$coord, 7)
fviz_pca_var(clim3.pca, col.var = "black")
# Positively correlated variables are grouped together.
# Negatively correlated variables are positioned on opposite sides
# of the plot origin (opposed quadrants).
# The distance between variables and the origin measures the quality of the variables on the
# factor map. Variables that are away from the origin are well represented on the factor map.
# In the case of clim.complete3 char, co2, and loi are posively correlated, while PC1.sh.O18,
# dca1.dob, PC1.ant.temp, temp.et are positively correlated.
# Char and co2 are better represented than loi (shorter arrow)
# dca1.dob, PC1.ant.temp, PC1.sh.O18 are better represented compared to temp.et1

# Quality of representation
# The quality of representation of the variables on factor map is called cos2 
# (square cosine, squared coordinates). You can access to the cos2 as follow:
head(var.clim3$cos2, 7)

library("corrplot")
dev.new()
corrplot(var.clim3$cos2, is.corr=FALSE)
# A high cos2 indicates a good representation of the variable on the principal component. 
# In this case the variable is positioned close to the circumference of the correlation circle.
# A low cos2 indicates that the variable is not perfectly represented by the PCs. In this case 
# the variable is close to the center of the circle.
# The cos2 values are used to estimate the quality of the representation
# The closer a variable is to the circle of correlations, the better its representation on the
# factor map (and the more important it is to interpret these components)
# Variables that are closed to the center of the plot are less important for the first components.

# Color by cos2 values: quality on the factor map
dev.new()
fviz_pca_var(clim3.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
# Contributions of variables to PCs
head(var.clim3$contrib, 7)
dev.new()
corrplot(var.clim3$contrib, is.corr =F) 

# Contributions of variables to PC1
dev.new()
fviz_contrib(clim3.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(clim3.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(clim3.pca, choice = "var", axes = 1:2, top = 10)

# clim.complete4 includes PC1.sh.temp, co21, loi1, char1, temp.et1, dca1.dob
# PC1.sh.temp includes antartic records + O18
head(clim.complete2)
clim.complete4 <- subset(clim.complete2, select = -c(1:2, 4:10, 16:18))
head(clim.complete4)
clim.complete4 <- as.data.frame(clim.complete4)
clim4.pca <- PCA(clim.complete4, scale.unit = TRUE, graph = FALSE)

eig.val <- get_eigenvalue(clim4.pca)
eig.val
dev.new()
fviz_eig(clim4.pca, addlabels = TRUE, ylim = c(0, 80))
var.clim4 <- get_pca_var(clim4.pca)
var.clim4
head(var.clim4$coord, 6)
fviz_pca_var(clim4.pca, col.var = "black")
# Positively correlated variables are grouped together.
# Negatively correlated variables are positioned on opposite sides
# of the plot origin (opposed quadrants).
# The distance between variables and the origin measures the quality of the variables on the
# factor map. Variables that are away from the origin are well represented on the factor map.
# In the case of clim.complete4 char, co2, and loi are posively correlated, while PC1.sh.O18,
# dca1.dob, PC1.sh.temp, temp.et are positively correlated.
# Char and co2 are better represented than loi (shorter arrow)
# dca1.dob, PC1.sh.temp, are better represented compared to temp.et1

# Quality of representation
# The quality of representation of the variables on factor map is called cos2 
# (square cosine, squared coordinates). You can access to the cos2 as follow:
head(var.clim4$cos2, 6)

dev.new()
corrplot(var.clim4$cos2, is.corr=FALSE)
# A high cos2 indicates a good representation of the variable on the principal component. 
# In this case the variable is positioned close to the circumference of the correlation circle.
# A low cos2 indicates that the variable is not perfectly represented by the PCs. In this case 
# the variable is close to the center of the circle.
# The cos2 values are used to estimate the quality of the representation
# The closer a variable is to the circle of correlations, the better its representation on the
# factor map (and the more important it is to interpret these components)
# Variables that are closed to the center of the plot are less important for the first components.

# Color by cos2 values: quality on the factor map
dev.new()
fviz_pca_var(clim4.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
# Contributions of variables to PCs
head(var.clim4$contrib, 6)
dev.new()
corrplot(var.clim4$contrib, is.corr =F) 

# Contributions of variables to PC1
dev.new()
fviz_contrib(clim4.pca, choice = "var", axes = 1, top = 6)
# Contributions of variables to PC2
fviz_contrib(clim4.pca, choice = "var", axes = 2, top = 6)
fviz_contrib(clim4.pca, choice = "var", axes = 1:2, top = 6)

# clim.complete5 includes PC1.sh.temp1, co21, loi1, char1, temp.et1, dca1.dob
# PC1.sh.temp1 is antartic and nzw complete (O18 & C13 isotopes)
head(clim.complete2)
clim.complete5 <- subset(clim.complete2, select = -c(1:3, 5:10, 16:18))
head(clim.complete5)
clim.complete5 <- as.data.frame(clim.complete5)
clim5.pca <- PCA(clim.complete5, scale.unit = TRUE, graph = FALSE)

eig.val <- get_eigenvalue(clim5.pca)
eig.val
dev.new()
fviz_eig(clim5.pca, addlabels = TRUE, ylim = c(0, 80))
var.clim5 <- get_pca_var(clim5.pca)
var.clim5
head(var.clim5$coord, 6)
fviz_pca_var(clim5.pca, col.var = "black")
# Positively correlated variables are grouped together.
# Negatively correlated variables are positioned on opposite sides
# of the plot origin (opposed quadrants).
# The distance between variables and the origin measures the quality of the variables on the
# factor map. Variables that are away from the origin are well represented on the factor map.
# In the case of clim.complete5 char, co2, and loi are posively correlated, while PC1.sh.O18,
# dca1.dob, PC1.sh.temp, temp.et are positively correlated.
# Char and co2 are better represented than loi (shorter arrow)
# dca1.dob, PC1.sh.temp, are better represented compared to temp.et1

# Quality of representation
# The quality of representation of the variables on factor map is called cos2 
# (square cosine, squared coordinates). You can access to the cos2 as follow:
head(var.clim5$cos2, 6)

dev.new()
corrplot(var.clim5$cos2, is.corr=FALSE)
# A high cos2 indicates a good representation of the variable on the principal component. 
# In this case the variable is positioned close to the circumference of the correlation circle.
# A low cos2 indicates that the variable is not perfectly represented by the PCs. In this case 
# the variable is close to the center of the circle.
# The cos2 values are used to estimate the quality of the representation
# The closer a variable is to the circle of correlations, the better its representation on the
# factor map (and the more important it is to interpret these components)
# Variables that are closed to the center of the plot are less important for the first components.

# Color by cos2 values: quality on the factor map
dev.new()
fviz_pca_var(clim5.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC5E07"), 
             repel = TRUE # Avoid text overlapping
)
# Temp1.et and dca1.dob perfectly overlap, but also temp1.et is the shorter env. vector 
# PC1.sh.temp1 and dca1.dob are highly correlated, although the local variable is closer to the  
# circle meaning is better represented in PC1 (0.93 vs. a.50)

# Contributions of variables to PCs
head(var.clim5$contrib, 6)
dev.new()
corrplot(var.clim5$contrib, is.corr =F) 

# Contributions of variables to PC1
dev.new()
fviz_contrib(clim5.pca, choice = "var", axes = 1, top = 6)
# Contributions of variables to PC2
fviz_contrib(clim5.pca, choice = "var", axes = 2, top = 6)
fviz_contrib(clim5.pca, choice = "var", axes = 1:2, top = 6)

# clim.complete6 includes PC1.sh.temp1, co21, loi1, char1, temp.et1, dca1.dob
# PC1.sh.temp 1 is antartic and nzw complete (O18 & C13 isotopes)
head(clim.complete2)
clim.complete6 <- subset(clim.complete2, select = -c(1:3, 5:9, 15:16))
head(clim.complete6)
clim.complete6 <- as.data.frame(clim.complete6)
clim6.pca <- PCA(clim.complete6, scale.unit = TRUE, graph = FALSE)

eig.val <- get_eigenvalue(clim6.pca)
eig.val
dev.new()
fviz_eig(clim6.pca, addlabels = TRUE, ylim = c(0, 80))
var.clim6 <- get_pca_var(clim6.pca)
var.clim6
head(var.clim6$coord, 6)
fviz_pca_var(clim6.pca, col.var = "black")

head(var.clim6$cos2, 6)
dev.new()
corrplot(var.clim6$cos2, is.corr = FALSE)

dev.new()
fviz_pca_var(clim6.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC5E07"), 
             repel = TRUE # Avoid text overlapping
)

# Contributions of variables to PCs
head(var.clim6$contrib, 6)
dev.new()
corrplot(var.clim6$contrib, is.corr =F) 

# Contributions of variables to PC1
dev.new()
fviz_contrib(clim6.pca, choice = "var", axes = 1, top = 6)
# Contributions of variables to PC2
fviz_contrib(clim6.pca, choice = "var", axes = 2, top = 6)
fviz_contrib(clim6.pca, choice = "var", axes = 1:2, top = 6)


# clim.complete7 includes PC1.ant.temp, PC1.nz.clim, PC2.nz.clim, co21, loi1, char1, 
# PC1.local.clim 
head(clim.complete2)
clim.complete7 <- subset(clim.complete2, select = -c(2:10, 14:16,20))
head(clim.complete7)
clim.complete7 <- as.data.frame(clim.complete7)
clim7.pca <- PCA(clim.complete7, scale.unit = TRUE, graph = FALSE)

# --------------------------------cca with local variables standadised-----------------------------#
local.env <- cbind(chiro.stand, loi.stand, et.stand, char.stand)
head(local.env)
local.env <- as.data.frame(local.env)

# Histograms environmental variables

dev.new()
par(mfrow=c(2, 2))
hist(local.env$chiro.stand)
hist(local.env$loi.stand)
hist(local.env$et.stand)
hist(local.env$char.stand)
skewness(local.env$chiro.stand)  # skewness = 1.2 and kurtosis = 2.9
kurtosis(local.env$chiro.stand)
skewness(local.env$loi.stand)  # skewness = -0.89 and kurtosis = 2.5
kurtosis(local.env$loi.stand)
skewness(local.env$et.stand)  # skewness = 0.0087 and kurtosis = 2.1
kurtosis(local.env$et.stand)
skewness(local.env$char.stand)# skewness = 2.6 and kurtosis = 11.
kurtosis(local.env$char.stand)

# Perhaps I should normalize before standardising

# Histogram from a single numeric vector 
# ggplot2.histogram(data=numVector)
# Basic histogram plot from the vector "weight"
(range(chiro.stand))

library(ggplot2)
dev.new()
par(mfrow=c(2, 2))
chiro <- ggplot2_histogram(data = local.env, xName = 'chiro.stand', scale = "density", binwidth = .7)
loi <- ggplot2.histogram(data = local.env, xName = 'loi.stand', scale = "density", binwidth = .7)
grid.arrange(dome, vos, nrow = 2) 
# Change the width of bars
# Change y axis values to density
dev.new()
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos2', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                         meanLineType = "dashed", meanLineSize = 2)
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome2', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                          meanLineType = "dashed", meanLineSize = 2)
grid.arrange(dome, vos, nrow = 2) 

# Add density curve
vos2 <- ggplot2.histogram(data = local.env, xName = 'temp.vos2',
                          fill = "white", color = "black",
                          addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)

dome2 <- ggplot2.histogram(data = antar.climate2, xName = 'temp.dome2',
                           fill = "white", color = "black",
                           addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)
set.seed(50)
mbig <- cca(macros.log ~  ., local.env)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, local.env)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + chiro.stand)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char.stand)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + loi.stand)
add1(m0, scope = formula(mbig), test = "permutation")

## -- included variables still significant?

local.mod1 <- cca(macros.log ~ chiro.stand + char.stand + loi.stand, 
                  local.env)
local.mod1
#               Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.4374     0.1623    3
# Unconstrained  2.2578     0.8377   51
# Inertia is scaled Chi-square 
summary(local.mod1)

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1    CCA2    CCA3
# Eigenvalue            0.3471 0.05422 0.03602
# Proportion Explained  0.7937 0.12396 0.08235
# Cumulative Proportion 0.7937 0.91765 1.00000

vif.cca(local.mod1)
# chiro.stand  char.stand   loi.stand 
# 2.491425    1.083385    2.366675 

local.cca <- cca(macros.log ~., local.env)
summary(local.cca)
RsquareAdj(local.cca)

set.seed(70)
cca.local.forward <- ordistep(cca(macros.log ~ 1, data = local.env),
                             scope = formula(local.cca),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
local.mod.pars <- cca(macros.log ~ chiro.stand + char.stand + loi.stand,
                     data = local.env)
summary(local.mod.pars)
# Partitioning of scaled Chi-square:
#                Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.4374     0.1623
# Unconstrained  2.2578     0.8377

# Accumulated constrained eigenvalues
# Importance of components:
#                       CCA1    CCA2    CCA3
# Eigenvalue            0.3471 0.05422 0.03602
# Proportion Explained  0.7937 0.12396 0.08235
# Cumulative Proportion 0.7937 0.91765 1.00000

control <- how(within = Within(type = "series", mirror =TRUE))
(check2 <- check(macros.log, control))
anova(local.mod.pars, permutations = control)
#          Df ChiSquare      F  Pr(>F)   
# Model     3   0.43737 4.9075 0.00625 **
# Residual 76   2.25781               

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(local.mod.pars, permutations = control, by = "axis") # CCA1 & CCA2 are important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.34714 11.6849 0.00625 **
# CCA2      1   0.05422  1.8250 0.12500   
# CCA3      1   0.03602  1.2124 0.30625   
# Residual 76   2.25781   

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(local.mod.pars)
#$r.squared
# [1] 0.1622792
# 
# $adj.r.squared
# [1] 0.1291338

vif.cca(local.mod.pars)
# chiro.stand  char.stand   loi.stand 
# 2.491425    1.083385    2.366675 

## CCA triplots (using lc site scores)
# Scaling 1: species scores scaled to relative eigenvalues, 
# sites are weighted averages of the species
quartz(title = "CCA triplot - scaling 1 - lc scores")
plot(local.mod.pars, scaling = 1, display = c("sp","lc","cn"), 
     main = "Triplot CCA spe ~ clim6 - scaling 1")

# Default scaling 2: site scores scaled to relative eigenvalues, 
# species are weighted averages of the sites
quartz(title="CCA triplot - scaling 2 - lc scores")
plot(local.mod.pars, display = c("sp","lc","cn"), 
     main="Triplot CCA spe ~ clim6 - scaling 2")

# CCA scaling 1 biplot without species (using lc site scores)
quartz(title = "CCA biplot - scaling 1")
plot(local.mod.pars, scaling = 1, display = c( "sp", "cn"), 
     type = "text",
     main="Biplot CCA spe ~ clim6 - scaling 1")


#-------------------------------------Defining final model-----------------------------------------#


# add1(m0, scope=formula(mbig), test = c("permutation"), 
#      permutations = how(within = Within(type = "series", mirror = TRUE)))

# clim.complete3 includes PC1.ant.temp, co21, loi1, char1, temp.et1, dca1.dob, PC1.sh.O18
# PC1.sh.O18 does include O18 isotope data from nzw and nza
head(clim.complete3)
set.seed(50)
mbig <- cca(macros.log ~  ., clim.complete3)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete3)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.sh.O18)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + co21)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")

## -- included variables still significant?

final.mod1 <- cca(macros.log ~ dca1.dob + PC1.sh.O18 + co21 + char1, 
                  clim.complete3)
final.mod1
#               Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.5697     0.2114    4
# Unconstrained  2.1255     0.7886   51
# Inertia is scaled Chi-square 

summary(final.mod1)
# Accumulated constrained eigenvalues
# Importance of components:
#                        CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3584 0.1361 0.04461 0.03058
# Proportion Explained  0.6291 0.2389 0.07831 0.05367
# Cumulative Proportion 0.6291 0.8680 0.94633 1.00000

vif.cca(final.mod1)
# dca1.dob       co21 PC1.sh.O18      char1 
# 4.326258   5.878850   3.452229   1.110397 

spe.cca <- cca(macros.log ~., clim.complete3)
summary(spe.cca)
RsquareAdj(spe.cca)

set.seed(200)
cca.step.forward <- ordistep(cca(macros.log ~ 1, data = clim.complete3),
                             scope = formula(spe.cca),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod1.pars.cca <- cca(macros.log ~ dca1.dob + PC1.sh.O18 + co21 + char1,
                     data = clim.complete3)
summary(mod1.pars.cca)
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.5697     0.2114
# Unconstrained  2.1255     0.7886

# Accumulated constrained eigenvalues
# Importance of components:
#                        CCA1   CCA2    CCA3    CCA4   
# Eigenvalue            0.3584 0.1361 0.04461 0.03058
# Proportion Explained  0.6291 0.2389 0.07831 0.05367
# Cumulative Proportion 0.6291 0.8680 0.94633 1.00000

control <- how(within = Within(type = "series", mirror =TRUE))
(check2 <- check(macros.log, control))
anova(mod1.pars.cca, permutations = control)
#          Df ChiSquare      F  Pr(>F)   
# Model     4   0.56971 5.0258 0.00625 **
# Residual 75   2.12547                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod1.pars.cca, permutations = control, by = "axis") # CCA1 & CCA2 are important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.35840 12.6466 0.00625 **
# CCA2      1   0.13613  4.8034 0.00625 **
# CCA3      1   0.04461  1.5742 0.20625   
# CCA4      1   0.03058  1.0789 0.47500   
# Residual 75   2.12547      

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod1.pars.cca)
#$r.squared
# [1] 0.2113822
# 
# $adj.r.squared
# [1] 0.1695182

vif.cca(mod1.pars.cca)
# dca1.dob PC1.sh.O18       co21      char1 
# 4.326258   3.452229   5.878850   1.110397 

head(clim.complete4)
# clim.complete4 includes PC1.sh.temp, co21, loi1, char1, temp.et1, dca1.dob
# PC1.sh.temp antartic records + O18nzw
set.seed(200)
mbig <- cca(macros.log ~  ., clim.complete4)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete4)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + co21)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.sh.temp)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + temp.et1)
add1(m0, scope = formula(mbig), test = "permutation")
## -- included variables still significant?

final.mod2 <- cca(macros.log ~  dca1.dob + co21 + PC1.sh.temp + temp.et1, clim.complete4)
final.mod2
#               Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.5323     0.1975    4
# Unconstrained  2.1629     0.8025   51
# Inertia is scaled Chi-square 

summary(final.mod2)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3634 0.1075 0.03344 0.02787
# Proportion Explained  0.6828 0.2020 0.06283 0.05236
# Cumulative Proportion 0.6828 0.8848 0.94764 1.000


vif.cca(final.mod2)
# dca1.dob      co21    PC1.sh.temp    temp.et1 
# 6.525658    4.917441    5.424722    1.496844

spe.cca1 <- cca(macros.log ~., clim.complete4)
summary(spe.cca1)


RsquareAdj(spe.cca1)
# $r.squared
# [1] 0.2290534
# 
# $adj.r.squared
# [1] 0.1662909
# 
set.seed(200)
cca.step.forward1 <- ordistep(cca(macros.log ~ 1, data = clim.complete4),
                             scope = formula(spe.cca1),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod2.pars.cca <- cca(macros.log ~ dca1.dob + co21 + PC1.sh.temp + temp.et1, data = clim.complete4)
summary(mod2.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.5323     0.1975
# Unconstrained  2.1628     0.8025

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4
# Eigenvalue            0.3632 0.1080 0.03298 0.02816
# Proportion Explained  0.6822 0.2029 0.06196 0.05290
# Cumulative Proportion 0.6822 0.8851 0.94710 1.00000

anova(mod2.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     4   0.55015 4.809 0.00625 **
# Residual 75   2.14503                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod2.pars.cca, permutations = control, by = "axis") # CCA1 important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36343 12.6023 0.00625 **
# CCA2      1   0.10753  3.7287 0.08125 . 
# CCA3      1   0.03344  1.1596 0.65625   
# CCA4      1   0.02787  0.9663 0.55625   
# Residual 75   2.16291

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod2.pars.cca)
# $r.squared
# [1] 0.1974913
# 
# $adj.r.squared
# [1] 0.1551383

vif.cca(mod2.pars.cca)
# dca1.dob        co21 PC1.sh.temp    temp.et1 
# 6.525658    4.917441    5.424722    1.496844 

site.scrs.mod2.final <-scores(mod2.pars.cca, display = "sites")  # site scores cca
sp.scrs.mod2.final <-scores(mod2.pars.cca, display = "species")  #species scores cca

DomeC <-as.data.frame(dome.interp.temp)
DomeC = DomeC$temp.dome1

Vostok <-as.data.frame(vos.interp.temp)
Vostok = Vostok$temp.vos1

mod2.scores <- cbind(PC1.sh.temp, co2.interp, O18.nzw, Vostok, DomeC, sit.sc, dobson)
head(mod2.scores)

macros.spp.rel <- decostand(macros5, method = "total")

spe.nmds <- metaMDS(macros.spp.rel, distance = "bray")
spe.nmds
spe.nmds$stress
dev.new(title = "NMDS on fish species - Percentage difference",
        noRStudioGD = TRUE
)
plot(
  spe.nmds,
  type = "t",
  main = paste(
    "NMDS/Percentage difference - Stress =",
    round(spe.nmds$stress, 3)
  )
)

# Shepard plot and goodness of fit
dev.new(title = "NMDS - Shepard plot",
        width = 12,
        height = 6,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300)

sit.sc <- scores(spe.nmds)



dev.new()
Stratiplot(age ~ . - depth, data = chooseTaxa(mod2.scores, n.occ = 5, min.abun = 1),type = c("h","l","g"))
Stratiplot(age ~ . - depth, data = chooseTaxa(mod2.scores, n.occ = 5, min.abun = 1),type = c("l"), sort = "wa", varTypes="absolute")
##------------------------------------main grou

head(clim.complete2)
clim.complete6
head(clim.complete5)
set.seed(200)
mbig <- cca(macros.log ~  ., clim.complete5)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete5)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + co21 )
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.sh.temp1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + temp.et1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + loi1)
## -- included variables still significant?

final.mod3 <- cca(macros.log ~ dca1.dob + co21 + char1 + PC1.sh.temp1 + 
                    temp.et1 + loi1, clim.complete5)
final.mod3
#             Inertia Proportion
# Total          2.6952     1.0000     
# Constrained    0.6163     0.2287    6
# Unconstrained  2.0788     0.7713   51

summary(final.mod3)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5  CCA6
# Eigenvalue            0.3663 0.1134 0.05667 0.03299 0.02750 0.01949
# Proportion Explained  0.5942 0.1840 0.09194 0.05353 0.04462 0.03162
# Cumulative Proportion 0.5942 0.7783 0.87022 0.92376 0.96838 1.00000

vif.cca(final.mod3)
# dca1.dob         co21        char1 PC1.sh.temp1     temp.et1         loi1 
# 6.920244     4.960198     1.156780     3.807416     1.481613     2.399837 

spe.cca2 <- cca(macros.log ~., clim.complete5)
summary(spe.cca2)

RsquareAdj(spe.cca2)
#$r.squared
# [1] 0.2286826
# 
# $adj.r.squared
# [1] 0.1655333

set.seed(200)
cca.step.forward2 <- ordistep(cca(macros.log ~ 1, data = clim.complete5),
                              scope = formula(spe.cca2),
                              direction = "forward",
                              permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod3.pars.cca <- cca(macros.log ~ dca1.dob + co21 + PC1.sh.temp1, data = clim.complete5)
summary(mod3.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.4874     0.1808
# Unconstrained  2.2078     0.8192

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3   
# Eigenvalue            0.3629 0.09666 0.02782
# Proportion Explained  0.7446 0.19833 0.05709
# Cumulative Proportion 0.7446 0.94291 1.00000

anova(mod3.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     3   0.48741 5.5929 0.00625 **
# Residual 76   2.20777                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod3.pars.cca, permutations = control, by = "axis") # CCA1 is important component
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36256 12.4807 0.00625 **
# CCA2      1   0.09751  3.3567 0.11250   
# CCA3      1   0.02735  0.9413 0.66250   
# Residual 76   2.20777                      

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod3.pars.cca)
# $r.squared
# [1] 0.1808346
# 
# $adj.r.squared
# [1] 0.1487863

vif.cca(mod3.pars.cca)
# dca1.dob         co21 PC1.sh.temp1 
# 5.254385     4.606373     3.716952 


head(clim.complete5)
set.seed(200)
mbig <- cca(macros.log ~  ., clim.complete5)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete5)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + co21 )
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.sh.temp1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + temp.et1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + loi1)
## -- included variables still significant?

final.mod3 <- cca(macros.log ~ dca1.dob + co21 + char1 + PC1.sh.temp1 + 
                    temp.et1 + loi1, clim.complete5)
final.mod3
#             Inertia Proportion
# Total          2.6952     1.0000     
# Constrained    0.6163     0.2287    6
# Unconstrained  2.0788     0.7713   51

summary(final.mod3)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5  CCA6
# Eigenvalue            0.3663 0.1134 0.05667 0.03299 0.02750 0.01949
# Proportion Explained  0.5942 0.1840 0.09194 0.05353 0.04462 0.03162
# Cumulative Proportion 0.5942 0.7783 0.87022 0.92376 0.96838 1.00000

vif.cca(final.mod3)
# dca1.dob         co21        char1 PC1.sh.temp1     temp.et1         loi1 
# 6.920244     4.960198     1.156780     3.807416     1.481613     2.399837 

spe.cca2 <- cca(macros.log ~., clim.complete5)
summary(spe.cca2)

RsquareAdj(spe.cca2)
#$r.squared
# [1] 0.2286826
# 
# $adj.r.squared
# [1] 0.1655333

set.seed(200)
cca.step.forward2 <- ordistep(cca(macros.log ~ 1, data = clim.complete5),
                              scope = formula(spe.cca2),
                              direction = "forward",
                              permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod3.pars.cca <- cca(macros.log ~ dca1.dob + co21 + PC1.sh.temp1, data = clim.complete5)
summary(mod3.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.4874     0.1808
# Unconstrained  2.2078     0.8192

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3   
# Eigenvalue            0.3629 0.09666 0.02782
# Proportion Explained  0.7446 0.19833 0.05709
# Cumulative Proportion 0.7446 0.94291 1.00000

anova(mod3.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     3   0.48741 5.5929 0.00625 **
# Residual 76   2.20777                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod3.pars.cca, permutations = control, by = "axis") # CCA1 is important component
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36256 12.4807 0.00625 **
# CCA2      1   0.09751  3.3567 0.11250   
# CCA3      1   0.02735  0.9413 0.66250   
# Residual 76   2.20777                      

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod3.pars.cca)
# $r.squared
# [1] 0.1808346
# 
# $adj.r.squared
# [1] 0.1487863

vif.cca(mod3.pars.cca)
# dca1.dob         co21 PC1.sh.temp1 
# 5.254385     4.606373     3.716952 


head(clim.complete7)
set.seed(200)
mbig <- cca(macros.log ~  ., clim.complete7)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete7)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + co21)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.ant.temp)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + PC1.local.clim)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")

## -- included variables still significant?

final.mod4 <- cca(macros.log ~  co21 + PC1.ant.temp + PC1.local.clim + char1, clim.complete7)
final.mod4
#               Inertia  Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.5372     0.1993    4
# Unconstrained  2.1580     0.8007   51

summary(final.mod4)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5  CCA6
# Eigenvalue            0.3600 0.1092 0.04906 0.01886
# Proportion Explained  0.6702 0.2033 0.09132 0.03512
# Cumulative Proportion 0.6702 0.8736 0.96488 1.00000
vif.cca(final.mod4)
# co21        PC1.ant.temp PC1.local.clim       char1 
# 3.921674       3.095369       2.648063       1.113807 

spe.cca3 <- cca(macros.log ~., clim.complete7)
summary(spe.cca3)

RsquareAdj(spe.cca3)
#$r.squared
# [1] 0.2468679
# 
# $adj.r.squared
# [1] 0.1737699

set.seed(100)
cca.step.forward3 <- ordistep(cca(macros.log ~ 1, data = clim.complete7),
                              scope = formula(spe.cca3),
                              direction = "forward",
                              permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod3.pars.cca <- cca(macros.log ~ co21 + PC1.ant.temp + PC1.local.clim + char1, data = clim.complete7)
summary(mod3.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.4874     0.1808
# Unconstrained  2.2078     0.8192

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3   
# Eigenvalue            0.3629 0.09666 0.02782
# Proportion Explained  0.7446 0.19833 0.05709
# Cumulative Proportion 0.7446 0.94291 1.00000

anova(mod3.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     3   0.48741 5.5929 0.00625 **
# Residual 76   2.20777                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod3.pars.cca, permutations = control, by = "axis") # CCA1 is important component
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36256 12.4807 0.00625 **
# CCA2      1   0.09751  3.3567 0.11250   
# CCA3      1   0.02735  0.9413 0.66250   
# Residual 76   2.20777                      

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod3.pars.cca)
# $r.squared
# [1] 0.1808346
# 
# $adj.r.squared
# [1] 0.1487863

vif.cca(mod3.pars.cca)
# dca1.dob         co21 PC1.sh.temp1 
# 5.254385     4.606373     3.716952 











# Scaling 1: species scores scaled to the relative eigenvalues,
# sites are weighted averages of the species
dev.new()
par(mfrow = c(2,2))

## CCA triplots (using lc site scores)
# Scaling 1: species scores scaled to relative eigenvalues, 
# sites are weighted averages of the species
quartz(title = "CCA triplot - scaling 1 - lc scores")
plot(mod3.pars.cca, scaling = 1, display = c("sp","lc","cn"), 
     main = "Triplot CCA spe ~ clim6 - scaling 1")

# Default scaling 2: site scores scaled to relative eigenvalues, 
# species are weighted averages of the sites
quartz(title="CCA triplot - scaling 2 - lc scores")
plot(mod3.pars.cca, display = c("sp","lc","cn"), 
     main="Triplot CCA spe ~ clim6 - scaling 2")

# CCA scaling 1 biplot without species (using lc site scores)
quartz(title = "CCA biplot - scaling 1")
plot(mod3.pars.cca, scaling = 1, display = c( "sp", "cn"), 
     type = "text",
     main="Biplot CCA spe ~ clim6 - scaling 1")




spe.sc1 <- scores(mod3.pars.cca,
                  choices = 1:2, scaling = 1,
                  display = "sp")
arrows(0,0,
       spe.sc1[,1]*0.92,
       spe.sc1[,2]*0.92,
       length = 0,
       lty = 1,
       col = "darkgrey"
)


source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/triplot.rda.R")
spe.good <- goodness(mod3.pars.cca)
sel.sp <- which(spe.good[,2] >= 0.6) 
triplot.rda(mod3.pars.cca,
            site.sc ="lc",
            scaling = 1,
            cex.char2 = 0.7,
            pos.env = 3,
            pos.centr = 1,
            mult.arrow = 1.1,
            mar.percent = 0.05,
            select.spe = sel.sp
            )
triplot.rda(mod3.pars.cca,
            site.sc = "lc",
            scaling = 2,
            ex.char2 = 0.7,
            pos.env = 3,
            pos.centr = 1,
            mult.arrow = 1.1,
            mar.percent = 0.05,
            select.spe = sel.sp
            )


# Manual model building 
# maximal model



# clim.complete explains about 24% of the total variance
# while cca2 and cca2 explain 57% and 22% of the total constrained variability
vif.cca(comp.model)
# vif values are all below 20
# PC2      char2    temp.et2 loi2     co22     O28.nza 
# 3.029563 2.222625 2.529972 2.920847 5.325257 2.520389 

control <- how(within = Within(type = "series", mirror = FALSE))
(check2 <- check(macros.log, control))
summary(check2)

(anova(comp.model,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0225
(anova(comp.model, by="term", permutations = control))
# significant env.var under restricted permutation PC2 composite temperature significant at 0.05 level (p-value = 0.0225) 
# & char marginally significant(0.0750) 
(anova(comp.model, by="axis", permutations = control))
# CCA2 & CCA2 are significant at 0.05 level, p-value = 0.0225 in both cases

# Running the model with only significant variables
set.seed(60)
final.mod <- cca(macros.log ~ Comp.2 + co22 + char2, data = clim.complete2)
final.mod
summary(final.mod)
# final.mod explains about 27% of the total variability and CCA2 explains about 76% of the constrained variability
#  CCA2 explains about % 27
# R2adjuted is around 24%
RsquareAdj(final.mod)

(anova(final.mod,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0225
(anova(final.mod, by = "term", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0225), char (0.05) and CO2 (0.025)
(anova(final.mod, by = "axis", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0225), char (0.05) and CO2 (0.025)
(anova(final.mod, by = "margin", permutations = control))
vif.cca(final.mod)

ource("/Users/giselleastorga/Documents/Thesis_GA/second_chapter/data_second/cleanplot.pca.R")
dev.new(width = 22,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(2, 2))
cleanplot.pca(temp.rda, scaling = 2, mar.percent = 0.08)
cleanplot.pca(temp.rda, scaling = 2, mar.percent = 0.04)


set.seed(200)
head(clim.complete3)
mod0 <- cca(macros.log ~ 2, data = clim.complete3)  # Model with intercept only
summary(mod0)
mod2 <- cca(macros.log ~ ., data = clim.complete3)  # Model with all explanatory variables
summary(mod2)
model.fit2 <- ordistep(mod0, scope=formula(mod2),direction = "forward")
summary(model.fit2)
vif.cca(model.fit2)

(anova(model.fit2,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0225
(anova(model.fit2, by = "term", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0225), char (0.05) and CO2 (0.025)
(anova(model.fit2, by = "axis", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0225), char (0.05) and CO2 (0.025)
(anova(model.fit2, by = "margin", permutations = control))
vif.cca(model.fit2)

# Install packages for ggplots
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# Libraries requiered for graphics
library("gridExtra")
library("ggpubr")
library("ggplot2")
library("devtools")
library("easyGgplot2")
library(moments)

# Histograms environmental variables
is.data.frame(antar.temp)
dev.new()
hist(antar.temp$temp.dome1)
hist(antar.temp$temp.vos1)

skewness(antar.temp$temp.dome1)# kutosis = -0.9 and kurtosis = 3.4
kurtosis(antar.temp$temp.dome1)

# Histogram from a single numeric vector 
# ggplot2.histogram(data=numVector)
# Basic histogram plot from the vector "weight"
(range(temp.dome1))

dev.new()
par(mfrow=c(2, 2))
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome1', scale = "density", binwidth = .7)
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos1', scale = "density", binwidth = .7)
grid.arrange(dome, vos, nrow = 2) 
# Change the width of bars
# Change y axis values to density
dev.new()
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos2', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                         meanLineType = "dashed", meanLineSize = 2)
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome2', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                          meanLineType = "dashed", meanLineSize = 2)
grid.arrange(dome, vos, nrow = 2) 

# Add density curve
vos2 <- ggplot2.histogram(data = antar.climate2, xName = 'temp.vos2',
                  fill = "white", color = "black",
                  addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)

dome2 <- ggplot2.histogram(data = antar.climate2, xName = 'temp.dome2',
                          fill = "white", color = "black",
                          addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)
dev.new()
grid.arrange(dome2, vos2, nrow = 2) 

NZ.O28 <- ggplot2.histogram(data = antar.climate2, xName = 'O28',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = 0.02)
NZ.C23 <- ggplot2.histogram(data = antar.climate2, xName = 'C23',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = 0.06)

dev.new()
grid.arrange(NZ.O28, NZ.C23, nrow = 2, ncol = 2) 
grid.arrange(vos2, dome2, NZ.O28, NZ.C23, nrow = 4, ncol = 2)

NZW.O28 <- ggplot2.histogram(data = antar.climate5, xName = 'O28',
                           fill = "white", color = "black",
                           addDensityCurve = TRUE, densityFill = '#FF6666', bindwidth = 0.05)
NZW.C23 <- ggplot2.histogram(data = antar.climate5, xName = 'C23',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', bindwidth = 0.05)
dev.new()
grid.arrange(NZW.O28, NZW.C23,vos2,dome2, nrow = 4, ncol = 2) 

skewness(antartic$temp.vos2)# kutosis = -0.9 and kurtosis = 3.4
kurtosis(antartic$temp.vos2)
head(antartic)


##****************CCA**********************************##
# CCA with PCAs of temp.pca and log2p(macros); 80 columns 52 variables
macros5 <- read.csv("combined_macros_shorted(5).csv")
str(macros5)
head(temp.pca)
macros.log <- log2p(macros5)

head(clim.complete2)
clim.complete3 <- subset(clim.complete2, select = -c(7))
head(clim.complete3)
set.seed(60)
mod0 <- cca(macros.log ~ 2, data = clim.complete3)  # Model with intercept only
summary(mod0)
mod2 <- cca(macros.log ~ ., data = clim.complete3)  # Model with all explanatory variables
summary(mod2)
model.fit <- ordistep(mod0, scope = formula(mod2), direction = "forward")
summary(model.fit)
# Automatic variable selection indicates both components are important
# PC2 and PC2 Antartic variation in deuterium explains about 22% of the total variability with
# CCA 2 explaining about 83% of the constrained varibility
# CCA 2 explains about 27% of the constrained variability
# Free permutation of our cca results is highly significant, but data is temporaly ordered
# so that restricted permutations are needed (e.g. dataset with 80 stratigraphic ordered samples,
# has just 80 posible permutations. I generated the complete set of permutations by suffleSet() instead of 
# a random set. This is also posible because dataset is relatively small.
(anova(model.fit))
control <- how(within = Within(type = "series", mirror = FALSE))
(check2 <- check(macros.log, control))
summary(check2)
(anova(model.fit,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0225
(anova(model.fit, by="term", permutations = control))
# only PC2 is significant at 0.05 (p-value = 0.0225), char (0.05) and CO2 (0.025)
(anova(model.fit, by="axis", permutations = control))
# CCA2 is significant at 0.05 level, p-value = 0.0225 

# antar.climate3 include Antartic tmeperature data + O.28 from Mt. Arthur
head(antar.climate3)
climate3.pca <- prcomp(antar.climate3, scale. = T, center = T)
climate3.pca 
summary(climate3.pca)

climate3.rda <- rda(antar.climate3, scale = T, center = T)
climate3.rda 
summary(climate3.rda)
# PC2 65%; PC2 26% 
# select the data frame with eigenvalues of particular axes:
(ev <- climate3.pca$CA$eig)#Eigenvalues
# Apply Kaiser-Guttman criterion to selected axes
ev[ev > mean(ev)]
# Broken stick model
source("evplot.R")

dev.new()
plot(climate3.pca, type = "l",main="")
ev <- climate3.pca$sdev^2
pdf("ev.pdf")
dev.new()
evplot(ev)
# only PC2 is significant
# To find the eigenvalues for each PC 
(ev <- climate3.pca$sdev^2)
climate3.pca$sdev^2 #eigen values, variance in each direction
climate3.pca$sdev^2/sum(climate3.pca$sdev^2)# = proportion of variance vector

pdf("ev.pdf")
dev.new()
evplot(ev)


scores.climate3.pca <- scores(climate3.pca, display = "sites")
write.csv(scores.climate3.pca, file="scores.climate3.pca.csv")
climate3 <- read.csv("scores.climate3.pca.csv")
head(climate3)
climate3.a <-subset(climate3, select = -c(2))# leaving only PC2-PC3

head(climate3.a)
climate3.b <-subset(climate3, select = -c(2,3,4))# leaving only PC2-PC2
head(climate3.b)

set.seed(600)
sh.temp.mod <- cca(macros.log ~ PC2, data = climate3.b)  # Model with intercept only
summary(sh.temp.mod)

# sh.temp explained about 22% of the total variability 
# using as input matrix only PC2

(anova(sh.temp.mod ,  permutations = control))
# only PC2 is important 
(anova(sh.temp.mod , by="term", permutations = control))
# only cca2 is important. 
(anova(sh.temp.mod , by="axis", permutations = control))
# cca2 & cca2 are significant after restricted permutations


# Last two lines show the comparison between real variation 
# represented by individual PCA axes, and relevant variation 
# calculated by broken-stick model (2 in the column with particular 
# axis means that this axis explains more than would explain the axis of the same order in a null model)
# ith a bit effort, the barplot visualizing relationship between variation represented by individual PCA 
# axis and variation explained by broken-stick model (the line % > bs% from the table above) can be drawn:
# dev.new()
# barplot (sig.rand[c('percentage of variance', 'broken-stick percentage'), ], beside = T, 
#          xlab = 'PCA axis', ylab = 'explained variation [%]', col = c('grey', 'black'), 
#          legend = TRUE)

#Exploring the plant macrofossil data with CA & DCA
#****************************************************************************************************
##Data used in paper: macros shorted(5).csv. 52 columns, 80 rows at species higher than 5% 

macros5<-read.csv("combined_macros_shorted(5).csv")
macros5
str(macros5)
colSums(macros5)# only two species have abundance fraction of 5
rowtotals <- rowSums(macros5)
coltotals <- colSums(macros5)
(meanrow <- mean(rowtotals))
(sdrow <- sd(rowtotals))
(CV <- (100*sdrow/meanrow))# CV = 92.6 meaning that relativisation would have a moderate effect on results

(meancol <- mean(coltotals))
(sdcol <- sd(coltotals))
(CV <- (100*sdcol/meancol))# CV = 297

macros.log <- log1p(macros5)


DCA.macros <- decorana(macros5, iweigh = 1, ira = 0) 
summary(DCA.macros)
DCA.macros$fraction #downweighting abundance fraction starts at 5 
cca(macros.log)
dca.samplescores <- scores(DCA.macros, display = c("sites"), choices = 1:2)
write.csv(dca.samplescores, file = "dca.samplescores.csv")

DCA_species <-scores(DCA.macros, display=c("species"))

# extracts only axis 2 scores for samples

#macros<-read.csv("combined_macros_shorted.csv")##macros shorted (47 columsn) at species higher than 5% 
#str(macros)
#macros.slog<-log2p(macros)
#DCA.macros2 <- decorana(macros.slog) 
#summary(DCA.macros2)
#cca(macros.slog)

dev.new()
pdf("DCA.sites.pdf")
ordiplot (DCA.macros, display = 'si', type = 'n')
points (DCA.macros, col = zones$group, pch = zones$group )
for (i in seq (2, 5)) ordihull (DCA.macros, groups = zones$group, show.groups = i, col = i, lty = 'dotted')
dev.off()

zones <- read.csv("zones.csv")
head(zones)
dev.new()
# pdf("DCA sites.pdf")
ordiplot (DCA.macros, display = 'si', type = 'n')
points (DCA.macros, col = zones$group, pch = zones$group )

dev.new()
# pdf("DCA.sites.pdf")
ordiplot (DCA.macros, display = 'si', type = 'n')
for (i in seq (2, 5)) ordispider (DCA.macros, groups = zones$group, show.groups = i, col = i, label = T)
for (i in seq (2, 5)) ordihull (DCA.macros, groups = zones$group, show.groups = i, col = i, lty = 'dotted')
dev.off()

source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')


#DCA
#****************************************************************************************************

DCA.macros <- decorana(macros.log) #performs DCA and downweights rare taxa
summary(DCA.macros)
macros_scores_dca<-scores(DCA.macros)
macros_scores_DCA<-DCA.macros$rproj
DCA_species<-scores(DCA.macros, display=c("species"))
DCA_sites<-scores(DCA.macros, display=c("sites"))
write.csv(macros_scores_DCA, file = "macros_scores_DCA.csv")
shnam <- make.cepnames(names(macros_log5))



dev.new()
plot(DCA.macros)

dev.new()
pdf("DCA species.pdf")

plot(DCA.macros, display=c("none"), cols=c(2,2))
# plot axis 2 & 2, but display no points or labels
text(DCA.macros, display="sites", pcol = "red", cex=0.7)
# plot samples as green crosses for axis 2 and 2
text(DCA.macros, display=c("species"), choices=2:2,
       cex=0.7)
# plot labels for taxa for axis 2 and 2, using cex to shrink
# size of labels. Larger plots may also be used to alleviate
# congestion of labels.
dev.off()



##to abbreviate Latin names
shnam <- make.cepnames(names(macros.log))
DCA.macros<-decorana(env2) 
###priority to the most abundant species
dev.new()
stems <- colSums(macros_log5)
plot(DCA.macros, dis="sites")
plot(DCA.macros, dis="species")
plot(DCA.macros, dis="species",   pcol = "red", pch="+")



plot(CA.macros, dis="sp", type="n")
sel <- orditorp(CA.macros, dis="species",   pcol = "red", pch="+")
plot(CA.macros,type="n",main="Correspondence Analysis") #creates empty plot(type=”n”)
text(CA.macros,display="species",col="red",cex=.6) #adds red labels for species
points(CA.macros,display="sites",pch=20)



#****************************************************************************************************
#  Correlation , scaling=2 
# (focus on correlation among species/variables, reflected in angle of particular vectors)
pdf("Figure_correlplot.pdf",  width = 20, height = 27, units = 'cm', res = 200)
dev.new()
ordiplot(model_fit2, scaling=2,main="Correlation",type="n")
segments(x0=0,y0=0,x2=scores(model_fit2, display="species", scaling=2)[,2],
y2=scores(model_fit2, display="species", scaling=2)[,2])
text(model_fit2, display="sp", scaling=2, col=2, pch=22, cex=.7)
text(model_fit2, display="bp", scaling=2,
row.names(scores(model_fit2, display="bp")), col=4)
text(model_fit2, display=c("sites"),pch=25,cex=.5, scaling=2,labels=rownames(macros_log5))
cor(macros.log,env2)
dev.off()

# Correlation , scaling=2
plot(model_fit2, scaling=2,main="Correlation",type="n")
segments(x0=0,y0=0,x2=scores(model_fit2, display="species", scaling=2)[,2],
         y2=scores(model_fit2, display="species", scaling=2)[,2])
text(model_fit2, display="sp", scaling=2, col=2)
text(model_fit2, display="bp", scaling=2,row.names(scores(model_fit2, display="bp")), col=4)
text(model_fit2, display=c("lc"), scaling=2,labels=rownames(macros_log5))

# Distance , scaling=2
eps("Figure_distplot.eps",  width = 20, height = 27, units = 'cm', res = 200)
x22()
plot(model_fit2, scaling=2,main="Distance triplot",type="n")
segments(x0=0,y0=0,x2=scores(model_fit2, display="species", scaling=2)[,2],
y2=scores(model_fit2, display="species", pch=2, cex=.2, scaling=2)[,2])
text(model_fit2, display="sp", scaling=2, col=2)
text(model_fit2, display="bp", scaling=2,row.names(scores(model_fit2, display="bp")), col=4)
text(model_fit2, display=c("sites"),pch=2,cex=.8, scaling=2,labels=rownames(macros_log5))
dev.off()




scores.macros_cca_model_fit2<-scores(model_fit2,display="sites")##scores cca
scores.macros_cca_spp_model_fit2<-scores(model_fit2,display="species")##scores cca
write.csv(scores.macros_cca_spp_model_fit2, file = "scores.macros_cca_spp_model_fit2.csv")



##ordistep(macros_cca_temp_stand_fire_3~2,direction="both",step=999, max.perm=999)
ordistep(cca(macros_log5 ~  2, env.var), reformulate(names(env.var)), direction="both",pstep=2000, step=2000,max.perm=2000)
ordistep(d, scope = formula(d),direction="both",step=2000, max.perm=2000)


## End(Not run)
## Manual model building
## -- define the maximal model for scope
mbig <- cca(macros_log5 ~  ., climatics_complete)
## -- define an empty model to start with
m0 <- cca(macros_log5 ~ 2, climatics_complete)
## -- manual selection and updating
add2(m0, scope=formula(mbig), test="perm")
m0 <- update(m0, . ~ . + temperature_all)
add2(m0, scope=formula(mbig), test="perm")
m0 <- update(m0, . ~ . + co2_dome)
add2(m0, scope=formula(mbig), test="perm")
## -- included variables still significant?
drop2(m0, test="perm")
add2(m0, scope=formula(mbig), test="perm")
ordistep(cca(macros_log5 ~  2, climatics_complete), reformulate(names(climatics_complete)), perm.max=2000)



#******************************RAREFACTION************************************************
macros<-read.csv("combined_rarefaction.csv")
modern<-read.csv("modern_sediment.csv")
raremax <- min(rowSums(macros))
raremax
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
out <- with(pars[2:47, ],
            rarecurve(macros, step = 20, sample = raremax, col = col,
                      lty = lty, label = FALSE))

mrar <- rarefy(macros, min(rowSums(macros)))

data(modern)
test<-rich(matrix=modern, nrandom=499,verbose=TRUE)
test$cr # observed cumulative species richness
[2] 25
test$mr # observed mean value of species richness over the n samples
[2] 8.25


test<-rich(matrix=macros, nrandom=499,verbose=TRUE)
test$cr 
[2] 47  # observed cumulative species richness
test$mr # observed mean value of species richness over the n samples
[2]  22.94
data(macros) # culture plot
rare_macros<-raref(matrix=macros, dens=sum(macros), nrandom=500)
rare_macros$Sinterp[2]


data(BCI, package = "vegan")
BCI2 <- BCI[2:26, ]
raremax <- min(rowSums(BCI2))
raremax




#######################Cluster analysis################################
# Bray curtis cluster using macros2 (80 samples, 52 species)
macros2<- read.csv("combined_macros_shorted(5).csv")
str(macros2)
macrosper2 <- percenta(macros2,first=2,last=52)
str(macrosper2)
write.csv(macrosper2, file = "macrosper2.csv")
macros2<-read.csv("macrosper2.csv")
## Bray-Curtis distances between samples, even if method is not especify
diss2 <- vegdist(macrosper2, method="bray")
clust2<- chclust(diss2, method="coniss")
dev.new()
bstick(clust2,20)

diss2 <- vegdist(sqrt(macrosper2/200), method = "bray")
clust2 <- chclust(diss2, method="coniss")
dev.new()
bstick(clust2)# points above the red line suggest 3 zones
spe.chclust2<-chclust(vegdist(macrosper2))
dev.new()
bstick(spe.chclust2,20)
k<-3
(gr4<-cutree(spe.chclust2, k = k))
dev.new()
plot(spe.chclust2, hang =-2, main = "CONISS clustering Lake Dobson")
rect.hclust(spe.chclust2, k = k)
## Bray curtis cluster using combined_macros_shorted at 5.csv (80 samples, 47 species)
macros2<- read.csv("combined_macros_shorted.csv")
str(macros2)
macrosper2 <- percenta(macros2,first=2,last=47)
str(macrosper2)

## Bray-Curtis distances between samples, even if method is not especify
diss2 <- vegdist(macrosper2, method="bray")
clust2<- chclust(diss2, method="coniss")
dev.new()
bstick(clust2,20)
diss2 <- vegdist(sqrt(macrosper2/200), method = "bray")
clust2 <- chclust(diss2, method="coniss")
dev.new()
bstick(clust2,20)# points above the red line suggest 3 zones

## Bray curtis cluster using macros_full_30_Dec_226 (80 samples, 64 species)
macros3<- read.csv("macros_full_30_Dec_226.csv")
str(macros3)
macrosper3 <- percenta(macros3,first=2,last=64)
str(macrosper3)
write.csv(macrosper3, file="macrosper3.csv")
macrosper3<-read.csv("macrosper3.csv")

diss3 <- vegdist(macrosper3, method="bray")
clust3<- chclust(diss3, method="coniss")
dev.new()
bstick(clust3,20)

diss3 <- vegdist(sqrt(macrosper3/200), method = "bray")
clust3 <- chclust(diss3, method="coniss")
dev.new()
bstick(clust3)# points above the red line suggest 3 zones

## Bray curtis cluster using macros_percenta_clu.csv (88 samples, 47 species)
macros4<- read.csv("macros_percenta_clu.csv")
str(macros4)
## Bray-Curtis distances between samples, even if method is not especify
diss4 <- vegdist(macros4, method="bray")
clust4<- chclust(diss4, method="coniss")
dev.new()
bstick(clust4)
diss4 <- vegdist(sqrt(macros4/200), method = "bray")
clust4 <- chclust(diss4, method="coniss")
dev.new()
bstick(clust4)# points above the red line suggest 3 zones

## Bray curtis cluster using macros_full_30_Dec_226 (80 samples, 62 species)
macros5<- read.csv("macros_full_27_may_2026.csv")
str(macros5)
macrosper5 <- percenta(macros5,first=2,last=62)
str(macrosper5 )
## Bray-Curtis distances between samples, even if method is not especify
diss5 <- vegdist(macrosper5, method="bray")
clust5<- chclust(diss5, method="coniss")
dev.new()
bstick(clust5)

diss5 <- vegdist(sqrt(macrosper5/200), method = "bray")
clust5 <- chclust(diss5, method="coniss")
dev.new()
bstick(clust5)# points above the red line suggest 3 zones

macros6<- read.csv("macros.full.percenta.csv")
str(macros6)
## Bray-Curtis distances between samples, even if method is not especify
diss6 <- vegdist(macros6, method="bray")
clust6<- chclust(diss6, method="coniss")
dev.new()
bstick(clust6)

diss6 <- vegdist(sqrt(macros6/200), method = "bray")
clust6 <- chclust(diss6, method="coniss")
dev.new()
bstick(clust6)# points above the red line suggest 3 zones


#88 columns, 47 species + depth+age
macros7<- read.csv("macrosper.csv")
str(macros7)
depth<-macros7[48:49]
macrosper7<-macros7[2:47]
## Bray-Curtis distances between samples, even if method is not especify
diss7 <- vegdist(macrosper7, method="bray")
clust7<- chclust(diss7, method="coniss")
dev.new()
bstick(clust7,20)

diss7 <- vegdist(sqrt(macrosper7/200), method = "bray")
clust7 <- chclust(diss7, method="coniss")
dev.new()
bstick(clust7)# points above the red line suggest 3 zones

#macros_percenta 
macros8<- read.csv("macrosp.csv")
str(macros8)
depth<-macros8[63:64]
str(depth)
macrosper8<-macros8[2:62]
str(macrosper8)
## Bray-Curtis distances between samples, even if method is not especify
diss8 <- vegdist(macrosper8, method="bray")
clust8<- chclust(diss8, method="coniss")
dev.new()
bstick(clust8,20)

diss8 <- vegdist(sqrt(macrosper8/200), method = "bray")
clust8 <- chclust(diss8, method="coniss")
dev.new()
bstick(clust8)


#macros_percenta 
macros9<- read.csv("macros_cluster.csv")
str(macros9)
depth<-macros9[48:49]
str(depth)
macros9<-macros9[2:47]
str(macrosper9)
macrosper9 <- percenta(macros9,first=2,last=47)
str(macrosper9 )

## Bray-Curtis distances between samples, even if method is not especify
diss9 <- vegdist(macrosper9, method="bray")
clust9<- chclust(diss9, method="coniss")
dev.new()
bstick(clust9,20)

diss9 <- vegdist(sqrt(macrosper9/200), method = "bray")
clust9 <- chclust(diss9, method="coniss")
dev.new()
bstick(clust9)


# Basic diagram
dev.new()
plot(clust_bray, hang=-2)
# Rotated through 90 degrees
dev.new()
plot(clust_bray, hang=-2, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
dev.new()
plot(clust_bray, xvar=macros_clu$age, hang=-2, horiz=TRUE, x.rev=TRUE)
plot(clust_bray, xvar=macros_clu$depth, hang=-2, horiz=FALSE, x.rev=FALSE)
bstick(clust_bray)
plot(clust_bray,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 6)


depth<-macros_cluster[48:49]
str(depth)
macrosclu <- macros_cluster[2:47]
str(macrosclu)
macroscluper <- percenta(macrosclu,first=2,last=47)
str(macroscluper)

## Bray-Curtis distances between samples, even if method id not especify
diss2 <- vegdist(macrospershort, method="bray")
clust2<- chclust(diss2, method="coniss")
dev.new()
bstick(clust2)
dev.new()
plot(clust2,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 4)

diss2 <- vegdist(sqrt(macrospershort/200), method = "bray")
clust2 <- chclust(diss2, method="coniss")
dev.new()
bstick(clust2)# points above the red line suggest 3 zones
dev.new()
plot(clust_bray,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 3)



# Basic diagram
dev.new()
plot(clust_bray, hang=-2)
# Rotated through 90 degrees
dev.new()
plot(clust_bray, hang=-2, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
dev.new()
plot(clust_bray, xvar=macros_clu$age, hang=-2, horiz=TRUE, x.rev=TRUE)
plot(clust_bray, xvar=macros_clu$depth, hang=-2, horiz=FALSE, x.rev=FALSE)
bstick(clust_bray)
plot(clust_bray,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 6)


#SQRT cluster wirth macros short (80 samples)
macros<-read.csv("macros_percenta5.csv")
diss <- dist(sqrt(macros/200))
clust_sqrt <- chclust(diss, method="coniss")
dev.new()
bstick(clust_sqrt)
dev.new()
plot(clust_sqrt,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
plot(clust_sqrt, xvar=macros$depth, hang=-2, horiz=FALSE, x.rev=FALSE)
rect.hclust(clust_sqrt, 5)
#whole data samples with non macros (22) were removed
macros_clu<-read.csv("macros_cluster.csv")
str(macros_clu)
macros_percenta_clu<-percenta(macros_clu,first=2,last=47)
write.csv(macros_percenta_clu, file = "macros_percenta_clu.csv")
macros_percenta_clu<-read.csv("macros_percenta_clu.csv")
diss_sqrt <- dist(sqrt(macros_percenta_clu/200))
clust_sqrt <- chclust(diss_sqrt, method="coniss")
dev.new()
bstick(clust_sqrt)
plot(clust_sqrt,hang=-2,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_sqrt, 5)

# Basic diagram
plot(clust, hang=-2)
# Rotated through 90 degrees
#plot(clust, hang=-2, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
#plot(clust, xvar=RLGH$depths$Depth, hang=-2, horiz=TRUE, x.rev=TRUE)

plot(clust,hang=-2)
plot(clust,hang=-2, horiz=TRUE ,x.rev=TRUE)
plot(clust,hang=-2, horiz=TRUE ,x.rev=TRUE,cex=.6,xvar=macros$depth,main="Cluster Analysis of Lake Dobson")

plot(clust,hang=-2,cex=.6,xvar=macros$Depth,main="Cluster Analysis of Lake Dobson") 
rect.hclust(dud.clust2, 6)

macros<-read.csv("macros.csv")
dud.dist2 <- vegdist(macros) 
dud.clust2 <- chclust(dud.dist2) 
par(mfrow=c(2,2)) 
plot(dud.clust2,hang=-2,cex=.6) 
plot(dud.clust2,hang=-2,cex=.6)
rect.hclust(clust, 6)


macros.norm<-decostand(macros,"normalize")
macros.ch<-vegdist(macros.norm,"euc")
macros.ch.coniss<-hclust(macros.ch,method="coniss")
# Rotated through 90 degrees
plot(macros.ch.single)
macros.ch.complete<-hclust(macros.ch,method="complete")
plot(macros.ch.complete)
macros.ch.complete<-hclust(macros.ch,method="complete")
plot(macros.ch.complete)
macros.ch.centroid<-hclust(macros.ch,method="centroid")
plot(macros.ch.centroid)
macros.ch.ward<-hclust(macros.ch,method="ward.D2")
plot(macros.ch.ward)

# Histogram plot with multiple groups
# Multiple histograms on the same plot
# Color the histogram plot by the groupName "sex"
# ggplot2.histogram(data=weight, xName='weight',
#                   groupName='sex', legendPosition="top")
# # Histogram plots with semi-transparent fill.
# # alpha is the transparency of the overlaid color
# ggplot2.histogram(data=weight, xName='weight',
#                   groupName='sex', legendPosition="top",
#                   alpha=0.5 )
# # Histogram plots with mean lines
# ggplot2.histogram(data=weight, xName='weight',
#                   groupName='sex', legendPosition="top",
#                   alpha=0.5, addDensity=TRUE,
#                   addMeanLine=TRUE, meanLineColor="white", meanLineSize=2.5)


