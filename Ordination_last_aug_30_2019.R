rm(list=ls())
# setwd("/Users/Gastorga/Dropbox/chapter2_15_may_2015/data_second")
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
str(macros5)
summary(macros5)
macros.log <- log1p(macros5)# computes log(1+x)
decorana(macros5, iweigh = 1, ira = 0)
# DCA performed on log(raw count data). 80 samples 52 species including individual conifer abundances
# argument ira = 0 for detrending and iweight = 1 for downweighting of rare species
# DCA1 2.9 DCA2 2.2
# If raw abundances are used DCA1 = 3.7; DCA2 = 2.1

spe.pca <- prcomp(macros.log, scale. = T, center = T) 
spe.pca
summary(spe.pca)
spe.rda <- rda(macros.log, scale = TRUE, center = TRUE)
spe.rda
summary(spe.rda)

#-------------------------------------Interpolate---------------------------------------------------#
# a try with antartic data with AICC2012
dobson <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/dobson.csv")
head(dobson)
dob.depth <- dobson$depth
dob.age <- dobson$age

# All the records are presented in the AICC2012 chronology
# Composite CO2 curve first row correspond to EDML core from Monin et al 2004; 
# Siegenthaier et al 2005 (in antarctica2015co2.xls) 
# The following ages are from EDC CO2 composite curve of Lüthi et al., 2008 transfered on 
# AICC2012 (icecores_data_on_AICC2012.xlsx). I could use the whole chronology proposed for EDC in 
# 2015 version of the AICC2012 chronology, but we will see the fit with the composite 
# chronology is better.

antartic <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Antartic_data.csv")
head(antartic)
age.dome <- antartic$age_dome
temp.dome1 <- antartic$temp_dome
temp.dome2 <- antartic$temp_dome


dome.temp <- cbind(temp.dome1, temp.dome2)
dome.interp.temp <-
  interp.dataset(
    y = dome.temp,
    x = age.dome,
    xout = dob.age,
    rep.negt = F
  )
head(dome.interp.temp)
dome.interp.temp <- subset(dome.interp.temp, select = -c(2))

age.vos <- antartic$age_vos

temp.vos1 <- antartic$deltaTS
temp.vos2 <- antartic$deltaTS
vos.temp <- cbind(temp.vos1, temp.vos2)
vos.interp.temp <-
  interp.dataset(
    y = vos.temp,
    x = age.vos,
    xout = dob.age,
    rep.negt = F
  )
head(vos.interp.temp)
vos.interp.temp <- subset(vos.interp.temp, select = -c(2))
head(vos.interp.temp)

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

#---------------------------------Eagle Tarn interpolation-----------------------------------------#

# interpolated to Dobson ages
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

#------------------------------Lake Dobson Chironomids interpolation-------------------------------#
# I should not interpolated to the LD age, because chironomids are from the same core as the macros
# In the original data from A. Rees depth 945 had an age of 14879, but I calculate my own 
# age-depth model. I could interpolated to the LD depths of plant macros
 
chiro_dob <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Chiro_Dob.csv")
head(chiro_dob)
age.chiro <- chiro_dob$age
depth.chiro <- chiro_dob$depth
dca1.dob <- chiro_dob$DCA1
dca1.dob2 <- chiro_dob$DCA1
dca.dob <- cbind(dca1.dob, dca1.dob2)
chiro.interp <- interp.dataset(y = dca.dob , x = depth.chiro, xout = dob.depth)
head(chiro.interp)
chiro.interp <- subset(chiro.interp, select = -c(2))
head(chiro.interp)

#---------------------------------------NZ isotopes interpolation----------------------------------#
# Speleothem data from NW South Island NZ Williams 2005
# Updated 10/03/2008

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
# write.csv(NZW.interp, file = "NZW.interp.csv")
# NZW.interp <- read.csv("NZW.interp.csv")
O18.nzw <- NZW.interp$O18.nzw
C13.nzw <- NZW.interp$C13.nzw

# NZ.stand <- decostand(NZ.interp, method = "standardize")
# head(NZ.stand)
# Oxygen isotope variation in Mt. Arthur speleothems primarily represents changes in 
# meteoric waters falling above the caves, possibly responding to latitudinal changes
# in the position of the Subtropical Front in the Tasman Sea. How is the STF in Tasmania compared
# to the NZ one? 
# Carbon isotope variations in the speleothems record, represent changes in forest productivity,
# closely matching existing paleovegetation records in NZ according to Hellstrom 1998.


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
# write.csv(NZA.interp, file = "NZA.interp.csv")
# NZA.interp <- read.csv("NZA.interp.csv")
O18.nza <- NZA.interp$O18.nza
C13.nza <- NZA.interp$C13.nza

#----------------------------------------charcoal interpolation------------------------------------#
# Charcoal data is from the same core as tha plant macros so the 2 records should have same age.
# However, counts of charcoal particles were made every 1 cm starting at 0.5 cm and plant macros 
# starting at 1 cm and counted every 10 cm. I interpoletated to LD depth.

char <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/charcoal.csv")
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
        x = depth.char,
        xout = dob.depth,
        rep.negt = F
    )
head(char.interp)

#-------------------------------Loss on Ignision interpolation-------------------------------------#
# loi data interpolated to plant macros depth
loi <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/loi.csv")
head(loi)
depth.loi <- loi$depth
age.loi <- loi$age
loi1 = loi$loi
loi2 = loi$loi
lois <- cbind(loi1, loi2)
head(lois)
lois <- subset(lois, select = -c(2))
loi.interp <-
    interp.dataset(
        y = lois,
        x = depth.loi,
        xout = dob.depth,
        rep.negt = F
    )
# loi_interpolated <- cbind(loi.interp, dob.age)
head(loi.interp)

#-------------------------------------Data exploration---------------------------------------------#
# assembling environmental variables based on Lake Dobson age scale
# clim.complete with antartic and NZ records (NZW and NZA) not joined in PCA

clim.complete <- cbind(dome.interp.temp, vos.interp.temp, 
                       char.interp, et.interp, chiro.interp, loi.interp, co2.interp, 
                       NZA.interp, NZW.interp, dob.age)
head(clim.complete)

# Patterns of correlations among env.var
(cor(clim.complete, use = "pairwise.complete.obs"))

antar.temp <- cbind(vos.interp.temp, dome.interp.temp)
head(antar.temp)
is.data.frame(antar.temp)
antar.temp <- as.data.frame(antar.temp)
is.data.frame(antar.temp)
# PCA based on a covariance matrix because variables are in the same scale 
# not sure if I should center
# Arguments scale = TRUE and centre = TRUE calls for a standardization of the variables
# antar.temp include the Vostok and DomeC variacion of temperature (delta TS)
head(antar.temp)
cor(antar.temp)
temp.pca <- prcomp(antar.temp) # temperature data of Antarctica 
summary(temp.pca) 
# values based on temperature difference PC1 = 90%; PC2 = 10%

temp.rda <- rda(antar.temp, scale = FALSE, center = TRUE)
temp.rda
summary(temp.rda)
names(temp.rda)
temp.rda$tot.chi# Total inertia or the sum of all eigenvalues.
temp.rda$CA

# calculate axis-importance and draw the barplots:
ev.antar.temp <- temp.rda$CA$eig
# source("/Users/Gastorga/Google Drive/Thesis_GA/MacrosLD/FunctionsLD/evplot.R")
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=1')

dev.new()
# pdf("ev.pdf")

evplot(ev.antar.temp)
# dev.off()
# only PC1 is important 

sites.temp.antar <- scores(temp.rda, display = "sites", choices = 1)
# write.csv(sites.temp.antar, file = "sites.temp.antar.csv")
# temp.comp <- read.csv("sites.temp.antar.csv")
head(sites.temp.antar)
temp.antar <- sites.temp.antar
head(temp.antar)#only PC1

source("/Users/Gastorga/Google Drive/Lake Dobson/Dead-pretty-plants/FunctionsLD/cleanplot.pca.R")
dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(temp.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(temp.rda, scaling = 2, mar.percent = 0.04)
# temp.comp  contains only PC1 accounting for 90% of the total variation in the temperature data
# although both variables are within the equilibrium circle
# I'll use PC1 as temperature in the complete matrix

# PCA based on a covariance matrix because variables are in the same scale 
# not sure if I should center
# Arguments scale=TRUE and centre = TRUE calls for a standardization of the variables
# antar.temp include the Vostok and DomeC variacion de deuterium data

# antartic and NZW.O18 all together in PCA like in the Thesis

southern.temp <- cbind(antar.temp, NZW.interp)
head(southern.temp)
southern.temp <- subset(southern.temp, select = -c(4))
head(southern.temp)
south.temp.pca <- prcomp(southern.temp, scale. = T, center = T) 
south.temp.pca
summary(south.temp.pca) 
# PC1 = 74.9%; PC2 = 18.4% based on temperature data

south.temp.rda <- rda(southern.temp, scale = TRUE, center = T)
south.temp.rda
summary(south.temp.rda)
# percentage of variance explained by each  component
(temp.eig <- south.temp.rda$CA$eig/south.temp.rda$tot.chi)

dev.new()
screeplot(south.temp.rda, bstick = TRUE)

ev <- south.temp.pca$sdev^2
# pdf("ev.pdf")

dev.new()
evplot(ev)

# only PC1 is significant
sites.south.temp.rda <- scores(south.temp.rda, choices = c(1), display = "sites")
# write.csv(sites.south.temp.rda, file = "sites.south.temp.rda.csv")
# sh.temp <- read.csv("sites.south.temp.rda.csv")
sh.temp <-sites.south.temp.rda
head(sh.temp)

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp.rda, scaling = 2, mar.percent = 0.04)

# same as above but with data from Port Arthur,
# I dont think this is correct, beciasue O18 from Port Arthur
# does not represent temerature, but instead meteoric water 
southern.temp1 <- cbind(antar.temp, NZA.interp)
head(southern.temp1)
southern.temp1 <- subset(southern.temp1, select = -c(4))
head(southern.temp1)
south.temp1.pca <- prcomp(southern.temp1, scale. = T, center = T) 
south.temp1.pca
south.temp1.rda <- rda(southern.temp1, scale = TRUE, center = T)
south.temp1.rda
summary(south.temp1.rda)

# percentage of variance explained by each  component
(temp.eig <- south.temp1.rda$CA$eig/south.temp1.rda$tot.chi)
summary(south.temp1.pca) 
# PC1 = 78%; PC2 = 15% based on temperature difference from antartic records

dev.new()
screeplot(south.temp1.rda, bstick = TRUE)

# select the data frame with eigenvalues of particular axes:
ev <- south.temp1.pca$sdev^2

# pdf("ev.pdf")
dev.new()
evplot(ev)
# only PC1 is significant

sites.south.temp1.rda <- scores(south.temp1.rda, display = "sites", choices = 1)
# write.csv(sites.south.temp1.rda, file="sites.south.temp1.rda.csv")
# sh.temp1 <- read.csv("sites.south.temp1.rda.csv")
sh.temp1 <- sites.south.temp1.rda
head(sh.temp1)
# temp.comp  contains only PC1

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp1.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp1.rda, scaling = 2, mar.percent = 0.04)

#  it is weird how NZW and NZA points towards the oldest sites from Lake Dobson
#  could it be because the 2 records represent precipitation values instead of temperature?

#--------------------------------O18 from nza and nzw------------------------------------- 
NZ.isotopes <- cbind(NZW.interp, NZA.interp)
head(NZ.isotopes)
O18.nz <- subset(NZ.isotopes, select = -c(2,4))
head(O18.nz)

O18.pca <- prcomp(O18.nz) 
O18.pca 
summary(O18.pca) # PC1 = 77%; PC2 = 14%
# percentage of variance explained by each  component

O18.rda <- rda(O18.nz, scale = F, centre = F)
O18.rda
summary(O18.rda)
(O18.eig <- O18.rda$CA$eig/O18.rda$tot.chi)
# PC1 = 0.7386376; PC2 = 0.2613624 

dev.new()
screeplot(O18.rda, bstick = TRUE)

# pdf("ev.pdf")
dev.new()
ev.O18 <- O18.pca$sdev^2
evplot(ev.O18)
# PC1 is signidicant

sites.O18.rda <- scores(O18.rda, display = "sites", choices = 1)
# write.csv(sites.south.temp1.rda, file="sites.south.temp1.rda.csv")
# sh.temp1 <- read.csv("sites.south.temp1.rda.csv")
O18.sh <- sites.O18.rda
head(O18.sh)
# temp.comp  contains only PC1

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(O18.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(O18.rda, scaling = 2, mar.percent = 0.04)


C13.nz <-subset(NZ.isotopes, select = -c(1,3))
C13.pca <- prcomp(C13.nz) 
C13.pca 
summary(C13.pca) 
# PC1 = 53.9%; PC2 = 46.01%
# percentage of variance explained by each  component

C13.rda <- rda(C13.nz)
C13.rda
summary(C13.rda)
(C13.eig <- C13.rda$CA$eig/C13.rda$tot.chi)


dev.new()
screeplot(C13.rda, bstick = TRUE)

# pdf("ev.pdf")
dev.new()
ev.C13 <- C13.pca$sdev^2
evplot(ev.C13)
# PC1 
sites.C13.rda <- scores(C13.rda, display = "sites", choices = 1)
# write.csv(sites.south.temp1.rda, file="sites.south.temp1.rda.csv")
# sh.temp1 <- read.csv("sites.south.temp1.rda.csv")
C13.sh <- sites.C13.rda
head(C13.sh)
# temp.comp  contains only PC1

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(C13.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(C13.rda, scaling = 2, mar.percent = 0.04)


southern.climate <- cbind(antar.temp, NZW.interp, co2.interp)
head(southern.climate)
southern.climate <- subset(southern.climate, select = -c(4,6))
head(southern.climate)
south.climate.pca <- prcomp(southern.climate, scale. = T, center = T) 
south.climate.pca
# percentage of variance explained by each  component
south.climate.rda <- rda(southern.climate, scale = TRUE, centre = T)
south.climate.rda
summary(south.climate.pca) 
summary(south.climate.rda)
(temp.eig <- south.climate.rda$CA$eig/south.climate.rda$tot.chi)

# PC1 = 77%; PC2 = 14% antar.temp, NZW.interp and CO2
# PC1 = 63.48%; PC2 = 20.044% antar.temp, NZW.interp and CO2

dev.new()
screeplot(south.climate.rda, bstick = TRUE)

# pdf("ev.pdf")
dev.new()
ev <- south.climate.pca$sdev^2
evplot(ev)
# only PC1 is significant in both cases

sites.south.climate.rda <- scores(south.climate.rda, choices = c(1), display = "sites")
# write.csv(sites.south.climate.rda, file="sites.south.climate.rda.csv")
# sh.climate <- read.csv("sites.south.climate.rda.csv")
sh.climate <- sites.south.climate.rda
head(sh.climate)

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate.rda, scaling = 2, mar.percent = 0.04)

southern.climate2 <- cbind(antar.temp, NZW.interp, NZA.interp, co2.interp)
head(southern.climate2)
southern.climate2 <- subset(southern.climate2, select = -c(4,6))
head(southern.climate2)
south.climate2.pca <- prcomp(southern.climate2, scale. = T, center = T) 
south.climate2.pca
summary(south.climate2.pca)

# percentage of variance explained by each  component
south.climate2.rda <- rda(southern.climate2, scale = TRUE, center = TRUE)
south.climate2.rda
summary(south.climate2.rda)

(temp.eig <- south.climate2.rda$CA$eig/south.climate2.rda$tot.chi)
# PC1 = 61.85%; PC2 = 18.75% with NZW both isotopes and O18.nza. First two axes are significant
# PC1 = 73.18%; PC2 = 11.49% with O18 and C13 from nzw and nza O18. Only PC1 is significant

dev.new()
screeplot(south.climate2.rda, bstick = TRUE)

# select the data frame with eigenvalues of particular axes:
ev <- south.climate2.rda$CA$eig

# pdf("ev.pdf")
dev.new()
evplot(ev)
# only PC1 is significant when running the analysis with O18 data from nzw and nza

sites.south.climate2.rda <- scores(south.climate2.rda, choices = c(1:2), display = "sites")
# write.csv(sites.south.climate2.rda, file = "sites.south.climate2.rda.csv")
# sh.climate2 <- read.csv("sites.south.climate2.rda.csv")
sh.climate2 <- sites.south.climate2.rda
head(sh.climate2)

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate2.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate2.rda, scaling = 2, mar.percent = 0.04)

# PCA of all env.variables keeping antartic and NZW both isotopes separated and not including
# loi.interp, dob.interp
head(temp.comp)
PC1.ant = temp.comp$PC1
O18.nzw = NZW.interp$O18.nzw
clim.complete1 <- cbind(PC1.ant, O18.nzw, co2.interp, loi.interp, char.interp, et.interp, chiro.interp)
head(clim.complete1)
clim.pca <- prcomp(clim.complete1, scale. = TRUE, center = T)
clim.pca
summary(clim.pca)
# PCA1 62%, PC2 15%, PC3 10%

clim.rda <- rda(clim.complete1, scale = TRUE, center = TRUE)
clim.rda
summary(clim.rda)

(loadings <- clim.pca$rotation)
is.data.frame(loadings)
as.data.frame(loadings)

# Principal components scores (pedicted data based on rotation)
# They are the rotated data (centred and scaled if requested) multiplied by the rotation matrix
scores.clim.pca <- clim.pca$x

# We can get the correlations among the PCs' scores and the initial data (correlation matrix)
(correlations <- t(loadings)*clim.pca$sdev)
# an correlation between variables
(cor(clim.complete1, use = "pairwise.complete.obs"))

dev.new()
screeplot(clim.rda, bstick = TRUE)
ev <- clim.pca$CA$eig
# calculate axis-importance and draw the barplots:

ev <- clim.pca$sdev^2
# pdf("ev.pdf")
dev.new()
evplot(ev)
# kaisser-Gunttman criterion reveals first PCs as important
# however broken stick only signals PC1
# Plot with cleanplot

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(clim.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(clim.rda, scaling = 2, mar.percent = 0.04)
# scaling 1 = distancia entre los sitios o muestras (scaling sites)
# scaling 1, incluye circulo de contribucion de equilibrio 
# variables mas largas que el circulo son importantes O18.nza, temp.vos, sst, temp.et
# However we can clearly see that 
# scaling 2 = correlacion entre variables representada por el angulo de los vectores (scaling environmental variables)
# I have the impression that I have multicollinearity between co2 and loi
# also perhaps collinearity (correlation) between antartic records  of temperature

# clim.complete2 include all the old variables with sh.temp (antartic and O18.nzw isotopes combined)
# non-including loi.interp and dob.interp
head(sh.temp)
PC1.sh = sh.temp$PC1
head(NZA.interp)
O18.nza <- subset(NZA.interp, select = -c(2))
head(O18.nza)

clim.complete2 <- cbind(PC1.ant, O18.nzw, co2.interp, char.interp, et.interp, chiro.interp, loi.interp)
head(clim.complete2)


clim.pca2 <- prcomp(clim.complete2, scale. = TRUE, center = T)
clim.pca2
summary(clim.pca2)

clim.rda2 <- rda(clim.complete2, scale = TRUE)
clim.rda2
summary(clim.rda2)
# PC1 explains about 51% ; PC2 20%; PC3 17%
(loadings <- clim.pca2$rotation)
is.data.frame(loadings)
as.data.frame(loadings)

# Principal components scores (pedicted data based on rotation)
# They are the rotated data (centred and scaled if requested) multiplied by the rotation matrix
scores.clim.pca2 <- clim.pca2$x

# We can get the correlations among the PCs' scores and the initial data (correlation matrix)
(correlations <- t(loadings)*clim.pca2$sdev)


dev.new()
screeplot(clim.rda2, bstick = TRUE)
ev <- clim.pca2$CA$eig
# kaisser-Gunttman & broken stick only signals PC1 
ev <- clim.pca2$sdev^2
# pdf("ev.pdf")
dev.new()
evplot(ev)

# Plot with cleanplot
source("cleanplot.pca.R")
dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(clim.rda2, scaling = 1, mar.percent = 0.08)
cleanplot.pca(clim.rda2, scaling = 2, mar.percent = 0.04)




# Scaling 1 distance biplot shows which variables contribute the most to the ordination in a few dimensions 
# It shows that sites on the right hand side of the biplot, sites 64-80 have the highest values of precipitation
# (O18.nza, dob.dca) and highest TWARM values.


#----------------------------------------individual cca--------------------------------------------#
# Individual CCA with each variable log1p(macros); 80 columns 52 variables
macros5 <- read.csv("combined_macros_shorted(5).csv")
str(macros5)
macros.log <- log1p(macros5)# computes log(1+x)

# dataset containing 47 species and 80 samples, alpine conifers are clumped together
macros.short <- read.csv("combined_macros_shorted.csv")
str(macros.short)
macros.short.log <- log1p(macros.short)

head(clim.complete)
# clim.complete <- subset(clim.complete, select = -c(1:2))
clim.complete <- cbind(clim.complete, PC1.ant, PC1.sh)
clim.complete <- as.data.frame(clim.complete)
head(clim.complete)
dome.mod <- cca (macros.log ~ temp.dome1, data = clim.complete)
dome.mod
summary(dome.mod)
# dome explains about 10% (9.9869%) of the total variability
#  and about 9.1% whwn using short macros 47 species
control <- how(within = Within(type = "series", mirror = TRUE))
(check2 <- check(macros.log, control))
summary(check2)
(anova(dome.mod,  permutations = control))
# 
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.26597 8.5403 0.00625 **
# Residual 78   2.42921                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
head(clim.complete)
vostok.mod <- cca (macros.log ~ temp.vos1, data = clim.complete)
vostok.mod
summary(vostok.mod)
# vostok temperature record explains about 11% (10.96%) of the total variability
# and about 9.72% when using short macros 47 species
(anova(vostok.mod,  permutations = control))
# model is significant at 0.01** (p-value = 0.00625) 
# although it explains fewer variability than dome the ststistical power is higher

# temp.comp = PC1comp.1 of combined temperature data from antarctica
antar.mod <- cca(macros.log ~ PC1.ant, data = clim.complete)  
antar.mod
summary(antar.mod)
# PC1 explains about ~12% (11.50%) of the constrained variability
# and 10.5% whwn using short macros 47 species
(anova(antar.mod,  permutations = control))
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.30994 10.135 0.00625 **
# Residual 78   2.38524                      
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete)  
sh.mod <- cca(macros.log ~ PC1.sh, data = clim.complete)  
sh.mod
summary(sh.mod)
# PC1 explains about ~12.35% of the constrained variability
# and 11.44% when using short macros 47 species
(anova(sh.mod,  permutations = control))

# Df ChiSquare      F  Pr(>F)   
# Model     1   0.33298 10.995 0.00625 **
# Residual 78   2.36220                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

co2.mod <- cca(macros.log ~ co21, data = clim.complete)  
co2.mod
summary(co2.mod)

# co2 explains about ~13% (12.51%) of the total variability
# and 11.85% when using short macros 47 species
(anova(co2.mod,  permutation = control))
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.33711 11.151  0.001 ***
# Residual 78   2.35807                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dca.mod <- cca(macros.log ~ dca1.dob, data = clim.complete)  
dca.mod
summary(dca.mod)
# Chironomids DCA1analysis explains about ~13% (12.68) of the total variability
#  and 11.55% when using short macros 47 species
(anova(dca.mod,  permutations = control))
#          Df ChiSquare      F Pr(>F)    
# Model     1    0.3417 11.325 0.00625 **
# Residual 78    2.3535     

loi.mod <- cca(macros.log ~ loi1, data = clim.complete)  
loi.mod
summary(loi.mod)
# loi explains about 9% (8.889) of the total variability
#  and 7.96% when using short macros 47 species
(anova(loi.mod,  permutations = control))
#          Df ChiSquare      F Pr(>F)    
# Model     1   0.23961 7.6112 0.0375 *
# Residual 78   2.45557         

# Model: cca(formula = macros.short.log ~ loi1, data = clim.complete)
#           Df ChiSquare      F  Pr(>F)  
# Model     1   0.18804 6.7482 0.04375 *
# Residual 78   2.17350                 


char.mod <- cca(macros.log ~ char1, data = clim.complete)  
char.mod
summary(char.mod)
# char explains just about 3 % (2.543%) of the total variability
# and 2.416% whwn using macros short
(anova(char.mod,  permutations = control))
# model is marginally significant at the 0.1 level (p-value = 0.08125) 

et.mod <- cca(macros.log ~ temp.et1, data = clim.complete)  
et.mod
summary(et.mod)
# Eagle Tarn temperature explains about 5% (5.244%) of total variability
#  and 4.529% when using macros short
(anova(et.mod,  permutations = control))
# model is significant at 0.01** level (p-value = 0.00625) 

O18.mod <- cca(macros.log ~ O18.nzw, data = clim.complete)  
O18.mod
summary(O18.mod)
# NZ.O18 explains about 8% (8.14%) of the total variability
#  and 7.702% when using macros short 47 species
(anova(O18.mod,  permutations =  control))
# model is highly significant 

#           Df ChiSquare      F Pr(>F)    
# Model     1    0.2194 6.9124 0.00625 **
# Residual 78    2.4758      
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

C13.mod <- cca(macros.log ~ C13.nzw, data = clim.complete)  
C13.mod
summary(C13.mod)
# NZ.C13 explains about 2 % (1.848%) of the total variability
#  and 1.696% and non significant when using macros short
(anova(C13.mod,  permutations =  control))
# marginally significant 

#           Df   ChiSquare      F Pr(>F)  
# Model     1    0.0498 1.4683  0.075 .
# Residual 78    2.6454                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


O18nza.mod <- cca(macros.log ~ O18.nza, data = clim.complete)  
O18nza.mod 
summary(O18nza.mod )
RsquareAdj(O18nza.mod)
# NZA.O18 explains about 9% (9.087%)of the total variability
(anova(O18nza.mod,  permutations = control))
# model is significant at 0.05 level(*) p-value = 0.0375

C13nza.mod <- cca(macros.log ~ C13.nza, data = clim.complete)  
C13nza.mod
summary(C13nza.mod)
# NZA.13C explains about 2.3% (2.298%)of the total variability
(anova(C13nza.mod,  permutations =  control))
# non-significant model

set.seed(60)
head(clim.complete)
cond.model <- cca(macros.log ~ dob.age, data = clim.complete)
cond.model
summary(cond.model)
(anova(cond.model,  permutations = control))
# age explains about 12% (12.28%) of the total variability 
# permutation test showed the sediment age was significantly 
# correlated with the plant macrofossil data at 0.01 level (p-value = 0.00625)
# therefore sediment age is used as covariable in subsequent ordination


# clim.complete include all variables with composite antartic temperature comp.1 
# and both NZW isotopes separated
set.seed(100)
head(clim.complete)
clim.complete <- as.data.frame(clim.complete)
cond.model.ant <- cca(macros.log ~ co21 + PC1.ant + O18.nza + temp.et1 + char1 + loi1 + dca1.dob + Condition(dob.age), data = clim.complete)
cond.model.ant
summary(cond.model.ant)
#                Inertia Proportion Rank
# Total          2.6952     1.0000     
# Conditional    0.3311     0.1228    1
# Constrained    0.3357     0.1245    6
# Unconstrained  2.0284     0.7526   51
# Inertia is scaled Chi-square 

# Accumulated constrained eigenvalues
# Importance of components:
#                       CCA1    CCA2    CCA3    CCA4    CCA5    CCA6
# Eigenvalue            0.1570 0.05693 0.04549 0.03175 0.02327 0.02125
# Proportion Explained  0.4677 0.16960 0.13552 0.09460 0.06932 0.06330
# Cumulative Proportion 0.4677 0.63727 0.77279 0.86738 0.93670 1.00000
# 
# Scaling 2 for species and site scores
# * Species are scaled proportional to eigenvalues
# * Sites are unscaled: weighted dispersion equal on all dimensions
control <- how(within = Within(type = "series", mirror = TRUE))
set.seed(60)
(anova(cond.model.ant, model = "full", permutations = control)) # Overall model is significant
(anova(cond.model.ant, by = "term", model = "full", permutations = control))
(anova(cond.model.ant, by = "axis", model = "full", permutations = control)) # only CCA1 is significant
(anova(cond.model.ant, by="margin", model = "full", permutations = control))# the effect of a particualr term when all other model terms are included
control <- how(within = Within(type = "series", mirror = TRUE))

set.seed(60)
head(clim.complete)
cond.model.sh <- cca(macros.log ~ co21 + PC1.sh + temp.et1 + char1 + C13.nzw + dca1.dob + Condition(dob.age), data = clim.complete)
cond.model.sh
summary(cond.model.sh)
#                Inertia Proportion Rank
# Total          2.6952     1.0000     
# Conditional    0.3311     0.1228    1
# Constrained    0.2988     0.1109    5
# Unconstrained  2.0653     0.7663   51
# Inertia is scaled Chi-square 

# Accumulated constrained eigenvalues
# Importance of components:
#                       CCA1    CCA2    CCA3    CCA4    CCA5    CCA6
# Eigenvalue            0.1620 0.04915 0.04042 0.02470 0.02251
# Proportion Explained  0.5423 0.16449 0.13526 0.08266 0.07532
# Cumulative Proportion 0.5423 0.70676 0.84201 0.92468 1.00000

# environmental constrains explain about 11% of the total variability and age about 12%
set.seed(60)
(anova(cond.model.sh, model = "full", permutations = control)) # Overall model is significant
(anova(cond.model.sh, by = "term", model = "full", permutations = control))
(anova(cond.model.sh, by = "axis", model = "full", permutations = control)) # only CCA1 is significant
(anova(cond.model.sh, by="margin", model = "full", permutations = control))

# inertcomp(cond.model.sh, prop = TRUE) # I dont really know horw to interpret this
# intersetcor(cond.model.sh)

# Manual model building
clim.complete0
head(clim.complete1)
head(clim.complete2)
head(clim.complete)

# clim.complete3 includes char1, temp.et1, dca1.dob, loi1, co21, C13.nzw, PC1.sh
clim.complete3 <- subset(clim.complete, select = -c(1:2, 9:12))
head(clim.complete3)
# clim.complete4 includes PC1.ant, CO2, O18.nzw, C13.nzw, cha1, temp.et, loi, dca1.dob
clim.complete4 <- subset(clim.complete, select = -c(1:2,6, 9:12, 14))
head(clim.complete4)
clim.complete4 <- as.data.frame(clim.complete4)
# clim.complete0 corresponds to the dataset used in my thesis
# the final mod included co2, PC1.sh and char explaining about 17.7% of variability instead of 16%
# and the two firs axis account for 92% of the constrained variability
# instead 89% declared in the thesis. This could be because in the thesis I used a different 
# data matrix 47 species, 80 samples. also because antartic records are in the AICC2012 chronology.

# clim.complete3 includes char1, temp.et1, dca1.dob, loi1, co21, C13.nzw, PC1.sh
# the final model includes variables dca1.dob + co21 + PC1.sh explaining about 18% of the total variability
# with CCA1 = 74% and CCA2 = 20%.
# CCA-based forward selection with 9999 permutations selected the following model:
# macros.log ~ dca1.dob + co21 + PC1.sh + temp.et1 
# parsimonious cca produce an overall significative model (0.00625 **) explaining 19.75% of the total variability
# with CCA1 = 68% and CCA2 = 20% of the constrained variability. CCA1 is significant
# Adj.r.squared = 0.1550955
# VIF = dca1.dob 6.497104, co21 4.947425, PC1.sh 5.371197 temp.et1 1.495513 

clim.complete1 <- as.data.frame(clim.complete1)
head(clim.complete1)

mbig <- cca(macros.log ~  ., clim.complete1)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete1)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope=formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))
m0 <- update(m0, . ~ . + co21)
add1(m0, scope=formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))
m0 <- update(m0, . ~ . + PC1.sh)
add1(m0, scope=formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))
m0 <- update(m0, . ~ . + char1)
add1(m0, scope=formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))
m0 <- update(m0, . ~ . + loi1)
add1(m0, scope=formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))

## -- included variables still significant?

final.mod1 <- cca(macros.log ~  dca1.dob +  co21 + PC1.sh , clim.complete1)
final.mod1
summary(final.mod1)

vif.cca(final.mod1)
# dca1.dob     co21   PC1.sh 
# 5.846661 4.907481 5.260713 
drop1(final.mod1, test = c("permutation"), 
      permutations = how(within = Within(type = "series", mirror = F)))
add1(m0, scope = formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))

spe.cca <- cca(macros.log ~., clim.complete1)
summary(spe.cca)
RsquareAdj(spe.cca)

set.seed(100)
cca.step.forward <- ordistep(cca(macros.log ~ 1, data = clim.complete1),
                             scope = formula(spe.cca),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  dca1.dob + co21 + PC1.sh + temp.et1 
spe.cca.pars <- cca(macros.log ~ dca1.dob +  co21 + PC1.sh + temp.et1 + char1, data = clim.complete1)
summary(spe.cca.pars)
anova(spe.cca.pars, permutations = control)
anova(spe.cca.pars, permutations = control, by = "axis")
RsquareAdj(spe.cca.pars)
vif.cca(spe.cca.pars)



set.seed(200)
head(clim.complete3)
mod0 <- cca(macros.log ~ 1, data = clim.complete3)  # Model with intercept only
summary(mod0)
mod1 <- cca(macros.log ~ ., data = clim.complete3)  # Model with all explanatory variables
summary(mod1)
RsquareAdj(mod1) 
# $r.squared
# [1] 0.2556905
# 
# $adj.r.squared
# [1] 0.1833542
# Scaling 1 : species scores scaled to the relative eigenvalues,
# sites are weighted averages of the species
dev.new()
par(mfrow = c(1,2))
plot(spe.cca.pars, scaling = 1, display = c("sp", "lc", "cn"),
     main ="Triplot CCA macros.log ~ clim.complete3 - Scaling 1")
# Default scaling 2: site scores sclaed to the relative eigenvalues,
# species are weigjted averages of the sites
plot(spe.cca.pars,
     display = c("lc", "cn"),
     main = "Triplot CCA macros.log ~ clim.complete3 - Scaling 2")



model.fit1 <- ordistep(mod0, scope = formula(mod1), direction = "forward")
summary(model.fit1)
vif.cca(model.fit1)

(anova(model.fit1,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0125
(anova(model.fit1, by = "term", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(model.fit1, by = "axis", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(model.fit1, by = "margin", permutations = control))
vif.cca(model.fit1)

cca.final.forward <- ordistep(cca(macros.log ~ 1, data = clim.complete3), scope = formula(mod1),
                      direction = "forward", permutations = how(nperm = 999))

summary(cca.final.forward)


## Automatic model building based on AIC but with permutation tests
step(cca(dune ~  1, dune.env), reformulate(names(dune.env)), test="perm")
## see ?ordistep to do the same, but based on permutation P-values
## Not run:
(ordistep(cca(dune ~  1, dune.env), reformulate(names(dune.env))))



# Manual model building 
# maximal model



# clim.complete explains about 24% of the total variance
# while cca1 and cca2 explain 57% and 21% of the total constrained variability
vif.cca(comp.model)
# vif values are all below 10
# PC1      char1    temp.et1 loi1     co21     O18.nza 
# 3.019563 1.211625 1.529971 1.910847 5.325157 2.520389 

control <- how(within = Within(type = "series", mirror = FALSE))
(check2 <- check(macros.log, control))
summary(check2)

(anova(comp.model,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0125
(anova(comp.model, by="term", permutations = control))
# significant env.var under restricted permutation PC1 composite temperature significant at 0.05 level (p-value = 0.0125) 
# & char marginally significant(0.0750) 
(anova(comp.model, by="axis", permutations = control))
# CCA1 & CCA2 are significant at 0.05 level, p-value = 0.0125 in both cases

# Running the model with only significant variables
set.seed(60)
final.mod <- cca(macros.log ~ Comp.1 + co21 + char1, data = clim.complete2)
final.mod
summary(final.mod)
# final.mod explains about 17% of the total variability and CCA1 explains about 76% of the constrained variability
#  CCA2 explains about % 17
# R2adjuted is around 14%
RsquareAdj(final.mod)

(anova(final.mod,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0125
(anova(final.mod, by = "term", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(final.mod, by = "axis", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(final.mod, by = "margin", permutations = control))
vif.cca(final.mod)

ource("/Users/giselleastorga/Documents/Thesis_GA/second_chapter/data_second/cleanplot.pca.R")
dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(temp.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(temp.rda, scaling = 2, mar.percent = 0.04)


set.seed(200)
head(clim.complete3)
mod0 <- cca(macros.log ~ 1, data = clim.complete3)  # Model with intercept only
summary(mod0)
mod1 <- cca(macros.log ~ ., data = clim.complete3)  # Model with all explanatory variables
summary(mod1)
model.fit1 <- ordistep(mod0, scope=formula(mod1),direction = "forward")
summary(model.fit1)
vif.cca(model.fit1)

(anova(model.fit1,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0125
(anova(model.fit1, by = "term", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(model.fit1, by = "axis", permutations = control))
# significant env.var under restricted permutation antartic temp (0.0125), char (0.05) and CO2 (0.025)
(anova(model.fit1, by = "margin", permutations = control))
vif.cca(model.fit1)

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
par(mfrow=c(1, 2))
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome1', scale = "density", binwidth = .7)
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos1', scale = "density", binwidth = .7)
grid.arrange(dome, vos, nrow = 1) 
# Change the width of bars
# Change y axis values to density
dev.new()
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos1', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                         meanLineType = "dashed", meanLineSize = 1)
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome1', scale = "density", binwidth = .7, addMeanLine = TRUE, meanLineColor = "red",
                          meanLineType = "dashed", meanLineSize = 1)
grid.arrange(dome, vos, nrow = 1) 

# Add density curve
vos2 <- ggplot2.histogram(data = antar.climate1, xName = 'temp.vos1',
                  fill = "white", color = "black",
                  addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)

dome2 <- ggplot2.histogram(data = antar.climate1, xName = 'temp.dome1',
                          fill = "white", color = "black",
                          addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = .4)
dev.new()
grid.arrange(dome2, vos2, nrow = 1) 

NZ.O18 <- ggplot2.histogram(data = antar.climate1, xName = 'O18',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = 0.02)
NZ.C13 <- ggplot2.histogram(data = antar.climate1, xName = 'C13',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', binwidth = 0.06)

dev.new()
grid.arrange(NZ.O18, NZ.C13, nrow = 1, ncol = 2) 
grid.arrange(vos2, dome2, NZ.O18, NZ.C13, nrow = 4, ncol = 1)

NZW.O18 <- ggplot2.histogram(data = antar.climate5, xName = 'O18',
                           fill = "white", color = "black",
                           addDensityCurve = TRUE, densityFill = '#FF6666', bindwidth = 0.05)
NZW.C13 <- ggplot2.histogram(data = antar.climate5, xName = 'C13',
                            fill = "white", color = "black",
                            addDensityCurve = TRUE, densityFill = '#FF6666', bindwidth = 0.05)
dev.new()
grid.arrange(NZW.O18, NZW.C13,vos2,dome2, nrow = 4, ncol = 1) 

skewness(antartic$temp.vos1)# kutosis = -0.9 and kurtosis = 3.4
kurtosis(antartic$temp.vos1)
head(antartic)


##****************CCA**********************************##
# CCA with PCAs of temp.pca and log1p(macros); 80 columns 52 variables
macros5 <- read.csv("combined_macros_shorted(5).csv")
str(macros5)
head(temp.pca)
macros.log <- log1p(macros5)

head(clim.complete2)
clim.complete3 <- subset(clim.complete2, select = -c(7))
head(clim.complete3)
set.seed(60)
mod0 <- cca(macros.log ~ 1, data = clim.complete3)  # Model with intercept only
summary(mod0)
mod1 <- cca(macros.log ~ ., data = clim.complete3)  # Model with all explanatory variables
summary(mod1)
model.fit <- ordistep(mod0, scope = formula(mod1), direction = "forward")
summary(model.fit)
# Automatic variable selection indicates both components are important
# PC1 and PC2 Antartic variation in deuterium explains about 11% of the total variability with
# CCA 1 explaining about 83% of the constrained varibility
# CCA 2 explains about 17% of the constrained variability
# Free permutation of our cca results is highly significant, but data is temporaly ordered
# so that restricted permutations are needed (e.g. dataset with 80 stratigraphic ordered samples,
# has just 80 posible permutations. I generated the complete set of permutations by suffleSet() instead of 
# a random set. This is also posible because dataset is relatively small.
(anova(model.fit))
control <- how(within = Within(type = "series", mirror = FALSE))
(check2 <- check(macros.log, control))
summary(check2)
(anova(model.fit,  permutations = control))
# overal model is significant at  0.05 level (*), p-value = 0.0125
(anova(model.fit, by="term", permutations = control))
# only PC1 is significant at 0.05 (p-value = 0.0125), char (0.05) and CO2 (0.025)
(anova(model.fit, by="axis", permutations = control))
# CCA1 is significant at 0.05 level, p-value = 0.0125 

# antar.climate3 include Antartic tmeperature data + O.18 from Mt. Arthur
head(antar.climate3)
climate3.pca <- prcomp(antar.climate3, scale. = T, center = T)
climate3.pca 
summary(climate3.pca)

climate3.rda <- rda(antar.climate3, scale = T, center = T)
climate3.rda 
summary(climate3.rda)
# PC1 65%; PC2 26% 
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
# only PC1 is significant
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
climate3.a <-subset(climate3, select = -c(1))# leaving only PC1-PC3
head(climate3.a)
climate3.b <-subset(climate3, select = -c(1,3,4))# leaving only PC1-PC2
head(climate3.b)

set.seed(600)
sh.temp.mod <- cca(macros.log ~ PC1, data = climate3.b)  # Model with intercept only
summary(sh.temp.mod)

# sh.temp explained about 11% of the total variability 
# using as input matrix only PC1

(anova(sh.temp.mod ,  permutations = control))
# only PC1 is important 
(anova(sh.temp.mod , by="term", permutations = control))
# only cca1 is important. 
(anova(sh.temp.mod , by="axis", permutations = control))
# cca1 & cca2 are significant after restricted permutations


# Last two lines show the comparison between real variation 
# represented by individual PCA axes, and relevant variation 
# calculated by broken-stick model (1 in the column with particular 
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
(CV <- (100*sdrow/meanrow))# CV = 91.6 meaning that relativisation would have a moderate effect on results

(meancol <- mean(coltotals))
(sdcol <- sd(coltotals))
(CV <- (100*sdcol/meancol))# CV = 197

macros.log <- log1p(macros5)


DCA.macros <- decorana(macros.log, iweigh = 1, ira = 0) 
summary(DCA.macros)
DCA.macros$fraction #downweighting abundance fraction starts at 5 
cca(macros.log)
dca.samplescores <- scores(DCA.macros,display = c("sites"), choices = 1)
write.csv(dca.samplescores, file = "dca.samplescores.csv")

# extracts only axis 1 scores for samples

#macros<-read.csv("combined_macros_shorted.csv")##macros shorted (47 columsn) at species higher than 5% 
#str(macros)
#macros.slog<-log1p(macros)
#DCA.macros1 <- decorana(macros.slog) 
#summary(DCA.macros1)
#cca(macros.slog)

dev.new()
pdf("DCA.sites.pdf")
ordiplot (DCA.macros, display = 'si', type = 'n')
points (DCA.macros, col = zones$group, pch = zones$group )
for (i in seq (1, 5)) ordihull (DCA.macros, groups = zones$group, show.groups = i, col = i, lty = 'dotted')
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
for (i in seq (1, 5)) ordispider (DCA.macros, groups = zones$group, show.groups = i, col = i, label = T)
for (i in seq (1, 5)) ordihull (DCA.macros, groups = zones$group, show.groups = i, col = i, lty = 'dotted')
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

plot(DCA.macros, display=c("none"), cols=c(1,2))
# plot axis 1 & 2, but display no points or labels
text(DCA.macros, display="sites", pcol = "red", cex=0.7)
# plot samples as green crosses for axis 1 and 2
text(DCA.macros, display=c("species"), choices=1:2,
       cex=0.7)
# plot labels for taxa for axis 1 and 2, using cex to shrink
# size of labels. Larger plots may also be used to alleviate
# congestion of labels.
dev.off()



##to abbreviate Latin names
shnam <- make.cepnames(names(macros.log))
DCA.macros<-decorana(env1) 
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
pdf("Figure_correlplot.pdf",  width = 20, height = 17, units = 'cm', res = 100)
dev.new()
ordiplot(model_fit1, scaling=2,main="Correlation",type="n")
segments(x0=0,y0=0,x1=scores(model_fit1, display="species", scaling=2)[,1],
y1=scores(model_fit1, display="species", scaling=2)[,2])
text(model_fit1, display="sp", scaling=2, col=2, pch=21, cex=.7)
text(model_fit1, display="bp", scaling=2,
row.names(scores(model_fit1, display="bp")), col=4)
text(model_fit1, display=c("sites"),pch=15,cex=.5, scaling=2,labels=rownames(macros_log5))
cor(macros.log,env1)
dev.off()

# Correlation , scaling=2
plot(model_fit1, scaling=2,main="Correlation",type="n")
segments(x0=0,y0=0,x1=scores(model_fit1, display="species", scaling=2)[,1],
         y1=scores(model_fit1, display="species", scaling=2)[,2])
text(model_fit1, display="sp", scaling=2, col=2)
text(model_fit1, display="bp", scaling=2,row.names(scores(model_fit1, display="bp")), col=4)
text(model_fit1, display=c("lc"), scaling=2,labels=rownames(macros_log5))

# Distance , scaling=1
eps("Figure_distplot.eps",  width = 20, height = 17, units = 'cm', res = 100)
x11()
plot(model_fit1, scaling=1,main="Distance triplot",type="n")
segments(x0=0,y0=0,x1=scores(model_fit1, display="species", scaling=1)[,1],
y1=scores(model_fit1, display="species", pch=2, cex=.2, scaling=1)[,2])
text(model_fit1, display="sp", scaling=1, col=2)
text(model_fit1, display="bp", scaling=1,row.names(scores(model_fit1, display="bp")), col=4)
text(model_fit1, display=c("sites"),pch=2,cex=.8, scaling=1,labels=rownames(macros_log5))
dev.off()




scores.macros_cca_model_fit1<-scores(model_fit1,display="sites")##scores cca
scores.macros_cca_spp_model_fit1<-scores(model_fit1,display="species")##scores cca
write.csv(scores.macros_cca_spp_model_fit1, file = "scores.macros_cca_spp_model_fit1.csv")



##ordistep(macros_cca_temp_stand_fire_3~1,direction="both",step=999, max.perm=999)
ordistep(cca(macros_log5 ~  1, env.var), reformulate(names(env.var)), direction="both",pstep=1000, step=1000,max.perm=1000)
ordistep(d, scope = formula(d),direction="both",step=1000, max.perm=1000)


## End(Not run)
## Manual model building
## -- define the maximal model for scope
mbig <- cca(macros_log5 ~  ., climatics_complete)
## -- define an empty model to start with
m0 <- cca(macros_log5 ~ 1, climatics_complete)
## -- manual selection and updating
add1(m0, scope=formula(mbig), test="perm")
m0 <- update(m0, . ~ . + temperature_all)
add1(m0, scope=formula(mbig), test="perm")
m0 <- update(m0, . ~ . + co2_dome)
add1(m0, scope=formula(mbig), test="perm")
## -- included variables still significant?
drop1(m0, test="perm")
add1(m0, scope=formula(mbig), test="perm")
ordistep(cca(macros_log5 ~  1, climatics_complete), reformulate(names(climatics_complete)), perm.max=1000)



#******************************RAREFACTION************************************************
macros<-read.csv("combined_rarefaction.csv")
modern<-read.csv("modern_sediment.csv")
raremax <- min(rowSums(macros))
raremax
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
out <- with(pars[1:47, ],
            rarecurve(macros, step = 20, sample = raremax, col = col,
                      lty = lty, label = FALSE))

mrar <- rarefy(macros, min(rowSums(macros)))

data(modern)
test<-rich(matrix=modern, nrandom=499,verbose=TRUE)
test$cr # observed cumulative species richness
[1] 25
test$mr # observed mean value of species richness over the n samples
[1] 8.15


test<-rich(matrix=macros, nrandom=499,verbose=TRUE)
test$cr 
[1] 47  # observed cumulative species richness
test$mr # observed mean value of species richness over the n samples
[1]  12.94
data(macros) # culture plot
rare_macros<-raref(matrix=macros, dens=sum(macros), nrandom=500)
rare_macros$Sinterp[2]


data(BCI, package = "vegan")
BCI2 <- BCI[1:26, ]
raremax <- min(rowSums(BCI2))
raremax




#######################Cluster analysis################################
# Bray curtis cluster using macros1 (80 samples, 52 species)
macros1<- read.csv("combined_macros_shorted(5).csv")
str(macros1)
macrosper1 <- percenta(macros1,first=1,last=52)
str(macrosper1)
write.csv(macrosper1, file = "macrosper1.csv")
macros1<-read.csv("macrosper1.csv")
## Bray-Curtis distances between samples, even if method is not especify
diss1 <- vegdist(macrosper1, method="bray")
clust1<- chclust(diss1, method="coniss")
dev.new()
bstick(clust1,10)

diss1 <- vegdist(sqrt(macrosper1/100), method = "bray")
clust1 <- chclust(diss1, method="coniss")
dev.new()
bstick(clust1)# points above the red line suggest 3 zones
spe.chclust1<-chclust(vegdist(macrosper1))
dev.new()
bstick(spe.chclust1,10)
k<-3
(gr4<-cutree(spe.chclust1, k = k))
dev.new()
plot(spe.chclust1, hang =-1, main = "CONISS clustering Lake Dobson")
rect.hclust(spe.chclust1, k = k)
## Bray curtis cluster using combined_macros_shorted at 5.csv (80 samples, 47 species)
macros2<- read.csv("combined_macros_shorted.csv")
str(macros2)
macrosper2 <- percenta(macros2,first=1,last=47)
str(macrosper2)

## Bray-Curtis distances between samples, even if method is not especify
diss2 <- vegdist(macrosper2, method="bray")
clust2<- chclust(diss2, method="coniss")
dev.new()
bstick(clust2,10)
diss2 <- vegdist(sqrt(macrosper2/100), method = "bray")
clust2 <- chclust(diss2, method="coniss")
dev.new()
bstick(clust2,10)# points above the red line suggest 3 zones

## Bray curtis cluster using macros_full_30_Dec_216 (80 samples, 64 species)
macros3<- read.csv("macros_full_30_Dec_216.csv")
str(macros3)
macrosper3 <- percenta(macros3,first=1,last=64)
str(macrosper3)
write.csv(macrosper3, file="macrosper3.csv")
macrosper3<-read.csv("macrosper3.csv")

diss3 <- vegdist(macrosper3, method="bray")
clust3<- chclust(diss3, method="coniss")
dev.new()
bstick(clust3,10)

diss3 <- vegdist(sqrt(macrosper3/100), method = "bray")
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
diss4 <- vegdist(sqrt(macros4/100), method = "bray")
clust4 <- chclust(diss4, method="coniss")
dev.new()
bstick(clust4)# points above the red line suggest 3 zones

## Bray curtis cluster using macros_full_30_Dec_216 (80 samples, 62 species)
macros5<- read.csv("macros_full_17_may_2016.csv")
str(macros5)
macrosper5 <- percenta(macros5,first=1,last=62)
str(macrosper5 )
## Bray-Curtis distances between samples, even if method is not especify
diss5 <- vegdist(macrosper5, method="bray")
clust5<- chclust(diss5, method="coniss")
dev.new()
bstick(clust5)

diss5 <- vegdist(sqrt(macrosper5/100), method = "bray")
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

diss6 <- vegdist(sqrt(macros6/100), method = "bray")
clust6 <- chclust(diss6, method="coniss")
dev.new()
bstick(clust6)# points above the red line suggest 3 zones


#88 columns, 47 species + depth+age
macros7<- read.csv("macrosper.csv")
str(macros7)
depth<-macros7[48:49]
macrosper7<-macros7[1:47]
## Bray-Curtis distances between samples, even if method is not especify
diss7 <- vegdist(macrosper7, method="bray")
clust7<- chclust(diss7, method="coniss")
dev.new()
bstick(clust7,10)

diss7 <- vegdist(sqrt(macrosper7/100), method = "bray")
clust7 <- chclust(diss7, method="coniss")
dev.new()
bstick(clust7)# points above the red line suggest 3 zones

#macros_percenta 
macros8<- read.csv("macrosp.csv")
str(macros8)
depth<-macros8[63:64]
str(depth)
macrosper8<-macros8[1:62]
str(macrosper8)
## Bray-Curtis distances between samples, even if method is not especify
diss8 <- vegdist(macrosper8, method="bray")
clust8<- chclust(diss8, method="coniss")
dev.new()
bstick(clust8,10)

diss8 <- vegdist(sqrt(macrosper8/100), method = "bray")
clust8 <- chclust(diss8, method="coniss")
dev.new()
bstick(clust8)


#macros_percenta 
macros9<- read.csv("macros_cluster.csv")
str(macros9)
depth<-macros9[48:49]
str(depth)
macros9<-macros9[1:47]
str(macrosper9)
macrosper9 <- percenta(macros9,first=1,last=47)
str(macrosper9 )

## Bray-Curtis distances between samples, even if method is not especify
diss9 <- vegdist(macrosper9, method="bray")
clust9<- chclust(diss9, method="coniss")
dev.new()
bstick(clust9,10)

diss9 <- vegdist(sqrt(macrosper9/100), method = "bray")
clust9 <- chclust(diss9, method="coniss")
dev.new()
bstick(clust9)


# Basic diagram
dev.new()
plot(clust_bray, hang=-1)
# Rotated through 90 degrees
dev.new()
plot(clust_bray, hang=-1, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
dev.new()
plot(clust_bray, xvar=macros_clu$age, hang=-1, horiz=TRUE, x.rev=TRUE)
plot(clust_bray, xvar=macros_clu$depth, hang=-1, horiz=FALSE, x.rev=FALSE)
bstick(clust_bray)
plot(clust_bray,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 6)


depth<-macros_cluster[48:49]
str(depth)
macrosclu <- macros_cluster[1:47]
str(macrosclu)
macroscluper <- percenta(macrosclu,first=1,last=47)
str(macroscluper)

## Bray-Curtis distances between samples, even if method id not especify
diss1 <- vegdist(macrospershort, method="bray")
clust1<- chclust(diss1, method="coniss")
dev.new()
bstick(clust1)
dev.new()
plot(clust1,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 4)

diss2 <- vegdist(sqrt(macrospershort/100), method = "bray")
clust2 <- chclust(diss2, method="coniss")
dev.new()
bstick(clust2)# points above the red line suggest 3 zones
dev.new()
plot(clust_bray,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 3)



# Basic diagram
dev.new()
plot(clust_bray, hang=-1)
# Rotated through 90 degrees
dev.new()
plot(clust_bray, hang=-1, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
dev.new()
plot(clust_bray, xvar=macros_clu$age, hang=-1, horiz=TRUE, x.rev=TRUE)
plot(clust_bray, xvar=macros_clu$depth, hang=-1, horiz=FALSE, x.rev=FALSE)
bstick(clust_bray)
plot(clust_bray,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_bray, 6)


#SQRT cluster wirth macros short (80 samples)
macros<-read.csv("macros_percenta5.csv")
diss <- dist(sqrt(macros/100))
clust_sqrt <- chclust(diss, method="coniss")
dev.new()
bstick(clust_sqrt)
dev.new()
plot(clust_sqrt,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
plot(clust_sqrt, xvar=macros$depth, hang=-1, horiz=FALSE, x.rev=FALSE)
rect.hclust(clust_sqrt, 5)
#whole data samples with non macros (11) were removed
macros_clu<-read.csv("macros_cluster.csv")
str(macros_clu)
macros_percenta_clu<-percenta(macros_clu,first=1,last=47)
write.csv(macros_percenta_clu, file = "macros_percenta_clu.csv")
macros_percenta_clu<-read.csv("macros_percenta_clu.csv")
diss_sqrt <- dist(sqrt(macros_percenta_clu/100))
clust_sqrt <- chclust(diss_sqrt, method="coniss")
dev.new()
bstick(clust_sqrt)
plot(clust_sqrt,hang=-1,cex=.6,main="Cluster Analysis of Lake Dobson CONISS") 
rect.hclust(clust_sqrt, 5)

# Basic diagram
plot(clust, hang=-1)
# Rotated through 90 degrees
#plot(clust, hang=-1, horiz=TRUE)
# Rotated and observations plotted according to sample depth.
#plot(clust, xvar=RLGH$depths$Depth, hang=-1, horiz=TRUE, x.rev=TRUE)

plot(clust,hang=-1)
plot(clust,hang=-1, horiz=TRUE ,x.rev=TRUE)
plot(clust,hang=-1, horiz=TRUE ,x.rev=TRUE,cex=.6,xvar=macros$depth,main="Cluster Analysis of Lake Dobson")

plot(clust,hang=-1,cex=.6,xvar=macros$Depth,main="Cluster Analysis of Lake Dobson") 
rect.hclust(dud.clust2, 6)

macros<-read.csv("macros.csv")
dud.dist2 <- vegdist(macros) 
dud.clust2 <- chclust(dud.dist2) 
par(mfrow=c(1,1)) 
plot(dud.clust2,hang=-1,cex=.6) 
plot(dud.clust2,hang=-1,cex=.6)
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
#                   addMeanLine=TRUE, meanLineColor="white", meanLineSize=1.5)



# setwd("/Users/Gastorga/Desktop/chapter_one")#modern sediments
#modern sediments
modern<-read.csv("modern_sediment.csv")
modern_percenta<-percenta(modern,first=1,last=25)
write.csv(modern_percenta, file = "modern_percenta.csv")
modern_per<-read.csv("modern_percenta.csv")
summary(modern_per)
modern_log<-log(modern_per+1)
#DCA modern
#****************************************************************************************************
DCA.modern<-decorana(modern_log,iweigh=1,ira=0) #performs DCA and downweights rare taxa
summary(DCA.modern)
modern_scores_DCA<-DCA.modern$rproj
write.csv(modern_scores_DCA, file = "modern_scores_DCA.csv")
shnam <- make.cepnames(names(macros_log5))
plot(DCA.modern)
names(modern_log)
##to abbreviate Latin names
shnam <- make.cepnames(names(modern_log))
###priority to the most abundant species
stems <- colSums(modern_log)
plot(DCA.modern, dis="sites")
sel <- orditorp(DCA.modern, dis="sites", lab=shnam,  pcol = "red", pch="+")

#CAmodern
#****************************************************************************************************
CA.modern<-decorana(modern_log,ira=1) #performs CA (actually reciprocal averaging)
summary(CA.modern)
modern_scores_CA<-CA.modern$rproj##sample scores
modern_scores_CA
write.csv(modern_scores_CA, file = "modern_scores_CA.csv")
plot(CA.modern, dis="sp", type="n")
sel <- orditorp(CA.modern, dis="species",   pcol = "red", pch="+")
plot(CA.modern,type="n",main="Correspondence Analysis") #creates empty plot(type=”n”)
text(CA.modern,display="species",col="red",cex=.6) #adds red labels for species
points(CA.modern,display="sites",pch=16)

#surveys
lakeside<-read.csv("lakeside.csv")
summary(lakeside)
lakeside_percenta<-percenta(lakeside,first=1,last=61)
write.csv(lakeside_percenta, file = "lakeside_percenta.csv")

#****************************************************************************************************

mod<-read.csv("modern.csv")
modern_stand<- decostand(mod2, "standardize")
write.csv(env_stand,file="env_stand.csv")
env_stand<-read.csv("env_stand.csv")
environmental <- cbind(temp.pca_all,env_stand)
write.csv(environmental, file = "environmental.csv")
env<-read.csv("environmental.csv")