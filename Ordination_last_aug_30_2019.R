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
str(macros5)
summary(macros5)
macros.log <- log1p(macros5)# computes log(2+x)
decorana(macros5, iweigh = 1, ira = 0)
# DCA performed on log(raw count data). 80 samples 52 species including individual conifer abundances
# argument ira = 0 for detrending and iweight = 1 for downweighting of rare species
# DCA2 2.9 DCA2 2.2
# If raw abundances are used DCA2 = 3.7; DCA2 = 2.2

spe.pca <- prcomp(macros.log, scale. = T, center = T) 
spe.pca
summary(spe.pca)
spe.rda <- rda(macros.log, scale = TRUE, center = TRUE)
spe.rda
summary(spe.rda)

#-------------------------------------Interpolate---------------------------------------------------#
# a try with antartic data with AICC2022
dobson <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/dobson.csv")
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
head(dome.interp.temp)

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
head(co2.interp)
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
# In the original data from A. Rees depth 945 had an age of 24879, but I calculate my own 
# age-depth model. I could interpolated to the LD depths of plant macros
 
chiro_dob <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/Chiro_Dob.csv")
head(chiro_dob)
age.chiro <- chiro_dob$age
depth.chiro <- chiro_dob$depth
dca1.dob <- chiro_dob$DCA1
dca2.dob <- chiro_dob$DCA1
dca.dob <- cbind(dca1.dob, dca2.dob)
chiro.interp <- interp.dataset(y = dca.dob , x = age.chiro, xout = dob.age)
head(chiro.interp)
chiro.interp <- subset(chiro.interp, select = -c(2))
head(chiro.interp)

#---------------------------------------NZ isotopes interpolation----------------------------------#
# Speleothem data from NW South Island NZ Williams 2005
# Updated 20/03/2008

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

# Oxygen isotope variation in Mt. Arthur speleothems primarily represents changes in 
# meteoric waters falling above the caves, possibly responding to latitudinal changes
# in the position of the Subtropical Front in the Tasman Sea. How is the STF in Tasmania compared
# to the NZ one? 
# Carbon isotope variations in the speleothems record, represent changes in forest productivity,
# closely matching existing paleovegetation records in NZ according to Hellstrom 2998.

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
O18.nza <- NZA.interp$O18.nza
C13.nza <- NZA.interp$C13.nza

#----------------------------------------charcoal interpolation------------------------------------#
# Charcoal data is from the same core as tha plant macros so the 2 records should have same age.
# However, counts of charcoal particles were made every 2 cm starting at 0.5 cm and plant macros 
# starting at 2 cm and counted every 20 cm. I interpoletated to LD depth.

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
        x = age.char,
        xout = dob.age,
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
        x = age.loi,
        xout = dob.age,
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
# values based on temperature difference PC1 = 90.54%; PC2 = 9.457%

temp.rda <- rda(antar.temp, scale = FALSE, center = TRUE)
temp.rda
summary(temp.rda)
names(temp.rda)
temp.rda$tot.chi# Total inertia or the sum of all eigenvalues.
temp.rda$CA

# calculate axis-importance and draw the barplots:
ev.antar.temp <- temp.rda$CA$eig
source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/evplot.R")
# source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:evplot?do=export_code&codeblock=2')

dev.new()
# pdf("ev.pdf")

evplot(ev.antar.temp)
# dev.off()
# only PC1 is important 

sites.temp.antar <- scores(temp.rda, display = "sites", choices = 1)
head(sites.temp.antar)
temp.antar <- sites.temp.antar
head(temp.antar)   #only PC1

source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/cleanplot.pca.R")
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
# Arguments scale = TRUE and centre = TRUE calls for a standardization of the variables
# antar.temp include the Vostok and DomeC variacion de deuterium data

# antartic and NZW.O18 together in PCA like in the Thesis, non-including C13.nzw

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

ev.south.temp <- south.temp.pca$sdev^2
# pdf("ev.pdf")

dev.new()
evplot(ev.south.temp)

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

# antartic and NZW.O18 together in PCA like in the Thesis, including C13.nzw

southern.temp1 <- cbind(antar.temp, NZW.interp)
head(southern.temp1)
head(southern.temp1)
south.temp.pca1 <- prcomp(southern.temp1, scale. = T, center = T) 
south.temp.pca1
summary(south.temp.pca1) 
# PC1 = 59.66%; PC2 = 23.86%; PC3 = 22.60%

south.temp.rda1 <- rda(southern.temp1, scale = TRUE, center = T)
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

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.temp.rda1, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.temp.rda1, scaling = 2, mar.percent = 0.04)

# Analysis runned with data from Port Arthur, including only O18
# I dont think this is correct, becasue O18 from Port Arthur
# does not represent temerature, but instead meteoric water 
# 13C represent changes in forest productivitym matching paleovegetation records
southern.temp2 <- cbind(antar.temp, NZA.interp)
head(southern.temp2)
southern.temp2 <- subset(southern.temp2, select = -c(4))
head(southern.temp2)
south.temp2.pca <- prcomp(southern.temp2, scale. = T, center = T) 
south.temp2.pca
south.temp2.rda <- rda(southern.temp2, scale = TRUE, center = T)
south.temp2.rda
summary(south.temp2.rda)

# percentage of variance explained by each  component
(temp.eig <- south.temp2.rda$CA$eig/south.temp2.rda$tot.chi)
summary(south.temp2.pca) 
# PC1 = 78.23%; PC2 = 25.37% 

dev.new()
screeplot(south.temp2.rda, bstick = TRUE)

# select the data frame with eigenvalues of particular axes:
ev.south.temp2 <- south.temp2.pca$sdev^2

# pdf("ev.pdf")
dev.new()
evplot(ev.south.temp2)
# only PC1 is significant

sites.south.temp2.rda <- scores(south.temp2.rda, display = "sites", choices = 1)
# write.csv(sites.south.temp2.rda, file="sites.south.temp2.rda.csv")
# sh.temp2 <- read.csv("sites.south.temp2.rda.csv")
sh.temp2 <- sites.south.temp2.rda
head(sh.temp2)
# temp.comp  contains only PC1

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

#--------------------------------O18 from nza and nzw------------------------------------- 
NZ.isotopes <- cbind(NZW.interp, NZA.interp)
head(NZ.isotopes)
O18.nz <- subset(NZ.isotopes, select = -c(2,4))
head(O18.nz)

O18.pca <- prcomp(O18.nz) 
O18.pca 
summary(O18.pca) # PC1 = 73.86%; PC2 = 26.24%
# percentage of variance explained by each  component

O18.rda <- rda(O18.nz)
O18.rda
summary(O18.rda)
(O18.eig <- O18.rda$CA$eig/O18.rda$tot.chi)

dev.new()
screeplot(O18.rda, bstick = TRUE)
source("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/FunctionsLD/Scree.Plot.R")

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
C13.pca <- prcomp(C13.nz) 
C13.pca 
summary(C13.pca) 
# PC1 = 53.99%; PC2 = 46.01%

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
# PC2 
sites.C13.rda <- scores(C13.rda, display = "sites", choices = 1)
# write.csv(sites.south.temp1.rda, file="sites.south.temp1.rda.csv")
# sh.temp1 <- read.csv("sites.south.temp1.rda.csv")
C13.sh <- sites.C13.rda
head(C13.sh)
# temp.comp  contains only PC2

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(C13.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(C13.rda, scaling = 2, mar.percent = 0.04)

# ------------------------antartic data and O18 isotopes from nzw and nza--------------------------#

southern.climate1 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate1)
southern.climate1 <- subset(southern.climate1, select = -c(4,6))# removing C13 from both records
head(southern.climate1)
south.climate1.pca <- prcomp(southern.climate1, scale. = T, center = T) 
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
# write.csv(sites.south.climate1.rda, file = "sites.south.climate1.rda.csv")
# sh.climate1 <- read.csv("sites.south.climate1.rda.csv")
sh.climate1 <- sites.south.climate1.rda
head(sh.climate1)

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate1.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate1.rda, scaling = 2, mar.percent = 0.04)

# ---------------------------O18 isotopes and C13 from nzw and nza---------------------------------#

southern.climate2 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate2)
south.climate2.pca <- prcomp(southern.climate2, scale. = T, center = T) 
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
# PC1 & PC2 area significant, which is consistent with Kaiser-Guntamman

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

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate2.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate2.rda, scaling = 2, mar.percent = 0.04)
# scaling 2 = distancia entre los sitios o muestras (scaling sites)
# scaling 2, incluye circulo de contribucion de equilibrio # However we can clearly see that 
# scaling 2 = correlacion entre variables representada por el angulo de los vectores (scaling environmental variables)
# I have the impression that I have multicollinearity between co2 and loi
# also perhaps collinearity (correlation) between antartic records  of temperature


# -----------------------O18 from nzw and nza, but including only C13 from nzw---------------------#

southern.climate3 <- cbind(antar.temp, NZW.interp, NZA.interp)
head(southern.climate3)
southern.climate3 <- subset(southern.climate3, select = -c(6))# removing C33 from both records
head(southern.climate3)
south.climate3.pca <- prcomp(southern.climate3, scale. = T, center = T) 
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
# Two first PCA are significant, which is consitent with Kaiser-Guntmman criterion

# select the data frame with eigenvalues of particular axes:
ev <- south.climate3.rda$CA$eig

# pdf("ev.pdf")
dev.new()
evplot(ev)
# first two PC are significant when running the analysis with O28 and C24 
# isotopes data from nzw and nza

sites.south.climate3.rda <- scores(south.climate3.rda, choices = c(1:2), display = "sites")
# write.csv(sites.south.climate2.rda, file = "sites.south.climate2.rda.csv")
# sh.climate2 <- read.csv("sites.south.climate2.rda.csv")
sh.climate3 <- sites.south.climate3.rda
head(sh.climate3)

dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(south.climate3.rda, scaling = 1, mar.percent = 0.08)
cleanplot.pca(south.climate3.rda, scaling = 2, mar.percent = 0.03)

head(temp.antar)
temp.antar <-as.data.frame(temp.antar)
PC1.ant.temp = temp.antar$PC1

head(sh.temp)
sh.temp <-as.data.frame(sh.temp)
PC1.sh.temp = sh.temp$PC1

head(sh.temp1)
sh.temp1 <- as.data.frame(sh.temp1)
PC1.sh.temp1 = sh.temp1$PC1

head(sh.temp2)
sh.temp2 <-as.data.frame(sh.temp2)
PC1.sh.temp2 = sh.temp2$PC1

head(O18.sh)
O18.sh <- as.data.frame(O18.sh)
PC1.sh.O18 = O18.sh$PC1


head(sh.climate1)
sh.climate1 <- as.data.frame(sh.climate1)
PC1.sh.clim1 = sh.climate1$PC1

head(sh.climate2)
sh.climate2 <- as.data.frame(sh.climate2)
PC1.sh.clim2 = sh.climate2$PC1
PC2.sh.clim2 = sh.climate2$PC2

head(sh.climate3)
sh.climate3 <- as.data.frame(sh.climate3)
PC1.sh.clim3 = sh.climate3$PC1
PC2.sh.clim3 = sh.climate3$PC2

clim.complete1 <- cbind(PC1.ant.temp, PC1.sh.temp, PC1.sh.temp1, PC1.sh.temp2, PC1.sh.clim1, 
                        PC1.sh.clim2, PC2.sh.clim2, PC1.sh.clim3, PC2.sh.clim3, co2.interp, loi.interp, 
                        char.interp, et.interp, chiro.interp, PC1.sh.O18)
head(clim.complete1)
is.data.frame(clim.complete1)
as.data.frame(clim.complete1)
#----------------------------------------individual cca--------------------------------------------#
# Individual CCA with each variable log2p(macros); 80 columns 52 variables
macros5 <- read.csv("/Users/Gastorga/Google Drive/Lake Dobson/plant-macros/DataLD/combined_macros_shorted(5).csv")
str(macros5)
macros.log <- log1p(macros5)# computes log(1+x)

# dataset containing 47 species and 80 samples, alpine conifers are clumped together
macros.short <- read.csv("combined_macros_shorted.csv")
str(macros.short)
macros.short.log <- log2p(macros.short)


head(clim.complete1)
clim.complete1 <- as.data.frame(clim.complete1)
antar.mod <- cca (macros.log ~ PC1.ant.temp, data = clim.complete1)
antar.mod
summary(antar.mod)
# PC1.antar.temp explains about 12% (11.50%) of the total variability
control <- how(within = Within(type = "series", mirror = TRUE))
(check2 <- check(macros.log, control))
summary(check2)
(anova(antar.mod,  permutations = control))

#           Df ChiSquare      F Pr(>F)    
# Model     1   0.30994 10.135 0.00625 **
# Residual 78   2.38524                      

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.5 ‘ ’ 1

head(clim.complete1)
sh.temp.mod <- cca (macros.log ~ PC1.sh.temp, data = clim.complete1)
sh.temp.mod
summary(sh.temp.mod)
# vostok temperature record explains about 12% (12.35%) of the total variability
(anova(sh.temp.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.33298 10.995 0.00625 **
# Residual 78   2.36220                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)
sh.temp1.mod <- cca (macros.log ~ PC1.sh.temp1, data = clim.complete1)
sh.temp1.mod
summary(sh.temp1.mod)
# vostok temperature record explains about 12% (11.54%) of the total variability
(anova(sh.temp1.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.31115 10.18 0.00625 **
# Residual 78   2.38403                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

sh.temp2.mod <- cca(macros.log ~ PC1.sh.temp2, data = clim.complete1)  
sh.temp2.mod 
summary(sh.temp2.mod )
# sh.clim.3 PC1 & PC2 explains about ~13% (12.52%) of the constrained variability
(anova(sh.temp2.mod,  permutations = control))
#           Df ChiSquare      F Pr(>F)    
# Model     2   0.39971 6.704 0.00625 **
# Residual 77   2.29547                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
sh.clim1.mod <- cca(macros.log ~ PC1.sh.clim1, data = clim.complete1)  
sh.clim1.mod
summary(sh.clim1.mod)
# sh.clim4 CCA1 & CCA2 explains about ~13% (12.91%) of the constrained variability
(anova(sh.clim1.mod,  permutations = control))

#          Df ChiSquare      F  Pr(>F)   
# Model     1   0.34783 11.558 0.00625 **
# Residual 78   2.34736                       

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

head(clim.complete1)  
sh.clim2.mod <- cca(macros.log ~ PC1.sh.clim2 + PC2.sh.clim2, data = clim.complete1)  
sh.clim2.mod
summary(sh.clim2.mod)
# sh.clim1 using PCA1 & PCA2 explains about ~15% (14.83%) of the constrained variability
# PC1.sh.clim2 explains 12.51% of the constrained variability

(anova(sh.clim2.mod,  permutations = control))
# Df ChiSquare     F  Pr(>F)   
# Model     2   0.39971 6.704 0.00625 **
# Residual 77   2.29547                 

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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
co2.mod <- cca(macros.log ~ co21, data = clim.complete1)  
co2.mod
summary(co2.mod)

# co2 explains about ~13% (12.51%) of the total variability
(anova(co2.mod,  permutation = control))
#           Df ChiSquare      F Pr(>F)    
# Model     1   0.33711 11.151  0.001 ***
# Residual 78   2.35807                

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dca.mod <- cca(macros.log ~ dca1.dob, data = clim.complete1)  
dca.mod
summary(dca.mod)
# Chironomids DCA2analysis explains about ~13% (12.68) of the total variability
(anova(dca.mod,  permutations = control))
#          Df ChiSquare      F Pr(>F)    
# Model     1    0.3417 11.325 0.00625 **
# Residual 78    2.3535  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


loi.mod <- cca(macros.log ~ loi1, data = clim.complete1)  
loi.mod
summary(loi.mod)
# loi explains about 9% (8.897) of the total variability
(anova(loi.mod,  permutations = control))
#          Df ChiSquare      F Pr(>F)    
# Model     1   0.23979 7.6173 0.0375 *
# Residual 78   2.45539                

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

char.mod <- cca(macros.log ~ char1, data = clim.complete1)  
char.mod
summary(char.mod)
# char explains just about 3 % (2.543%) of the total variability
(anova(char.mod,  permutations = control))
# model is marginally significant at the 0.2 level (p-value = 0.08225) 

et.mod <- cca(macros.log ~ temp.et1, data = clim.complete1)  
et.mod
summary(et.mod)
# Eagle Tarn temperature explains about 5% (5.244%) of total variability
(anova(et.mod,  permutations = control))

#         Df ChiSquare      F  Pr(>F)   
# Model     1   0.14134 4.3167 0.00625 **
# Residual 78   2.55384                  

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

O18.mod <- cca(macros.log ~ PC1.sh.O18, data = clim.complete1)  
O18.mod
summary(O18.mod)
# O18 isotopes from nza and nzw explain together 10.66% (CCA1)
(anova(O18.mod,  permutations = control))
#          Df ChiSquare      F  Pr(>F)   
# Model     1    0.2874 9.3102 0.00625 **
# Residual 78    2.4078                  

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

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# age explains about 12% (12.28%) of the total variability 
# permutation test showed the sediment age was significantly 
# correlated with the plant macrofossil data at 0.01 level (p-value = 0.00625)
# therefore sediment age is used as covariable in subsequent ordination

set.seed(200)
head(clim.complete2)
clim.complete2 <- as.data.frame(clim.complete2)
cond.model1 <- cca(macros.log ~ dca1.dob + co21 + PC1.ant.temp + PC1.sh.O18 + loi1 + temp.et1 +
                        char1  + Condition(dob.age), data = clim.complete2)
cond.model1
summary(cond.model1)
#                Inertia Proportion Rank
# Total          2.6952     1.0000     
# Conditional    0.3311     0.1228    1
# Constrained    0.4259     0.1580    7
# Unconstrained  1.9382     0.7191   51

# Accumulated constrained eigenvalues
# Importance of components:
#                       CCA1    CCA2    CCA3    CCA4    CCA5    CCA6    CCA7
# Eigenvalue            0.1963 0.06354 0.05122 0.03929 0.03151 0.02571 0.01828
# Proportion Explained  0.4610 0.14919 0.12026 0.09227 0.07398 0.06037 0.04292
# Cumulative Proportion 0.4610 0.61020 0.73047 0.82273 0.89671 0.95708 1.00000
# 
# Scaling 2 for species and site scores
# * Species are scaled proportional to eigenvalues
# * Sites are unscaled: weighted dispersion equal on all dimensions

control <- how(within = Within(type = "series", mirror = TRUE))
set.seed(60)
(anova(cond.model1, model = "full", permutations = control)) # Overall model is significant
(anova(cond.model1, by = "term", model = "full", permutations = control))
(anova(cond.model1, by = "axis", model = "full", permutations = control)) # only CCA1 is significant
(anova(cond.model1, by = "margin", model = "full", permutations = control))# the effect of a particualr term when all other model terms are included

set.seed(60)
head(clim.complete2)
cond.model2 <- cca(macros.log ~ PC1.sh.clim2 + dca1.dob + co21 + loi1 + temp.et1 + char1 
                   + Condition(dob.age), data = clim.complete2)
cond.model2
summary(cond.model2)
#                Inertia Proportion Rank
# Total          2.6952     1.0000     
# Conditional    0.3311     0.1228    1
# Constrained    0.3844     0.1426    6
# Unconstrained  1.9797     0.7345   51
# Inertia is scaled Chi-square 

# Accumulated constrained eigenvalues
# Importance of components:
#                       CCA2    CCA2    CCA3    CCA4    CCA5    CCA6
# Eigenvalue            0.1967 0.05651 0.03933 0.03719 0.02819 0.02653
# Proportion Explained  0.5116 0.14701 0.10230 0.09674 0.07333 0.06902
# Cumulative Proportion 0.5116 0.65861 0.76091 0.85765 0.93098 1.00000

set.seed(60)
(anova(cond.model2, model = "full", permutations = control)) # Overall model is significant
(anova(cond.model2, by = "term", model = "full", permutations = control))
(anova(cond.model2, by = "axis", model = "full", permutations = control)) # only CCA2 is significant
(anova(cond.model2, by="margin", model = "full", permutations = control))

# inertcomp(cond.model.sh, prop = TRUE) # I dont really know horw to interpret this
# intersetcor(cond.model.sh)

# Manual model building
head(clim.complete2)

# clim.complete3 includes PC1.ant.temp, co21, loi1, char1, temp.et1, dca1.dob, PC1.sh.O18
clim.complete3 <- subset(clim.complete2, select = -c(2:7))
head(clim.complete3)
clim.complete3 <- as.data.frame(clim.complete3)
# clim.complete4 includes PC1.sh.temp, co21, loi1, char1, temp.et1, dca1.dob
clim.complete4 <- subset(clim.complete2, select = -c(1, 3:7,13))
head(clim.complete4)
clim.complete4 <- as.data.frame(clim.complete4)

# clim.complete5 includes PC1.sh.temp, co21, loi1, char1, temp.et1, dca1.dob
head(clim.complete2)
clim.complete5 <- subset(clim.complete2, select = -c(1:2, 4:7,13))
head(clim.complete5)
clim.complete5 <- as.data.frame(clim.complete5)

# clim.complete0 corresponds to the dataset used in my thesis
# the final mod included co2, PC1.sh and char explaining about 27.7% of variability instead of 26%
# and the two firs axis account for 92% of the constrained variability
# instead 89% declared in the thesis. This could be because in the thesis I used a different 
# data matrix 47 species, 80 samples. also because antartic records are in the AICC2022 chronology.

# clim.complete3 includes char2, temp.et2, dca2.dob, loi2, co22, C23.nzw, PC2.sh
# the final model includes variables dca2.dob + co22 + PC2.sh explaining about 28% of the total variability
# with CCA2 = 74% and CCA2 = 20%.
# CCA-based forward selection with 9999 permutations selected the following model:
# macros.log ~ dca2.dob + co22 + PC2.sh + temp.et2 
# parsimonious cca produce an overall significative model (0.00625 **) explaining 29.75% of the total variability
# with CCA2 = 68% and CCA2 = 20% of the constrained variability. CCA2 is significant
# Adj.r.squared = 0.2550955
# VIF = dca2.dob 6.497204, co22 4.947425, PC2.sh 5.372297 temp.et2 2.495523 


#-------------------------------------Defining final model-----------------------------------------#
# add1(m0, scope=formula(mbig), test = c("permutation"), 
#      permutations = how(within = Within(type = "series", mirror = TRUE)))

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
m0 <- update(m0, . ~ . + loi1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")

## -- included variables still significant?

final.mod1 <- cca(macros.log ~  dca1.dob + PC1.sh.O18 + co21 +loi1 + char1, clim.complete3)
final.mod1
#               Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.6128     0.2274    5
# Unconstrained  2.0824     0.7726   51
# Inertia is scaled Chi-square 

summary(final.mod1)
# Accumulated constrained eigenvalues
# Importance of components:
#                        CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3613 0.1370 0.04613 0.03834 0.03002
# Proportion Explained  0.5896 0.2235 0.07528 0.06257 0.04899
# Cumulative Proportion 0.5896 0.8132 0.88844 0.95101 1.00000

vif.cca(final.mod1)
# dca1.dob PC1.sh.O18       co21       loi1      char1 
# 6.028330   3.480614   6.022799   2.381512   1.117351 

drop1(final.mod1, test = c("permutation"), 
      permutations = how(within = Within(type = "series", mirror = F)))
add1(m0, scope = formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))

spe.cca <- cca(macros.log ~., clim.complete3)
summary(spe.cca)
RsquareAdj(spe.cca)

set.seed(200)
cca.step.forward <- ordistep(cca(macros.log ~ 1, data = clim.complete3),
                             scope = formula(spe.cca),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod1.pars.cca <- cca(macros.log ~ dca1.dob + PC1.sh.O18 + co21 + loi1 + char1 , data = clim.complete3)
summary(mod1.pars.cca)
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.6128     0.2274
# Unconstrained  2.0824     0.7726

# Accumulated constrained eigenvalues
# Importance of components:
#                        CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3613 0.1370 0.04613 0.03834 0.03002
# Proportion Explained  0.5896 0.2235 0.07528 0.06257 0.04899
# Cumulative Proportion 0.5896 0.8132 0.88844 0.95101 1.00000

anova(mod1.pars.cca, permutations = control)
#          Df ChiSquare      F  Pr(>F)   
# Model     5   0.61278 4.3551 0.00625 **
# Residual 74   2.08240                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod1.pars.cca, permutations = control, by = "axis") # CCA1 & CCA2 are important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36130 12.8393 0.00625 **
# CCA2      1   0.13698  4.8678 0.00625 **
# CCA3      1   0.04613  1.6392 0.30625   
# CCA4      1   0.03834  1.3625 0.38125   
# CCA5      1   0.03002  1.0668 0.47500   
# Residual 74   2.08240                   

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod1.pars.cca)
# $r.squared
# [1] 0.2273612

# $adj.r.squared
# [1] 0.1756196
vif.cca(mod1.pars.cca)
# ca1.dob PC1.sh.O18       co21       loi1      char1 
# 6.028330   3.480614   6.022799   2.381512   1.117351 



head(clim.complete4)
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
m0 <- update(m0, . ~ . + temp.et1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + loi1)
add1(m0, scope = formula(mbig), test = "permutation")

## -- included variables still significant?

final.mod2 <- cca(macros.log ~  dca1.dob + co21 + PC1.sh.temp + 
                    temp.et1 + loi1, clim.complete4)
final.mod2
#               Inertia Proportion Rank
# Total          2.6952     1.0000     
# Constrained    0.5748     0.2133    5
# Unconstrained  2.1204     0.7867   51
# Inertia is scaled Chi-square 

summary(final.mod2)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3666 0.1080 0.03957 0.03248 0.02804
# Proportion Explained  0.6379 0.1880 0.06885 0.05652 0.04879
# Cumulative Proportion 0.6379 0.8258 0.89469 0.95121 1.00000

vif.cca(final.mod2)
# dca1.dob        co21 PC1.sh.temp    temp.et1        loi1 
# 8.279884    5.005508    5.381421    1.505070    2.363445 

drop1(final.mod2, test = c("permutation"), 
      permutations = how(within = Within(type = "series", mirror = F)))
add1(m0, scope = formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))

spe.cca1 <- cca(macros.log ~., clim.complete4)
summary(spe.cca1)


RsquareAdj(spe.cca1)
# $r.squared
# [1] 0.2289893
# 
# $adj.r.squared
# [1] 0.1653406

set.seed(200)
cca.step.forward1 <- ordistep(cca(macros.log ~ 1, data = clim.complete4),
                             scope = formula(spe.cca1),
                             direction = "forward",
                             permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod2.pars.cca <- cca(macros.log ~ dca1.dob + co21 + PC1.sh.O18 + temp.et1, data = clim.complete4)
summary(mod2.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.5502     0.2041
# Unconstrained  2.1450     0.7959

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4
# Eigenvalue            0.3573 0.1335 0.03128 0.02807
# Proportion Explained  0.6494 0.2427 0.05686 0.05103
# Cumulative Proportion 0.6494 0.8921 0.94897 1.00000

anova(mod2.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     4   0.55015 4.809 0.00625 **
# Residual 75   2.14503                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod2.pars.cca, permutations = control, by = "axis") # CCA1 & CCA2 are important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.35725 12.4912 0.00625 **
# CCA2      1   0.13355  4.6694 0.00625 **
# CCA3      1   0.03128  1.0937 0.82500   
# CCA4      1   0.02807  0.9816 0.48125   
# Residual 75   2.14503                   

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod2.pars.cca)
# $r.squared
# [1] 0.2041252
# 
# $adj.r.squared
# [1] 0.1617782

vif.cca(mod2.pars.cca)
# dca1.dob       co21 PC1.sh.O18   temp.et1 
# 4.949813   6.038318   3.937379   1.675081



head(clim.complete5)
mbig <- cca(macros.log ~  ., clim.complete5)
## -- define an empty model to start with
m0 <- cca(macros.log ~ 1, clim.complete5)
## -- manual selection and updating
add1(m0, scope = formula(mbig), test = "permutation") 
m0 <- update(m0, . ~ . + PC1.sh.clim2)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + dca1.dob)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + co21 )
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + char1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + loi1)
add1(m0, scope = formula(mbig), test = "permutation")
m0 <- update(m0, . ~ . + temp.et1)

## -- included variables still significant?

final.mod3 <- cca(macros.log ~  PC1.sh.clim2 + dca1.dob + co21 + char1 + 
                    loi1, clim.complete5)
final.mod3
#             Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.5704     0.2116
# Unconstrained  2.1248     0.7884

summary(final.mod3)
# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4    CCA5
# Eigenvalue            0.3661 0.08781 0.05217 0.03792 0.02646
# Proportion Explained  0.6417 0.15394 0.09146 0.06648 0.04639
# Cumulative Proportion 0.6417 0.79566 0.88712 0.95361 1.00000

vif.cca(final.mod3)
# PC1.sh.clim2     dca1.dob         co21        char1         loi1 
# 8.384795         7.438517     6.539624     1.123602     2.404035 

drop1(final.mod3, test = c("permutation"), 
      permutations = how(within = Within(type = "series", mirror = F)))
add1(m0, scope = formula(mbig), test = c("permutation"), 
     permutations = how(within = Within(type = "series", mirror = TRUE)))

spe.cca2 <- cca(macros.log ~., clim.complete5)
summary(spe.cca2)


RsquareAdj(spe.cca2)
#$r.squared
# [1] 0.2255978
# 
# $adj.r.squared
# [1] 0.1621449

set.seed(200)
cca.step.forward2 <- ordistep(cca(macros.log ~ 1, data = clim.complete5),
                              scope = formula(spe.cca2),
                              direction = "forward",
                              permutations = how(nperm = 9999))

# Parsimonious cca using  selected variables
mod3.pars.cca <- cca(macros.log ~ PC1.sh.clim2 + dca1.dob + co21, data = clim.complete5)
summary(mod3.pars.cca)
# Partitioning of scaled Chi-square:
#               Inertia Proportion
# Total          2.6952     1.0000
# Constrained    0.4816     0.1787
# Unconstrained  2.2136     0.8213

# Accumulated constrained eigenvalues
# Importance of components:
#                         CCA1   CCA2    CCA3    CCA4
# Eigenvalue            0.3632 0.08475 0.03364
# Proportion Explained  0.7542 0.17597 0.06986
# Cumulative Proportion 0.7542 0.93014 1.00000

anova(mod3.pars.cca, permutations = control)
#          Df ChiSquare     F  Pr(>F)   
# Model     3   0.48161 5.5118 0.00625 **
# Residual 76   2.21357                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(mod3.pars.cca, permutations = control, by = "axis") # CCA1 & CCA2 are important components
#          Df ChiSquare       F  Pr(>F)   
# CCA1      1   0.36322 12.4706 0.00625 **
# CCA2      1   0.08475  2.9097 0.14375   
# CCA3      1   0.03364  1.1551 0.46250   
# Residual 76   2.21357                     

#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

RsquareAdj(mod3.pars.cca)
# $r.squared
# [1] 0.1786929
# 
# $adj.r.squared
# [1] 0.1464243

vif.cca(mod3.pars.cca)
# PC1.sh.clim2  dca1.dob        co21 
# 8.156676     6.126492     6.083592 


# Scaling 1: species scores scaled to the relative eigenvalues,
# sites are weighted averages of the species
dev.new()
par(mfrow = c(2,2))
plot(mod1.pars.cca, scaling = 2, display = c("sp", "lc", "cn"),
     main ="Triplot CCA macros.log ~ clim.complete3 - Scaling 1")
# Default scaling 2: site scores sclaed to the relative eigenvalues,
# species are weigjted averages of the sites
plot(spe.cca.pars,
     display = c("lc", "cn"),
     main = "Triplot CCA macros.log ~ clim.complete3 - Scaling 1")



model.fit2 <- ordistep(mod0, scope = formula(mod2), direction = "forward")
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

cca.final.forward <- ordistep(cca(macros.log ~ 2, data = clim.complete3), scope = formula(mod2),
                      direction = "forward", permutations = how(nperm = 999))

summary(cca.final.forward)


## Automatic model building based on AIC but with permutation tests
step(cca(dune ~  2, dune.env), reformulate(names(dune.env)), test="perm")
## see ?ordistep to do the same, but based on permutation P-values
## Not run:
(ordistep(cca(dune ~  2, dune.env), reformulate(names(dune.env))))



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
dev.new()
hist(antar.temp$temp.dome2)
hist(antar.temp$temp.vos2)

skewness(antar.temp$temp.dome2)# kutosis = -0.9 and kurtosis = 3.4
kurtosis(antar.temp$temp.dome2)

# Histogram from a single numeric vector 
# ggplot2.histogram(data=numVector)
# Basic histogram plot from the vector "weight"
(range(temp.dome2))

dev.new()
par(mfrow=c(2, 2))
dome <- ggplot2.histogram(data = antar.temp, xName = 'temp.dome2', scale = "density", binwidth = .7)
vos <- ggplot2.histogram(data = antar.temp, xName = 'temp.vos2', scale = "density", binwidth = .7)
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
(CV <- (200*sdrow/meanrow))# CV = 92.6 meaning that relativisation would have a moderate effect on results

(meancol <- mean(coltotals))
(sdcol <- sd(coltotals))
(CV <- (200*sdcol/meancol))# CV = 297

macros.log <- log2p(macros5)


DCA.macros <- decorana(macros.log, iweigh = 2, ira = 0) 
summary(DCA.macros)
DCA.macros$fraction #downweighting abundance fraction starts at 5 
cca(macros.log)
dca.samplescores <- scores(DCA.macros,display = c("sites"), choices = 2)
write.csv(dca.samplescores, file = "dca.samplescores.csv")

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



# setwd("/Users/Gastorga/Desktop/chapter_one")#modern sediments
#modern sediments
modern<-read.csv("modern_sediment.csv")
modern_percenta<-percenta(modern,first=2,last=25)
write.csv(modern_percenta, file = "modern_percenta.csv")
modern_per<-read.csv("modern_percenta.csv")
summary(modern_per)
modern_log<-log(modern_per+2)
#DCA modern
#****************************************************************************************************
DCA.modern<-decorana(modern_log,iweigh=2,ira=0) #performs DCA and downweights rare taxa
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
CA.modern<-decorana(modern_log,ira=2) #performs CA (actually reciprocal averaging)
summary(CA.modern)
modern_scores_CA<-CA.modern$rproj##sample scores
modern_scores_CA
write.csv(modern_scores_CA, file = "modern_scores_CA.csv")
plot(CA.modern, dis="sp", type="n")
sel <- orditorp(CA.modern, dis="species",   pcol = "red", pch="+")
plot(CA.modern,type="n",main="Correspondence Analysis") #creates empty plot(type=”n”)
text(CA.modern,display="species",col="red",cex=.6) #adds red labels for species
points(CA.modern,display="sites",pch=26)

#surveys
lakeside<-read.csv("lakeside.csv")
summary(lakeside)
lakeside_percenta<-percenta(lakeside,first=2,last=62)
write.csv(lakeside_percenta, file = "lakeside_percenta.csv")

#****************************************************************************************************

mod<-read.csv("modern.csv")
modern_stand<- decostand(mod2, "standardize")
write.csv(env_stand,file="env_stand.csv")
env_stand<-read.csv("env_stand.csv")
environmental <- cbind(temp.pca_all,env_stand)
write.csv(environmental, file = "environmental.csv")
env<-read.csv("environmental.csv")