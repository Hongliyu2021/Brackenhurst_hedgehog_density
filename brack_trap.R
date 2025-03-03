
# Documentation of coding for Brackenhurst hedgehog density project, for paper "Spatiotemporal density of hedgehogs in a changing agroecosystem".
# Hongli Yu, Nottingham Trent University, hongli.yu2020@my.ntu.ac.uk, hongli.yu2021@gmail.com

# = -----------------------------------------------------------------------

library(tidyverse)
library(gdata)
library(oSCR)
library(viridisLite)
library(viridis)

# = -----------------------------------------------------------------------

## step 1: read in edf, tdf

# edf, captures,
hog_edf = read.csv("hog.edf.csv")
# note: sex binary: Female = 0, Male = 1

# tdf, trap array, search effort
hog.tdf1 <- read.csv("hog.2009.tdf.csv")
hog.tdf2 <- read.csv("hog.2010.tdf.csv")
hog.tdf3 <- read.csv("hog.2011.tdf.csv")
hog.tdf4 <- read.csv("hog.2012.tdf.csv")
hog.tdf5 <- read.csv("hog.2013.tdf.csv")
hog.tdf6 <- read.csv("hog.2014.tdf.csv")
hog.tdf7 <- read.csv("hog.2015.tdf.csv")
hog.tdf8 <- read.csv("hog.2017.tdf.csv")
hog.tdf9 <- read.csv("hog.2018.tdf.csv")
hog.tdf10 <- read.csv("hog.2021.tdf.csv")
hog.tdf11 <- read.csv("hog.2022.tdf.csv")
      
# = -----------------------------------------------------------------------

## step 2: format data; make scrFrame; make state space ssDF

hog.data <- 
  data2oscr(edf = hog_edf,
            tdf = list(hog.tdf1,hog.tdf2,hog.tdf3,hog.tdf4, hog.tdf5, hog.tdf6,hog.tdf7,hog.tdf8,hog.tdf9,hog.tdf10,hog.tdf11), # 11 tdf files - one for each session
            sess.col = which(colnames(hog_edf)%in%"SESSION"), 
            id.col = which(colnames(hog_edf)%in%"ANIMAL_ID"),
            occ.col = which(colnames(hog_edf)%in%"SO"), 
            trap.col = which(colnames(hog_edf)%in%"grid_50m_ID"),
            sex.col = which(colnames(hog_edf)%in%"SEX"),
            sex.nacode = "NA",
            K = c(42,60,52,60,60,50,14,10,9,60,23), # no of sampling occ. per session
            ntraps = rep(828,11), # no traps - same across sessions
            trapcov.names = c("hab"),
            tdf.sep = "/")

save(hog.data, file = "hog.data.RData")

# make scrFrame
hog.sf <- hog.data$scrFrame

# get summary of annual capture history and plot
hog.sf

par(mfrow=c(2,6),mar=c(1,1,1,1),oma=c(0,0,0,0))
plot(hog.sf, ax = F) 

# make make state space ssDF, called habitat.ssDF here
# since a new badger sett was identified in year 2015, the covariate 'Badger sett' (distance to the nearest badger sett) has changed the values
# read in habitat covariate values
habitat_before = read.csv("habitat_2009_2014_before_new_badger_sett.csv") # before the new badger sett identified, apply to 2009 to 2014
habitat_before$Land_3 = as.character(habitat_before$Land_3)

habitat_after = read.csv("habitat_2015_2022_after_new_badger_sett.csv") # after the new badger sett identified, apply to 2015 to 2022
habitat_after$Land_3 = as.character(habitat_after$Land_3)

# make habitat.ssDF, which is a list

habitat2009= habitat_before%>%mutate(Tr=1)
habitat2010 = habitat_before%>%mutate(Tr=2)
habitat2011 = habitat_before%>%mutate(Tr=3)
habitat2012 = habitat_before%>%mutate(Tr=4)
habitat2013 = habitat_before%>%mutate(Tr=5)
habitat2014 = habitat_before%>%mutate(Tr=6)
habitat2015 = habitat_after %>%mutate(Tr=7)
habitat2017 = habitat_after %>%mutate(Tr=8)
habitat2018 = habitat_after %>%mutate(Tr=9)
habitat2021 = habitat_after %>%mutate(Tr=10)
habitat2022 = habitat_after %>%mutate(Tr=11)

habitat.ssDF = list(list(habitat2009), list(habitat2010), list(habitat2011), 
                    list(habitat2012),list(habitat2013), list(habitat2014),
                    list(habitat2015), list(habitat2017), list(habitat2018), 
                    list(habitat2021), list(habitat2022))

# = -----------------------------------------------------------------------

## step 3: modelling

m_null1 <- oSCR.fit(list(D~1,p0~1,sig~1), hog.sf, habitat.ssDF, trimS=1)

m_null2 <- oSCR.fit(list(D~1,p0~sex,sig~1), hog.sf, habitat.ssDF, trimS=1)

m0 <- oSCR.fit(list(D~1,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m1 <- oSCR.fit(list(D~Building,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m2 <- oSCR.fit(list(D~session,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m3 <- oSCR.fit(list(D~Badger sett,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m4 <- oSCR.fit(list(D~Edge,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m5 <- oSCR.fit(list(D~Soil,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m6 <- oSCR.fit(list(D~Soil + Edge,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

# as above, to run all combinations of covariates for D (density) as referred in main text of the paper.

# = -----------------------------------------------------------------------

## step 4: model comparison 

f1 <- fitList.oSCR(list(m_null1, m_null2, m0,m1,m2, m3,m4,m5,m6),rename=TRUE) # to add all the relevant models

ms1 <- modSel.oSCR(f1)

ms1$aic.tab

# = -----------------------------------------------------------------------
