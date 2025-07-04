
# documentation of coding for Brackenhurst hedgehog density project
# for paper "Habitat and predator heterogeneity influence density of a declining mammal"
# Hongli Yu, Nottingham Trent University, hongli.yu2020@my.ntu.ac.uk, hongli.yu2021@gmail.com

# = -----------------------------------------------------------------------

library(tidyverse)
library(gdata)
library(oSCR)
library(boot)
# = -----------------------------------------------------------------------

## step 1: read in edf, tdf

# edf, captures,
hog_edf = read.csv("hog.edf.csv")
# note: sex binary: Female = 0, Male = 1

# tdf, trap array, search effort
hog.tdf1 = read.csv("hog.2009.tdf.csv")
hog.tdf2 = read.csv("hog.2010.tdf.csv")
hog.tdf3 = read.csv("hog.2011.tdf.csv")
hog.tdf4 = read.csv("hog.2012.tdf.csv")
hog.tdf5 = read.csv("hog.2013.tdf.csv")
hog.tdf6 = read.csv("hog.2014.tdf.csv")
hog.tdf7 = read.csv("hog.2015.tdf.csv")
hog.tdf8 = read.csv("hog.2017.tdf.csv")
hog.tdf9 = read.csv("hog.2018.tdf.csv")
hog.tdf10 = read.csv("hog.2021.tdf.csv")
hog.tdf11 = read.csv("hog.2022.tdf.csv")
      
# = -----------------------------------------------------------------------

## step 2: format data; make scrFrame; make state space ssDF

hog.data = 
  data2oscr(edf = hog_edf,
            tdf = list(hog.tdf1,hog.tdf2,hog.tdf3,hog.tdf4, hog.tdf5, hog.tdf6,hog.tdf7,hog.tdf8,hog.tdf9,hog.tdf10,hog.tdf11), # 11 tdf files - one for each session
            sess.col = which(colnames(hog_edf)%in%"SESSION"), 
            id.col = which(colnames(hog_edf)%in%"ANIMAL_ID"),
            occ.col = which(colnames(hog_edf)%in%"SO"), 
            trap.col = which(colnames(hog_edf)%in%"grid_50m_ID"),
            sex.col = which(colnames(hog_edf)%in%"SEX"),
            sex.nacode = "NA",
            K = c(42,60,52,60,60,50,14,10,9,60,23), # no. sampling occasion per session
            ntraps = rep(828,11), # no. traps
            trapcov.names = c("hab"),
            tdf.sep = "/")

# make scrFrame
hog.sf = hog.data$scrFrame

# get summary of annual capture history and plot
hog.sf

par(mfrow=c(2,6),mar=c(1,1,1,1),oma=c(0,0,0,0))
plot(hog.sf, ax = F) 

# make make state space ssDF, called habitat.ssDF
# read in habitat covariate values
habitat_before = read.csv("habitat_2009_2014_before_new_badger_sett.csv") # before the new badger sett identified, apply to 2009 to 2014
habitat_after = read.csv("habitat_2015_2022_after_new_badger_sett.csv") # after the new badger sett identified, apply to 2015 to 2022

# make habitat.ssDF (a list file)

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

# example code

m_null1 = oSCR.fit(list(D~1,p0~1,sig~1), hog.sf, habitat.ssDF, trimS=1)

m_null2 = oSCR.fit(list(D~1,p0~sex,sig~1), hog.sf, habitat.ssDF, trimS=1)

m0 = oSCR.fit(list(D~1,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m1 = oSCR.fit(list(D~BUILDING,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m2 = oSCR.fit(list(D~session,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m3 = oSCR.fit(list(D~BADGER,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m4 = oSCR.fit(list(D~EDGE,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m5 = oSCR.fit(list(D~SOIL,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m6 = oSCR.fit(list(D~SOIL + EDGE,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

m17 <- oSCR.fit(list(D~SESSION + SOIL + EDGE + BADGER,p0~sex,sig~sex), hog.sf, habitat.ssDF, trimS=1)

# to run all combinations of covariates for density (D) modelling

# = -----------------------------------------------------------------------

## step 4: model comparison 

f1 = fitList.oSCR(list(m_null1, m_null2, m0, m1, m2, m3, m4, m5, m6, m17),rename=TRUE) # to add all the relevant models

ms1 = modSel.oSCR(f1)

ms1$aic.tab

top1 = m17 # top model (the one with lowest AIC)

m17 # get maximum likelihood estimates and standard errors 

plot(pred_top1)

# = -----------------------------------------------------------------------

## step 5: calculate detection, sigma, and probability of being a male

# detection

pred.df = data.frame(session = factor(c(1)), sex = factor(c(0,1)))

pred.det.top1 = get.real(model = top1, type = "det", newdata = pred.df)

pred.det.top1$sex = factor(pred.det.top1$sex)

pred.det.top1

# sigma

pred_sig = data.frame(session = factor(c(1)), sex = factor(c(0,1)))

pred.sig.top1 = get.real(model = top1, type = "sig", newdata = pred_sig)

pred.sig.top1$sex = factor(pred.sig.top1$sex)

pred.sig.top1

# probability of being a male

inv.logit(-0.394) # psi.constant estimate (in model output) based on top1 m17

# = -----------------------------------------------------------------------

## step 6: calculate estimated density

# example code

# define the values to make prediction for

ED.sex.top1 = get.real(top1, type = "dens") # session-specific density model
ED.sex = ED.sex.top1

# estimated density per session

ED.f.1 = ED.sex[[1]][which(ED.sex[[1]]$sex=="f"),]
ED.m.1 = ED.sex[[1]][which(ED.sex[[1]]$sex=="m"),]

# session 1
sum(ED.f.1$estimate) # total number of females in session 1
sum(ED.m.1$estimate) # total number of males in session 1
sum(ED.f.1$estimate) + sum(ED.m.1$estimate) # total number of females + males in session 1

# session 2
ED.f.2 = ED.sex[[2]][which(ED.sex[[2]]$sex=="f"),]
ED.m.2 = ED.sex[[2]][which(ED.sex[[2]]$sex=="m"),]

sum(ED.f.2$estimate) # total number of females in session 2
sum(ED.m.2$estimate) # total number of males in session 2
sum(ED.f.2$estimate) + sum(ED.m.2$estimate) # total number of females + males in session 2

# to run for all sessions as above

# = -----------------------------------------------------------------------

## step 7: realised density

# realised density per session

pred_top1 = predict.oSCR(top1, hog.sf, habitat.ssDF, override.trim=TRUE)

# mean realised density across sessions
pred_top1_mean = ((pred_top1$r[[1]] + pred_top1$r[[2]] + pred_top1$r[[3]] + pred_top1$r[[4]] + pred_top1$r[[5]] + 
                     pred_top1$r[[6]] + pred_top1$r[[7]]+ pred_top1$r[[8]] + pred_top1$r[[9]] + pred_top1$r[[10]] + pred_top1$r[[11]])/11)

pred_top1$r[[12]] = pred_top1_mean

plot(pred_top1$r[[12]])

# change after the new badger sett appeared in 2015 (session 7), i.e., those after minus those before

pred_top1_change = (pred_top1$r[[7]]+ pred_top1$r[[8]] + pred_top1$r[[9]] + pred_top1$r[[10]] + 
                      pred_top1$r[[11]])/5 - (pred_top1$r[[1]] + pred_top1$r[[2]] + pred_top1$r[[3]] + pred_top1$r[[4]] + pred_top1$r[[5]] + pred_top1$r[[6]])/6

pred_top1$r[[13]] = pred_top1_change

plot(pred_top1$r[[13]])

# = -----------------------------------------------------------------------















