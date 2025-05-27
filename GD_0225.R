#### install and library the pacakge ####
install.packages("GD")
library("GD")

#### set work path ####
setwd("C:/Users/27519/Desktop/")

#### try NDVI ####
df <- read.csv("MODIS_Topt.csv")

#### discretation
## set optional parameters of optimal discretization
## optional methods: equal, natural, quantile, geometric, sd and manual
discmethod <- c("equal","quantile","natural","geometric","sd")
discitv <- c(6:10)
# optimal discretization
odc1 <- optidisc(Topt ~ Tgs, data = df, discmethod, discitv)
odc1
plot(odc1)

#### one-step: disc + factor detector
## "gdm" function
for (idisc in discmethod)
{
  for (itv in discitv){
      print(idisc, itv)
      ndvigdm <- gdm(Topt ~ Tgs + PFT,
                 continuous_variable = c("Tgs"),
                 data = df,
                 discmethod = idisc, discitv = itv)
      print(ndvigdm$Factor.detector$Factor$qv)
      print(ndvigdm$Factor.detector$Factor$sig)
  }

}
#ndvigdm <- gdm(Topt ~ Tgs + PFT,
 #              continuous_variable = c("Tgs"),
  #             data = df,
   #            discmethod = discmethod, discitv = discitv)
#ndvigdm$Factor.detector$Factor$qv
# ndvigdm$Factor.detector$Factor$sig
# plot(ndvigdm)

#### test robustness ####


#### try ndvi example
## Example 1
## NDVI: ndvi_40
## set optional parameters of optimal discretization
## optional methods: equal, natural, quantile, geometric, sd and manual
discmethod <- c("equal","natural","quantile")
discitv <- c(4:6)
## "gdm" function
## In this case, Climatezone and Mining are categorical variables,
## and Tempchange and GDP are continuous variables.
ndvigdm <- gdm(NDVIchange ~ Climatezone + Mining + Tempchange + GDP,
               continuous_variable = c("Tempchange", "GDP"),
               data = ndvi_40,
               discmethod = discmethod, discitv = discitv) # ~3s
ndvigdm
plot(ndvigdm)

## Example 2
## H1N1: h1n1_100
## set optional parameters of optimal discretization
discmethod <- c("equal","natural","quantile","geometric","sd")
discitv <- c(3:7)
continuous_variable <- colnames(h1n1_100)[-c(1,11)]
## "gdm" function
h1n1gdm <- gdm(H1N1 ~ .,
               continuous_variable = continuous_variable,
               data = h1n1_100,
               discmethod = discmethod, discitv = discitv)
h1n1gdm
plot(h1n1gdm)
