library(tidyverse)

#define ages/years for assessment
age <- 1:15
catch <- read.table(file='data-raw/catchage.txt',header=F,col.names= c('Year',paste('Age',5:max(age),sep="")))[,2:12]
land <- read.table("data-raw/landings.txt", header = F,col.names = c("Year", "Landings"))

#crl transformation
cA <- ncol(catch)
catcht <- rowSums(catch)

p_ya <- t(catch/catcht)

pya_sum <- apply(p_ya, 2, FSA::rcumsum)
pA <- nrow(pya_sum)
pya_sum <- pya_sum[1:(pA-1),]

pi_ya <- p_ya[1:(pA-1),]/(pya_sum)

x_ya <- log(pi_ya/(1-pi_ya))
crl.mat  <- t(x_ya)

load("R/sysdata.rda")

usethis::use_data(matur.mat,cw.mat,crl.mat, land, internal=T,overwrite = T)
