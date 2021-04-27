library(usethis)
library(dplyr)
library(reshape2)
library(tidyverse)

#define ages/years for assessment
age <- 1:15
end.year <- 2017
assess.year <- 1960:end.year

#to read in data
catch <- read.table(file='data-raw/catchage.txt',header=F,col.names= c('Year',paste('Age',5:max(age),sep="")))
landings <- read.table("data-raw/landings.txt", header = F,col.names = c("Year", "Landings"))
SRV.matrix <-  read.table("data-raw/SpringSurvey.txt", header = F,col.names =  c('Year',paste('Age',1:20,sep="")))
FRV.matrix <- read.table(file='data-raw/FallSurvey.txt',header=F,col.names=c('Year',paste('Age',0:20,sep="")))
mat.all <- read.table(file='data-raw/mat.txt',header=F,col.names=c('Year',paste('Age',5:15,sep="")))
catch_wt <- read.table(file='data-raw/catch_wts.txt',header=F,col.names =c('Year',paste('age',5:15,sep="")))
weights <- read.csv("data-raw/stock_weights.csv", header = TRUE)

#################################Survey data
#Fall survey; to get subset of ages for assessment/convert to long format
FRV.temp <- FRV.matrix %>%
  pivot_longer(!Year, names_to = "Age",values_to = "index") %>%
  mutate(temp = as.numeric(factor(Age, levels = paste('Age',0:20,sep=""))))

FRV.plus <- FRV.temp  %>% filter(temp>15) %>%
  group_by(Year) %>%
  mutate(Age_15 = sum(index)) %>%
  distinct(Age_15)

FRV.vec <- FRV.temp %>% filter(temp<16,temp>1) %>%
  select(-temp) %>%
  pivot_wider(names_from = "Age", values_from="index") %>%
  mutate(Age15 = FRV.plus$Age_15) %>%
  pivot_longer(!Year, names_to = "Age",values_to = "index") %>%
  mutate(survey = "Fall")

#Spring survey; to get subset of ages for assessment/convert to long format
SRV.temp <- SRV.matrix %>%
  pivot_longer(!Year, names_to = "Age",values_to = "index") %>%
  mutate(temp = as.numeric(factor(Age, levels = paste('Age',1:20,sep=""))))

SRV.plus <- SRV.temp  %>% filter(temp>14) %>%
  group_by(Year) %>%
  mutate(Age_15 = sum(index)) %>%
  distinct(Age_15)

SRV.vec <- SRV.temp %>% filter(temp<15) %>%
  select(-temp) %>%
  pivot_wider(names_from = "Age", values_from="index") %>%
  mutate(Age15 = SRV.plus$Age_15) %>%
  pivot_longer(!Year, names_to = "Age",values_to = "index") %>%
  mutate(survey = "Spring")

#Spanish survey (removed due to confidentiality issues)
SS.matrix <- NULL
SS.vec<- NULL

#################################Stock weights
#filter and rescale weights and get subset for assessment
sw.matrix <- weights %>% filter(year<=end.year,year>=min(assess.year)) %>%
  mutate(index = springWt_3LNO/1000) %>%
  select(year, age, index) %>%
  pivot_wider(names_from = "age", values_from="index") %>%
  select(-year)

sw.temp <- as.matrix(sw.matrix)

#to get mean value for first NAs
for(i in length(unique(weights$age)):1){
  for(j in length(assess.year):1){
    if(is.na(sw.temp[j,i])){sw.temp[j,i] <- mean(sw.temp[(j+1):(j+3),i])
    break}
  }
}

#to set rest of NAs to mean value
for(i in length(unique(weights$age)):1){
  for(j in length(assess.year):1){
    if(is.na(sw.temp[j,i])){sw.temp[j,i] <- sw.temp[j+1,i]}
  }
}

mean.sw.mat <- as.data.frame(sw.temp[,age])
# mean.sw.mat$Year<-assess.year
# colnames(mean.sw.mat) <- c(paste('Age',age,sep=""),"Year")

#################################Catch numbers/weights
# cwt.vec <- catch_wt %>%
#   melt(id = "Year", variable.name = "Age",value.name = "index")
#
# catch.vec <- catch %>%
#   melt(id = "Year", variable.name = "Age",value.name = "catch") %>%
#   mutate(weight = cwt.vec$index)

#################################maturity
mat.all<- mat.all[,2:12]
catch.wt<-catch_wt[,2:12]

#####crl transformation
cA <- ncol(catch)
catcht <- rowSums(catch[,2:cA])

p_ya <- t(catch[,2:cA]/catcht)

pya_sum <- apply(p_ya, 2, FSA::rcumsum)
pA <- nrow(pya_sum)
pya_sum <- pya_sum[1:(pA-1),]

pi_ya <- p_ya[1:(pA-1),]/(pya_sum)

x_ya <- log(pi_ya/(1-pi_ya))
crl  <- t(x_ya)

iyear <- as.numeric(factor(assess.year))-1
iage <- as.numeric(factor(age))-1

indices <- dplyr::bind_rows(FRV.vec,SRV.vec,SS.vec)%>%
  dplyr::mutate(fs = if_else(survey=="Fall",(10.5/12),(5.5/12))) %>%
  mutate(Age= as.numeric(factor(Age, levels = paste('Age',1:15,sep=""))))

indices <- indices %>%
  mutate(iyear = iyear[match(indices$Year,assess.year)],
                iage = iage[match(indices$Age,age)],
                isurvey = as.numeric(factor(indices$survey))-1,
                surv_age = as.factor(paste(indices$survey,"_",str_pad(indices$Age , 2, pad = "0"),sep='')))

indices <- indices %>%
  mutate(isd = as.numeric(indices$surv_age)-1,
                surv_year = as.factor(paste(indices$survey,"_",indices$Year,sep='')))

indices <- indices %>% mutate(is_year = as.numeric(indices$surv_year)-1)

usethis::use_data(indices, landings, mean.sw.mat,
                  catch.wt,  mat.all, crl, overwrite = TRUE)
