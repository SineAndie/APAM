library(usethis)
library(tidyverse)

#define ages/years for assessment
age <- 1:15
end.year <- 2017
assess.year <- 1960:end.year

#to read in data
SRV.matrix <-  read.table("data-raw/SpringSurvey.txt", header = F,col.names =  c('Year',paste('Age',1:20,sep="")))
FRV.matrix <- read.table(file='data-raw/FallSurvey.txt',header=F,col.names=c('Year',paste('Age',0:20,sep="")))

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

#age/year indicator
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

usethis::use_data(indices, overwrite = TRUE)
