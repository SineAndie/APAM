library(tidyverse)

#define ages/years for assessment
age <- 1:15
end.year <- 2017
assess.year <- 1960:end.year
weights <- read.csv("data-raw/stock_weights.csv", header = TRUE)

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

sw.mat <- as.data.frame(sw.temp[,age])

usethis::use_data(sw.mat, overwrite = TRUE)

