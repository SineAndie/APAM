
mat.all <- read.table(file='data-raw/mat.txt',header=F,col.names=c('Year',paste('Age',5:15,sep="")))
catch_wt <- read.table(file='data-raw/catch_wts.txt',header=F,col.names =c('Year',paste('age',5:15,sep="")))

matur.mat <- mat.all[,2:12]
cw.mat <-catch_wt[,2:12]

usethis::use_data(matur.mat,cw.mat,internal=T,overwrite = T)

