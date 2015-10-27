library(intercure)
set.seed(30)


dataset <- read.delim("~/intercure_pkg/smoking_data.txt")
cov_theta = cov_beta = c("SexF","Duration","SI.UC","F10Cigs")
dataset <- plyr::arrange(dataset,Zip)
set.seed(5151)
dataset <- dataset[sample(1:nrow(dataset), 50, replace = TRUE),]
L <- dataset$Timept1
R <- dataset$Timept2
R <- ifelse(is.na(R),Inf,R)
grupo=dataset$Zip
delta<-ifelse(R==Inf,0,1)
M=10
b=0.001
n_int_max=100
BASE=dataset
NAME_DIF=""
ncores=1


dataset

set.seed(5151)
inter_frailty_cl(dataset,L,R,delta,cov_theta,cov_beta,grupo,10,0.0001,par_cl = NULL, outputFiles = TRUE)



library(MLEcens)
data(cosmesis)
cosmesis2 <- data.frame(cosmesis)
cosmesis2$x2 <- ifelse(cosmesis2$x2 == 100, Inf, cosmesis2$x2)
cosmesis2$delta <- ifelse(cosmesis2$x2 == Inf, 0, 1)
set.seed(131)
cosmesis2$group <- rbinom(nrow(cosmesis2), 1, 0.5)


library(doParallel)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
stopCluster(cl)

cl = makeCluster(4,type="SOCK")
registerDoSNOW(cl)


