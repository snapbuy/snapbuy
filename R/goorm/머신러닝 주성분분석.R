getwd()
setwd("C:/Users/snapb/Dropbox/R")
coin<-read.csv("C:/Users/snapb/Dropbox/구름/코드 및 데이터_구름/coin.csv",header=TRUE)
head(coin,10)
tail(coin,10)
dim(coin)
str(coin)
summary(coin)

na.fail(coin)
coin[!complete.cases(coin),]

cor(coin[2:9])
coin_2<-(coin[,-8])
cor(coin_2[2:8])

C.pca<-prcomp(coin_2[2:8],scale. = T)
C.pca
summary(C.pca)

plot(C.pca,type='l')
screeplot(C.pca)

W.pca<-as.matrix(coin_2[,2:8])%*%C.pca$rotation
head(W.pca)

coin.pc<-cbind(as.data.frame(W.pca),coin_2$BITCOIN)
head(coin.pc)
