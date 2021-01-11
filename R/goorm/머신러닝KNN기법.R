getwd()
setwd("E:/data/COIN")
coin<-read.csv("coin.csv",header = TRUE)
View(coin)

summary(coin)
dim(coin)
str(coin)
head(coin)
tail(coin)

N=nrow(coin)
tr.index<- sample(1:N,size = N*0.7,replace = FALSE)
tr.index

train<-coin[tr.index,]
test<-coin[-tr.index,]

x.train<-as.matrix(train[1:8])
x.test<-as.matrix(test[1:8])

y.train<-train$UPDOWN

head(x.test)

install.packages("class")
library(class)

model_1<-knn(x.train,x.test,y.train,k = 5)

table<-table(real = test$UPDOWN,predict=model_1)
table[1,1]+table[2,2]/sum(table)

library(e1071)
tune.out<-tune.knn(x = x.train,y=as.factor(y.train),k=1:100)
tune.out
plot(tune.out)

model_1<-knn(x.train,x.test,y.train,k=30)
table<- table(real = test$UPDOWN,predict = model_1)
(table[1,1]+ table[2,2])/sum(table)

#install.packages("Epi")
library(Epi)

ROC(test = model_1, stat = test$UPDOWN, plot = "ROC", AUC = T,main= "knn")

new.data<- data.frame(Ethereum = 170.91, XRP=0.321024, BITCOIN.CASH = 290.81, BIBOX.TOKEN = 22.90, Nebulas=1.21,
                      EOS = 5.17, BINANCE.COIN = 22.9, CARDANO = 0.078)

summary(new.data)

knn(x.train,new.data,y.train, k = 30)
