setwd("C:/PostDoc/PhDExtension/archive")

# adani = read.csv("ADANIPORTS.csv")
# tm = as.numeric(as.POSIXct(adani$Date, origin = "1970-01-01"))
# v = as.numeric(adani$VWAP)


temp = list.files(pattern="\\.csv$")
myfiles = lapply(temp, read.csv)

dat = NULL
b = NULL
for(i in 1:length(myfiles)){
  d = myfiles[[i]]
  
  if(ncol(d)>10){
    d1 = d[(nrow(d)-999):nrow(d),c("Date", "Symbol", "VWAP")]
    d1$VWAPscaled = as.numeric(scale(d1$VWAP))
    
    # d1$Date = as.numeric(as.POSIXct(d1$Date, origin = "1970-01-01"))
    # dat$Date = (dat$Date - min(dat$Date))/86400
    
    d1 = d1[order(d1$Date),]
    dat = rbind(dat, d1)
  }else{
    b = c(b,i)
  }
  
}
dat$Date = as.numeric(as.POSIXct(dat$Date, origin = "1970-01-01"))
dat$Date = dat$Date - min(dat$Date)
dat$Date = dat$Date/86400
# ggplot(data = dat, aes(x = Date,  y = VWAPscaled, group = Symbol, colour = Symbol))+geom_line()


w = which(dat$Symbol == "ZEEL")

#Spectral 
sgfun1 = function(lambda){
  ifelse(abs(lambda)<pi/3, 1, 0)
}
  
sgfun2 = function(lambda){
  ifelse(abs(lambda)>pi/3 & abs(lambda)<(2*pi)/3 , 1, 0)
}

sgfun3 = function(lambda){
  ifelse(abs(lambda)>(2*pi)/3 , 1, 0)
}

sgfun4 = function(lambda){
  3/4*(1-abs(lambda))
}

bgfun1 = function(lambda){
  ifelse(lambda[1]^2 + lambda[2]^2 < 0.5^2, 1, 0)
}
  

bgfun2 = function(lambda){
  sqrt(lambda[1]^2 + lambda[2]^2)
}

bgfun3 = function(lambda){
  cos(2*lambda[1])*sin(3*lambda[2])
}


bgfun4 = function(lambda){
  cos(lambda[1] + lambda[2])
}
# Bispectral

feature_mat = NULL

stock = unique(dat$Symbol)
# pdf("Stock2.pdf")
for(i in 1:length(stock)){
  d = dat[which(dat$Symbol == stock[i]),]
  if(i == 33){d = d[1:1000,]}
  if(nrow(d)>800){
     p= periodogram(diff(d$VWAPscaled))
     period = 1/p$freq[which(p$spec == max(p$spec))]
     spectralFeature  = c(spectral.Mean.Estimate(diff(d$VWAPscaled), sgfun1, 20),
                          spectral.Mean.Estimate(diff(d$VWAPscaled), sgfun2, 20),
                          spectral.Mean.Estimate(diff(d$VWAPscaled), sgfun3, 20),
                          spectral.Mean.Estimate(diff(d$VWAPscaled), sgfun4, 20))
     BispectraFeature = c(bispectral.Mean.Estimate(diff(d$VWAPscaled), bgfun1, 20),
                          bispectral.Mean.Estimate(diff(d$VWAPscaled), bgfun2, 20),
                          bispectral.Mean.Estimate(diff(d$VWAPscaled), bgfun3, 20),
                          bispectral.Mean.Estimate(diff(d$VWAPscaled), bgfun4, 20))
     diffStartEnd = d$VWAPscaled[nrow(d)] - d$VWAPscaled[1]
     meanDiff = mean(diff(d$VWAPscaled))
     maxDiff = max(d$VWAPscaled) - min(d$VWAPscaled)
     feature = data.frame(t(c(period, spectralFeature, BispectraFeature, 
                           diffStartEnd, meanDiff, maxDiff)))
     names(feature) = c("Period", "Spectral.Mean.1", "Spectral.Mean.2",
                        "Spectral.Mean.3", "Spectra.Mean.4",
                        "Bispectral.Mean.1", "Bispectral.Mean.2",
                        "Bispectral.Mean.3", "Bispectral.Mean.4",
                        "DiffEndStart", "MeanDiff", "MaxDiff")
     feature_mat = rbind(feature_mat, feature)
  }
  print(i)
  print(feature_mat)
  # plot(d$Date, d$VWAPscaled, ty="l", main = paste("Stock Name: ", stock[i], sep=""))
  
}
# dev.off()

feature_mat = read.csv("C:/PostDoc/PhDExtension/archive/Output/featureOct11.csv")
feature_mat = feature_mat[,-1]
fm = apply(feature_mat,2,scale)
row.names(fm) = stock

library(factoextra)
library(cluster)
library(fpc)

fviz_nbclust(fm, kmeans, method = "silhouette")
gap_stat <- clusGap(fm, FUN = kmeans, nstart = 25, K.max = 10, B = 10)
fviz_gap_stat(gap_stat)

km.res = kmeans(fm, 3, nstart=25)

grp1 = names(which(km.res$cluster==1))
grp2 = names(which(km.res$cluster==2))
grp3 = names(which(km.res$cluster==3))
pdf("Group3.pdf")
for(i in 1:length(grp3)){
  dg = dat[which(dat$Symbol == grp3[i]),]
  plot(as.numeric(dg$VWAPscaled), ty="l", main = grp3[i])
}
dev.off()


### Random Forest

library(randomForest)
rt.fit = randomForest(x = fm, y = NULL, ntree=10000, proximity = TRUE,
                      oob.prox = TRUE)
hclust.rf <- hclust(as.dist(1-rt.fit$proximity), method = "ward.D2")
rf.cluster = cutree(hclust.rf, k=3)


library(metricsgraphics)
library(readxl)
library(dplyr)
library(ggplot2)
library(randomForest)
library(varImp)
library (NbClust)
library (cluster)
library (clustertend)
library (factoextra)
library (fpc)
library (clValid)
library(hopkins)

hopkins = hopkins(fm, n = nrow(fm)-1)
pam.res3 <- pam(fm, 3,  metric = "euclidean", stand = FALSE)


library(clusterSim)
index.DB(fm, pam.res3$clustering)


fviz_ch <- function(data) {
  ch <- c()
  for (i in 2:10) {
    km <- kmeans(data, i) # perform clustering
    ch[i] <- calinhara(data, # data
                       km$cluster, # cluster assignments
                       cn=max(km$cluster) # total cluster number
    )
  }
  ch <-ch[2:10]
  k <- 2:10
  plot(k, ch,xlab =  "Cluster number k",
       ylab = "Caliński - Harabasz Score",
       main = "Caliński - Harabasz Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
}

library(NbClust)

fviz_db <- function(data) {
  k <- c(2:10)
  nb <- NbClust(data, min.nc = 2, max.nc = 10, index = "db", method = "kmeans")
  db <- as.vector(nb$All.index)
  plot(k, db,xlab =  "Cluster number k",
       ylab = "Davies-Bouldin Score",
       main = "Davies-Bouldin Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(db==min(db)) + 1, lwd=1, col="red", lty="dashed")
}


library(clValid)

fviz_dunn <- function(data) {
  k <- c(2:10)
  dunnin <- c()
  for (i in 2:10) {
    dunnin[i] <- dunn(distance = dist(data), clusters = kmeans(data, i)$cluster)
  }
  dunnin <- dunnin[2:10]
  plot(k, dunnin, xlab =  "Cluster number k",
       ylab = "Dunn Index",
       main = "Dunn Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(dunnin==max(dunnin)) + 1, lwd=1, col="red", lty="dashed")
}

fviz_dunn(df)



library(FeatureImpCluster)
library(flexclust)
library(tibble)

fm2 = data.frame(fm)
row.names(fm2) = as.numeric(as.factor(row.names(fm2)))
final = kmeans(fm2, 5)
m_df <- as.data.frame(fm2)
m_df$cluster <- part1$cluster
m_df$cluster <- as.character(m_df$cluster)

ggpairs(m_df, 1:5, mapping = ggplot2::aes(color = cluster, alpha = 0.5), 
        diag = list(continuous = wrap("densityDiag")), 
        lower=list(continuous = wrap("points", alpha=0.9)))


km2 <- kmeans(fm2, 2)
km3 <- kmeans(fm2, 3)
km4 <- kmeans(fm2, 4)
km5 <- kmeans(fm2, 5)
km6 <- kmeans(fm2, 6)
km7 <- kmeans(fm2, 7)
km8 <- kmeans(fm2, 8)
km9 <- kmeans(fm2, 9)
km10 <- kmeans(fm2, 10)
p1 <- fviz_cluster(km2, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 2") 
p2 <- fviz_cluster(km3, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 3")
p3 <- fviz_cluster(km4, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 4")
p4 <- fviz_cluster(km5, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 5")
p5 <- fviz_cluster(km6, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 6")
p6 <- fviz_cluster(km7, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 7")
p7 <- fviz_cluster(km5, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 8")
p8 <- fviz_cluster(km6, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 9")
p9 <- fviz_cluster(km7, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 10")

library(cowplot)
plot_grid(p1, p2, p3, p4, p5, p6,p7,p8,p9, labels = c("k2", "k3", "k4", "k5", "k6", "k7","k8","k9","k10"))

plot_grid(p1, p2, p3, p4, labels = c("k2", "k3", "k4", "k5"))
plot_grid(p5, p6, p7, p8, labels = c("k6", "k7", "k8", "k9"))


p9 <- fviz_cluster(final, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("k = 5")


part1 = kmeans(fm2, 5)
part2 = pam(fm2, 5)
part3 = clara(fm2, 5)
part4 = fanny(fm2, 5)


p1 <- fviz_cluster(part1, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("") 
p2 <- fviz_cluster(part2,  data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("")
p3 <- fviz_cluster(part3, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("")
p4 <- fviz_cluster(part4, data = fm2, ellipse.type = "convex") + theme_minimal() + ggtitle("")


library(cowplot)
plot_grid(p1, p2, p3, p4, labels = c("kmeans", "PAM", "CLARA","FANNY"))


fviz_cluster(final, data = fm2) + theme_minimal() + ggtitle("k = 5")