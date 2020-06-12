## ----knitr_settings, echo = FALSE----------------------------------------
options(width = 55, 
        repos = c(CRAN = "https://cran.at.r-project.org"))
set.seed(201702)

## ----s3d, out.width=".7\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
require(scatterplot3d)
#gute Werte
x<-c(0.84, -0.17,  1.48, -0.79,  0.545, -1.681, -0.47, -0.38,  0.03,  0.877)
y<-c(-0.40, -0.00, -0.20, -0.07, -0.98, 0.47, -1.200, -1.01, -0.537,  0.189)
z<- c(-2.02,  2.21,  1.46, -3.297, -1.63, -1.64,-0.33,  1.08,  0.97,  1.601)
data <- cbind(x,y,z)
m<-lm(z~x+y)
#3d plot
s3d<-scatterplot3d(x,y,z, type="h", pch=16,xlab="x1",ylab="x2",zlab = "y",xlim=c(-3,3),ylim=c(-3,3),zlim=c(-5,5), angle=40)
s3d$plane3d(m)

## ----s3d1, out.width=".5\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
require(scatterplot3d)
#gute Werte
x1<-c(0.84, -0.17,  1.48, -0.79,  0.545, -1.681, -0.47, -0.38,  0.03,  0.877)
x2<-c(-0.40, -0.00, -0.20, -0.07, -0.98, 0.47, -1.200, -1.01, -0.537,  0.189)
x3<- c(-2.02,  2.21,  1.46, -3.297, -1.63, -1.64,-0.33,  1.08,  0.97,  1.601)
data <- data.frame(cbind(x1,x2,x3))
#3d plot
s3d<-scatterplot3d(x1,x2,x3, type="h", pch=16, xlim=c(-3,3),ylim=c(-3,3),zlim=c(-5,5), angle=40)

## ------------------------------------------------------------------------
library(flexclust)
data(milk)
summary(milk)

## ----milk, out.width=".6\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
plot(milk)

## ----milk2, out.width=".6\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
milkpca <- prcomp(milk)
biplot(milkpca, cex=c(0.8,1))

## ----s3d2, out.width=".9\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
s3d<-scatterplot3d(x1,x2,x3, type="h", pch=16, xlim=c(-3,3),ylim=c(-3,3),zlim=c(-5,5), angle=40)
pcData <- prcomp(data, scale=TRUE) #, center=F, scale=F)
rotationVecs <- rbind(c(0,0,0),pcData$rotation[,1],c(0,0,0),pcData$rotation[,2],c(0,0,0),pcData$rotation[,3])
s3d$points3d(rotationVecs*3,  col="red", lwd=2, type = 'l')

## ----s3d3, out.width=".9\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
plot(pcData$x[,1],pcData$x[,2], xlab="PC1", ylab="PC2")

## ----pcanv, out.width=".8\\linewidth",fig.width=10,fig.height=6, echo=FALSE----
library(mvtnorm)
sigma <- matrix(c(1,0.7,0.7,1), ncol=2)
x <- rmvnorm(n=5000, mean=c(0,0), sigma=sigma)
plot(x,main="Comparison MND Data and PCs")
#eigen(sigma)

## ----pcanv1, out.width=".8\\linewidth",fig.width=10,fig.height=6, echo=FALSE----
plot(x,main="Comparison MND Data and PCs")
abline(h=0,col="red")
abline(v=0,col="red")

## ----pcanv2, out.width=".8\\linewidth",fig.width=10,fig.height=6, echo=FALSE----
plot(x,main="Comparison MND Data and PCs")
y<-prcomp(x)
abline(h=0,col="red")
abline(v=0,col="red")
points(y$x,col="red")
legend(-4,3,c("MND","PC"),col=c(1,2),pch = c(1, 1))

## ----pca-----------------------------------------------------------------
head(USArrests)
apply(USArrests, 2, var)

## ------------------------------------------------------------------------
prcomp(USArrests)  # inappropriate

## ----pcaplot,out.width=".55\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
plot(prcomp(USArrests))

## ------------------------------------------------------------------------
pr.out <- prcomp(USArrests, scale=TRUE)
pr.out

## ------------------------------------------------------------------------
summary(pr.out)

## ----pcaplotscale,out.width=".55\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
plot(pr.out)

## ----biplot,out.width=".65\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
biplot(pr.out, cex=c(0.8,1.2), scale=0)

## ------------------------------------------------------------------------
names(pr.out)
pr.out$rotation

## ------------------------------------------------------------------------
head(pr.out$x)

## ----cluster-dichte, echo=FALSE, out.width=".7\\linewidth",fig.width=6,fig.height=6----
## The following function is copied from package mlbench
mlbench.smiley <- 
function (n = 500, sd1 = 0.1, sd2 = 0.05)
{
    n1 <- round(n/6)
    n2 <- round(n/4)
    n3 <- n - 2 * n1 - n2
    x1 <- cbind(rnorm(n1, -0.8, sd1), rnorm(n1, 1, sd1))
    x2 <- cbind(rnorm(n1, 0.8, sd1), rnorm(n1, 1, sd1))
    x3 <- cbind(runif(n2, -0.2, 0.2), runif(n2, 0, 0.75))
    x3[, 1] <- x3[, 1] * (1 - x3[, 2])
    x4 <- runif(n3, -1, 1)
    x4 <- cbind(x4, x4^2 - 1 + rnorm(n3, 0, sd2))
    x <- retval <- list(x = rbind(x1, x2, x3, x4), classes = factor(c(rep(1,
        n1), rep(2, n1), rep(3, n2), rep(4, n3))))
    class(retval) <- c("mlbench.smiley", "mlbench")
    retval
}
x <- mlbench.smiley()
smiley <- x$x
colnames(smiley) <- NULL
plot(smiley, col=as.integer(x$classes))

## ----cluster-quantisierung, echo=FALSE, out.width=".7\\linewidth",fig.width=6,fig.height=6----
classes <- cutree(hclust(dist(x$x), "average"), 9)
plot(smiley, col=classes, xlab="")

## ----smiley, echo=TRUE, out.width=".8\\linewidth",fig.width=8,fig.height=6----
par(mfrow = c(1,2))
plot(hclust(dist(smiley), "single"), labels = FALSE)
plot(hclust(dist(smiley), "average"), labels = FALSE)

## ----dentitio------------------------------------------------------------
dentitio <- read.table("dentitio.txt")
dim(dentitio)

## ------------------------------------------------------------------------
dentitio[1:5,]

## ----hclust, echo=TRUE, out.width=".99\\linewidth",fig.width=8,fig.height=6----
dent.dist <- dist(dentitio, method="manhattan")
dent.hclust <- hclust(dent.dist, method="average")
dent.hclust
plot(dent.hclust)
table(cutree(dent.hclust, 3))

## ----heatmap, echo=TRUE, out.width=".7\\linewidth",fig.width=6,fig.height=6----
table(as.matrix(dentitio))
breaks <- c(-1:3, 8)
heatmap(as.matrix(dentitio), 
        hclustfun=function(x) hclust(x, "average"), 
        distfun= function(x) dist(x, "manhattan"),
        scale="none", breaks=breaks,
        col = heat.colors(5))
legend(0,1, legend = c(0:3, "4+"), fill = heat.colors(5), 
       bg = "white")

## ----k-means-var, echo=TRUE, out.width=".99\\linewidth",fig.width=10,fig.height=6----
smiley.km1 <- kmeans(smiley, centers=9)
smiley.km2 <- kmeans(smiley, centers=9)
sum(smiley.km1$withinss)
sum(smiley.km2$withinss)
par(mfrow = c(1,2))
plot(smiley, col = smiley.km1$cluster)
plot(smiley, col = smiley.km2$cluster)

## ----k-means-anz, echo=TRUE, out.width=".7\\linewidth",fig.width=6,fig.height=6----
smiley.kmeans <- lapply(2:10, function(i) 
   kmeans(smiley, i, nstart = 20))
barplot(sapply(smiley.kmeans, function(x) 
   sum(x$withinss)), names.arg = 2:10)

## ----k-means, echo=TRUE, out.width=".99\\linewidth",fig.width=10,fig.height=6----
par(mfrow = c(1,2))
plot(smiley, col = smiley.kmeans[[3]]$cluster)
plot(smiley, col = smiley.kmeans[[5]]$cluster)

