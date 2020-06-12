## ----knitr_settings, echo = FALSE----------------------------------------
options(width = 55, 
        repos = c(CRAN = "https://cran.wu.ac.at/"))
set.seed(201702)
rm(list=ls())

## ------------------------------------------------------------------------
jim <- list(height = 1.8, weight = 82, name = "Jim")
class(jim) <- "lecturer"
class(jim)
print(jim)

## ------------------------------------------------------------------------
print.lecturer <- function(x, ...) {
  cat("name:", x$name, "\n")
  cat("height:", x$height, "meters", "\n")
  cat("weight:", x$weight, "kilograms", "\n")
}
print(jim)

## ------------------------------------------------------------------------
myRep <- representation(height = "numeric", weight = "numeric",
name = "character")
setClass("lecturerS4", representation = myRep)
getClass("lecturerS4")
jimS4 <- new("lecturerS4")
jimS4

## ------------------------------------------------------------------------
jimS4 <- new("lecturerS4", height = 1.8, weight = 82, name = "Jim")
jimS4

## ------------------------------------------------------------------------
jimS4@name
validObject(jimS4)
jimS4@height <- "2"

## ------------------------------------------------------------------------
setMethod("show", signature("lecturerS4"),
function(object) {
cat("name:", object@name, "\n")
cat("height:", object@height, "meters", "\n")
cat("weight:", object@weight, "kilograms", "\n")
})
jimS4

## ------------------------------------------------------------------------
getMethod("show", signature("lecturerS4"))

## ------------------------------------------------------------------------
setGeneric("BMI", function(object) standardGeneric("BMI"))
setMethod("BMI", "lecturerS4", function(object) {
object@weight / object@height^2
})
BMI(jimS4)

## ----eval=FALSE----------------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite()

## ----biobase, out.width=".8\\linewidth",fig.width=10,fig.height=6, echo=FALSE----
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0,100),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")
polygon(c(45,80,80,45),c(10,10,70,70),col=rgb(1,0,0,.5),border=NA)
polygon(c(45,80,80,45),c(68,68,70,70),col=rgb(1,0,0,.5),border=NA)
text(62.5,40,"assay(s)", cex = 1)
text(62.5,30,"e.g. 'exprs'", cex = 1)
polygon(c(20,40,40,20),c(10,10,70,70),col=rgb(0,0,1,.5),border=NA)
polygon(c(20,40,40,20),c(68,68,70,70),col=rgb(0,0,1,.5),border=NA)
text(30,40,"featureData", cex = 1)
polygon(c(45,80,80,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
polygon(c(45,47,47,45),c(75,75,90,90),col=rgb(.5,0,.5,.5),border=NA)
text(62.5,82.5,"phenoData", cex = 1)

## ----affy, eval=FALSE----------------------------------------------------
## library(affy)

## ----affy2, echo=FALSE, warning=FALSE, results="hide", message=FALSE-----
library(affy)

## ----data----------------------------------------------------------------
Data <- ReadAffy()

## ----pData, warning=FALSE------------------------------------------------
pData(Data)
Data

## ----exprs---------------------------------------------------------------
summary(exprs(Data))

## ----logexprs------------------------------------------------------------
summary(log2(exprs(Data)))

## ----pm-mm---------------------------------------------------------------
lrp <- grep("lrp",featureNames(Data))
lrp
pm(Data, "lrp_b0889_st")[1:5,]

## ----pmmm, out.width=".7\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
par(mfrow = c(2,2))
matplot(pm(Data, "lrp_b0889_st"), type = "l", xlab = 
        "Probe No.", ylab = "PM Probe Intensity")
matplot(t(pm(Data, "lrp_b0889_st")), type = "l", xlab =
        "Array No.", ylab = "PM Probe Intensity")
matplot(mm(Data, "lrp_b0889_st"), type = "l", xlab = 
        "Probe No.", ylab = "MM Probe Intensity")
matplot(t(mm(Data, "lrp_b0889_st")), type = "l", 
        xlab = "Array No.", ylab = "MM Probe Intensity")

## ----distr, out.width=".5\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
par(mfrow = c(1,1), las = 2, mar = c(6, 4, 2, 2) + 0.1)
boxplot(Data, col = cols, 
        ylab = "unprocessed log scale probe-level data")

## ----density, out.width=".5\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
hist(Data, col = cols, lty = 1, 
     xlab = "Log (base 2) intensities")
legend(14, 1.2, c("nolrp_1", "nolrp_2", "nolrp_3",
       "nolrp_4", "wt_1", "wt_2", "wt_3", "wt_4"), 
       lty=1, col=cols)

## ----rma, out.width=".7\\linewidth",fig.width=6,fig.height=6-------------
eset.rma <- rma(Data)
par(mfrow = c(1,1), las = 2, mar = c(6, 4, 2, 2) + 0.1)
boxplot(exprs(eset.rma), col = cols, main = "RMA")

## ----gcrma, out.width=".7\\linewidth",fig.width=6,fig.height=6, message=FALSE----
library(gcrma)
eset <- gcrma(Data)
par(mfrow = c(1,1), las = 2, mar = c(6, 4, 2, 2) + 0.1)
boxplot(exprs(eset), col = cols, main = "GCRMA")

## ----genefilter, message=FALSE-------------------------------------------
library(genefilter)
## intensities above 200 in at least 25% of the samples
##just so we have fewer genes and it runs in finite time
f1 <- pOverA(.25, log2(200))
f2 <- function(x)(IQR(x)>0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
esetSub <- eset[selected,]

## ----limma, message=FALSE------------------------------------------------
library(limma)
strain <- c("lrp-", "lrp-", "lrp-", "lrp-", "lrp+",
            "lrp+", "lrp+", "lrp+")
design <- model.matrix(~ factor(strain))
colnames(design) <- c("lrp-", "lrp+vs-")

## ------------------------------------------------------------------------
design

## ----toptable------------------------------------------------------------
fit <- lmFit(esetSub, design)
fit <- eBayes(fit)
options(digits = 2)
tt <- topTable(fit, coef = 2, n = 50, adjust = "fdr")

## ------------------------------------------------------------------------
tt[1:10,]

## ----volcano, out.width=".7\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
volcanoplot(fit, highlight = 40, coef = 2, 
            names = rownames(fit))

## ----overview------------------------------------------------------------
apply(fit$coef, 2, function(x)
      table(abs(x)>1))
write.fit(fit, file = "Ecoli_fit.txt", 
          sep = "\t", adjust = "fdr")

## ----heatmap, out.width=".5\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
nn <- exprs(esetSub)[rownames(tt),]
heatmap(nn)

## ----dendr1, out.width=".6\\linewidth",fig.width=10,fig.height=6, echo=FALSE----
genes.distance <- as.dist(1 - cor(t(nn)))
genes.clusters <- hclust(genes.distance, method = "ward.D") 
plot(genes.clusters) 
genes.dendrogr <- as.dendrogram(genes.clusters)
genes.members <- cutree(genes.clusters, k = 2) 

## ----dendr2, out.width=".6\\linewidth",fig.width=6,fig.height=6, echo=FALSE----
samples.distance <- as.dist(1-cor(nn))
samples.clusters <- hclust(samples.distance, 
                           method = "ward.D") 
plot(samples.clusters) 
samples.dendrogr <- as.dendrogram(samples.clusters)
samples.members <- cutree(samples.clusters, k = 2)

## ------------------------------------------------------------------------
rowclust <- genes.members 
rowclust[rowclust == 1] <- "black" 
rowclust[rowclust == 2] <- "red" 

colclust <- samples.members 
colclust[colclust == 1] <- "black" 
colclust[colclust == 2] <-  "red" 

## ----heatmap2, out.width=".6\\linewidth",fig.width=6,fig.height=6, echo=FALSE, message=FALSE----
library(geneplotter)
hv <- heatmap(as.matrix(nn), Rowv = genes.dendrogr, 
              labRow = "",  Colv = samples.dendrogr, 
              RowSideColors = rowclust, 
              ColSideColors = colclust, 
              main = "Heatmap", margin = c(10,5), 
              col = rev(dChip.colors(50)))  

