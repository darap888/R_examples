library(tidyverse)
#Here is just a working R script to look at data 
#from LC-MS/MS analysis

#There is a commercial pipeline which produce some common statistical analysis
#like PCA, t-test, volcano-plots. 
#However, this statistics are based on assumption of data normality 
#which I suspected to be violated and decided to repeat data transformation and normalization
#but on raw data outputted from MZMine tool (so I controlled data processing)


#Table of LC-MS/MS features abundancies derived from MzMine
all_old_feattab = read.csv("all_old_feattab.txt")

#loook for outliers, assess quality of samples based on total abundancies
peaksum=apply(all_old_feattab[,-c(1:3)],2,function(x) sum(x,na.rm = TRUE))
summary(peaksum)
names(peaksum)
#only for selected samples
summary(peaksum[1:135])
boxplot(peaksum[1:135]) 
peaksum[1:135][peaksum[1:135]<quantile(peaksum[1:135])[2]]

sort(peaksum[1:137])[1:10]

sort(peaksum[1:137]) # 34 = check it? K1_K5 ok, but K3 not?
##########
#Here we change the column names to more appropriate ones ################
feattab.col0=str_remove(colnames(all_old_feattab),"^X")
feattab.col0=str_remove(feattab.col0,".Peak.area")
#read metadata
meta=read.table("meta3_lcms_inter28days_zk.txt",header = T)


match(feattab.col0,meta$filename)
feattab.col0[is.na(match(feattab.col0,meta$filename))]

feattab.meta=all_old_feattab[,!(is.na(match(feattab.col0,meta$filename)))]
feattab.col.meta=feattab.col0[!(is.na(match(feattab.col0,meta$filename)))]
order=match(feattab.col.meta,meta$filename)

colnames(feattab.meta)=meta$ATTRIBUTE_Plates[order]
rownames(feattab.meta)=all_old_feattab$row.ID

###############

#####################Transformation and normalization #####################################
###########instread of autoscaling (each variable is mean centered and divided by SD) 
#[bad=inflating of small variables which may contain measurement errors]
### use Pareto scaling (uses sqrt(SD) and covariance)

replacezero = function(x) "[<-" (x, !x | is.na(x), min(x[x>0],na.rm = TRUE)/2)


###############################  
##########################


paretoscale = function(rawtab) {
  rowmean = apply(rawtab,1,mean)
  rowsd.sqrt = sqrt(apply(rawtab,1,sd))
  scaledtab = sweep(rawtab,1,rowmean,"-") #mean centered
  scaledtab = sweep(scaledtab,1,rowsd.sqrt,"/") #divided
  return(scaledtab)
}

########### Transformation and normalization

#first delete rows (features) not detected in our selected samples:
sum(rowSums(feattab.meta)==0)
feattab.meta.nonzero=feattab.meta[rowSums(feattab.meta)!=0,]

feattab.meta.nonzero2=apply(feattab.meta.nonzero,2,replacezero)

logdata=log(feattab.meta.nonzero2,2)

pareto.logdata=paretoscale(t(logdata))


########check normality = not-normal:

#on random feature
par(mfrow=c(3,2))    # set the plotting area into a 1*2 array

hist(t(feattab.meta.nonzero2)[,1145])
qqnorm(t(feattab.meta.nonzero2)[,1145])
qqline(t(feattab.meta.nonzero2)[,1145], col="red")

hist(t(logdata)[,1145])
qqnorm(t(logdata)[,1145])
qqline(t(logdata)[,1145], col="red")

hist(pareto.logdata[,1145])

qqnorm(pareto.logdata[,1145])
qqline(pareto.logdata[,1145], col="red")  #looks more normal than raw data
shapiro.test(pareto.logdata[,1145])  #However, not normal based on the test

#

##########
sum(is.na(feattab.meta.nonzero2))
peaksums=colSums(feattab.meta.nonzero2) # need to be sample-normalized

# Sum normalization function - apply to each column
sumnorm <- function(z) {
  colsum <- apply(z, 2, sum)
  rv <- sweep(z, 2, colsum, "/")
  return(rv)
}
# Sum normalize data  = so same sum of intensities in each sample
set1.norm <- as.data.frame(sumnorm(feattab.meta.nonzero2))

# Log transform data
logd1.norm = log2(set1.norm)
plates=colnames(logd1.norm)


sc.col=str_detect(plates,"sc")
iSA.col=str_detect(plates,"i..SA")

# Perform t-test to get p-values = based on selected groups
pvalue.sc_iSA <- apply(logd1.norm, 1, function(x) { t.test(x[str_detect(plates,"sc")], x[str_detect(plates,"i..SA")])$p.value } )

pvalue.sc_iSA.wilcox <- apply(logd1.norm, 1, function(x) { wilcox.test(x[str_detect(plates,"sc")], x[str_detect(plates,"i..SA")])$p.value } )
sum(pvalue.sc_iSA.wilcox>0.05)

sum(pvalue.sc_iSA.wilcox>0.05)

# Apply Benjamini-Hochberg correction (FDR)
p.BHcorr.sc_iSA.wilcox <- p.adjust(pvalue.sc_iSA.wilcox, method = "BH")
sum(p.BHcorr.sc_iSA.wilcox<0.05)

# Apply Benjamini-Hochberg correction (FDR)
p.BHcorr.sc_iSA <- p.adjust(pvalue.sc_iSA, method = "BH")

# Apply Benjamini-Hochberg correction (FDR)
p.BHcorr.sc_iSA.wilcox <- p.adjust(pvalue.sc_iSA.wilcox, method = "BH")

#how many features are significantly different
sum(p.BHcorr.sc_iSA<0.05)
sum(pvalue.sc_iSA<0.05)


# Calculate negative log of adjusted p value
p.BH.nlog.sc_iSA <- -log10(p.BHcorr.sc_iSA)
p.BH.nlog.sc_iSA.wilcox <- -log10(p.BHcorr.sc_iSA.wilcox)

#############Now I want to calculate Fold Changes to later on compare volcano plots
### based on t-statistics vs non-parametric wilcox statistic

# Calculate row-wise group means for sum normalized data
sc <- apply(set1.norm[str_detect(plates,"i..SA")], 1, FUN=mean)
iSA <- apply(set1.norm[sc.col], 1, FUN=mean)

# Calculate log2 Fold Change from group means
FC.iSA <- iSA/sc
log2FC.iSA <- log(FC.iSA,2)

# Make new data frame with volcano plot data
volcano.dat1 <- data.frame(sc, iSA, pvalue.sc_iSA, p.BHcorr.sc_iSA, 
                           p.BH.nlog.sc_iSA, FC.iSA, log2FC.iSA,
                           p.BH.nlog.sc_iSA.wilcox)

# Sort in ascending order by BH-p
volcano.dat1 <- volcano.dat1[order(p.BHcorr.sc_iSA.wilcox),]

# make a simple volcano plot
plot(volcano.dat1$log2FC.iSA, volcano.dat1$p.BH.nlog.sc_iSA.wilcox)


# add cutoff lines to show significant variables
abline(v=2, col="red")
abline(v=-2, col="red")
abline(h=2, col="red")

#similarly simple plot for t-test derived p-value ==> regardless the less power of non-parametric test, 
#it detected the same features as t-test 