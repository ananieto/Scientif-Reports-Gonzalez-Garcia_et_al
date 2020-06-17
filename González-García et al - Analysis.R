
#Bioconductor

source("https://bioconductor.org/biocLite.R")
biocLite()

#Packages:

BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")
BiocManager::install("geneplotter")
BiocManager::install("sva")
BiocManager::install("limma")
install.packages("readxl")
install.packages("scatterplot3d")
install.packages("NMF")
install.packages("MASS")
install.packages("ggplot2")
if (!require(devtools)) install.packages("devtools")
library(devtools)
install_github("villardon/MultBiplotR", dependencies = TRUE)


library(limma)
library(affy)
library(AnnotationDbi)
library(hgu133plus2.db)
library(geneplotter)
library(sva)
library(readxl)
library(scatterplot3d)
library(NMF)
library(MultBiplotR)
library(rCUR)
library(MASS)
library(ggplot2)


############################################################
###############     1 - GEO data reading...      ###########

Data268<-ReadAffy()#Set the directory path containing CEL files
Samples_data <- read_excel("Samples data.xlsx")#Read Samples information

#Pre-normalization
hist(Data268)
boxplot(Data268, col=Samples_data$SerieCod)

Dataeset268 <- rma(Data268) #RMA normalization

data268 <- exprs(Dataeset268) #Gene expression matrix

#Post-normalization
boxplot(data268, col=Samples_data$SerieCod)
multidensity(data268, legend=F)


data268<-data268[-c(54614:54675),]#Affymetrix control probes

ind<-scan() #Selecting genes with Annotation information
#colnames(data268)<-scan(what=character())
data.44723<-data268[ind,] 


#PCA PRE-COMBAT
pca_pre<-prcomp(t(data.44723))

colors<-Samples_data$SerieCod
colors[which(colors==10)]<-"brown1"
colors[which(colors==1)]<-"cornflowerblue"
colors[which(colors==2)]<-"navy"
colors[which(colors==3)]<-"magenta"
colors[which(colors==4)]<-"orange"
colors[which(colors==5)]<-"lightseagreen"
colors[which(colors==6)]<-"black"
colors[which(colors==7)]<-"mediumpurple1"
colors[which(colors==8)]<-"limegreen"
colors[which(colors==9)]<-"hotpink"

shape<-Samples_data$SerieCod
shape[which(shape==10)]<-16
shape[which(shape==1)]<-16
shape[which(shape==2)]<-16
shape[which(shape==3)]<-17
shape[which(shape==4)]<-17
shape[which(shape==5)]<-17
shape[which(shape==6)]<-16
shape[which(shape==7)]<-17
shape[which(shape==8)]<-16
shape[which(shape==9)]<-17

scatterplot3d(pca_pre$x[,1:3],pch=shape,color = colors, col.grid = "gainsboro", lty.grid=par("lty"), 
              angle=70)

#COMBAT series normalization

data<-t(data.44723)
batch = Samples_data$Serie

#Since we are not adjusting for any other variables in this analysis, only an intercept is included in
#the model.

modcombat = model.matrix(~1,data=as.data.frame(data))
dim(modcombat)
combat_edata = ComBat(dat=t(data), batch=batch, mod=modcombat)

data268.no_batch<-combat_edata

#PCA POST-COMBAT
pca_pos<-prcomp(t(data268.no_batch))
scatterplot3d(pca_pos$x[,1:3],pch=shape,color = colors, col.grid = "gainsboro", lty.grid=par("lty"), 
              angle=60)
legend("right", legend = c("GSE43289","GSE29796","GSE43378","GSE45921","GSE62802","GSE15824","GSE19728","GSE2817","GSE33331","GSE4290"),
       col =  c("cornflowerblue","navy","magenta","orange","lightseagreen","black","mediumpurple1","limegreen","hotpink","brown1"), 
       pch = c(16,16,17,17,17,16,17,16,17,16))


#Discovery cohort
disc.data<-data268.no_batch[,which(Samples_data$Database=="Discovery")]
disc_info<-Samples_data[which(Samples_data$Database=="Discovery"),]

#Validation cohort
val.data<-data268.no_batch[,which(Samples_data$Database=="Validation")]
val_info<-Samples_data[which(Samples_data$Database=="Validation"),]

#Global cohort
global.data<-data268.no_batch

#No gender - WHO grade differences: 
table(Samples_data$Gender, Samples_data$Grado)
fisher.test(table(Samples_data$Gender, Samples_data$Database), alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95)

#No gender differences between discovery and validation cohort
table(Samples_data$Gender, Samples_data$Database)
fisher.test(table(Samples_data$Gender, Samples_data$Database), alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95)

#No age differences between discovery and validation cohort
library(nortest)
lillie.test(Samples_data$Age[which(Samples_data$Database=="Discovery")]) #Non-normal distribution
lillie.test(Samples_data$Age[which(Samples_data$Database=="Validation")]) #Normal distribution

wilcox.test(Samples_data$Age ~ Samples_data$Database, paired=FALSE) #Mann-Whitney test

######################################################################################
###########################  DISCOVERY    ############################################

#Number of components function for CUR decomposition (Mahoney & Drineas, 2009)
number_factors<-function(X)
{
  sv = svd(X)
  cs = cumsum(sv$d)
  k=length(cs[cs < cs[length(cs)] * 1])
  if(k==0)
  {
    k=2
  }
  if(k==1)
  {
    k=2
  }
  return(k)
}

#CUR decomposition

t.data<-t(disc.data)
cur1<-CUR(t.data,k=number_factors(t.data),c=(dim(t.data)[2]-1),method = "top.scores",weighted = FALSE)
plotLeverage(cur1,top.n=1000)

e.1000<-t(t.data[,cur1@C.index[1:1000]]) #Expression data of 1000 genes seleted by CUR

#Fold Change values: II vs IV grade comparison

e1.2<-e.1000[,c(which(disc_info$Grado==2),which(disc_info$Grado==4))]
f<-as.factor(c(rep(1,19),rep(2,108))) # 1=DA, 2=GBM

labels<-c("Healthy", "affected")
design<- model.matrix(~ 0 + f)
colnames(design)<-c("Healthy","affected")
fit<- lmFit(e1.2, design)
contrast.matrix<-makeContrasts(affected-Healthy,levels = design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)
q<-dim(e1.2)[1]
topTable1<-topTable(fit2,number=q,p=0.05, adjust="fdr")

#Fold Change values: III vs IV grade comparison

e1.2<-e.1000[,c(which(disc_info$Grado==3),which(disc_info$Grado==4))]
f<-c(rep(1,28),rep(2,108)) #1=AA, 2=GBM
f<-as.factor(f)
labels<-c("Healthy", "affected")
design<- model.matrix(~ 0 + f)
colnames(design)<-c("Healthy","affected")

fit<- lmFit(e1.2, design)
contrast.matrix<-makeContrasts(affected-Healthy,levels = design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)
q<-dim(e1.2)[1]
topTable1<-topTable(fit2,number=q,p=0.05, adjust="fdr")

#Fold Change values: II vs III grade comparison
e1.2<-e.1000[,c(which(disc_info$Grado==2),which(disc_info$Grado==3))]
f<-as.factor(c(rep(1,19),rep(2,28))) #1=DA, 2=AA
labels<-c("Healthy", "affected")
design<- model.matrix(~ 0 + f)
colnames(design)<-c("Healthy","affected")
fit<- lmFit(e1.2, design)
contrast.matrix<-makeContrasts(affected-Healthy,levels = design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)
q<-dim(e1.2)[1]
topTable1<-topTable(fit2,number=q,p=0.05, adjust="fdr")


#NMF
data.nmf<-read.table("clipboard", header=T)#Expression matrix 27 genes; 155 samples
res5<-nmf(data.nmf,2,seed="nndsvd")
w5<-res5@fit@W
h5<-res5@fit@H

basismap(res5)
coefmap(res5,main="Scores")

barplot(w5[,1],main="Factor 1",col="black")
barplot(w5[,2],main="Factor 2",col="black")

colors.ord<-c(rep("green",19),rep("blue",28),rep("red",108))
grade.ord<-c(rep(2,19),rep(3,28),rep(4,108))
plot(h5[1,],type="p",col=colors.ord,ylab="F1", pch=19)
plot(h5[2,],type="p",col=colors.ord,ylab="F2", pch=19)

#LDA


lda.disc <- lda(formula = grade.ord ~ ., 
                  data = as.data.frame(t(data.nmf)))

prop.lda = lda.disc$svd^2/sum(lda.disc$svd^2)
prop.lda

plda.disc <- predict(object = lda.disc,
                       newdata = as.data.frame(t(data.nmf)))

table(orig=grade.ord,lda=plda.disc$class)

dataset = data.frame(species = as.factor(grade.global),
                     lda = plda.disc$x)

#Canonical Biplot

res.biplot<-CanonicalBiplot(t(data.nmf), as.factor(grade.ord), InitialTransform = 5)  #3

plot(res.biplot, PlotGroups = F, PlotVars = TRUE, PlotInd = TRUE, LabelInd =
       F, CexGroup = 1, PchGroup = 16, pch=19,
     AddLegend = T, ShowAxes = T, LabelAxes = TRUE,
     LabelGroups = TRUE, PlotCircle = F, ConvexHulls =
       FALSE, LegendPos = "topright", ColorInd = NULL,
     voronoi = TRUE, mode = "a", TypeScale = "Complete",
     ValuesScale = "Original", MinQualityVars = 0, dpg = 0,
     dpi = 0, PredPoints = 0, PlotAxis = T, CexInd =
       NULL, CexVar = NULL, PchInd = NULL, PchVar = NULL,
     ColorVar = NULL, ShowAxis = T, VoronoiColor =
       "black", ShowBox=F)

######################### VALIDATION #################################################
######################################################################################

val.fc4<-val.data[row.names(data.fc4),c(which(val_info$Grado==2),which(val_info$Grado==3),which(val_info$Grado==4))]

val.nmf<-read.table("clipboard", header=T)
grade.val<-c(rep(2,14),rep(3,24),rep(4,75))
      
        #NMF
        res.val<-nmf(val.nmf,2,seed="nndsvd")
        w.val<-res.val@fit@W
        h.val<-res.val@fit@H
        
        basismap(res.val)
        coefmap(res.val,main="Scores",cexCol=5,cexRow=3)

        colors.val<-c(rep("green",14),rep("blue",24),rep("red",75))
        plot(h.val[1,],type="p",col=colors.val,ylab="F1", pch=19)
        plot(h.val[2,],type="p",col=colors.val,ylab="F2", pch=19)

        
        #LDA

        lda.val <- lda(formula = grade.val ~ ., 
                   data = as.data.frame(t(val.nmf)))
        
        prop.lda = lda.val$svd^2/sum(lda.val$svd^2)
        prop.lda
        
        plda <- predict(object = lda.val,
                        newdata = as.data.frame(t(val.nmf)))
        
        table(orig=grade.val,lda=plda$class)
        
        dataset = data.frame(species = as.factor(grade.val),
                             lda = plda$x)
        
        theme_set(theme_bw())
        p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = species, shape = species), size = 2.5) + 
          labs(x = paste("LD1"),
               y = paste("LD2"))
        p1
        
        #BIPLOT
        res.biplot3<-CanonicalBiplot(t(val.nmf), as.factor(grade.val), InitialTransform = 5)  #3
        plot(res.biplot3, PlotGroups = F, PlotVars = TRUE, PlotInd = TRUE, LabelInd =
               F, CexGroup = 1, PchGroup = 16, margin = 0.01,pch=19,
             AddLegend = F, ShowAxes = T, LabelAxes = TRUE,
             LabelGroups = TRUE, PlotCircle = F, ConvexHulls =
               FALSE, LegendPos = "topright", ColorInd = NULL,
             voronoi = TRUE, mode = "a", TypeScale = "Complete",
             ValuesScale = "Original", MinQualityVars = 0, dpg = 0,
             dpi = 0, PredPoints = 0, PlotAxis = T, CexInd =
               NULL, CexVar = NULL, PchInd = NULL, PchVar = NULL,
             ColorVar = NULL, ShowAxis = T, VoronoiColor =
               "black", ShowBox=F)

#GLOBAL

global.nmf<-read.table("clipboard", header=T)
grade.global<-c(rep(2,33),rep(3,52),rep(4,183))

        
        #NMF
        res.global<-nmf(global.nmf,2,seed="nndsvd")
        w.global<-res.global@fit@W
        h.global<-res.global@fit@H

        basismap(res.global)
        coefmap(res.global,main="Scores")
        
        colors.global<-c(rep("green",33),rep("blue",52),rep("red",183))
        plot(h.global[1,],type="p",col=colors.global,ylab="F1", pch=19)
        plot(h.global[2,],type="p",col=colors.global,ylab="F2", pch=19)
        
        
        
       
