library(devtools)
install_github("guokai8/o2plsda")


setwd("C:/Users/theer/OneDrive - Khon Kaen University/Desktop/Demodata")

library(o2plsda)


library(readr)
metabolite_table <- read_csv("16S/metabolite_table.csv")

OTU = read.delim("C:/Users/theer/OneDrive - Khon Kaen University/Desktop/Demodata/16S/otu_table.txt", header = TRUE, sep = "\t", dec = ".")

SampleInfo <- read_csv("16S/SampleInfo.csv")

group = SampleInfo$ClassNote
met_name = metabolite_table$Name
Denovo = OTU$OTU_ID
Sample_name = SampleInfo$MetaID

metabolite_table$Name

X = t(metabolite_table[, 4:33])
colnames(X) = metabolite_table$Name

dim(OTU)
Y = t(OTU[,2:31])
colnames(Y) = OTU$OTU_ID


dim(Y)
dim(X)

#set.seed(123)
# sample * values
#X = matrix(rnorm(5000),50,100)
# sample * values
#Y = matrix(rnorm(5000),50,100)
#rownames(X) <- paste("S",1:50,sep="")
#rownames(Y) <- paste("S",1:50,sep="")
#colnames(X) <- paste("Gene",1:100,sep="")
#colnames(Y) <- paste("Lipid",1:100,sep="")
#X = scale(X, scale=T)
#Y = scale(Y, scale=T)
## group factor could be omitted if you don't have any group
#group <- rep(c("Control","Treat"),each = 25)

X
dim(Y)

set.seed(123)
## nr_folds : cross validation k-fold (suggest 10)
## ncores : parallel paramaters for large datasets
cv <- o2cv(X,Y,1:5,1:3,1:3,group=group,nr_folds = 10)
?o2cv
#####################################
# The best paramaters are nc =  5 , nx =  3 , ny =  3
#####################################
# The Qxy is  0.08222935  and the RMSE is:  2.030108
#####################################

fit <- o2pls(X,Y,5,3,3)
summary(fit)



Xl <- loadings(fit,loading="Xjoint")
Xs <- scores(fit,score="Xjoint")
plot(fit,type="score",var="Xjoint", group=group)
plot(fit,type="loading",var="Xjoint", group=group,repel=F,rotation=TRUE)



res <- oplsda(fit,group, nc=5)
plot(res,type="score", group=group)
vip <- vip(res)
plot(res,type="vip", group = group, repel = FALSE,order=TRUE)

