install.packages("packMBPLSDA")
install.packages("o2plsda")


library(packMBPLSDA)


data(status)
data(medical)
data(omics)
data(nutrition)
ktabX <- ktab.list.df(list(medical = medical, nutrition = nutrition, omics = omics))
disjonctif <- (disjunctive(status))
dudiY   <- dudi.pca(disjonctif , center = FALSE, scale = FALSE, scannf = FALSE)
modelembplsQ <- mbplsda(dudiY, ktabX, scale = TRUE, option = "uniform", scannf = FALSE, nf = 2)


bloYobs <- 2
resdim <- testdim_mbplsda(object = modelembplsQ, nrepet = 30, threshold = 0.5,
                          bloY = bloYobs, cpus = 1, algo = c("max"), outputs = c("ER"))

ncpopt <- 1

rtsPermut <- permut_mbplsda(modelembplsQ, nrepet = 30, npermut = 100, optdim = ncpopt,
                            outputs = c("ER"), bloY = bloYobs, nbObsPermut = 10, cpus=1, algo = c("max"))
plot_permut_mbplsda(rtsPermut,"plotPermut_nf1_30rep_100perm")


predictions <- pred_mbplsda(modelembplsQ,
                            optdim = ncpopt,
                            threshold = 0.5,
                            bloY=bloYobs,
                            algo = c("max", "gravity", "threshold"))

plot_pred_mbplsda(predictions,"plotPred_nf1", propbestvar=0.20)

CVpred <- cvpred_mbplsda(modelembplsQ, nrepet = 90, threshold = 0.5, bloY = bloYobs,
                         optdim = ncpopt, cpus = 1, algo = c("max"))
plot_cvpred_mbplsda(CVpred,"plotCVPred_nf1_90rep")

resboot <- boot_mbplsda(modelembplsQ, optdim = ncpopt, nrepet = 30, cpus=1)
plot_boot_mbplsda(resboot,"plotBoot_nf1_30rep", propbestvar=0.20)
