################################################################################
# 1. Setup
################################################################################
pkgs <- c("devtools", "randomForest", "dplyr", "rjags", "arm", "xtable","caret",
          "ggplot2","knitr")
cran_no <- pkgs[!pkgs %in% installed.packages()[,1]]
for(i in cran_no){
  install.packages(i)
}
gh_pkgs <- c("usepa/LakeTrophicModelling")
gh_no <- gh_pkgs[!basename(gh_pkgs) %in% installed.packages()[,1]]
for(i in gh_no){
  devtools::install_github(i)
}
all_pkgs <- c(basename(gh_pkgs),pkgs)
lapply(all_pkgs, library, character.only = T)

################################################################################
# 2. Data
################################################################################

data(LakeTrophicModelling)

##################################################
#All Variables
#Clean Up Data - Complete Cases
predictors_all <- predictors_all[predictors_all!="DATE_COL"]
all_dat <- data.frame(ltmData[predictors_all],LogCHLA=log10(ltmData$CHLA))
row.names(all_dat)<-ltmData$NLA_ID
all_dat <- all_dat[complete.cases(all_dat),]  

##################################################
#GIS Variables
#Clean Up Data - Complete Cases
gis_dat_NTL <- data.frame(ltmData[predictors_gis],LogNTL=log10(ltmData$NTL))
gis_dat_PTL <- data.frame(ltmData[predictors_gis],LogPTL=log10(ltmData$PTL))

# row.names(gis_dat)<-ltmData$NLA_ID
row.names(gis_dat_NTL)<-ltmData$NLA_ID
row.names(gis_dat_PTL)<-ltmData$NLA_ID

# gis_dat <- gis_dat[complete.cases(gis_dat),]
gis_dat_NTL <- gis_dat_NTL[complete.cases(gis_dat_NTL),]
gis_dat_PTL <- gis_dat_PTL[complete.cases(gis_dat_PTL),]

################################################################################
# 3. Random Forest for Variable Selection
################################################################################

##################################################
#Model 1: All (GIS + WQ) Variables
#Variable Selection

all_vs <- varsel_regression_rf(all_dat$LogCHLA,all_dat[,names(all_dat)!="LogCHLA"]
                               , ntree=5000,prog=T)  

# Used plot to identify approxiamte minimum
all_vs_plot <- varsel_plot(all_vs)
# top 20 variables
all_vars_select <- unlist(all_vs$vars[19])

# Final Model - all variables
all_rf<-randomForest(y=all_dat$LogCHLA,x=all_dat[,all_vars_select]
                     , ntree=5000, importance=TRUE, proximity=TRUE
                     , keep.forest=TRUE,keep.inbag=TRUE)

##################################################
#Model 2: GIS Only Variables
#Variable Selection
gis_vs_NTL <- varsel_regression_rf(gis_dat_NTL$LogNTL
                                   , gis_dat_NTL[,names(gis_dat_NTL)!="LogNTL"]
                                   , ntree=5000,prog=T)
gis_vs_PTL <- varsel_regression_rf(gis_dat_PTL$LogPTL
                                   , gis_dat_PTL[,names(gis_dat_PTL)!="LogPTL"]
                                   , ntree=5000,prog=T)
#Used plot to identify approxiamte minimum
gis_vs_plot_NTL <- varsel_plot(gis_vs_NTL) 
gis_vs_plot_PTL <- varsel_plot(gis_vs_PTL)

gis_vars_select_NTL <- unlist(gis_vs_NTL$vars[9]) #top 10 variables
gis_vars_select_PTL <- unlist(gis_vs_PTL$vars[9]) #top 10 variables

#Final Model - gis variables
gis_rf_NTL<-randomForest(y=gis_dat_NTL$LogNTL
                         , x=gis_dat_NTL[,gis_vars_select_NTL]
                         , ntree=5000,importance=TRUE,proximity=TRUE
                         , keep.forest=TRUE,keep.inbag=TRUE)
gis_rf_PTL<-randomForest(y=gis_dat_PTL$LogPTL
                         , x=gis_dat_PTL[,gis_vars_select_PTL]
                         , ntree=5000,importance=TRUE,proximity=TRUE
                         , keep.forest=TRUE,keep.inbag=TRUE)

################################################################################
# 4. Random Forest for Variable Selection - Evaluation
################################################################################

data_def <- read.csv("data_def.csv", stringsAsFactors = FALSE)
all_imp <- importance(all_rf)
gis_imp_NTL <- importance(gis_rf_NTL)
gis_imp_PTL <- importance(gis_rf_PTL)

var_importance <- varImportance(all_imp,"All variables")
var_importance_NTL <- varImportance(gis_imp_NTL,"GIS variables")
var_importance_PTL <- varImportance(gis_imp_PTL,"GIS variables")

dplyr::arrange(var_importance,desc(mean_decrease_acc))
dplyr::arrange(var_importance_PTL,desc(mean_decrease_acc))
dplyr::arrange(var_importance_NTL,desc(mean_decrease_acc))

importancePlot(all_rf, data_def=data_def,type='acc',size=3)
importancePlot(gis_rf_NTL, data_def=data_def,type='acc',size=3)
importancePlot(gis_rf_PTL, data_def=data_def,type='acc',size=3)


dplyr::arrange(var_importance,desc(mean_decrease_gini))
dplyr::arrange(var_importance_PTL,desc(mean_decrease_gini))
dplyr::arrange(var_importance_NTL,desc(mean_decrease_gini))

importancePlot(all_rf, data_def=data_def,type='gini',size=3)
importancePlot(gis_rf_NTL, data_def=data_def,type='gini',size=3)
importancePlot(gis_rf_PTL, data_def=data_def,type='gini',size=3)

##################################################
#All variables partial dependence
co <- 100
all_rf_turb_pd <- partialPlot(all_rf,all_dat,"TURB",n.pt=co,plot=FALSE)
all_rf_ntl_pd <- partialPlot(all_rf,all_dat,"NTL",n.pt=co,plot=FALSE)
all_rf_ptl_pd <- partialPlot(all_rf,all_dat,"PTL",n.pt=co,plot=FALSE)
all_rf_elev_pd <- partialPlot(all_rf,all_dat,"ELEV_PT",n.pt=co,plot=FALSE)
all_rf_AlberY_pd <- partialPlot(all_rf,all_dat,"AlbersY",n.pt=co,plot=FALSE)
partial_plot(all_rf_turb_pd, x="Turbidity (NTU)"
             , y=expression(paste('Log10 Chl ', italic("a"),' (',mu,'g/L)')))
partial_plot(all_rf_ntl_pd, x=expression(paste('Total Nitrogen',' (',mu,'g/L)'))
             , y=expression(paste('Log10 Chl ', italic("a"),' (',mu,'g/L)')))
partial_plot(all_rf_ptl_pd, x=expression(paste('Total Phosporus',' (',mu,'g/L)'))
             , y=expression(paste('Log10 Chl ', italic("a"),' (',mu,'g/L)')))
partial_plot(all_rf_elev_pd, x="Elevation (m)"
             , y=expression(paste('Log10 Chl ', italic("a"),' (',mu,'g/L)')))
partial_plot(all_rf_AlberY_pd, x="AlbersY"
             , y=expression(paste('Log10 Chl ', italic("a"),' (',mu,'g/L)')))
gis_rf_N_Lat_pd <- partialPlot(gis_rf_NTL,gis_dat_NTL,"AlbersY",n.pt=co,plot=FALSE)
gis_rf_N_Evergreen_pd <- partialPlot(gis_rf_NTL,gis_dat_NTL,"EvergreenPer_3000m",n.pt=co,plot=FALSE)
gis_rf_N_Eco_pd <- partialPlot(gis_rf_NTL,gis_dat_NTL,"WSA_ECO9",n.pt=co,plot=FALSE)
gis_rf_N_Eco_pd$x <- factor(gis_rf_N_Eco_pd$x,levels=c("NAP","CPL","SAP","TPL","UMW","SPL","NPL","WMT","XER"),ordered=TRUE)
partial_plot(gis_rf_N_Eco_pd, x="Ecoregion", 
             y=expression(paste('Log10 N ',' (',mu,'g/L)')))
partial_plot(gis_rf_N_Lat_pd,
             x="Latitude",
             y=expression(paste('Log10 NTL ',' (',mu,'g/L)')))
partial_plot(gis_rf_N_Evergreen_pd,
             x="%Evergreen",
             y=expression(paste('Log10 NTL ',' (',mu,'g/L)')))
gis_rf_P_Lat_pd<- partialPlot(gis_rf_PTL,gis_dat_PTL,"AlbersY",n.pt=co,plot=FALSE)
gis_rf_P_Evergreen_pd<- partialPlot(gis_rf_PTL,gis_dat_PTL,"EvergreenPer_3000m",n.pt=co,plot=FALSE)
gis_rf_P_Eco_pd<- partialPlot(gis_rf_NTL,gis_dat_PTL,"WSA_ECO9",n.pt=co,plot=FALSE)
gis_rf_P_Eco_pd$x <- factor(gis_rf_P_Eco_pd$x,levels=c("NAP","CPL","SAP","TPL","UMW","SPL","NPL","WMT","XER"),ordered=TRUE)
partial_plot(gis_rf_P_Eco_pd, x="Ecoregion", 
             y=expression(paste('Log10 P',' (',mu,'g/L)')))
partial_plot(gis_rf_P_Lat_pd,
             x="Latitude",
             y=expression(paste('Log10 PTL ',' (',mu,'g/L)')))
partial_plot(gis_rf_P_Evergreen_pd,
             x="%Evergreen",
             y=expression(paste('Log10 PTL ',' (',mu,'g/L)')))

################################################################################
# 5.TSI POLR model
################################################################################

#################################################
# Functions
expected <- function(x, c1.5, c2.5, c3.5, sigma){
  p1.5 <- invlogit((x-c1.5)/sigma)
  p2.5 <- invlogit((x-c2.5)/sigma)
  p3.5 <- invlogit((x-c3.5)/sigma)
  return((1*(1-p1.5)+2*(p1.5-p2.5)+3*(p2.5-p3.5)+4*p3.5))
}

# for plotting logistic regression model
jitter.binary <- function(a, jitt=.05, up=1){
  up*(a + (1-2*a)*runif(length(a),0,jitt))
}

logit <- function(x) return(log(x/(1-x)))

invlogit <- function(x) return(1/(1+exp(-x)))

#################################################
# ConsistentTrophic State Classification using 
# PTL, NTL, and SDD

NLA2007 <- ltmData
Consistent <- ifelse((NLA2007$TS_CHLA_4==NLA2007$TS_PTL 
                      & NLA2007$TS_CHLA_4==NLA2007$TS_NTL), 1, 0)
NLA2007$Consistent <- Consistent

set.seed(100)
Sample <- sample(nrow(NLA2007[NLA2007$Consistent==0,]),size= round(0.1*dim(NLA2007[NLA2007$Consistent==0,])[1]),replace=FALSE)

Evaluation <- NLA2007[Sample,]
Model <- NLA2007[-Sample,]
#################################################
TSI.polrAllVar <- bayespolr(factor(Model[,"TS_CHLA_4"])
                            ~ log(Model[,'SECMEAN']) 
                            + log(Model[,"NTL"])
                            + log(Model[,"PTL"])
                            + Model[,"ELEV_PT"])
xtable(summary(TSI.polrAllVar)$coefficients[,1:2])

#################################################
# TSI.polrAllVar
extractAIC(TSI.polrAllVar)

beta_AllVar <- coef(TSI.polrAllVar)
kappa_AllVar <- TSI.polrAllVar$zeta

c1.5_AllVar <- kappa_AllVar[1]
c2.5_AllVar <- kappa_AllVar[2]
c3.5_AllVar <- kappa_AllVar[3]
sigma_AllVar <- 1/1

# Figure
# Graphical presentation of the POLR model.
# The x-axis is the trophic state index, the y-axis is each lake's trophic state, vertical lines show estimated cutpoints, and curve shows expected trophic state as estimated using ordered logistic regression.
par(mar=c(3,3,0.25,0.25), mgp=c(1.5,0.25,0), tck=-0.005)
plot(0, 0, xlim=c(-10,15), ylim=c(1,4), xlab="Trophic State Index", ylab="",
     type="n", axes=F)
axis(1)
axis(2, at=1:4, labels=c("Oligo","Meso","Eutro","Hyper"), las=1)
lines(rep(c1.5_AllVar, 2), c(1,2))
lines(rep(c2.5_AllVar, 2), c(2,3))
lines(rep(c3.5_AllVar, 2), c(3,4))
curve(expected(x, c1.5_AllVar, c2.5_AllVar, c3.5_AllVar, sigma_AllVar), add=TRUE)
with(NLA2007,
     points(cbind(log(Model$SECMEAN), log(Model$NTL), log(Model$PTL), Model$ELEV_PT)%*%beta_AllVar,
            jitter.binary(as.numeric(ordered(Model$TS_CHLA_4))), col="cyan4"))

#################################################
summary(TSI.polrAllVar)
beta_AllVar <- coef(TSI.polrAllVar)
kappa_AllVar <- TSI.polrAllVar$zeta

c1.5_AllVar <- kappa_AllVar[1]
c2.5_AllVar <- kappa_AllVar[2]
c3.5_AllVar <- kappa_AllVar[3]
sigma_AllVar <- 1/1

X <- cbind(log(Model$SECMEAN), log(Model$NTL), log(Model$PTL), Model$ELEV_PT)
TSI <- X%*% beta_AllVar
TSI <- TSI[!is.na(TSI)]
# se of kappas
se.c <- summary(TSI.polrAllVar)$coef[5:7,2] 
Ibcg <- seq(range(TSI)[1],range(TSI)[2], length.out = 50)
pA <- invlogit(kappa_AllVar[1] - Ibcg)
pB <- invlogit(kappa_AllVar[2] - Ibcg) -  invlogit(kappa_AllVar[1] - Ibcg)
pC <- invlogit(kappa_AllVar[3] - Ibcg) -  invlogit(kappa_AllVar[2] - Ibcg)
pNA <- 1.0 - invlogit(kappa_AllVar[3] - Ibcg)

# Figure
# Graphical presentation of the POLR model. 
# The x-axis is the trophic state index, the y-axis is the probability of being classified into one of the 4 trophic state classes, and the vertical lines and blue bars are the cutpoints $\pm$ one standard error.
par(mar=c(3,3,2,0.25), mgp=c(1.5,0.5,0), tck=-0.01)
plot(range(Ibcg), c(0,1), type="n",     xlab="Tropic State Index", ylab="Prob")
polygon(x=c(c1.5_AllVar-se.c[1], c1.5_AllVar+se.c[1], c1.5_AllVar+se.c[1],c1.5_AllVar-se.c[1]),
        y=c(0,0,1,1), col="cyan4", density=-1, border=NA)
polygon(x=c(c2.5_AllVar-se.c[2], c2.5_AllVar+se.c[2], c2.5_AllVar+se.c[2],c2.5_AllVar-se.c[2]),
        y=c(0,0,1,1), col="cyan4", density=-1, border=NA)
polygon(x=c(c3.5_AllVar-se.c[3], c3.5_AllVar+se.c[3], c3.5_AllVar+se.c[3],c3.5_AllVar-se.c[3]),
        y=c(0,0,1,1), col="cyan4", density=-1, border=NA)
segments(x0=c(c1.5_AllVar,c2.5_AllVar,c3.5_AllVar), y0=rep(0,3),
         x1=c(c1.5_AllVar,c2.5_AllVar,c3.5_AllVar), y1=rep(1,3),col=grey(0.3))
axis(3, at=c(c1.5_AllVar,c2.5_AllVar,c3.5_AllVar), labels=c("Oligo|Meso","Meso|Eu","Eu|Hyper"))
lines(Ibcg, pA)
lines(Ibcg, pB, lty=2)
lines(Ibcg, pC, lty=3)
lines(Ibcg, pNA, lty=4)
legend(-6, 0.5, legend=c("Oligo", "Meso","Eu", "Hyper"),
       lty=1:4, cex=0.75, bty="n")

################################################################################
# 6. POLR Evaluation
################################################################################
predict.EvaluationAllVar <- (cbind(log(Evaluation[,"SECMEAN"])
                                   , log(Evaluation[,"NTL"])
                                   , log(Evaluation[,"PTL"])
                                   , Evaluation[,"ELEV_PT"]))%*%beta_AllVar

predict.EvaluationAllVar <- predict.EvaluationAllVar[!is.na(predict.EvaluationAllVar)]

Predict.CatAllVar <-  vector(length = length(predict.EvaluationAllVar))

for (i in 1:length(predict.EvaluationAllVar)){
  if (predict.EvaluationAllVar[i]< kappa_AllVar[1]) Predict.CatAllVar[i] <- "Oligo"
  if (predict.EvaluationAllVar[i]< kappa_AllVar[2] && predict.EvaluationAllVar[i]> kappa_AllVar[1]) Predict.CatAllVar[i] <- "Meso"
  if (predict.EvaluationAllVar[i]< kappa_AllVar[3] && predict.EvaluationAllVar[i]> kappa_AllVar[2]) Predict.CatAllVar[i] <- "Eu"
  if (predict.EvaluationAllVar[i]> kappa_AllVar[3]) Predict.CatAllVar[i] <- "Hyper"
}

Pred.CatAllVar <- factor(Predict.CatAllVar, levels=c("Oligo", "Meso", "Eu", "Hyper"), ordered=TRUE)

True.CatAllVar <- Evaluation[, "TS_CHLA_4"]

True.CatAllVar <- True.CatAllVar[!is.na(log(Evaluation$SECMEAN))]

CM.AllVar <- confusionMatrix(Pred.CatAllVar, True.CatAllVar)
CM.AllVar
CM.AllVar$byClass[,"Balanced Accuracy"] 
CM.AllVar$overall["Accuracy"]
CMTable <- CM.AllVar$table
xtable(CMTable)

################################################################################
# 7. JAGS Model
################################################################################
set.seed(100)

#set up the initializations 
cutpt.inits <- array(dim= c(3))

for (k in 1:3){
  cutpt.inits[k] <- rnorm(1)
}

inits <- function () {list("cutpt_raw" = cutpt.inits)}

TSLogit_Data <- NLA2007[!is.na(NLA2007[,"SECMEAN"]) & !is.na(NLA2007[,"EvergreenPer_3000m"]),]

# Removing the missing values:
# TSLogit_Data <- Model[!is.na(Model[,"SECMEAN"]) & !is.na(Model[,"EvergreenPer_3000m"]),]
# replacing 0's to avoid Inf when logit transformed setting epsilon to half of the smallest non-zero value and replacing all 0 values with epsilon and all 1 values with 1-epsilon. Then apply the logit transformation.
# EvergreenPer_3000m <- TSLogit_Data[,"EvergreenPer_3000m"]/100
# min( TSLogit_Data[,"EvergreenPer_3000m"][TSLogit_Data[,"EvergreenPer_3000m"]!=min(TSLogit_Data[,"EvergreenPer_3000m"])] )/100
# EvergreenPer_3000m[EvergreenPer_3000m ==0] <-0.0001/2
# logit(EvergreenPer_3000m)

SDD.C <- as.numeric(scale(log(TSLogit_Data$SECMEAN), center = TRUE, scale = TRUE))
NTL.C <- as.numeric(scale(log(TSLogit_Data$NTL), center = TRUE, scale = TRUE))
PTL.C <- as.numeric(scale(log(TSLogit_Data$PTL), center = TRUE, scale = TRUE))
Latitude.C <- as.numeric(scale((TSLogit_Data$AlbersY), center = TRUE, scale = TRUE))
Elevation.C <- as.numeric(scale((TSLogit_Data$ELEV_PT), center = TRUE, scale = TRUE))

DataList = list('TS' = factor(TSLogit_Data[,"TS_CHLA_4"])
                ,'SD' = SDD.C
                ,'Nitrogen' = NTL.C
                ,'Phosphorus' = PTL.C
                ,'Elevation' = Elevation.C
                ,'Evergreen' =  TSLogit_Data[,"EvergreenPer_3000m"]/100
                ,'Eco_Region' = factor(TSLogit_Data[,"WSA_ECO9"])
                ,'Latitude' = Latitude.C)


#The parameter(s) to be monitored
parameters = c('alpha_SD', 'alpha_N', 'alpha_P', 'alpha_E'
               , 'beta_L','beta_E', 'beta_ER'
               , 'gamma_L','gamma_E', 'gamma_ER'
               , 's'
               , 'C')

# Number of steps to "tune" the samplers.
adaptSteps = 3000          

# Number of steps to "burn-in" the samplers.
#changes from 1000 the initail setting
burnInSteps = 5000    

# Number of chains to run.       
nChains = 1

# Total number of steps in chains to save.     
numSavedSteps=100000      

# Number of steps to "thin" (1=keep every step).
thinSteps= 5

# Steps per chain.
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

# Start the clock!
ptm <- proc.time()
JAGS.TSLogit <- jags.model('TSLogit.R',data = DataList 
                        , inits, n.chains = nChains, n.adapt = adaptSteps)
# Stop the clock
proc.time() - ptm
################################################################################
#8. JAGS Model Diagnostics
################################################################################
# Start the clock!
ptm <- proc.time()
Coda.TS.N <- coda.samples(JAGS.TSLogit, parameters, n.iter=100000)
# Stop the clock
proc.time() - ptm
#################################################
plot(Coda.TS.N[,1:3])
plot(Coda.TS.N[,4:6])
plot(Coda.TS.N[,7:9])
plot(Coda.TS.N[,10:12])
plot(Coda.TS.N[,13:15])
plot(Coda.TS.N[,16:18])
plot(Coda.TS.N[,19:21])
plot(Coda.TS.N[,22:24])
plot(Coda.TS.N[,24:27])
plot(Coda.TS.N[,29:30])

print(xtable(cbind(summary(Coda.TS.N)$quantiles, summary(Coda.TS.N)$statistics[,2])), floating=FALSE)
################################################################################
#9. JAGS Model Evaluation
################################################################################

#################################################
# 3 Chains Combined
simCodaOne.Coda.TS.N <- NULL
for (i in 1:1) simCodaOne.Coda.TS.N <- rbind(simCodaOne.Coda.TS.N, Coda.TS.N[[i]])

MCMC.TS.N <- as.mcmc(simCodaOne.Coda.TS.N)

print(summary(simCodaOne.Coda.TS.N[,1:30]))
colnames(MCMC.TS.N)
dim(MCMC.TS.N)


Coeff.TS.N.Summary <- matrix(NA, 30, 4)
for (i in 1:30){ Coeff.TS.N.Summary[i,] <- cbind(mean(simCodaOne.Coda.TS.N[,i])
                                                 , sd(simCodaOne.Coda.TS.N[,i])
                                                 , quantile(simCodaOne.Coda.TS.N[,i], c(0.025), type = 1)
                                                 , quantile(simCodaOne.Coda.TS.N[,i], c(0.975), type = 1))}
colnames(Coeff.TS.N.Summary) <- cbind("mean", "sd", "2.5%", "97.5%")

rownames(Coeff.TS.N.Summary ) <-colnames(MCMC.TS.N)
print(xtable(Coeff.TS.N.Summary, floating=FALSE))



Alpha <- rbind(Coeff.TS.N.Summary["alpha_SD",],  Coeff.TS.N.Summary["alpha_E",],  Coeff.TS.N.Summary["alpha_N",],  Coeff.TS.N.Summary["alpha_P",])

Eval.SDD.C <- as.numeric(scale(log(Evaluation$SECMEAN), center = TRUE, scale = TRUE))
Eval.NTL.C <- as.numeric(scale(log(Evaluation$NTL), center = TRUE, scale = TRUE))
Eval.PTL.C <- as.numeric(scale(log(Evaluation$PTL), center = TRUE, scale = TRUE))
Eval.Latitude.C <- as.numeric(scale((Evaluation$AlbersY), center = TRUE, scale = TRUE))
Eval.Elevation.C <- as.numeric(scale((Evaluation$ELEV_PT), center = TRUE, scale = TRUE))

predict.EvaluationAll <- (cbind(Eval.SDD.C, Eval.Elevation.C, Eval.NTL.C, Eval.PTL.C)) %*% Alpha[,"mean"]
predict.EvaluationAll <- predict.EvaluationAll[!is.na(predict.EvaluationAll)]

Predict.CatAll <-  vector(length = length(predict.EvaluationAll))
C <- rbind(Coeff.TS.N.Summary["C[1]",],  Coeff.TS.N.Summary["C[2]",],  Coeff.TS.N.Summary["C[3]",])

for (i in 1:length(predict.EvaluationAll)){
  if (predict.EvaluationAll[i]< C[1]) Predict.CatAll[i] <- "Oligo"
  if (predict.EvaluationAll[i]< C[2] && predict.EvaluationAll[i]> C[1]) Predict.CatAll[i] <- "Meso"
  if (predict.EvaluationAll[i]< C[3] && predict.EvaluationAll[i]> C[2]) Predict.CatAll[i] <- "Eu"
  if (predict.EvaluationAll[i]> C[3]) Predict.CatAll[i] <- "Hyper"
}

Pred.CatAll <- factor(Predict.CatAll, levels=c("Oligo", "Meso", "Eu", "Hyper"), ordered=TRUE)

True.CatAll <- Evaluation[, "TS_CHLA_4"]

True.CatAll <- True.CatAll[!is.na(log(Evaluation$SECMEAN))]

CM.TS.Multilevel <- confusionMatrix(Pred.CatAll, True.CatAll)
CM.TS.Multilevel
xtable(CM.TS.Multilevel$table)
CM.TS.Multilevel$overall["Accuracy"]
CM.TS.Multilevel$byClass[,"Balanced Accuracy"]

