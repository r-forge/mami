pkgname <- "MAMI"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('MAMI')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("HIV")
### * HIV

flush(stderr()); flush(stdout())

### Name: HIV
### Title: HIV treatment data
### Aliases: HIV
### Keywords: Datasets HIV

### ** Examples

data(HIV)
str(HIV)



cleanEx()
nameEx("lae")
### * lae

flush(stderr()); flush(stdout())

### Name: lae
### Title: Lasso Averaging Estimation
### Aliases: lae
### Keywords: Model Averaging Shrinkage

### ** Examples

library(lasso2)
data(Prostate)
lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)



cleanEx()
nameEx("mami")
### * mami

flush(stderr()); flush(stdout())

### Name: mami
### Title: Model averaging (and model selection) after multiple imputation
### Aliases: mami
### Keywords: Multiple Imputation Model Averaging Model Selection
###   Bootstrapping

### ** Examples

####################################################
# Example 1: Freetrade example from Amelia package #
#            Cross-Section-Time-Series Data        #
#            Linear and linear mixed model         #
####################################################
set.seed(24121980)
library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, ts = "year", cs = "country", noms="signed", polytime = 2, 
                        intercs = TRUE, empri = 2)

# AIC based model averaging and model selection in a linear model after MI 
# (with and without bootstrapping)
mami(freetrade_imp, method="criterion.selection", ycol="tariff",add.factor=c("country"))
mami(freetrade_imp, method="criterion.average", ycol="tariff",add.factor=c("country"))
m1 <- mami(freetrade_imp, method="criterion.selection", ycol="tariff",add.factor=c("country"),
           inference="+bootstrapping",B=25,X.org=freetrade,print.time=TRUE)
m1         # be patient with bootstrapping, increase B>=200 for better results
plot(m1, plots.p.page="4")

# For comparison: Mallow's model averaging (MMA) and complete case analysis 
mami(freetrade_imp, method="MMA", ycol="tariff",add.factor=c("country"))
mami(freetrade, method="criterion.selection", missing.data="CC", ycol="tariff",
     add.factor=c("country")) #Note the difference to imputed analysis  (e.g. usheg)
     

# Use linear mixed model with random intercept for country as alternative, just specify "id" option
# Note: same procedure for random intercepts (frailty) in logistic, Poisson, or Cox Model 
mami(freetrade_imp, ycol="tariff", id="country")

####################################################################
# Example 2: HIV treatment data, linear model and Cox model        #
####################################################################

# Impute with Amelia
data(HIV)
HIV_imp <-  amelia(HIV, m=5, idvars="patient",noms=c("hospital","sex","dead","tb","cm"),
                   ords=c("period","stage"),logs=c("futime","cd4"),
                   bounds=matrix(c(3,7,9,11,0,0,0,0,3000,5000,200,150),ncol=3,nrow=4))

# i)  Cox PH model
# Model selection (with AIC) to select risk factors for the hazard of death, 
# reported as hazard ratios
# Also: add transformations and interaction terms to candidate models 
mami(HIV_imp, method="criterion.selection",model.family="coxph", ycol=c("futime","dead"), 
     add.factor=c("hospital","stage","period"), add.transformation=c("cd4^2","age^2"), 
     add.interaction=list(c("cd4","age")), report.exp=TRUE, var.remove=c("patient","cd4slope6")) 
# Similar as above (= same but no hazard ratios reported, no interaction, hospitals as stratum),
# but with boostrap CI and visualization of bootstrap distribution (be patient...it's worth it)
m2 <- mami(HIV_imp, method="criterion.selection",model.family="coxph", inference="+bootstrapping",
           X.org=HIV, ycol=c("futime","dead"), add.factor=c("stage","period"),
           add.transformation=c("cd4^2","age^2"), use.stratum="hospital", B=25, 
           var.remove=c("patient","cd4slope6"),print.time=TRUE,print.warnings=FALSE) 
m2
plot(m2)

# ii) Linear model
# Model selection and averaging to identify predictors for immune recovery 6 months
# after starting antiretroviral therapy, presented as CD4 slope which is the average
# change in number of CD4 cells per week (deaths are ignored for this example)

# AIC based model selection (stepAIC) after multiple imputation 
mami(HIV_imp, method="criterion.selection", ycol="cd4slope6", 
     add.factor=c("hospital","stage","period"), var.remove=c("patient","dead","futime"))
# Model averaging (AIC weights) for variables typically captured
mami(HIV_imp,ycol="cd4slope6", add.factor=c("hospital","stage","period"),
     var.remove=c("patient","dead","futime","tb","cm","haem"))
# Mallow's model averaging
mami(HIV_imp, method="MMA", ycol="cd4slope6", add.factor=c("hospital","stage","period"),
     var.remove=c("patient","dead","futime"))


#########################################################################################
# Example 3:   Model selection/averaging with no missing data, using shrinkage          #
# Example from Tibshirani, R. (1996) Regression shrinkage and selection via the lasso,  #
# Journal of the Royal Statistical Society, Series B 58(1), 267-288.                    #
# Useful to use mami to obtain Bootstrap CI after model selection/averaging             #
#########################################################################################

library(lasso2)
data(Prostate)
mami(Prostate,method="LASSO/LAE",missing.data="none",ycol="lpsa"
     ,kfold=10) # LASSO (selection/averaging) based on 10-fold CV
mami(Prostate,missing.data="none",ycol="lpsa")  # AIC based averaging
m3 <- mami(Prostate,missing.data="none",ycol="lpsa", inference="+bootstrapping",
           B=50,print.time=TRUE) # with Boostrap CI
m3
plot(m3) # a few bimodal distributions: effect or not?



cleanEx()
nameEx("mma")
### * mma

flush(stderr()); flush(stdout())

### Name: mma
### Title: MMA: Mallow's Model Averaging
### Aliases: mma
### Keywords: Model Averaging

### ** Examples

library(lasso2)
data(Prostate)
mma(Prostate,modelformula=lpsa~.,ycol="lpsa")



cleanEx()
nameEx("plot.lae")
### * plot.lae

flush(stderr()); flush(stdout())

### Name: plot.lae
### Title: Visualize estimates from Lasso Averaging Estimation
### Aliases: plot.lae
### Keywords: Shrinkage Model Averaging

### ** Examples

library(lasso2)
data(Prostate)
l1 <- lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
plot(l1,xaxis="realnumbers")



cleanEx()
nameEx("plot.mami")
### * plot.mami

flush(stderr()); flush(stdout())

### Name: plot.mami
### Title: Plot bootstrap results when utilizing model selection/averaging
###   after multiple imputation
### Aliases: plot.mami
### Keywords: Model Averaging Model Selection Multiple Imputation

### ** Examples

library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, ts = "year", cs = "country", noms="signed", polytime = 2,
                        intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="criterion.selection", ycol="tariff",add.factor=c("country"),
           inference="+bootstrapping", B=25,X.org=freetrade,print.time=TRUE)
m1 # be patient with bootstrapping, increase B>=200 for better results
plot(m1) 



cleanEx()
nameEx("print.lae")
### * print.lae

flush(stderr()); flush(stdout())

### Name: print.lae
### Title: Printing results of Lasso Averaging Estimation
### Aliases: print.lae
### Keywords: Model Averaging Shrinkage

### ** Examples

library(lasso2)
data(Prostate)
l1<-lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
print(l1)



cleanEx()
nameEx("print.mami")
### * print.mami

flush(stderr()); flush(stdout())

### Name: print.mami
### Title: Printing results of model selection and model averaging after
###   multiple imputation
### Aliases: print.mami
### Keywords: Model Averaging Model Selection Multiple Imputation

### ** Examples

library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, p2s=0, ts = "year", cs = "country", noms="signed",
                        polytime = 2, intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="criterion.average", ycol="tariff",add.factor=c("country"))
print(m1)
m1



cleanEx()
nameEx("print.mma")
### * print.mma

flush(stderr()); flush(stdout())

### Name: print.mma
### Title: Printing results of Mallow's Model Averaging
### Aliases: print.mma
### Keywords: Model Averaging

### ** Examples

library(lasso2)
data(Prostate)
m1 <- mma(Prostate,modelformula=lpsa~.,ycol=9)
m1



cleanEx()
nameEx("summary.lae")
### * summary.lae

flush(stderr()); flush(stdout())

### Name: summary.lae
### Title: Summarizes results from Lasso Averaging Estimation
### Aliases: summary.lae
### Keywords: Shrinkage Model Averaging

### ** Examples

library(lasso2)
data(Prostate)
l1<-lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
summary(l1)



cleanEx()
nameEx("summary.mami")
### * summary.mami

flush(stderr()); flush(stdout())

### Name: summary.mami
### Title: Summarizes results of model selection and model averaging after
###   multiple imputation
### Aliases: summary.mami
### Keywords: Model Averaging Model Selection Multiple Imputation

### ** Examples

library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, p2s=0, ts = "year", cs = "country", noms="signed",
                        polytime = 2, intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="criterion.average", ycol="tariff",add.factor=c("country"))
summary(m1)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
