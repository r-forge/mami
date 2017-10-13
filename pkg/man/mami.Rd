\name{mami}
\alias{mami}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model averaging (and model selection) after multiple imputation
}
\description{                                       
Performs model selection/averaging on multiply imputed data and combines the resulting estimates. The package also provides access to less frequently used 
model averaging techniques and offers integrated bootstrap estimation.
}
\usage{
mami(X, missing.data=c("imputed","none","CC"),  
  model=c("gaussian","binomial","poisson","cox"), outcome=1, id=NULL, 
  method=c("MA.criterion","MS.criterion","LASSO","LAE","MMA"), 
  criterion=c("AIC","BIC","BIC+","CV","GCV"), kfold=5, cvr=F,
  inference=c("standard","+boot"), B=20, X.org=NULL, CI=0.95,
  var.remove=NULL, aweights=NULL, add.factor=NULL, add.interaction=NULL, 
  add.transformation=NULL, add.stratum=NULL,    
  candidate.models=c("all","restricted","very restricted"), screen=0,
  report.exp=FALSE, print.time=FALSE, print.warnings=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Either a list of multiply imputed datasets (each of them a dataframe), or an object of class `amelia' created by \pkg{Amelia II}, or an object of class `mids' created by \pkg{mice}, or a single dataframe. 
}
  \item{missing.data}{
A character string, typically \code{"imputed"} when multiply imputed data are provided under \code{X}, or \code{"CC"} if a complete case analysis
is desired, or \code{"none"} if there is no missing data.
}
  \item{model}{
A character string specifying the model family, either \code{"gaussian"} for linear regression models, \code{"binomial"} for logistic regression models,
\code{"poisson"} for Poisson regression models, or \code{"coxph"} for Cox's proportional hazards models.
}
  \item{outcome}{
A character vector or integer specifying the variable(s) or column(s) which should be treated as outcome variable(s).
}
  \item{id}{
A character vector or integer specifying the variable or column to be used for a random intercept in the analysis model.
}
  \item{method}{
A character string specifying the model selection or model averaging technique: \code{"MA.criterion"} for model averaging based on exponential AIC or BIC weights, 
\code{"MS.criterion"} for stepwise variable selection based on AIC or BIC; or (all subset) model selection based on (generalized) cross validation, \code{"MMA"} for Mallow's model averaging 
(only linear model, see \code{\link{mma}}), \code{"LASSO"} for model selection with LASSO and additional lasso averaging (see \code{\link{lae}}) if \code{method="LAE"} 
("LAE" not available for Cox model).
}
  \item{criterion}{
A character string specifying the model selection criterion used for criterion-based model selection and averaging; currently either \code{"AIC"}
for Akaike's Information Criterion, \code{"BIC"} for the Bayes criterion of Schwarz, \code{"BIC+"} for the Bayes criterion of Schwarz with quicker model averaging based on the leaps algorithm of the 
\pkg{BMA} package. Model selection can also be done based on the cross validation error (\code{CV}, based on the squared loss function), or \code{GCV} for generalized cross validation.  
}
  \item{kfold}{
An integer specifying \code{kfold} cross validation; to be used when applying shrinkage estimation (\code{method="LASSO"}, \code{method="LAE"}) or criterion \code{CV} is used for model selection.
}
 \item{cvr}{
A logical value specifying whether the subsets used for cross validation should be shuffled or not.
}
  \item{inference}{
A character string, either \code{"standard"} for applying multiple imputation combining rules upon the model selection/averaging estimates or
\code{"+boot"} if additional bootstrap based inference is required. See reference section for more information.
}
  \item{B}{
An integer indicating the number of bootstrap replications to use (when \code{inference} \code{=} \code{"+boot"} is chosen).
}
  \item{X.org}{
A dataframe consisting of the original unimputed data (which needs to be specified when \code{inference} \code{=} \code{"+boot"} is chosen).
}
  \item{CI}{
A scalar greater than 0 and less than 1 specifying the confidence of the confidence interval.
}
  \item{var.remove}{
Either a vector of character strings or integers, specifying the variables or columns which are part of the data but not to be considered 
in the model selection/averaging procedure.
}
  \item{aweights}{
A weight vector that is relevant to the analysis model.
}
  \item{add.factor}{
Either a vector of character strings or integers, specifying the variables which should be treated as categorical/factors in the analysis. Variables which are already defined to be factors in the data are detected automatically and do not necessarily need to be specified with this option.
}
  \item{add.interaction}{
A list of either character strings or integers, specifying the variables which should be added as interactions in the analysis model.
}
  \item{add.transformation}{
A vector of character strings, specifying transformations of variables which should be added to the analysis models.
}
  \item{add.stratum}{
A character vector or integer specifying the variable used as a stratum in Cox regression analysis.
}
  \item{candidate.models}{
A character string specifying whether for criterion based model selection/averaging all possible candidate models should be considered (\code{"all"}),
or only candidate models with a limited amount of variables (\code{"restricted","very restricted"}).
}
 \item{screen}{
Number of variables which should be removed in an initial screening step with LASSO.
}
  \item{report.exp}{
A logical value specifying whether exponentiated coefficients should be reported or not.
}
  \item{print.time}{
A logical value specifying whether analysis time and anticipated estimation time for bootstrap estimation should be printed.
}
  \item{print.warnings}{
A logical value specifying whether warnings and any other output from the function should be printed or not.
}
  \item{\dots}{
Further arguments to be passed, i.e. for functions \code{\link{lae}}; \code{dredge} from the \pkg{MuMIn} package; \code{bic.glm} and \code{bic.surv} from the \pkg{BMA} package; or \code{cv.glmnet} from package \pkg{glmnet}.
}
}
\details{
Model selection/averaging will be performed on each imputed dataset. The results will  be combined according to formulae (7)-(10) in 
Schomaker and Heumann (CSDA, 2014), see \emph{References} below for more details. If \code{inference="+boot"} is chosen, then the procedure
described in Table 1 will be performed in addition to standard MI inference. For longitudinal data (specified via \code{id}) the bootstrap is
based on the subject/person/id level. To obtain insightful results from bootstrap estimation \code{B} should be large, at least \emph{B>200} 
and \code{\link{plot.mami}} may be used. 

Note that a variable will be formally selected if it is selected (by means of either model selection or averaging) in at least one imputed
set of data, but its overall impact will depend on how often it is  chosen. As a result, effects of variables which are not supported 
throughout imputed datasets and candidate models will simply be less pronounced. Variable importance measures based on model averaging weights
are calculated for each imputed dataset and will then be averaged.

If \code{method="MA.criterion"} is chosen and the number of variables is large, then computation time might be a burden and obtaining
results can even become difficult. The package provides several options to decrease computation time, see Section 6.1 of the manual at 
\url{http://mami.r-forge.r-project.org/}.

The function provides access to linear, logistic, Poisson and Cox proportional hazard models; one may add a random intercept to each of these models with the
\code{id} option. Other models are not supported yet. Variables used for the imputation model but not needed for the analysis model can be removed with 
option \code{var.remove}.
}
\value{
Returns an object of \code{class} `mami': 
\item{coefficients.ma}{A matrix of coefficients, standard errors and confidence intervals for model \emph{averaging} estimators.}
\item{coefficients.ma.boot}{A matrix of coefficients and bootstrap results (confidence intervals, mean, standard error) for model \emph{averaging}.}
\item{coefficients.s}{A matrix of coefficients, standard errors and confidence intervals for model \emph{selection} estimators.}
\item{coefficients.boot.s}{A matrix of coefficients and bootstrap results (confidence intervals, mean, standard error) for model \emph{selection}.}
\item{variable.importance}{A vector containing the variable importance for each variable based on model averaging weights.}
\item{boot.results}{A list of detailed estimation results for each bootstrap sample. The first list element refers to the results from model
selection, the second entry the results from model averaging.}
}

\references{
Schomaker, M., Heumann, C. (2014) \emph{Model Selection and Model Averaging after Multiple Imputation},
Computational Statistics & Data Analysis, 71:758-770
}
\author{
Michael Schomaker
}
\seealso{
\code{\link{plot.mami}} to visualize bootstrap results, \code{\link{lae}} and \code{\link{mma}} for model averaging techniques.
}
\examples{
citation("MAMI")

####################################################
# Example 1: Freetrade example from Amelia package #
#            Cross-Section-Time-Series Data        #
#            Linear and linear mixed model         #
####################################################
set.seed(24121980)
library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, ts = "year", cs = "country", noms="signed", 
                        polytime = 2, intercs = TRUE, empri = 2)

# AIC based model averaging and model selection in a linear model after MI 
# (with and without bootstrapping)
mami(freetrade_imp, method="MS.criterion", outcome="tariff",add.factor=c("country"))
mami(freetrade_imp, method="MA.criterion", outcome="tariff",add.factor=c("country"))
m1 <- mami(freetrade_imp, method="MS.criterion", outcome="tariff",add.factor=c("country"),
           inference="+boot",B=25,X.org=freetrade,print.time=TRUE)

m1            # be patient with bootstrapping, increase B>=200 for better results
summary(m1)   # easier to read
plot(m1, plots.p.page="4")

# For comparison: Mallow's model averaging (MMA) and complete case analysis 
mami(freetrade_imp, method="MMA", outcome="tariff",add.factor=c("country"))
mami(freetrade, method="MS.criterion", missing.data="CC", outcome="tariff",
     add.factor=c("country")) #Note the difference to imputed analysis (e.g. usheg)     

# Use linear mixed model with random intercept for country using "id"
mami(freetrade_imp, outcome="tariff", id="country")

####################################################################
# Example 2: HIV treatment data, linear model and Cox model        #
####################################################################

# Impute with Amelia
data(HIV)
HIV_imp <- amelia(HIV, m=5, idvars="patient",noms=c("hospital","sex","dead","tb","cm"),
                  ords=c("period","stage"),logs=c("futime","cd4"),
                  bounds=matrix(c(3,7,9,11,0,0,0,0,3000,5000,200,150),ncol=3,nrow=4))

# i)  Cox PH model
# Model selection (with AIC) to select risk factors for the hazard of death, 
# reported as hazard ratios
# Also: add transformations and interaction terms to candidate models 
mami(HIV_imp, method="MS.criterion",model="cox", outcome=c("futime","dead"), 
     add.factor=c("hospital","stage","period"), add.transformation=c("cd4^2","age^2"), 
     add.interaction=list(c("cd4","age")), report.exp=TRUE, var.remove=c("patient","cd4slope6"))
      
# Similar as above (no hazard ratios reported, no interaction, hospitals as stratum),
# but with boostrap CI and visualization of bootstrap distribution 
m2 <- mami(HIV_imp, method="MS.criterion",model="cox", inference="+boot",
           X.org=HIV, outcome=c("futime","dead"), add.factor=c("stage","period"),
           add.transformation=c("cd4^2","age^2"), add.stratum="hospital", B=25, 
           var.remove=c("patient","cd4slope6"),report.exp=TRUE,
           print.time=TRUE,print.warnings=FALSE) 
summary(m2)
plot(m2)

# ii) Linear model
# Model selection and averaging to identify predictors for immune recovery 6 months
# after starting antiretroviral therapy, presented as CD4 slope which is the average
# change in number of CD4 cells per week (deaths are ignored for this example)

# AIC based model selection (stepAIC) after multiple imputation 
mami(HIV_imp, method="MS.criterion", outcome="cd4slope6", 
     add.factor=c("hospital","stage","period"), var.remove=c("patient","dead","futime"))

# Model averaging (BIC weights + Occams Window) 
m3 <- mami(HIV_imp,outcome="cd4slope6", add.factor=c("hospital","stage","period"),
           criterion="BIC+",var.remove=c("patient","dead","futime"))
summary(m3)

#########################################################################################
# Example 3:   Model selection/averaging with no missing data, using shrinkage          #
# Example from Tibshirani, R. (1996) Regression shrinkage and selection via the lasso,  #
# Journal of the Royal Statistical Society, Series B 58(1), 267-288.                    #
# Useful to use mami to obtain Bootstrap CI after model selection/averaging             #
#########################################################################################

data(Prostate)
mami(Prostate,method="LASSO",missing.data="none",outcome="lpsa" 
     ,kfold=10) # LASSO (selection) based on 10-fold CV
mami(Prostate,method="LAE",missing.data="none",outcome="lpsa" 
     ,kfold=10) # LASSO averaging based on 10-fold CV
m4 <- mami(Prostate,missing.data="none",outcome="lpsa", inference="+boot",
           B=50,print.time=TRUE) # AIC based averaging with Boostrap CI
summary(m4)
plot(m4) # a few bimodal distributions: effect or not?


###################################################
# Example 4: use utilities from other packages    #
###################################################

# Examples (without missing data for simplicity)
# Model Averaging with AIC: use restrictions as done in "dredge"
# Example: Candidate models cannot contain "svi" and "lcp" at the same time
mami(Prostate,missing.data="none",outcome="lpsa", subset = !(svi && lcp))   
# Example: Make Occam's Window smaller when using BMA package
mami(Prostate,missing.data="none",outcome="lpsa",criterion="BIC+", OR=5)  
# Example: Use Elastic Net instead of Lasso (passing on alpha to glmnet())
mami(Prostate,missing.data="none",outcome="lpsa",method="LASSO",alpha=0.5)

}

\keyword{Multiple Imputation}
\keyword{Model Averaging}
\keyword{Model Selection}
\keyword{Bootstrapping}
