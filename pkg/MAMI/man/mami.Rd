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
mami(X, method = c("criterion.average", "criterion.selection", "MMA", "LASSO/LAE", "WAML"), criterion = c("AIC", "BIC", "BIC+", "CV"),
 B = 20, X.org = NULL, inference = c("standard", "+bootstrapping"), missing.data = c("imputed", "none", "CC"), var.remove = NULL, 
 user.weights = NULL, candidate.models = c("all", "restricted", "very restricted"), model.family = c("gaussian", "binomial", "poisson",
 "coxph"), add.factor = NULL, add.interaction = NULL, add.transformation = NULL, ycol = 1, CI = 0.95, kfold = 5, id = NULL,
 use.stratum = NULL, report.exp = FALSE, print.time = FALSE, print.warnings = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Either a list of multiply imputed datasets (each of them a dataframe), or an object of class `amelia' created by \pkg{Amelia II}, or a single dataframe. 
}
  \item{method}{
A character string specifying the model selection or model averaging technique: \code{"criterion.average"} for model averaging based on exponential AIC or BIC weights, 
\code{"criterion.selection"} for stepwise variable selection based on AIC or BIC, \code{"MMA"} for Mallow's model averaging (only linear model), 
\code{"LASSO/LAE"} for shrinkage estimation and averaging with LASSO (only linear model), and \code{"WALS"} for weighted-average least squares (linear and logistic model, available soon).
}
  \item{criterion}{
A character string specifying the model selection criterion used for criterion-based model selection/averaging; currently either \code{"AIC"}
for Akaike's Information Criterion, \code{"BIC"} for the Bayes criterion of Schwarz, \code{"BIC+"} for the Bayes criterion of Schwarz with quicker model averaging based on the leaps algorithm of the 
\code{BMA} package, or \code{CV} for the cross validation error. \code{CV} can currently only be used 
for GLMs and may be more useful for model selection than model averaging. 
}
  \item{B}{
An integer indicating the number of bootstrap replications to use (when \code{inference} \code{=} \code{"+bootstrapping"} is chosen).
}
  \item{X.org}{
A dataframe consisting of the original unimputed data (which needs to be specified when \code{inference} \code{=} \code{"+bootstrapping"} is chosen).
}
  \item{inference}{
A character string, either \code{"standard"} for applying multiple imputation combining rules upon the model selection/averaging estimates or
\code{"+bootstrapping"} if additional bootstrap based inference is required. See reference section for more information.
}
  \item{missing.data}{
A character string, typically \code{"imputed"} when multiply imputed data are provided under \code{X}, or \code{"CC"} if a complete case analysis
is desired, or \code{"none"} if a there is no missing data.
}
  \item{var.remove}{
Either a vector of character strings or integers, specifying the variables or columns which are part of the data but not to be considered 
in the model selection/averaging procedure.
}
  \item{user.weights}{
A weight vector that is relevant to the analysis model.
}
  \item{candidate.models}{
A character string specifying whether for criterion based model selection/averaging all possible candidate models should be considered (\code{"all"}),
or only candidate models with a limited amount of variables (\code{"restricted","very restricted"}).
}
  \item{model.family}{
A character string specifying the model family, either \code{"gaussian"} for linear regression models, \code{"binomial"} for logistic regression models,
\code{"poisson"} for Poisson regression models, or \code{"coxph"} for Cox's proportional hazards models.

}
  \item{add.factor}{
Either a vector of character strings or integers, specifying the variables which should be treated as categorical/factors in the analysis. Variables which are already defined to be factors in the data are detected automatically and do not necessarily need to be specified with this option.
}
  \item{add.interaction}{
A list of vectors of either character strings or integers, specifying the variables which should be added as interactions in the analysis model.
}
  \item{add.transformation}{
A vector of character strings, specifying transformations of variables which should be added to the analysis models.
}
  \item{ycol}{
A vector or integer specifying the variable(s) or column(s) which should be treated as outcome variable.
}
  \item{CI}{
A value greater than 0 and less than 1 specifying the confidence of the confidence interval.
}
  \item{kfold}{
An integer specifying \code{kfold} cross validation; to be used when applying shrinkage estimation (\code{method="LASSO/LAE"}) or criterion \code{CV}.
}
  \item{id}{
A character vector or integer specifying the variable or column to be used for a random intercept in the analysis model.
}
  \item{use.stratum}{
A character vector or integer specifying the variable used as a stratum in Cox regression analysis.
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
Further arguments to be passed, i.e. for functions \code{\link{lae}}, \code{dredge} from the MuMIn package or \code{bic.glm} and \code{bic.surv} from the \code{BMA} package.
}
}
\details{
Model selection/averaging will be performed on each imputed dataset. The results will  be combined according to formulae (7)-(10) in 
Schomaker and Heumann (CSDA, 2014), see \emph{References} below for more details. If \code{inference="+bootstrapping"} is chosen, then the procedure
described in Table 1 will be performed in addition to standard MI inference. For longitudinal data (specified via \code{id}) the bootstrap is
based on the subject/person/id level. To obtain insightful results from bootstrap estimation \code{B} should be large, at least \emph{B>200} 
and \code{\link{plot.mami}} may be used. 

Note that a variable will be formally selected if it is selected (by means of either model selection or averaging) in at least one imputed
set of data, but its overall impact will depend on how often it is  chosen. As a result, effects of variables which are not supported 
throughout imputed datasets and candidate models will simply be less pronounced. Variable importance measures based on model averaging weights
are calculated for each imputed dataset and will then be averaged.

If \code{method="criterion.average"} is chosen and the number of variables is large, then computation time might be a burden and obtaining
results can even become unfeasible. The reason for this is that for criterion based model \emph{averaging} the implementation of package \pkg{MuMIn} 
is used, which considers all possible candidate models, that is \emph{2^p} different candidate models for p parameters to estimate. If it is clear that
only a subset of variables are relevant then the options \code{candidate.models="restricted/very restricted"} may be useful which 
essentially specifies that only up to a half/fourth of the provided variables can be added to a single candidate model. However, this option
should be used with caution. Alternatively \code{criterion="BIC+"} can be used which utilizes efficient Bayesian Model Averaging based on the leaps 
algorithm of package \code{"BMA"}. Also, one may consider a model selection or averaging strategy not implemented here and combine estimates "by hand"
according to formulae (7)-(10) in the cited reference below.

The function provides access to linear, logistic, Poisson and Cox proportional hazard models; one may add a random intercept to each of these models with the
\code{id} option. Other models are not supported yet. Variables used for the imputation model but not needed for the analysis model can be removed with 
option \code{var.remove}.
}
\value{
Returns an object of \code{class} `mami': 
\item{coefficients}{A matrix of coefficients, standard errors and confidence intervals for model \emph{averaging} estimators.}
\item{coefficients.boot}{A matrix of coefficients and bootstrap results (confidence intervals, mean, standard error) for model \emph{averaging}.}
\item{coefficients.s}{A matrix of coefficients, standard errors and confidence intervals for model \emph{selection} estimators.}
\item{coefficients.boot.s}{A matrix of coefficients and bootstrap results (confidence intervals, mean, standard error) for model \emph{selection}.}
\item{variable.importance}{A vector containing the variable importance for each variable based on model averaging weights.}
\item{boot.results}{A list of detailed estimation results for each bootstrap sample. The first list element refers to the results from model
selection the second entry the results from model averaging.}
}

\references{
Schomaker, M., Heumann, C. (2014) \emph{Model Selection and Model Averaging after Multiple Imputation},
Computational Statistics & Data Analysis, 71:758-770
}
\author{
Michael Schomaker
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.mami}} to visualize bootstrap results, \code{\link{lae}} and \code{\link{mma}} for model averaging techniques.
}
\examples{
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
}

\keyword{Multiple Imputation}
\keyword{Model Averaging}
\keyword{Model Selection}
\keyword{Bootstrapping}
