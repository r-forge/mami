\name{HIV}
\alias{HIV}
\docType{data}
\title{
HIV treatment data
}
\description{
A hypothetical, simulated dataset which provides typical data from HIV treatment research. 
}
\usage{data(HIV)}
\format{
  A data frame with 5000 observations on the following 14 variables:
  \describe{
    \item{\code{patient}}{The patient's id number}
    \item{\code{hospital}}{The treatment facility}
    \item{\code{futime}}{Follow-up time of the respective patient}
    \item{\code{dead}}{Outcome at the end of follow-up: (1=dead, 0=alive)}
    \item{\code{sex}}{The patients's sex (1=female, 0=male)}
    \item{\code{age}}{The patient's age}
    \item{\code{cd4}}{The patient's baseline CD4 cell count}
    \item{\code{cd4slope6}}{The patients CD4 slope (average CD4 gain/loss per week) at 6 month of follow up}
    \item{\code{weight}}{The patient's baseline weight (in kg)}
    \item{\code{period}}{Categorical variable specifying the period of treatment initiation: before 2004, between 2004 and 2007, or after 2007}
    \item{\code{haem}}{The patient's baseline haemoglobin level}
    \item{\code{stage}}{The patient's baseline WHO stage indicating disease severity, ordinal variable}
    \item{\code{tb}}{The patient's baseline tuberculosis status (1=yes,0=no)}
    \item{\code{cm}}{The patient's baseline cryptococcal meningitis status (1=yes,0=no)}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%
}
\references{
Schomaker, M. et al. (2013) \emph{Immune Recovery After Starting ART in HIV-infected Patients Presenting and not Presenting with Tuberculosis in South Africa}, Journal of Acquired Immune Deficiency Syndroms, 63:142-145

Hoffmann, C. et al. (2013)  \emph{CD4 Count Slope and Mortality in HIV-infected Patients on Antiretroviral Therapy: Multi-cohort Analysis from South Africa}, Journal of Acquired Immune Deficiency Syndroms, 63:34-41
}
\examples{
data(HIV)
str(HIV)
}
\keyword{Datasets}
\keyword{HIV}