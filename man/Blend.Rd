\name{Blend}

\alias{Blend}

\title{Joint analysis of experiments with mixtures and random effects.}

\description{Joint analysis of experiments with mixtures and random effects, assuming a process variable, represented by a covariate.}

\usage{Blend(Exp, X, Y, Conc = NULL, Effects = NULL)}
\arguments{
  \item{Exp}{Vector with the names of the experiments.}
  \item{X}{Regressor variables, not including the covariate.}
  \item{Y}{Dependent variable}
  \item{Conc}{Vector with concentrations (covariate) of the experiments.}
  \item{Effects}{Vector of the effects of the blends in a reference mixture (example: centroid).}
}

\value{
\item{MPred}{Matrix with predicted and observed values.}
\item{MCPred}{Matrix with the values predicted by the components.}
\item{MExp}{Desing of the experiments }
\item{Theta}{Vector with Theta estimates.}
}

\references{
Kalirajan, K. P.. On the estimation of a regression model with fixed and random coefficients. \emph{Journal of Applied Statistics}, 17: 237-244, 1990. doi:10.1080/757582835
}

\author{
Marcelo Angelo Cirillo

Paulo Cesar Ossani
}

\seealso{\code{\link{Plot.Blend}}}

\examples{
data(DataCD) # dataset  
#data(DataNAT) # dataset 

Dados <- DataCD

Exp <- Dados[,2]   # identification of experiments

X   <- Dados[,3:5] # independent variable
Y   <- Dados[,11]  # dependent variable 
 
# effects of the blends in a reference mixture 
Effects <- rep(c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),4) 

Conc <- as.matrix(Dados[,7]) # covariate (process variable)

Res <- Blend(Exp = Exp, X = X, Y = Y, Conc = Conc, Effects = Effects)

print("Predicted values:"); Res$MPred
print("Values predicted by the components:"); Res$MCPred
print("Design matrix:"); Res$MExp
print("Estimates of linear model parameters:"); Res$Theta 

Tit <- c("Covariate (process variable)")
Plot.Blend(Res, Title = Tit, PosLeg = 2, BoxLeg = TRUE, Color = TRUE)

}

\keyword{Joint analysis}