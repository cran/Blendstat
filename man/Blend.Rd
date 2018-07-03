\name{Blend}

\alias{Blend}

\title{Joint analysis of experiments with mixtures and random effects.}

\description{Joint analysis of experiments with mixtures and random effects, taking on a process variable represented by a covariable.}

\usage{Blend(Exp, X, Y, Conc = NULL, Effects = NULL)}
\arguments{
  \item{Exp}{Vector with the names of the experiments.}
  \item{X}{Mixture variables (components), without the vector of the concentrations (covariable).}
  \item{Y}{Response variable.}
  \item{Conc}{Vector with the concentrations (covariable) of the experiments.}
  \item{Effects}{Vector of the effects of the mixtures in a reference mixture (example: centroid).}
}

\value{
\item{MPred}{Matrix with the predicted and observed values.}
\item{MCPred}{Matrix with the values predicted by components.}
\item{MExp}{Matrix with the design of the experiments.}
\item{Theta}{Vector with the Theta estimates.}
}

\references{
Kalirajan, K. P.. On the estimation of a regression model with fixed and random coefficients. \emph{Journal of Applied Statistics}, 17(2): 237-244, 1990. doi:10.1080/757582835

Swany, P. A. V. B. \emph{Statistical Inference in Random Coefficient Regression Models.} Amsterdam: Springer Science & Business Media, 1971. 209 p.
}

\author{
Marcelo Angelo Cirillo

Paulo Cesar Ossani
}

\seealso{\code{\link{Plot.Blend}}}

\examples{
data(DataNAT) # dataset 

Exp <- DataNAT[,2]   # identification of experiments

X   <- DataNAT[,3:6] # independent variable
Y   <- DataNAT[,11]  # dependent variable 
 
# effects of the blends in a reference mixture 
Effects <- rep(c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),4) 

Conc <- as.matrix(DataNAT[,7]) # covariate (process variable)

Res <- Blend(Exp = Exp, X = X, Y = Y, Conc = Conc, Effects = Effects)

print("Predicted and observed values"); Res$MPred
print("Values predicted by components:"); Res$MCPred
print("Design of the experiments:"); Res$MExp
print("Estimates of the linear model parameters:"); Res$Theta 

Tit <- c("Covariate (process variable)","Variable")

Xlab = "Effects"  # label of the X axis
Ylab = "Predicted values"  # label of the Y axis

Plot.Blend(Res, Titles = Tit, PosLeg = 2, xlabel = Xlab,
           ylabel = Ylab, BoxLeg = TRUE, Color = TRUE,
           Casc = TRUE)

}

\keyword{Joint analysis}