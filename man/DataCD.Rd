\name{DataCD}

\alias{DataCD}

\docType{data}

\title{Dataset, peeled cherry coffee.}

\description{Database of coffee blends of different varieties processed via wet (peeled cherry).}

\usage{data(DataCD)}

\format{Database of coffee blends of different varieties processed via wet (peeled cherry). Formed by the variables: Exp (code of the experiments); CEB (specialty Bourbon Yellow coffee produced at an altitude above 1,200m); CT (roasted commercial coffee); CC (Conillon coffee); CEA (Acaia specialty coffee produced at altitude below 1,100m); Conc (concentrations at 7\% and 10\% (m/v) of roasted and ground coffee beans in 100 ml of water). Response variables defined by the sensorial attributes: Body, Taste, Acidity, Bitterness, Score.}

\references{
Project yield and research  entitled  by "Quality of blends of specialty and non-specialty coffees of the region of the Mantiqueira Mountains - treatment of discrepant scores in tests with consumers". CNPq for their aid via grant number 304974/2015-3.
}

\examples{
data(DataCD) # dataset

Exp <- DataCD[,2]   # identification of the experiments

X   <- DataCD[,3:6] # independent variables (components)
Y   <- DataCD[,11]  # dependent variable (response Bitterness)
 
# effects o the mixtures in the reference mixture
Effects <- rep(c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),4) 

Conc <- as.matrix(DataCD[,7]) # covariable (process variable)

Res <- Blend(Exp = Exp, X = X, Y = Y, Conc = Conc, Effects = Effects)

print("Predicted and observed values"); Res$MPred
print("Values predicted by components:"); Res$MCPred
print("Design of the experiments:"); Res$MExp
print("Estimates of the linear model parameters:"); Res$Theta 

}

\keyword{Dataset}