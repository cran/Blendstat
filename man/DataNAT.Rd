\name{DataNAT}

\alias{DataNAT}

\docType{data}

\title{Dataset, natural cherry coffee.}

\description{Database of coffee blends of different varieties dry processed.}

\usage{data(DataNAT)}

\format{Database formed by variables: Exp (identification of experiments); CEB (specialty Arabica Genotype Yellow Bourbon); CT (Commercial coffee) ; CC (Canephora coffee); CEA (specialty Arabica Genotype Acaia); Conc (Concentration: 7\% m/v and 10\% in 35 g/500 ml) concentrations were used in the experiments dry Process Coffee. Dependent variables defined by sensory attributes: Body, Flavour, Acidity, Bitterness, Score.}

\references{
CNPq for their aid via grant number 304974/2015-3.
}

\examples{
data(DataNAT) # dataset

Dados <- DataNAT

Exp <- Dados[,2]   # identification of experiments

X   <- Dados[,3:6] # independent variable
Y   <- Dados[,11]  # dependent variable 
 
# effects of the blends in a reference mixture 
Effects <- rep(c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),4) 

Conc <- as.matrix(Dados[,7]) # covariate (process variable)

Res <- Blend(Exp = Exp, X = X, Y = Y, Conc = Conc, Effects = Effects)

print("Predicted values:"); Res$MPred
print("Values predicted by the components:"); Res$MCPred
print("Design matrix:"); Res$MExp
print("Estimates of linear model parameters:"); Res$Theta 

}

\keyword{Dataset}