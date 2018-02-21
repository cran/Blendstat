\name{Plot.Blend}

\alias{Plot.Blend}

\title{Optimization graphs.}

\description{Trace plot of the experiments as a function of the components.}

\usage{Plot.Blend(BL, Title = NULL, PosLeg = 2, BoxLeg = FALSE, Color = TRUE)}
 
\arguments{
  \item{BL}{Function data Blend.}
  \item{Title}{Title for the graph of the effects of concentrations. If not set, assume the default text.}
  \item{PosLeg}{1 - caption left,\cr
                2 - caption right (default),\cr
                3 - caption above,\cr
                4 - caption below.}
  \item{BoxLeg}{Place frame in legend (default = FALSE).}
  \item{Color}{Colorful graphics (default = TRUE)}

}

\value{Return several graphics.}

\author{
Marcelo Angelo Cirillo

Paulo Cesar Ossani
}

\seealso{\code{\link{Blend}}}

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

\keyword{Optimization}
