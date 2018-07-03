\name{Plot.Blend}

\alias{Plot.Blend}

\title{Plots of the results.}

\description{Plots of the results of the joint analysis of the experiments.}

\usage{Plot.Blend(BL, Titles = c(NA,NA), PosLeg = 2, xlabel = NA, 
           ylabel = NA, BoxLeg = FALSE, Color = TRUE, Casc = TRUE)}
 
\arguments{
  \item{BL}{Data of the Blend function.}
  \item{Titles}{Titles for the plot of the effects of the concentrations and components. If it is not defined, it assumes the default text.}
  \item{PosLeg}{1 for caption in the left upper corner,\cr
                2 for caption in the right upper corner (default),\cr
                3 for caption in the right lower corner,\cr
                4 for caption in the left lower corner.}              
  \item{xlabel}{Names the X axis, if not set, assumes the default text.}
  \item{ylabel}{Names the Y axis, if not set, assumes the default text.}
  \item{BoxLeg}{Puts frame on the caption (default = TRUE).}
  \item{Color}{Coloridos plots (default = TRUE).}
  \item{Casc}{Cascade effect in the presentation of the plots (default = TRUE).}

}

\value{Return several plots.}

\author{
Marcelo Angelo Cirillo

Paulo Cesar Ossani
}

\seealso{\code{\link{Blend}}}

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

Tit <- c("Covariable (process variable)","Variable")

Xlab = "Effects"  # label of the X axis
Ylab = "Predicted values"  # label of the Y axis

Plot.Blend(Res, Titles = Tit, PosLeg = 2, xlabel = Xlab,
           ylabel = Ylab, BoxLeg = TRUE, Color = TRUE,
           Casc = TRUE)

}

\keyword{Plots of the results}
