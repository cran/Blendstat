\name{Plot.Blend}

\alias{Plot.Blend}

\title{Plots of the results.}

\description{Plots of the results of the joint analysis of the experiments.}

\usage{Plot.Blend(BL, titles = c(NA,NA), posleg = 2, xlabel = NA, 
           ylabel = NA, boxleg = FALSE, color = TRUE, expcolor = NA, 
           casc = TRUE)}
 
\arguments{
  \item{BL}{Data of the Blend function.}
  \item{titles}{Titles for the plot of the effects of the concentrations and components. If it is not defined, it assumes the default text.}
  \item{posleg}{1 for caption in the left upper corner,\cr
                2 for caption in the right upper corner (default),\cr
                3 for caption in the right lower corner,\cr
                4 for caption in the left lower corner.}              
  \item{xlabel}{Names the X axis, if not set, assumes the default text.}
  \item{ylabel}{Names the Y axis, if not set, assumes the default text.}
  \item{boxleg}{Puts frame on the caption (default = TRUE).}
  \item{color}{Colorful plots (default = TRUE).}
  \item{expcolor}{Vector with the colors of the experiments.}
  \item{casc}{Cascade effect in the presentation of the plots (default = TRUE).}

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

X <- DataCD[,3:6] # independent variables (components)
Y <- DataCD[,11]  # dependent variable (response Bitterness)
 
# effects o the mixtures in the reference mixture
Effects <- rep(c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),4) 

Conc <- as.matrix(DataCD[,7]) # covariable (process variable)

Res <- Blend(exp = Exp, X = X, Y = Y, conc = Conc, effects = Effects)

print("Predicted and observed values"); Res$MPred
print("Values predicted by components:"); Res$MCPred
print("Design of the experiments:"); Res$MExp
print("Estimates of the linear model parameters:"); Res$Theta 

Tit <- c("Covariable (process variable)","Variable")

Xlab = "Effects"  # label of the X axis
Ylab = "Predicted values"  # label of the Y axis

Plot.Blend(Res, titles = Tit, posleg = 2, xlabel = Xlab,
           ylabel = Ylab, boxleg = TRUE, color = TRUE,
           expcolor = c("goldenrod3","gray53","red2", "blue2"), 
           casc = TRUE)

}

\keyword{Plots of the results}
