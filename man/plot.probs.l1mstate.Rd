\name{plot.probs.l1mstate}
\alias{plot.probs.l1mstate}
\title{
Plot the transition probabilities
}
\description{
Plot the transition probabilities produced by probs.l1mstate
}
\usage{
\method{plot}{probs.l1mstate}(x,from,type=c("stacked","filled","single","separate"),
     ord,cols,xlab="Years since transplant",ylab="Probability",
     xlim,ylim,lwd,lty,cex,legend,legend.pos,bty="o",...)
}
\arguments{
  \item{x}{object produced from probs.l1mstate}
  \item{from}{the starting state from which the probabilities are used to plot}
  \item{type}{with "stacked" type, the transition probabilities are stacked and the distance between two adjacent curves indicates the probability; with "filled" type, it is the same but the space between adjacent curves are filled; with "single" type, the probabilities are shown as different curves in a single plot; with "separate" type, separate plots are shown for the estimated transition probabilities}
  \item{ord}{A vector indicates the order of plotting in case "stacked" or "filled"}
  \item{cols}{A vector specifying colors for the different transitions}
  \item{xlab}{A title for the x-axis; default is "Years since transplant"}
  \item{ylab}{A title for the y-axis; default is "Probability"}
  \item{xlim}{The x limits of the plot; default is range of time}
  \item{ylim}{The y limits of the plot}
  \item{lwd}{The line width; default is 3}
  \item{lty}{The line type}
  \item{cex}{Character size}
  \item{legend}{The transition numbers; if missing, these will be taken from the transition matrix contained in cumhaz.l1mstate object.}
  \item{legend.pos}{The position of the legend; default is "topleft"}
  \item{bty}{The box type of the legend}
  \item{\dots}{Further arguments to plot}
}
\value{No return value.}
\details{
A plot is produced, and nothing is returned.
}
\author{
Xuan Dang \email{xuandang11289@gmail.com}
}

\examples{
library(L1mstate)
set.seed(1001)

N <- 200

x1 <- rbinom(N, 1, .5)
x2 <- rbinom(N, 1, .5)
x3 <- rbinom(N, 1, .5)
x4 <- rbinom(N, 1, .5)
x5 <- rbinom(N, 1, .5)
x6 <- rbinom(N, 1, .5)
x <-  data.frame(x1,x2,x3,x4,x5,x6)
col_names  <-  c(sprintf("X\%d", seq(1,6)))
colnames(x)  <-  col_names

beta12 <- c(-.65,-.65,-.65,0,0,0)
beta13 <- c(-.65,0,0,0,-.65,0) 
beta23 <- c(0,-.65,-.65,0,0,-.65)

N12  <-  N-2*N\%/\%5
N13 <- N-N12
hx12 <- exp(as.matrix(x[1:N12,]) \%*\% beta12) 
ty12 <- rexp(N12,hx12) 
hx13 <- exp(as.matrix(x[(N12+1):N,]) \%*\% beta13) 
ty13 <- rexp(N13,hx13)
hx23 <- exp(as.matrix(x[1:N12,]) \%*\% beta23) 
ty23 <- rexp(N12,hx23) 
t2 <- rep(0,N)
t2[1:N12]  <-  ty12
t2[(N12+1):N]  <-  ty13
s2  <-  rep(0,N)
s2[1:N12] <-  rep(1,N12)
t3 <- rep(0,N)
t3[1:N12]  <-  ty12+ty23
t3[(N12+1):N]  <-  ty13
#censor case for state 3
tcens3 <- rbinom(n=N, prob = 0.3, size = 1)
s3 <- 1-tcens3

# pre-process data into long-format using msprep
dt <- data.frame(illt=t2, ills=s2, dt=t3, ds=s3, x)
tmat <- matrix(c(NA,NA,NA,1,NA,NA,2,3,NA), nrow = 3)
longdt <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 keep = col_names, data=dt,trans=tmat)

out <- cv.l1mstateR(longdt, nlambda = 100, nfolds = 10)

ptsA <- longdt[which(longdt$X1==1 & longdt$X2==0 & longdt$X3==1 & longdt$X4==0 
               & longdt$X5==1 & longdt$X6==1),]
## observed transitions (ground truth)
# predicted time = 0
events(ptsA)
# use models to predict the transition probabilities
# L1MSTATE
ptA <- ptsA[which(ptsA$id == unique(ptsA$id)[1]),]
ptA <- ptA[,c(4,9:14)]

cumhazA <- cumhaz.l1mstate(object=out, longdt=longdt, newdata=ptA, cv.method="pcvl")
probA_0 <- probs.l1mstate(cumhazA, longdt = longdt, tmat, predt = 0, direction = "forward")
statecols <- heat_hcl(6, c = c(90, 10), l = c(20, 80), power = c(1/5, 2))[c(6, 5, 3, 4, 2, 1)]
ord <- c(1, 2, 3, 4, 5, 6)
plot.probs.l1mstate(probA_0, main = "L1MSTATE", from = 1, ord = ord,
                    las=1, cex=0, type = "filled", col = statecols[ord])
}
\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
