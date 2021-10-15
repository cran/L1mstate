\name{plot.cumhaz.l1mstate}
\alias{plot.cumhaz.l1mstate}
\title{
 Plot the estimated cumulative hazard rates of the multi-state model. 
}
\description{
 Plot the estimated cumulative hazard rates of the multi-state model.
}
\usage{
\method{plot}{cumhaz.l1mstate}(x,type=c("single","separate"),cols,
     xlab="Years since transplant",ylab="Cumulative hazard",
     ylim,lwd=3,lty,legend,legend.pos,bty="o",...)
}
\arguments{
  \item{x}{fitted cumhaz.l1mstate object}
  \item{type}{One of "single" or "separate"; in case of "single", all estimated cumulative hazards are drawn in a single plot, in case of "separate", separate plots are shown for the estimated cumulative hazards.}
  \item{cols}{A vector specifying colors for the different transitions}
  \item{xlab}{A title for the x-axis; default is "Years since transplant"}
  \item{ylab}{A title for the y-axis; default is "Cumulative hazard"}
  \item{ylim}{The y limits of the plot}
  \item{lwd}{The line width; default is 3}
  \item{lty}{The line type}
  \item{legend}{The transition numbers; if missing, these will be taken from the transition matrix contained in cumhaz.l1mstate object.}
  \item{legend.pos}{The position of the legend; default is "topleft"}
  \item{bty}{The box type of the legend}
  \item{\dots}{Further arguments to plot}
}
\value{No return value}
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

ptsA <- longdt[which(longdt$X1==1 & longdt$X2==0 & longdt$X3==1 & longdt$X4==0 & 
                     longdt$X5==1 & longdt$X6==1),]
## observed transitions (ground truth)
# predicted time = 0
events(ptsA)
# use models to predict the transition probabilities
# L1MSTATE
ptA <- ptsA[which(ptsA$id == unique(ptsA$id)[1]),]
ptA <- ptA[,c(4,9:14)]

cumhazA <- cumhaz.l1mstate(object=out, longdt=longdt, newdata=ptA, cv.method="pcvl")

plot.cumhaz.l1mstate(cumhazA, main = "L1MSTATE", type = "single", lwd=2,las=1,ylim = c(0,4))

legend("topleft",inset=.01,legend=c("1","2","3"),col=c("black","red","green"),
       lty=1,cex=0.8,title="Transitions",text.font=2,bg='white')
}

\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
