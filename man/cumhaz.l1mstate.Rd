\name{cumhaz.l1mstate}
\alias{cumhaz.l1mstate}

\title{
Compute subject-specific transition hazard rates
}
\description{
This function computes subject-specific cumulative transition hazard rates for each of the possible transitions in the multi-state model.
}
\usage{
cumhaz.l1mstate(object, longdt, newdata, cv.method = c("pcvl", "min"))
}

\arguments{
  \item{object}{fitted cv.l1mstateR object.}
  \item{longdt}{long-format data input}
  \item{newdata}{a new data with the same long-format with the same covariate names as longdt.}
  \item{cv.method}{the cross-validation method used to select the optimal result.}
}

\value{
 \item{time}{A list of all time points for each of the transitions in the multi-state model.}
 \item{baseHaz}{A list of the estimated subject-specific baseline hazards for each of the transitions in the multi-state model.}
 \item{var}{A list of the variance of the estimated subject-specific hazards for each of the transitions in the multi-state model.}
 \item{cumbaseHaz}{A list of the estimated subject-specific cumulative hazards for each of the transitions in the multi-state model.}
  \item{Haz}{A list of the estimated subject-specific hazards for each of the transitions in the multi-state model.}
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
}
\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
