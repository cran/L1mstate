\name{probs.l1mstate}
\alias{probs.l1mstate}
\title{
Compute subject-specific or overall transition probabilities
}
\description{
Compute subject-specific or overall transition probabilities.
}
\usage{
probs.l1mstate(object, longdt, tmat, predt, direction=c("forward","fixedhorizon"))
}
\arguments{
  \item{object}{fitted cumhaz.l1mstate object}
  \item{longdt}{long-format data input}
  \item{tmat}{a transition matrix for multi-state model}
  \item{predt}{a prediction time}
  \item{direction}{indicates whether prediction is forward or for a fixed horizon}
}

\value{
   A list with each element [[s]] containing a data frame with the estimated transition probabilities from state s. It also includes transition and tmat information for plotting purpose
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
}

\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
