# SpLoc

Spatially Localized Signals

R code to apply SpLoc to longitudinal cortical thickness data. The current version supports parallel computing using the *doParallel* package.

* Park, J.Y., Fiecas, M. "Permutation-based inference for spatially localized signals in longitudinal MRI data." Preprint.


## Installation
To install the latest development builds directly from GitHub, run this:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("junjypark/SpLoc")
```

## Usage
Fitting SpLoc for longitudinal data requires a number of components:

* **data**: a MxS matrix, where S is the number of spatial locations and M is the sum of the number of visits for all subjects. Each row represents an image.

* **X**: a NxQ matrix, where N is the number of subjects and Q is the number of covariates, including intercepts.

* **NNmatrix**: a QxS binary sparse matrix (1 or 0), where Q is the number of neighbors pre-specified by user.

* **group**: a binary vector (1 or -1) of length N that specifies the group.

* **time** (for longitudinal data only): a vector of length M that stores time of the visit for the image corresponding to each row of the data matrix.

* **n.visits** (for longitudinal data only): a vector of length N that contains the number of visits (scans) for each subject.

The SpLoc can be performed using the followings:

```R
getResid=getSummaryMatrix(data, X, mask=rep(1, ncol(data)),longitudinal=T, n.visits, randomslope=T,  NNmatrix, time.var)
fit=SpLoc(ymat, NNmatrix, group=group, nperm=1000, alpha=0.05, seed=1234)
select=ClusterSearch(fit$Tstat, fit$threshold, NNmatrix)
```


## Example

```R
#Load R Package
library(SpLoc)

#Specify parameters for LME
tau2=0.5                          #Variance of noise
b.cov=matrix(c(3,0.5,0.5,0.2),2)  #Covariance of random intercept and slope

alpha0=1
alpha1=c(1,-1,0.5)
beta0=0.5
beta1=1
gamma=0.5                         #Parameter of interest

n.subj=50                         #Number of subjects

#Load nearest neighbor information (sparse matrix format in R)
NNmatLH=readRDS("NNmat_LH.rds") #Left hemisphere
NNmatRH=readRDS("NNmat_RH.rds") #Right hemisphere

# Define signals in Design 1
ind.signal.lh=which(NNmatLH[,28000]==1)
ind.signal.rh=which(NNmatRH[,28000]==1)

# Define signals in Design 2
# ind.signal.lh=c(which(NNmatLH[,14000]==1),which(NNmatLH[,15000]==1), which(NNmatLH[,16000]==1))
# ind.signal.rh=c(which(NNmatRH[,14000]==1),which(NNmatRH[,15000]==1), which(NNmatRH[,16000]==1))

# Define signals in Design 3
# ind.signal.lh=c(which(NNmatLH[,6500]==1),which(NNmatLH[,8500]==1), which(NNmatLH[,10500]==1),
#                 which(NNmatLH[,12500]==1),which(NNmatLH[,14500]==1))
# ind.signal.rh=c(which(NNmatRH[,6500]==1),which(NNmatRH[,8500]==1), which(NNmatRH[,10500]==1),
#                 which(NNmatRH[,12500]==1),which(NNmatRH[,14500]==1))

###Step 1 : generate simulation data
n.visits=sample(3:4, n.subj, replace = T)
time=NULL 
for (i in 1:n.subj){
  time=c(time, (0:(n.visits[i]-1))/2 ) 
}

X=matrix(rnorm(n.subj*3), n.subj) 
X.expand=X[rep(1:n.subj, n.visits),]
dx.status=c(rep(-1, n.subj/2),rep(1,n.subj/2)) 
dx.expand=rep(dx.status, n.visits)

pred=c(alpha0+X.expand%*%alpha1+dx.expand*beta0+time*beta1)
ymat.lh=matrix(NA,sum(n.visits),2329) 
ymat.rh=matrix(NA,sum(n.visits),2332) 

#Generate data for the left and right hemispheres
for (j in 1:2329){ 
  b=mvrnorm(n.subj,c(0,0), b.cov) 
  b.expand=b[rep(1:50, n.visits),]
  epsilon=rnorm(sum(n.visits),0,sqrt(tau2))
  ymat.lh[,j]=pred+b.expand[,1]+time*b.expand[,2]+epsilon 
}

for (j in 1:2332){
  b=mvrnorm(n.subj,c(0,0), b.cov)
  b.expand=b[rep(1:50, n.visits),]
  epsilon=rnorm(sum(n.visits),0,sqrt(tau2))
  ymat.rh[,j]=pred+b.expand[,1]+time*b.expand[,2]+epsilon
}

#Add signals
ymat.lh[,ind.signal.lh]=ymat.lh[,ind.signal.lh]+matrix(rep(gamma*dx.expand*time, length(ind.signal.lh)), sum(n.visits) )
ymat.rh[,ind.signal.rh]=ymat.rh[,ind.signal.rh]+matrix(rep(gamma*dx.expand*time, length(ind.signal.rh)), sum(n.visits) )
```
