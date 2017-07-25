library(quantmod)
library(corrplot) 
library(QRM)
library(tseries)

#Set the directory using 'setwd()'
#Import 10 year monthly  historical data of the stocks and market
SP500 <- read.csv("SP500TR.csv",header=T)
MMM <- read.csv("MMM.csv",header=T)
JNJ <- read.csv("JNJ.csv",header=T)
PG <- read.csv("PG.csv",header=T)
WMT <- read.csv("WMT.csv",header=T)

#Plot closing prices for each stock
plot(MMM[,"Close"], main = "MMM", xlab="Months", ylab="Price", col="dark red")
plot(JNJ[,"Close"], main = "JNJ", xlab="Months", ylab="Price", col="dark red")
plot(PG[,"Close"], main = "PG", xlab="Months", ylab="Price", col="dark red")
plot(WMT[,"Close"], main = "WMT", xlab="Months", ylab="Price", col="dark red")
plot(SP500[,"Close"], main = "SP500", xlab="Months", ylab="Price", col="dark red")


#Combine price columns for all 4 stocks
data <- cbind.data.frame(SP500$Close, MMM$Close, JNJ$Close, PG$Close, WMT$Close)

#Calculate Mothly returns
#Denote n the number of time periods:
n <- nrow(data)
SP500_returns <- ((data[2:n, 1] - data[1:(n-1), 1])/data[1:(n-1), 1])
MMM_returns <- ((data[2:n, 2] - data[1:(n-1), 2])/data[1:(n-1), 2])
JNJ_returns <- ((data[2:n, 3] - data[1:(n-1), 3])/data[1:(n-1), 3])
PG_returns <- ((data[2:n, 4] - data[1:(n-1), 4])/data[1:(n-1), 4])
WMT_returns <- ((data[2:n, 5] - data[1:(n-1), 5])/data[1:(n-1), 5])

#Combine returns
data_returns <- cbind.data.frame(SP500_returns, MMM_returns, JNJ_returns, PG_returns, WMT_returns)

#Calculate mean and standard deviation
Mean <- sapply(data_returns,mean)
SD=sapply(data_returns,sd)

#Combine calculated mean and standard-deviation
cbind(Mean,SD)

#Convert numeric matrix into data.frame
data_returns <- as.data.frame(data_returns)

#Run the regression to get the beta of stocks seperately
lm.MMM <- lm(MMM_returns~SP500_returns,data_returns)
summary(lm.MMM)     # we get Beta of 'MMM'

lm.JNJ <- lm(JNJ_returns~SP500_returns,data_returns)
summary(lm.JNJ)    #we get beta of 'JNJ'

lm.PG <- lm(PG_returns~SP500_returns,data_returns)
summary(lm.PG)    #we get beta of 'PG'

lm.WMT <- lm(WMT_returns~SP500_returns,data)
summary(lm.WMT)   #we get beta of 'WMT

#Assign Beta values
Beta_MMM <- 0.1694
Beta_JNJ <- 0.0520
Beta_PG  <- 0.0619
Beta_WMT <- 0.1271

#Calculate Expected Return using CAPM "ER=Rf-BETA*(E(Rm)-Rf)
#Expected return of market equals to its mean value
ER_MMM     <- 0.1694*0.007143058
ER_JNJ    <- 0.0520*0.007143058
ER_PG <- 0.0619*0.007143058
ER_WMT   <- 0.1271*0.007143058


# Calculate Variance "var(stock)= Beta*Market var+ residual error var"
Var_MMM <- Beta_MMM^2*sd(SP500_returns)^2+0.05392^2
Var_JNJ <- Beta_JNJ^2*sd(SP500_returns)^2+0.05392^2
Var_PG <- Beta_PG^2*sd(SP500_returns)^2+0.05392^2
Var_WMT <- Beta_WMT^2*sd(SP500_returns)^2+0.05392^2

#Calculate covariation matrix (exclude market data)
covmat <- cov(data_returns[,-1])

#Build expected return matrix of the parameters
ER_matrix <- matrix(c(ER_MMM,ER_JNJ,ER_PG,ER_WMT), nrow=1)

#Get the optimal weights
weights <- portfolio.optim(ER_matrix,covmat = covmat,shorts=F)

#MMM=0.175; JNJ=0.207; PG=0.251; WMT=0.367

 
#Calculate and plot a correlation matrix 
#library(corrplot)
cormat <- cor(data_returns[,-1])
colnames(cormat) <- c("MMM","JNJ","PG","WMT")
rownames(cormat) <- c("MMM","JNJ","PG","WMT")
corrplot(cormat, method = "circle",tl.col = "black",main = "correlation matrix") #plot matrix

#############################################################################################

#Start to define variables for simulation
#Maturity in years
maturity <- 10

#Use monthly steps
nsteps <- maturity*12
dt <-  maturity / nsteps

#Number of assets
nAssets = 4

#Number of simulations
nTrails = 10000

#Expected return p.a. for each asset, stored in vector 
E.R <- rep(NA,nAssets)
E.R[1] <- ER_MMM   
E.R[2] <- ER_JNJ
E.R[3] <- ER_PG
E.R[4] <- ER_WMT

#Define variable size
simulated.Returns   <- array(NA,    dim = c(nsteps+1, nTrails, nAssets))
cumulative.PortReturns  <- matrix(rep(NA,nsteps*nTrails), nrow = nsteps, ncol = nTrails)
ES  <- rep(NA, maturity)

#Define portfolio weights
port.weights <- weights$pw
  

#Perform cholesky decomposition
Chol <- (chol(cormat))

#Generate standard-normal, random variables
x <- array(rnorm(nsteps*nTrails*nAssets), c(nsteps*nTrails,nAssets))

#Generate correlated standard-normal, random variables
ep <- x %*% Chol

#Define the drift
drift <- E.R/12 - 0.5 * diag(covmat)

#Generate asset paths
temp = array(exp(as.vector(drift %*% t(dt)) + t(ep *sqrt(diag(covmat)))), c(nAssets,nsteps,nTrails))

for(i in 2:nsteps) temp[,i,] = temp[,i,] * temp[,(i-1),] 

#Change dimension of the array temp from dim(nAssets, nsteps, nTrails) to dim(nsteps, nAssets, nTrails)
simulated.Returns <- aperm(temp,c(2,1,3))

#Compute portfolio returns for each simulation (nTrails). To do this, each step is weighted with "port.weights"
#Since I generate continuous returns, I first transform them into discrete, multiply with weights and then transform back into continuous.
for (z in 1:nTrails) {
  for (i in 1:nsteps) cumulative.PortReturns[i,z] = log(1+((exp(simulated.Returns[i,,z]-1)-1) %*% port.weights))
}

#Finally, compute the monthly expected shortfall (5%-level) by taking the average of the 5% worst portfolio yields
#Do steps of 12 as  the expected shortfall is calculated at the end of each year
z = 0
for (i in seq(12, nsteps, by = 12 )) {
  z = z + 1 
  ES[z]   <- mean(sort(cumulative.PortReturns[i,]) [1:(0.05*nTrails)])
}

#plot a sample of simulated portfolio returns (P.S every time we run different simulation,
#we get a new plot)
#library(QRM)
par(mar=c(2,2,1,1)) #Set margins
plot(as.timeSeries(cumulative.PortReturns[,1:100]), plot.type = 'single')

















