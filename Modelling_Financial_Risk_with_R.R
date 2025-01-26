#Part 1: Set up the data
rm(list=ls())
require(zoo)

setwd("C:/Users/Admin/Downloads/Semester 5/Financial analysis with R") # set working directory

fund     <- "ALL004.mst"
raw.data <- read.csv2(file = fund, head=TRUE, sep=",", dec=".")
dateForm <- "%Y%m%d"
prices <- raw.data$"X.OPEN."
dates  <- as.Date(as.character(raw.data$"X.DTYYYYMMDD."), format=dateForm)

P      <- zoo(prices , order.by = dates)
r      <- diff(log(P))

#Part 2: Descriptive Statistics, check for normality and autocorrelation
require(zoo)
require(ggplot2)
require(rugarch)

# Panel of figures for full sample 
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Level", xlab="", ylab="")
plot(r, main="returns", xlab="", ylab="")
# Settings for analysis                                  
T  <- length(r)                        # full sample, nobs
N  <- 1250    # we shorten the sample to the last 5 years
P_shorten <- tail(P,N)
r  <- tail(r,N)                        # log-returns
R  <- coredata(r)
# Descriptive statistics
plot(P_shorten, main="Level", xlab="", ylab="")
plot(r, main="returns", xlab="", ylab="") 

mu = sum(R)/N # mean

R0    <- R - mu
M2    <- sum(R0^2)/N
M3    <- sum(R0^3)/N
M4    <- sum(R0^4)/N

sig <- sqrt(M2)       ## volatility
S   <- M3/(sig^3)     ## skewness
K   <- M4/(sig^4)     ## kurtosis
require(moments)
M1   = moment(R, order = 1, central = FALSE, na.rm = TRUE)
M2   = moment(R, order = 2, central = TRUE, na.rm = TRUE)
M3   = moment(R, order = 3, central = TRUE, na.rm = TRUE)
M4   = moment(R, order = 4, central = TRUE, na.rm = TRUE)
mu0  = M1
sig0 = sqrt(M2)
S0   = M3/(sig0^3)
K0   = M4/(sig0^4)

# Annualization
Nyear <- 250
muA   <- mu*Nyear
sigA  <- sig*sqrt(Nyear)
# Normality test  
#D'Agostino test for skewness
agostino.test(R)
#Anscombe-Glynn test of kurtosis for normal samples
anscombe.test(R)
#Jarque-Bera test of normality
jarque.test(R)

# Empirical density  
R0        <- (R-mu)/sig #z-score
bwdth     <- 0.1

ggplot(data.frame(R0), aes(x = R0)) +
  theme_bw() +
  geom_histogram(binwidth = bwdth, colour = "white", fill = "yellow4", size = 0.1) +
  #stat_function(fun = function(x) dnorm(x)*N*bwdth, color = "red", size = 1)                        # normal distribution
  stat_function(fun = function(x) ddist("std",x,shape=5)*N*bwdth, color = "black", size = 1)    # t-Student distribution
+ xlim(-7,-2)

# Estimation of t-Studenta parameters        

# Method of moments (K - kurtosis)
v0 <- 4 + 6/(K-3)
# Max likelihood method 
require(MASS)
d0 <- fitdistr(R0, "normal") # fit normal distribution to the standardized return
d0$loglik #Extracts the log-likelihood value for the normal distribution fit. 
#A higher log-likelihood indicates a better fit.
d1 <- fitdistr(R0, "t", m = 0, start = list(s=sqrt((v0-2)/v0), df=v0), lower=c(0.001,3))
d1$loglik # fit into the t-distribution
v=d1$estimate[[2]]
# Comparison
(d1$loglik-d0$loglik)/N #-> T-student model is a better fit
# QQ plot
q        <- seq(0.001, 0.999, 0.001)
Qemp     <- quantile(R0,q)                     
Qteo     <- qt(q,v)*sqrt((v-2)/v)              
lim0    <- c(-5,5)                             
par(mfrow=c(1,1), cex = 0.7, bty="l")
plot(Qemp,Qteo, main="QQplot", col="red", xlim = lim0, ylim = lim0,
     xlab="Qemp", ylab="Qteo") 
abline(a=0,b=1, lwd=2)
# Autocorrelations                     
require(tseries)
# ACF for returns
z0 <- acf(R,20, plot=TRUE)
plot(z0$acf[2:20], type="h", main="ACF for returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
# ACF for squarred returns 
z1 <- acf(R^2,20, plot=FALSE)
plot(z1$acf[2:20], type="h", main="ACF for sq returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
# Cross-correlations
z2 <- ccf(R^2, R, 20, type="correlation", plot=TRUE)
plot(z2$acf[22:40], type="h", main="correlation for sq. ret and ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[20:1], type="h", main="correlation for ret and sq. ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
# panel of plots
par(mfrow=c(2,2), cex = 0.7, bty="l")
plot(z0$acf[2:20], type="h", main="ACF for returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z1$acf[2:20], type="h", main="ACF for sq returns", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[22:40], type="h", main="correlation for sq. ret and ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[20:1], type="h", main="correlation for ret and sq. ret lags", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
# Ljunga-Box test
Box.test(R, lag = 20, type = c("Ljung-Box"))
Box.test(R^2, lag = 20, type = c("Ljung-Box"))

# Part 3: Calculate VaR and ES for Parametric Method, Historical Simulation, Monte Carlo and CF
# Parametric Method
require(PerformanceAnalytics)
# Define parameters
confidence_level <- 0.95
z_alpha <- qnorm(confidence_level)  # Z-score for the confidence level
mu <- mean(R)
sigma <- sd(R)
# VaR calculation
VaR_parametric <- -(mu + z_alpha * sigma)
# ES calculation
ES_parametric <- -(mu + sigma * dnorm(z_alpha) / (1 - confidence_level))

# Historical Simulation
# Historical VaR
VaR_historical <- -quantile(R, 1 - confidence_level)
# Historical ES
ES_historical <- -mean(R[R <= VaR_historical])
print(VaR_historical)
print(ES_historical)

# Monte Carlo Simulation
set.seed(123)  # For reproducibility
num_simulations <- 100000
simulated_returns <- rnorm(num_simulations, mean = mu, sd = sigma)
# Monte Carlo VaR
VaR_montecarlo <- -quantile(simulated_returns, 1 - confidence_level)
# Monte Carlo ES
ES_montecarlo <- -mean(simulated_returns[simulated_returns <= VaR_montecarlo])
print(VaR_montecarlo)
print(ES_montecarlo)

# Cornish-Fisher Adjustment
# Calculate skewness and kurtosis
skewness <- PerformanceAnalytics::skewness(R)
kurtosis <- PerformanceAnalytics::kurtosis(R)
# Adjusted z-score for Cornish-Fisher
z_adjusted <- z_alpha + (1/6) * (z_alpha^2 - 1) * skewness +
  (1/24) * (z_alpha^3 - 3*z_alpha) * kurtosis -
  (1/36) * (2*z_alpha^3 - 5*z_alpha) * skewness^2
# Cornish-Fisher VaR
VaR_cornish <- -(mu + z_adjusted * sigma)
# Cornish-Fisher ES
ES_cornish <- -(mu + sigma * dnorm(z_adjusted) / (1 - confidence_level))
print(VaR_cornish)
print(ES_cornish)


# Part 4: Fitting models: MA, EWMA and GARCH family
# Load return data
# MA (Rolling Window)
rolling_window <- 20 # Example: 20-day rolling window
ma_volatility <- rollapply(R, width = rolling_window, FUN = sd, fill = NA)
# EWMA
lambda <- 0.94 # Decay factor (commonly used value)
ewma_volatility <- sqrt(filter(lambda * R^2, 1 - lambda, method = "recursive"))
# GARCH(1,1) Specification
spec_garch <- ugarchspec(variance.model = list(model = "sGARCH"), 
                         mean.model = list(armaOrder = c(0, 0)),
                         distribution.model = "std") # Student-t for fat tails

# Fit GARCH(1,1)
garch_fit <- ugarchfit(spec = spec_garch, data = R)

# EGARCH Specification
spec_egarch <- ugarchspec(variance.model = list(model = "eGARCH"), 
                          mean.model = list(armaOrder = c(0, 0)),
                          distribution.model = "std")

# Fit EGARCH
egarch_fit <- ugarchfit(spec = spec_egarch, data = R)

# GJR-GARCH Specification
spec_gjr <- ugarchspec(variance.model = list(model = "gjrGARCH"), 
                       mean.model = list(armaOrder = c(0, 0)),
                       distribution.model = "std")

# Fit GJR-GARCH
gjr_fit <- ugarchfit(spec = spec_gjr, data = R)
# Calculating VaR for GARCH(1,1)
garch_volatility <- sigma(garch_fit)# Extract conditional volatility from the GARCH model
confidence_level <- 0.05# Define confidence level
z_alpha <- qnorm(confidence_level)
mu <- fitted(garch_fit) # Fitted mean
var_garch <- mu + z_alpha * garch_volatility # Calculate VaR (VaR = μ + z_alpha * σ)
# Calculating VaR for EGARCH
garch_volatility <- sigma(egarch_fit)# Extract conditional volatility from the GARCH model
confidence_level <- 0.05# Define confidence level
z_alpha <- qnorm(confidence_level)
mu <- fitted(egarch_fit) # Fitted mean
var_garch <- mu + z_alpha * garch_volatility # Calculate VaR (VaR = μ + z_alpha * σ)
# Calculating VaR for GJR-GARCH
garch_volatility <- sigma(gjr_fit)# Extract conditional volatility from the GARCH model
confidence_level <- 0.05# Define confidence level
z_alpha <- qnorm(confidence_level)
mu <- fitted(gjr_fit) # Fitted mean
var_garch <- mu + z_alpha * garch_volatility # Calculate VaR (VaR = μ + z_alpha * σ)

#Part 5: 
# Backtest GARCH VaR
violations <- sum(R < var_garch) / length(R)
# Expected violation rate
expected_violations <- confidence_level
cat("Actual Violation Rate: ", violations, "\n")
cat("Expected Violation Rate: ", expected_violations, "\n")
# Backtest using PerformanceAnalytics
VaRTest(confidence_level, actual = R, VaR = var_garch)


