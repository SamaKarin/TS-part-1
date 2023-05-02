# Load necessary libraries


install.packages('car')
install.packages("lmtest")
install.packages("forecast")
install.packages("AICcmodavg")
library(forecast)
library(tseries)
library(AICcmodavg)
library(lmtest)
library(car)
library(FinTS)
library(fma)

# Get the database

mydata <- read.csv("C:/Users/flavi/OneDrive/Bureau/L3/R/Question 1 Project TS.csv", sep = ";")

# Convert the data to a time series object

ts_data <- ts(mydata$Romania.Total.Population, start = c(1960,frequency = 1))
ts_data
summary(ts_data) # Some information about the Time Series

# Plot the Time Series

plot(ts_data,
     type='l',
     ylab="Total Population in Romania",
     xlab="Year",
     col="blue",
     main="Total Population of Romania (1960 to 2021)",
     xlim=c(1960,2021),
     ylim=c(18406905,23201835))
grid()
legend("bottomright", legend = "Population in Millions", col="blue",lty = 1,bty = "n")


# Correlogram

sacf <- acf(ts_data, lag.max = 60, main='SACF of the Time Series') 
spacf <-pacf(ts_data, lag.max = 61, main='SPACF of the Time Series ')

# Q stat and associated P-values from lag 1 up to lag 61

q_stat <- numeric(61)
p_value <- numeric(61)

for (i in 1:61) {
  q_stat[i] <- Box.test(ts_data, lag = i, type = "Ljung-Box")$statistic
  p_value[i] <- Box.test(ts_data, lag = i, type = "Ljung-Box")$p.value
}

# Table with all SAcF, SPACF, Q-stat and Prob

table <- data.frame(SACF = sacf$acf, SPACF=spacf$acf, Q=q_stat, Prob=p_value)
print(table)

# Perform unit root tests to check for stationarity (Augmented Dickey Fuller, Philips Peron, KPSS)

pp_test <- pp.test(ts_data) # Philips Peron, H0 <- stationary TS ; H1 <- nonstationary TS
adf_test <- adf.test(ts_data) # Augmented Dickey Fuller, H0 <- stationary TS ; H1 <- nonstationary TS
kpss_test <- kpss.test(ts_data) # KPSS , H0 <- nonstationary TS, H1 <- stationary TS

data.frame(ProbADF=adf_test$p.value,
           ProbPP=pp_test$p.value,
           ProbKPSS=kpss_test$p.value)

# If the time series is non-stationary, perform transformations to achieve stationarity

if (adf_test$p.value > 0.05 | kpss_test$p.value < 0.05) {
  diff_ts_data <- diff(ts_data)
}

# Plot the Population(-1)

plot(diff_ts_data,
     ylab="Population(-1)",
     xlab="Year",
     col="blue",
     main="Differenced time series of the total Romanian Population (1960-2021)")
grid()
legend("bottomleft",legend = "D(population)", col = "blue", lty = 1, bty = "n")

# Check for Unit Root of the differenced Time series (Augmented Dickey Fuller, Philips Peron, KPSS)

pp_test2 <- pp.test(diff_ts_data)
adf_test2 <- adf.test(diff_ts_data)
kpss_test2 <- kpss.test(diff_ts_data)

data.frame(ProbPP=pp_test2$p.value,
           ProbADF=adf_test2$p.value,
           ProbKPSS=kpss_test2$p.value)

# If the differenced time series is non-stationary, perform transformations to achieve stationarity

if(adf_test2$p.value > 0.05 | kpss_test2$p.value < 0.05) {diff_ts_data2 <- diff(diff_ts_data)}

# Plot ∇^2 Pop

plot(diff_ts_data2,
     ylab="Population(-2)",
     xlab="Years",
     col="blue",
     main="The Time series is integrated of order 2")

grid()
legend("bottomleft", legend = "D(D(Population))", col = "blue", lty = 1, bty="n")

# Check for Unit root of the twice differenced time series
# TRY LOG
pp_test3 <- pp.test(diff_ts_data2)
adf_test3 <- adf.test(diff_ts_data2)
kpss_test3 <- kpss.test(diff_ts_data2)

data.frame(ProbPP=pp_test3$p.value,ProbADF=adf_test3$p.value,ProbKPSS=kpss_test3$p.value)

# The Time Series is Integrated of order 2 (The 2nd difference is stationnary)

# Step 1 : Identify p and q for ARIMA(p,d=2,q)

# Correlogram of the SPAC and SPACF of D(D(Population)) to find ARIMA(p,2,q)

sacfDD <- acf(diff_ts_data2, lag.max = 10, main='SACF Correlogram of D(D(Population))', )
spacfDD <-pacf(ts_data, lag.max = 10, main='PACF of the Time Series ')

par(mfrow = c(1, 2))

plot(sacfDD, main = "SACF Correlogram of D(D(Pop))", ylab = "SACF", xlab = "Lag", 
     xlim = c(0,10), ylim = c(-0.5, 1))
axis(1, at = 1:10, labels = 1:10)
axis(2, at = seq(-1, 1, 0.2))

grid(lty = "dotted")

plot(spacfDD, main = "SPACF Correlogram of D(D(Pop))", ylab = "SPACF", xlab = "Lag", 
     xlim = c(0,10), ylim = c(-0.5, 1))
axis(1, at = 1:10, labels = 1:10)
axis(2, at = seq(-1, 1, 0.2))
grid(lty = "dotted")


# max p = 1
# max q = 2
# Possible choices :
#  
# p/q            0                         1                       2
# 0          ARMA(0,0)            ARMA(0,1) <=> MA(1)      ARMA(0,2) <=> MA(2)
# 1       ARMA(1,0) <=> AR(1)         ARMA(1,1)               ARMA(1,2)

# Identify the best model using the Box-Jenkins methodology
auto_arima_model <- auto.arima(ts_data,
                               d=2,
                               max.p = 1,
                               max.q = 2,
                               trace = TRUE,
                               ic = c("aicc", "aic", "bic"))


# Check the validity of the identified model

checkresiduals(auto_arima_model)  # Serial Correlation, maybe heteroskedasticity...

# Verify if the function auto.arima gave us a correct output.

arima021 <- arima(ts_data, c(0,2,1))
arima022 <- arima(ts_data, c(0,2,2))
arima120 <- arima(ts_data, c(1,2,0))
arima121 <- arima(ts_data, c(1,2,1))
arima122 <- arima(ts_data, c(1,2,2))


# Create a 5x5 data frame to compare results

my_df <- data.frame(matrix(nrow=3, ncol=3))

# Set custom column names

row.names(my_df) <- c("AIC", "BIC", "LL")
# Set custom row names

colnames(my_df) <- c("ARMA(1,0)", "ARMA(0,1)", "ARMA(1,1)", "ARMA(0,2)", "ARMA(1,2)")

summary(arima122)

#Input values into table

my_df[,2] <- c(arima021$aic,BIC(arima021),arima021$loglik)
my_df[,1] <- c(arima120$aic,BIC(arima120),arima120$loglik)
my_df[,3] <- c(arima121$aic,BIC(arima121),arima121$loglik)
my_df[,4] <- c(arima022$aic,BIC(arima022),arima022$loglik)
my_df[,5] <- c(arima122$aic,BIC(arima122),arima122$loglik)

my_df

testarima021 <- coeftest(arima021)
testarima022 <- coeftest(arima022)
testarima120 <- coeftest(arima120)
testarima121 <- coeftest(arima121)
testarima122 <- coeftest(arima122)

print(testarima021) #MA 1 valid in ARIMA(0,2,1)
print(testarima022) # MA1not valid in ARIMA(0,2,2)
print(testarima120) # AR1 not valid in ARIMZ(1,2,0)
print(testarima121) #AR1 valid in ARIMA(1,2,1)
print(testarima122) #MA 1 valid at 10% and AR1 not valid in ARIMA(1,2,2) # In general the model is valid


# Test ARIMA(1,2,1)

par(mfrow = c(1, 2))
acf(ts(arima121$residuals))
pacf(ts(arima121$residuals))
Box.test(arima121$residuals, 
         lag = 9,
         type = "Ljung-Box")

# Test ARIMA(1,2,2) for serial correlation

par(mfrow = c(1, 2))
acf(ts(arima122$residuals))
pacf(ts(arima122$residuals))
Box.test(arima122$residuals, 
         lag = 10,
         type = "Ljung-Box") # Absence of autocorrelation

# Are the residuals white noise ?

hist(arima122$residuals)
jarque.bera.test(arima122$residuals) # Not Normally distributed

# Homoskedasticity ? 

ArchTest(arima122$residuals, lags = 11) # Yes

#Go for ARIMA(1,2,2)

modelPop <- arima122

# Forecast

forecast_Pop <- forecast(modelPop,
                        level = c(95),
                        h=15)

summary(forecast_Pop)

forecasttable <- data.frame(forecast_Pop)
colnames(forecasttable) <- c("Population Forecast", "Lower 95%", "Upper 95%")

forecasttable 

# Plot the forecast

plot(forecast_Pop,
     col = "darkred",
     type = "l",
     xlim=c(1960,2036))
grid()

# Test the forecast

Box.test(forecast_Pop$residuals,
         lag = 5,
         type = "Ljung-Box") # No significant autocorrelation

hist(forecast_Pop$residuals)

jarque.bera.test(forecast_Pop$residuals) # Not Normal residuals

ArchTest(forecast_Pop$residuals,
         lags=8) # Constant variance

# Inverse Roots of polynomial

# Compute AR roots
arroots <- function(object)
{
  if(class(object) != "Arima" & class(object) != "ar")
    stop("object must be of class Arima or ar")
  if(class(object) == "Arima")
    parvec <- object$model$phi
  else
    parvec <- object$ar
  if(length(parvec) > 0)
  {
    last.nonzero <- max(which(abs(parvec) > 1e-08))
    if (last.nonzero > 0)
      return(structure(list(roots=polyroot(c(1,-parvec[1:last.nonzero])),
                            type="AR"), class='armaroots'))
  }
  return(structure(list(roots=numeric(0),type="AR"),class='armaroots'))
}

# Compute MA roots
maroots <- function(object)
{
  if(class(object) != "Arima")
    stop("object must be of class Arima")
  parvec <- object$model$theta
  if(length(parvec) > 0)
  {
    last.nonzero <- max(which(abs(parvec) > 1e-08))
    if (last.nonzero > 0)
      return(structure(list(roots=polyroot(c(1,parvec[1:last.nonzero])),
                            type="MA"), class='armaroots'))
  }
  return(structure(list(roots=numeric(0),type="MA"),class='armaroots'))
}

plot.armaroots <- function(x, xlab="Real",ylab="Imaginary",
                           main=paste("Inverse roots of",x$type,"characteristic polynomial"),
                           ...)
{
  oldpar <- par(pty='s')
  on.exit(par(oldpar))
  plot(c(-1,1),c(-1,1),xlab=xlab,ylab=ylab,
       type="n",bty="n",xaxt="n",yaxt="n", main=main, ...)
  axis(1,at=c(-1,0,1),line=0.5,tck=-0.025)
  axis(2,at=c(-1,0,1),label=c("-i","0","i"),line=0.5,tck=-0.025)
  circx <- seq(-1,1,l=501)
  circy <- sqrt(1-circx^2)
  lines(c(circx,circx),c(circy,-circy),col='gray')
  lines(c(-2,2),c(0,0),col='gray') 
  lines(c(0,0),c(-2,2),col='gray')
  if(length(x$roots) > 0) {
    inside <- abs(x$roots) > 1
    points(1/x$roots[inside],pch=19,col='black')
    if(sum(!inside) > 0)
      points(1/x$roots[!inside],pch=19,col='red')
  }
}

plot(arroots(modelPop),
     main = "Inverse roots of the AR polynomial") # Lie inside the unit circle

# Compute MA roots

plot(maroots(modelPop),
     main = "Inverse roots of the MA polynomial") # Inside the unit circle 

# The process is stationary and invertible


