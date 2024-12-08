library(maxLik)
library(purrr)
library(resampledata)

#### ROZKŁAD NORMALNY ####
# estymowane parametry - mu, sigma

NormLogLik <- function(x, param){
  mu <- param[1]
  sigma <- param[2]
  -0.5*length(x)*log(2*pi) - length(x)*log(sigma) - 0.5*sum((x-mu)^2/sigma^2)
}

NormGrad <- function(x, param){
  mu <- param[1]
  sigma <- param[2]
  result <- numeric(2)
  result[1] <- sum((x-mu)/sigma^2)
  result[2] <- -length(x)/sigma + sum((x-mu)^2/sigma^3)
  return(result)
}

NormHess <- function(x, param){
  mu <- param[1]
  sigma <- param[2]
  result <- matrix(nrow = 2, ncol = 2)
  result[1,1] <- -length(x)/sigma^2
  result[1,2] <- -2*sum((x-mu)/sigma^3)
  result[2,1] <- -2*sum((x-mu)/sigma^3)
  result[2,2] <- length(x)/sigma^3 - 3*sum((x-mu)^2/sigma^4)
  return(result)
}

#### ROZKŁAD GEOMETRYCZNY ####
# estymowany parametr: p

GeomLogL <- function(k, p){
  log(1 - p) * (sum(k) - length(k)) + length(k) * log(p)
}

GeomGrad <- function(k, p){
  (length(k) - sum(k)) / (1 - p) + length(k) / p
}

GeomHess <- function(k, p){
  (length(k) - sum(k)) / (1 - p)^2 - length(k) / p^2
}


#### Zadanie 1 ####

x <- c(-1.13,-0.62,0.06,-1.82,-0.27,-0.74,0.49,-0.35,-0.41,-0.05,-0.20,0.06,
       -1.06,0.58,1.20, 1.59,0.20,0.65,-0.53,-0.73,-1.16,1.34,1.77,0.70,-0.10,
       -0.30,0.60,0.49,-0.45,-1.15,-0.31,1.30,0.91,0.47,-0.44,1.23,-0.15,1.22,
       -1.26,-0.18)

maxLik(logLik = function(mu) {NormLogLik(x = x, c(mu, 1))}, 
       start = 0)

mu <- 0.03625 # wynik estymacji MNW

# Test istotności - model I (wariancja znana)
# Statystyka testowa: U
# hipoteza zerowa: mu = 0

U <- (mu - 0) / 1 * sqrt(length(x))
pnorm(U)
p.value <- 2*(1 - pnorm(U)) 
# p-value = 0.8187 zatem nie ma podstaw do odrzucenia H0


#### Zadanie 2 ####

k <- c(0,0,3,1,3,11,0,7,2,3,4,2,0,1,2,1,3,2,4,2,1,4,6,4,0,6,0,0,4,2)

maxLik(logLik = partial(GeomLogL, k = k), 
       grad = partial(GeomGrad, k = k), 
       hess = partial(GeomHess, k = k), 
       start = 0.1)

p <- 0.3846154

# porównanie z estymatorem metodą momentów: 

p.mm <- 1 / (mean(k) + 1)

# test istotności
# Do testowania istotności wykorzystam fakt, że estymator MNW ma rozkład normalny
# Statystyka testowa: U
# hipoteza zerowa: p = 0.25

U <- (p - 0.25) / sqrt(0.25*(1 - 0.25) / length(k))
p.value <- 2*(1 - pnorm(U)) 

# p-value = 0.08861204
# NIe można odrzucić H0 (p=0.25) na poziomie istotności 5%


#### Zadanie 3 ####

ratings.df <- read.csv2("ratings_Musical_Instruments.csv")

PoissLogLik <- function(k, lambda) {
  N <- length(k)
  L <- -N*lambda - N*log(1-exp(-lambda)) + log(lambda) * sum(k)
  return(L)
}

maxLik(logLik = partial(PoissLogLik, k = ratings.df$Ratings), 
       start = 6)

lambda.est <- 6.00807

# Do testowania istotności wykorzystam fakt, że estymator MNW ma rozkład normalny
# odchylenie standardowe policzyłem jako odwrotność hesjanu. 
# U - statystyka testowa

n <- length(ratings.df$Ratings)
sd <- n * (1 / lambda.est - exp(lambda.est) / (exp(lambda.est) - 1)^2)
U <- (lambda.est - 6) / sd 
p.value <- 2*(1 - pnorm(U)) 
# p.value = 0.9999995
# nie ma podstaw do odrzucenia hipotezy zerowej

# Test Likelihood Ratio
LR <- 2*(PoissLogLik(ratings.df$Ratings, lambda.est) - 
           PoissLogLik(ratings.df$Ratings, 6))

critical.value <- qchisq(0.95, 1)
LR > critical.value
# nie ma podstaw do odrzucenia hipotezy zerowej


#### Zadanie 5 ####

Service

GammaLogLik <- function(x, param){
  k <- param[[1]]
  theta <- param[[2]]
  N <- length(x)
  (k-1)*sum(log(x)) - sum(x)/theta - N*k*log(theta) - N*log(gamma(k))
}

maxLik(partial(GammaLogLik, x = Service$Times), start = c(1,1))

k.est <- 2.813057
theta.est <- 0.247035

# Wizualizacja rzeczywistego rozkładu (histogramu) oraz 
# dopasowanej funkcji gęstości (czerwona linia) z parametrami MNW.
# Zielona linia przedstawia funkcję gęstości z parametrami hipotezy zerowej

hist(Service$Times, freq = FALSE)
lines(x = seq(0,2.5,0.01), 
      y = dgamma(seq(0,2.5,0.01), shape = k.est, scale = theta.est), 
      col = "red", 
      lwd = 2)
lines(x = seq(0,2.5,0.01), 
      y = dgamma(seq(0,2.5,0.01), shape = 1, scale = 2), 
      col = "green", 
      lwd = 2)

# Test wykonam metodą Likelihood Ratio
LR <- 2*(GammaLogLik(Service$Times, c(k.est, theta.est)) - 
           GammaLogLik(Service$Times, c(1, 2)))

critical.value <- qchisq(0.95, 2)
LR > critical.value
# Możemy odrzucić hipotezę zerową



#### Zadanie 8 ####

contest.df <- read.csv("contest_data.csv")

TruncNormLogDensity <- function(x, mu, sigma){
  if (x < 0) {
    return(0)
  } else if (x == 0) {
    return(log(1 - pnorm(mu/sigma, mean = 0, sd = 1)))
  } else {
    return(-0.5*log(2*pi) - log(sigma) - 0.5*(x-mu)^2/sigma^2)
  }
}

TruncNormLogLik <- function(x, param){
  mu <- param[1]
  sigma <- param[2]
  return(sum(sapply(x, TruncNormLogDensity, mu = mu, sigma = sigma)))
}

maxLik()
