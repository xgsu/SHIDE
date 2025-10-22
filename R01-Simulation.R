
# rm(list = ls(all.names = TRUE))
source("Functions-SHIDE.R")

# ============
# I. NORMAL 
# ============

set.seed(666)
N <- c(50, 500)
nrun <- 300
mu0 <- 0; sigma0 <- 1   #  normal
m <- 10; c0 <- 1; nclass <- NULL  # SHIDE related 
MISE <- NULL
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    # generate data from GMM
    x <- rnorm(n, mean=mu0, sd=sigma0)
    # hist(x, nclass=20); 
    
    # 1. KDE
    d.KDE <- density(x, bw = "SJ", n = 1024)
    x0 <- d.KDE$x
    # 2(a) SHIDE - opt
    fit.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "opt", 
                       m=m, c0=c0, nclass=nclass) 
    d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
    # 2(b) SHIDE - per
    fit1.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "per", 
                       m=m, c0=c0, nclass=nclass) 
    d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)    

    # True Density
    d.true <- list(x=x0, 
                   y=dnorm(x=x0, mean=mu0, sd=sigma0))
    
    MISE <- rbind(MISE, c(model="normal", n=n, run=i, 
                          mise.KDE=compute_mise(d.true, d.KDE),
                          mise.SHIDE.opt=compute_mise(d.true, d.SHIDE),
                          mise.SHIDE.per=compute_mise(d.true, d1.SHIDE)
    ))
  }
}
MISE <- as.data.frame(MISE); MISE
Summarize.MISE(MISE)
MISE.normal <- MISE





# ==================================
# II. GAUSSIAN MIXTURE MODEL (GMM)
# ==================================

set.seed(321)
N <- c(50, 500)
nrun <- 300
pi=0.35; mu1=-1;  sd1=1; mu2=2; sd2=2; # GMM
m <- 10; c0 <- 1; nclass <- NULL  # SHIDE related 
MISE <- NULL
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    # generate data from GMM
    x <- rmixnorm2(n=n, pi=pi, mu1=mu1, sd1=sd1, mu2=mu2, sd2=sd2)
    # hist(x); hist(x, nclass=20); 
    
    # 1. KDE
    d.KDE <- density(x, bw = "SJ", n = 1024)
    x0 <- d.KDE$x
    # 2(a) SHIDE - opt
    fit.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "opt", 
                       m=m, c0=c0, nclass=nclass) 
    d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
    # 2(b) SHIDE - per
    fit1.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "per", 
                        m=m, c0=c0, nclass=nclass) 
    d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)  
    # True Density
    d.true <- list(x=x0, 
        y=dmixnorm2(x=x0, pi=pi, mu1=mu1, sd1=sd1, mu2=mu2, sd2=sd2))
    
    MISE <- rbind(MISE, c(model="GMM", n=n, run=i, 
                          mise.KDE=compute_mise(d.true, d.KDE),
                          mise.SHIDE.opt=compute_mise(d.true, d.SHIDE),
                          mise.SHIDE.per=compute_mise(d.true, d1.SHIDE)
                          ))
  }
}
MISE <- as.data.frame(MISE)
Summarize.MISE(MISE)
MISE.GMM <- MISE



# ==========================
# III. Cauchy Distribution 
# ==========================

set.seed(789)
N <- c(50, 500)
nrun <- 300
location = 0; scale = 1  #  Cauchy
m <- 10; c0 <- 1; nclass <- NULL  # SHIDE related 
MISE <- NULL
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    # generate data from GMM
    x <- rcauchy(n, location=location, scale=scale)
    # hist(x, nclass=20); 
    
    # 1. KDE
    d.KDE <- density(x, bw = "SJ", n = 1024)
    x0 <- d.KDE$x
    # 2(a) SHIDE - opt
    fit.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "opt", 
                       m=m, c0=c0, nclass=nclass) 
    d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
    # 2(b) SHIDE - per
    fit1.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "per", 
                        m=m, c0=c0, nclass=nclass) 
    d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)  
    # True Density
    d.true <- list(x=x0, 
                   y=dcauchy(x=x0, location=location, scale=scale))
    
    MISE <- rbind(MISE, c(model="Cauchy", n=n, run=i, 
                          mise.KDE=compute_mise(d.true, d.KDE),
                          mise.SHIDE.opt=compute_mise(d.true, d.SHIDE),
                          mise.SHIDE.per=compute_mise(d.true, d1.SHIDE)
    ))
  }
}
MISE <- as.data.frame(MISE); MISE
Summarize.MISE(MISE)
MISE.Cauchy <- MISE




# ===============================================
# IV. Exponential Distribution (positive case)
# ===============================================

set.seed(123)
N <- c(50, 500)
nrun <- 300
rate0=1; # Exponential rate
m <- 10; c0 <- 1; nclass <- NULL  # SHIDE related 
MISE <- NULL
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    # generate data from GMM
    x <- rexp(n=n, rate=rate0)
    # hist(x); hist(x, nclass=20); 
    
    # 1. KDE
    d.KDE <- density(x, bw = "SJ", n = 1024)
    x0 <- d.KDE$x
    d.KDE$y[x0 <0] <- 0
    # 2(a) SHIDE - opt
    fit.shide <- shide(x, LB=0, UB=NULL, method.bound = "opt", 
                       m=m, c0=c0, nclass=nclass) 
    d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
    # 2(b) SHIDE - per
    fit1.shide <- shide(x, LB=0, UB=NULL, method.bound = "per", 
                        m=m, c0=c0, nclass=nclass) 
    d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)  
    # True Density
    d.true <- list(x=x0, y=dexp(x=x0, rate=rate0))
    
    MISE <- rbind(MISE, c(model="exp", n=n, run=i, 
              mise.KDE=compute_mise(d.true, d.KDE, lower_lim=0),
                mise.SHIDE.opt=compute_mise(d.true, d.SHIDE, lower_lim=0),
                mise.SHIDE.per=compute_mise(d.true, d1.SHIDE, lower_lim=0)
    ))
  }
}
MISE <- as.data.frame(MISE)
Summarize.MISE(MISE)
MISE.exp <- MISE





# =============================================
# V. Bounded Distribution - Truncated Normal
# =============================================

set.seed(555)
N <- c(50, 500)
nrun <- 300
a=-1; b=0.5; mu0=0; sigma0=3  # truncated normal
m <- 10; c0 <- 1; nclass <- NULL  # SHIDE related 
MISE <- NULL
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    # generate data from GMM
    x <- rtruncnorm(n, a = a, b = b, mean = mu0, sd = sigma0)
    # hist(x); hist(x, nclass=20); 
    
    # 1. KDE
    d.KDE <- density(x, bw = "SJ", n = 1024)
    x0 <- d.KDE$x
    d.KDE$y[x0 <a] <- 0;  d.KDE$y[x0 > b] <- 0
    
    # 2(a) SHIDE - opt
    fit.shide <- shide(x, LB=a, UB=b, method.bound = "opt", 
                       m=m, c0=c0, nclass=nclass) 
    d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
    # 2(b) SHIDE - per
    fit1.shide <- shide(x, LB=a, UB=b, method.bound = "per", 
                        m=m, c0=c0, nclass=nclass) 
    d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)  
    # True Density
    d.true <- list(x=x0, 
            y=dtruncnorm(x0, a = a, b = b, mean = mu0, sd = sigma0))
    
    MISE <- rbind(MISE, c(model="truncated.Normal", n=n, run=i, 
        mise.KDE=compute_mise(d.true, d.KDE, lower_lim=a, upper_lim =b),
        mise.SHIDE.opt=compute_mise(d.true, d.SHIDE, lower_lim=a, upper_lim =b),
        mise.SHIDE.per=compute_mise(d.true, d1.SHIDE, lower_lim=a, upper_lim =b)
    ))
  }
}
MISE <- as.data.frame(MISE)
Summarize.MISE(MISE)
MISE.bounded <- MISE





# ====================
# SUMMARIZE RESULTS
# ====================

MISE <- rbind(MISE.normal, MISE.GMM, MISE.Cauchy, MISE.exp, MISE.bounded)

MISE_summary <- Summarize.MISE(MISE)

library(knitr)

kable(MISE_summary, 
      format = "latex",
      caption = "MISE Summary Statistics",
      digits = 6,
      booktabs = TRUE,
      col.names = c("Model", "method", "median.n50",
                    "mad.n50", "median.n500", "mad.n500"))







# --------------------------------------
# FURTHER EXPLORATION OF CAUCHY DIST
# --------------------------------------

set.seed(7)
x <- rcauchy(500, location=0, scale=1)
# 1. KDE
d.KDE <- density(x, bw = "SJ", n = 1024)
x0 <- d.KDE$x
# 2(a) SHIDE - opt
fit.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "opt", nclass=100) 
d.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit.shide)
# 2(b) SHIDE - per
fit1.shide <- shide(x, LB=NULL, UB=NULL, method.bound = "per", nclass=100) 
d1.SHIDE <- predict.SHIDE(x.new=x0, fit.shide=fit1.shide)  
# True Density
d.true <- list(x=x0, 
               y=dcauchy(x=x0, location=location, scale=scale))

c(mise.KDE=compute_mise(d.true, d.KDE),
  mise.SHIDE.opt=compute_mise(d.true, d.SHIDE),
  mise.SHIDE.per=compute_mise(d.true, d1.SHIDE))
    
hist(x, plot=TRUE, probability =TRUE, col=gray(level = 0.5, alpha = 0.3), 
     ylim=c(0, 3.0), xlim=c(-80, 80), nclass=200)
lines(d.true, col="green4", lty=1, lwd=2)
lines(d.KDE, col="red", lty=1)
lines(d.SHIDE, col="blue", lty=1, lwd=1)
# lines(d1.SHIDE, col="blue", lty=1, lwd=1)







#