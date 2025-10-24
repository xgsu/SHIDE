
#############################################
#############################################
#  FUNCTIONS IN SHIDE 
# (SIMUALTION AND HISTOGRAM INTERPOLATION
# FOR DENSITY ESTIMATION)
#############################################
#############################################


require(tidyverse)
library(data.table)
library(truncnorm)   # Truncated normal dist
library(pracma)   # For numerical integration

expit <- function(x) (tanh(x/2)+1)/2  	# VERIFIED
# expit0 <- function(x) exp(x)/(1+exp(x))
logit <- function(x) 2*atanh(2*x-1) # VERIFIED
# logit0 <- function(x) {log(x) - log1p(-x)} 
# x0 <- 1:10/100
# cbind(x0, logit(x0), logit0(x0))


# ===================================================================
# FUNCTION bound() COMPUTES BOUND OR BANDWIDTH FOR FUNCTION rNoise()
# ===================================================================
bound <- function(x, alpha=0.20) {
  d0 <- sort(diff(sort(x)))
  b <- quantile(d0, probs=alpha)
  return(b)
}


# ===================================================================
# Function bound.opt() computes OPTIMAL BOUND OR BANDWIDTH IN SHIDE
# ===================================================================
# Two approaches implemented; 
# k >= 3 is one necessary condition for theoretical analysis

bound.opt <- function(x, k=3, c=1, 
                      alpha = 0.5,            # using median spacing
                      method.Psi=c("pilot.density", "normal.approx"))
{
  n <- length(x)
  # get pilot density estimates
  pilot <- density(x, bw="nrd0")  # Gaussian kernel pilot with normal reference bw
  
  # h.opt: MINIMIZING MISE
  # -------------------------------
  R.K <- (k/2)*choose(2*k, k) / (4^k) 
  sigmaK <- 1/sqrt(3*k) 
  # normal approaximation for Psi.hat
  method.Psi <- ifelse(is.null(method.Psi), "pilot.density", method.Psi)
  if (method.Psi=="normal.approx") {
    # Silverman-type normal-reference estimate of Psi(f'')
    # s <- sd(x)                  
    s <- IQR(x)/1.349 # more robust
    Psi_hat <- (3 / (8 * sqrt(pi))) * s^(-5)
  } else if (method.Psi=="pilot.density") {
    # Estimate Psi via numeric differentiation of pilot density
    x_grid <- seq(min(x), max(x), length.out = 1000)
    # second derivative by finite difference on pilot$y or using smooth.spline on pilot for higher accuracy
    spline_fit <- smooth.spline(pilot$x, pilot$y, spar=0.6)    # fit a spline through pilot density for smooth derivative
    f2_vals <- predict(spline_fit, x_grid, deriv = 2)$y
    Psi_hat <- trapz(x_grid, f2_vals^2)        # integrate using trapezoidal rule (requires library(pracma) or implement manually)
  } else stop("Wrong entry for method.Psi.")
  const <- ((R.K * (1 + 1/c)) / (sigmaK^4 * Psi_hat))^(1/5)  
  h.opt <- const * n^(-1/5) 
  
  # h_perc: NEAREST NEIGHBOR SPACING PERCENTILE
  # --------------------------------------------- 
  spacings <- diff(sort(x))
  d_alpha <- quantile(spacings, probs=alpha)
  f_hat_xalpha <- approx(pilot$x, pilot$y, xout = quantile(x, probs=alpha))$y 
  q_alpha <- -log(1 - alpha)
  lambda_n <- n^(4/5) * ((R.K*(1+1/c))/(sigmaK^4 * Psi_hat))^(1/5) * (f_hat_xalpha / q_alpha)
  h.perc <- lambda_n * d_alpha
  details <- list(R.K=R.K, Psi.hat=Psi_hat, const=const, q_alpha=q_alpha, lambda_n=lambda_n) 
  return(list(h.opt=h.opt, h.perc=h.perc, details=details))
}




# ===========================================================
# FUNCTION rNoise() GENERATES NOISES WITH POLYNOMIAL KERNEL 
# DENSITY AND BOUNDED SUPPORT
# ===========================================================
rNoise <- function(n, k=2, bound=1){
  U <- matrix(runif(n*k, min = -1/2, max=1/2), nrow =n, ncol=k) 
  U <- rowSums(U)   # BOUND [-k/2, k/2]
  U <- U*bound/(k/2) 
  return(U)
}


# ===============================================================================
# FUNCTION shide() ESTIMATES DENSITY BASED ON SIMULATION AND HISTOGRAM SMOOTHING
# ===============================================================================
# k >= 3 is one necessary condition for theoretical analysis

shide <- function(x, LB=NULL, UB=NULL, 
                  m=10, k=3,  epsilon=1e-10, nonnegativity.correction=TRUE,  
                  method.bound = c("opt", "perc"), c0=1, alpha=0.5,
                  plot.it=FALSE, nclass=NULL,
                  bw = "SJ", kernel = "epanechnikov",    # KDE OPTIONS
                  ...) {
  if (!is.vector(x)) stop("x must be a vector.")
  # A QUICK LOOK AT THE HISTOGRAM AS WELL AS KDE 
  if (plot.it) {
    hist(x, prob = TRUE, col = gray(level=0.75, alpha=0.2),
         nclass=nclass,
         main = "Histogram with KDE density curve")
    lines(density(x, bw = bw, adjust = 1, kernel = kernel), 
          col = "red", lwd = 2) # KDE
  }
    
  # SHIDE
  LB <- ifelse(is.null(LB), -Inf, LB)
  UB <- ifelse(is.null(UB), Inf, UB)
  if (LB > UB) stop("LB should be no greater than UB.")
  # x0 <- x
  n <- length(x)
  if (is.infinite(LB) && is.infinite(UB)) {
    if (!is.null(method.bound)) {
      Bound <- bound.opt(x, k=k, c=c0, alpha = alpha)
      bd <- ifelse(method.bound=="perc", Bound$h.perc, Bound$h.opt)
    } else bd <- bound(x, alpha=alpha) 
    x0 <- rep(x, m) + rNoise(n*m, k=k, bound=bd)
  } else if (is.finite(LB) && is.infinite(UB)) {
    x1 <- log(x - LB + epsilon) 
    if (!is.null(method.bound)) {
      Bound <- bound.opt(x1, k=k, c=c0, alpha = alpha)
      bd <- ifelse(method.bound=="perc", Bound$h.perc, Bound$h.opt)
    } else bd <- bound(x1, alpha=alpha) 
    x2 <- rep(x1, m) + rNoise(n*m, k=k, bound=bd)
    x0 <- exp(x2) + LB - epsilon 
  } else if (is.infinite(LB) && is.finite(UB)) {
    x1 <- log(UB-x + epsilon) 
    if (!is.null(method.bound)) {
      Bound <- bound.opt(x1, k=k, c=c0, alpha = alpha)
      bd <- ifelse(method.bound=="perc", Bound$h.perc, Bound$h.opt)
    } else bd <- bound(x1, alpha=alpha) 
    x2 <- rep(x1, m) + rNoise(n*m, k=k, bound=bd)
    x0 <- UB - exp(x2) + epsilon 
  } else if (is.finite(LB) && is.finite(UB)) {
    x1 <- (x-LB+epsilon) / (UB-LB)
    x1 <- logit(x1) 
    if (!is.null(method.bound)) {
      Bound <- bound.opt(x1, k=k, c=c0, alpha = alpha)
      bd <- ifelse(method.bound=="perc", Bound$h.perc, Bound$h.opt)
    } else bd <- bound(x1, alpha=alpha) 
    x2 <- rep(x1, m) + rNoise(n*m, k=k, bound=bd)
    x0 <- expit(x2)*(UB-LB) + LB - epsilon 
  } else stop("Hmmm. How did you get here?")
  
  # PARALLEL BOXPLOTS OF OBSERVED AND SIMULATED DATA
  fig.boxplot <- NULL
  if (plot.it) {
    df0 <- data.frame(y=c(x, x0), 
                      group=rep(c("Obsered", "Simulated"), c(n, n*m)))
    fig.boxplot <- ggplot(df0, aes(x = group, y = y,
                  colour = group,
                  shape = group)) + 
    geom_boxplot(outlier.shape = NA, notch=TRUE, varwidth=TRUE) +
    geom_jitter(cex = 0.8, width=0.2)
  }

  # HISTOGRAM OF SIMULATED DATA
  hist.out <- suppressWarnings(hist(x0, plot=plot.it, probability=TRUE, nclass=nclass, 
                   col=gray(level = 0.5, alpha = 0.3)))
  
  # names(hist.out)
  midpoints <- hist.out$mids
  # counts <- hist.out$counts
  d.hist <- hist.out$density
  if (nonnegativity.correction) d.hist <- sqrt(hist.out$density)
  # print(d.hist)
  # sum(density); cbind(counts/length(x0), density)

  # SMOOTH INTERPOLATION WITH HISTOGRAM

  # lines(midpoints, density, lty=1, col="blue", lwd=0.5) # LINEAR INTERPOLATION
  fit.spline <- spline(x=midpoints, y=d.hist, method = "natural") # n=100) 
  density.function <- splinefun(x=midpoints, y=d.hist, method = "natural")
  # names(fit.spline); cbind(midpoints, fit.spline$x) 
  
  # Adjust density estimates based on specified bounds
  y <- fit.spline$y
  if (nonnegativity.correction) y <- (fit.spline$y)^2
  if (is.finite(LB)) y[fit.spline$x < LB] <- 0; 
  if (is.finite(UB)) y[fit.spline$x > UB] <- 0
  
  if (plot.it) {
    points(midpoints, density, pch=19, cex=0.6, col="blue")
    lines(fit.spline$x, y=y, lty=1, col="green4", lwd=2)
  }

  return(list(data.observed=x, data.simulated=x0, bd=bd, 
              nonnegativity.correction=nonnegativity.correction, 
              fit.spline=fit.spline, x=fit.spline$x, y=y, 
              density.function=density.function, LB=LB, UB=UB, 
              fig.boxplot=fig.boxplot))
}


# ================================================================
# Function predict.SHIDE() estiamtes density at given new x values
# ================================================================
# fit.shide - A SHIDE object
# x.new - new x values

predict.SHIDE <- function(x.new, fit.shide) {
  f.den <- fit.shide$density.function
  nonnegativity.correction <- fit.shide$nonnegativity.correction
  LB <- fit.shide$LB; UB <- fit.shide$UB
  y <- f.den(x.new)
  if (nonnegativity.correction) y <- (f.den(x.new))^2 #不要忘了平方! 
  # bound adjustment
  if (is.finite(LB)) y[x.new < LB] <- 0; 
  if (is.finite(UB)) y[x.new > UB] <- 0
  return(list(x=x.new, y=y))   
}



# ====================================================================
# Function compute_mise() computes MISE between two density estimates
# ====================================================================

compute_mise <- function(true_density, estimated_density, 
                         lower_lim = -Inf, upper_lim = Inf, 
                         n_points = 1000) {
  # Create evaluation grid covering both density supports
  x_min <- min(true_density$x, estimated_density$x)
  x_max <- max(true_density$x, estimated_density$x)
  x_seq <- seq(x_min, x_max, length.out = n_points)
  
  # Interpolate densities on common grid
  f_true <- approx(true_density$x, true_density$y, xout = x_seq)$y
  f_est <- approx(estimated_density$x, estimated_density$y, xout = x_seq)$y
  
  # Handle NA values (outside kernel support)
  f_true[is.na(f_true)] <- 0
  f_est[is.na(f_est)] <- 0
  
  # Apply integration limits
  if (is.finite(lower_lim) || is.finite(upper_lim)) {
    in_range <- x_seq >= lower_lim & x_seq <= upper_lim
    x_seq <- x_seq[in_range]
    f_true <- f_true[in_range]
    f_est <- f_est[in_range]
  }
  
  # Compute squared differences and integrate
  squared_diff <- (f_true - f_est)^2
  mise <- trapz(x_seq, squared_diff)  # Numerical integration
  
  return(mise)
}



# ======================================================
# Functions related to Gaussian mixture models (GMM)
# ======================================================

## Simulator: draw from a 2-component normal mixture 
rmixnorm2 <- function(n, pi, mu1, sd1, mu2, sd2) {
  z <- rbinom(n, size = 1, prob = pi)         # 1 => comp 1, 0 => comp 2
  x <- numeric(n)
  n1 <- sum(z == 1); n2 <- n - n1
  x[z == 1] <- rnorm(n1, mean = mu1, sd = sd1)
  x[z == 0] <- rnorm(n2, mean = mu2, sd = sd2)
  x
}

## True density of the 2-component mixture at points x
dmixnorm2 <- function(x, pi, mu1, sd1, mu2, sd2, log = FALSE) {
  f <- pi * dnorm(x, mean = mu1, sd = sd1) +
    (1-pi)  * dnorm(x, mean = mu2, sd = sd2)
  if (log) return(log(f)) else return(f)
}





# =======================
# SUMMARY OF RESULTS
# ======================

Summarize.MISE <- function(MISE) {
    MISE %>%
      mutate(mise.KDE = as.numeric(mise.KDE), 
             mise.SHIDE.opt = as.numeric(mise.SHIDE.opt),
             mise.SHIDE.per = as.numeric(mise.SHIDE.per)) %>% 
      pivot_longer(
        cols = c(mise.KDE, mise.SHIDE.opt, mise.SHIDE.per),
        names_to = "method",
        values_to = "mise_value"
      ) %>%
      mutate(method = gsub("mise\\.", "", method)) %>%
      group_by(model, n, method) %>%
      summarise(
        median = median(mise_value, na.rm = TRUE),
        mad = mad(mise_value, na.rm = TRUE),
        .groups = 'drop'
      ) %>% 
    pivot_wider(
      id_cols = c(model, method),
      names_from = n,
      values_from = c(median, mad),
      names_glue = "{.value}.n{n}"
    ) %>%
    # Reorder columns for better presentation
    select(model, method, 
           median.n50, mad.n50, 
           median.n500, mad.n500)
}






















#