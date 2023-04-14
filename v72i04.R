library("tseries")
library("midasr")
library("ggplot2")
data("oos_prec", package = "midasr")
qplot(x = n, y = value, data = oos_prec, geom = "line", colour = Constraint, 
  ylab = "") + facet_wrap(~ Type, scales = "free_y") + xlab("Sample size") + 
  theme_bw()

x <- 1:12

fmls(x, k = 2, m = 3)

plot(x = 0:16, nealmon(p = c(2, 0.5, -0.1), d = 17), type = "l", 
  xlab = "High frequency lag", ylab = "Weights", col = 4)
lines(x = 0:7, nealmon(p = c(1, -0.5), d = 8), col = 2)

set.seed(1001)
n <- 250
trend <- c(1:n)
x <- rnorm(4 * n)
z <- rnorm(12 * n)
fn_x <- nealmon(p = c(1, -0.5), d = 8)
fn_z <- nealmon(p = c(2, 0.5, -0.1), d = 17)
y <- 2 + 0.1 * trend + mls(x, 0:7, 4) %*% fn_x + mls(z, 0:16, 
  12) %*% fn_z + rnorm(n)

eq_u <- lm(y ~ trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, 
  m = 12))

eq_u <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), 
  start = NULL)

eq_u <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), 
  start = NULL, data = list(y = y, trend = trend, x = x, z = z))

eq_u <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12), 
  start = NULL, data = list(data.frame(y = y, trend = trend), 
    x = x, z = z))

eq_r <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 
  0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 
  0.5, -0.1)))
summary(eq_r)

## eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + 
##   mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), 
##   z = c(2, 0.5, -0.1)), Ofunction = "optim", method = "Nelder-Mead")

## eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + 
##   mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), 
##   z = c(2, 0.5, -0.1)), Ofunction = "nls", method = "plinear")

## eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + 
##   mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), 
##   z = c(2, 0.5, -0.1)), Ofunction = "optim", method = "Nelder-Mead")
## eq_r2 <- update(eq_r2, Ofunction = "nls")

eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 
  0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 
  0.5, -0.1)), Ofunction = "optim", method = "Nelder-Mead")
eq_r2$opt
eq_r2$convergence

nealmon_gradient <- function(p, d, m) {
  i <- 1:d
  pl <- poly(i, degree = length(p) - 1, raw = TRUE)
  eplc <- exp(pl %*% p[-1])[, , drop = TRUE]
  ds <- colSums(pl * eplc)
  s <- sum(eplc)
  cbind(eplc/s, p[1] * (pl * eplc/s - eplc %*% t(ds)/s^2))
}

eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 
  0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 
  0.5, -0.1)), weight_gradients = list(nealmon = nealmon_gradient))

deriv_tests(eq_r, tol = 1e-06)

coef(eq_r, midas = TRUE)

amweights(p = c(1, -0.5), d = 8, m = 4, weight = nealmon, type = "C")

nealmon(p = c(1, -0.5), d = 4)

eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4, amweights, nealmon, 
  "C") + mls(z, 0:16, 12, nealmon), start = list(x = c(1, -0.5), 
  z = c(2, 0.5, -0.1)))

hAh_test(eq_r)
hAhr_test(eq_r)

eq_rb <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 
  0:12, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 
  -0.1)))
hAh_test(eq_rb)
hAhr_test(eq_rb)

set_x <- expand_weights_lags(weights = c("nealmon", "almonp"), 
  from = 0, to = c(5, 10), m = 1, start = list(nealmon = c(1, 
    -1), almonp = c(1, 0, 0)))
set_z <- expand_weights_lags(c("nealmon", "nealmon"), 0, c(10, 
  20), 1, start = list(nealmon = c(1, -1), nealmon = c(1, -1, 
  0)))

expand_weights_lags(weights = c("nealmon", "nbeta"), from = 1, 
  to = c(2, 3), m = 1, start = list(nealmon = c(1, -1), nbeta = rep(0.5, 
    3)))

## eqs_ic <- midas_r_ic_table(y ~ trend + mls(x, 0, m = 4) + fmls(z, 
##   0, m = 12), table = list(z = set_z, x = set_x))
## eqs_ic$candlist[[5]]  <- 
##   update(eqs_ic$candlist[[5]], 0function = "nls")
## eqs_ic <- update(eqs_ic)
## modsel(eqs_ic, IC = "AIC", type = "restricted")

newx <- rnorm(4)
newz <- rnorm(12)
forecast(eq_rb, newdata = list(x = newx, z = newz, trend = 251))

eq_f <- midas_r(y ~ trend + mls(x, 4 + 0:7, 4, nealmon) + mls(z, 
  12 + 0:16, 12, nealmon), start = list(x = c(1, -0.5), z = c(2, 
  0.5, -0.1)))
forecast(eq_f, newdata = list(x = rep(NA, 4), z = rep(NA, 12), 
  trend = 251))

## cbfc <- select_and_forecast( y ~ trend + mls(x, 0, 4) + 
##   mls(z, 0, 12), from = list(x = c(4, 8, 12), z = c(12, 24, 36)), 
##   to = list(x  = rbind(c(14, 19), c(18, 23), c(22, 27)), 
##   z = rbind(c(22, 27), c(34, 39), c(46, 51))), insample = 1:200, 
##   outsample = 201:250, weights  = list(x = c("nealmon", "almonp"), 
##   z = c("nealmon", "almonp")), wstart = list(nealmon = rep(1, 3),
##   almonp = rep(1, 3)), IC = "AIC", seltype = "restricted",
##   ftype = "fixed", measures = c("MSE", "MAPE", "MASE"),
##   fweights = c("EW", "BICW", "MSFE", "DMSFE"))
## cbfc$accuracy$individual cbfc$accuracy$average
## mod1 <- midas_r(y ~ trend + mls(x, 4:14, 4, nealmon) + 
##   mls(z, 12:22, 12, nealmon), start = list(x = c(10, 1, -0.1), 
##   z = c(2, -0.1))) avgf <- average_forecast(list(mod1), 
##   data = list(y = y, x  = x, z = z, trend = trend),
##   insample = 1:200, outsample = 201:250, type = "fixed",
##   measures = c("MSE", "MAPE", "MASE"),
##   fweights = c("EW", "BICW", "MSFE", "DMSFE"))

data("USqgdp", package = "midasr")
data("USpayems", package = "midasr")
y <- window(USqgdp, end = c(2011, 2))
x <- window(USpayems, end = c(2011, 7))

yg <- diff(log(y)) * 100
xg <- diff(log(x)) * 100

nx <- ts(c(NA, xg, NA, NA), start = start(x), frequency = 12)
ny <- ts(c(rep(NA, 33), yg, NA), start = start(x), frequency = 4)

plot.ts(nx, xlab = "Time", ylab = "Percentages", col = 4, ylim = c(-5, 
  6))
lines(ny, col = 2)

xx <- window(nx, start = c(1985, 1), end = c(2009, 3))
yy <- window(ny, start = c(1985, 1), end = c(2009, 1))

beta0 <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbeta), 
  start = list(xx = c(1.7, 1, 5)))

coef(beta0)

betan <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbetaMT), 
  start = list(xx = c(2, 1, 5, 0)))
coef(betan)

um <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3), start = NULL)
coef(um)

fulldata <- list(xx = window(nx, start = c(1985, 1), end = c(2011, 
  6)), yy = window(ny, start = c(1985, 1), end = c(2011, 2)))
insample <- 1:length(yy)
outsample <- (1:length(fulldata$yy))[-insample]

avgf <- average_forecast(list(beta0, betan, um), data = fulldata, 
  insample = insample, outsample = outsample)
sqrt(avgf$accuracy$individual$MSE.out.of.sample)

harstep <- function(p, d, m) {
  if (d != 20) 
    stop("HAR(3)-RV process requires 20 lags")
  out <- rep(0, 20)
  out[1] <- p[1] + p[2]/5 + p[3]/20
  out[2:5] <- p[2]/5 + p[3]/20
  out[6:20] <- p[3]/20
  out
}

data("rvsp500", package = "midasr")
spx2_rvol <- 100 * sqrt(252 * as.numeric(rvsp500[, "SPX2.rv"]))
mh <- midas_r(rv ~ mls(rv, 1:20, 1, harstep), data = list(rv = spx2_rvol), 
  start = list(rv = c(1, 1, 1)))
summary(mh)

mr <- midas_r(rv ~ mls(rv, 1:20, 1, nealmon), data = list(rv = spx2_rvol), 
  start = list(rv = c(0, 0, 0)), weight_gradients = list())
summary(mr)

hAhr_test(mh)
hAhr_test(mr)

plot_midas_coef(mh, title = "")
lines(0:19, coef(mr, midas = TRUE, term_names = "rv"), col = 3)

tb <- expand_weights_lags("nealmon", from = 1, to = c(5, 15), 
  start = list(nealmon = c(0, 0, 0)))
mtb <- midas_r_ic_table(rv ~ mls(rv, 1:20, 1, nealmon), data = list(rv = spx2_rvol), 
  table = list(rv = tb), test = "hAh_test", weight_gradients = list(), 
  show_progress = FALSE)
mtb$candlist <- lapply(mtb$candlist, update, Ofunction = "nls")
mtb$test <- "hAhr_test"
mtb <- update(mtb)

bm <- modsel(mtb)
bm <- update(bm, Ofunction = "optim", method = "BFGS")

ar20 <- midas_r(rv ~ mls(rv, 1:20, 1), data = list(rv = spx2_rvol), 
  start = NULL)
forc <- average_forecast(list(ar20, mh, bm), data = list(rv = spx2_rvol), 
  insample = 1:1000, outsample = 1001:1100, type = "rolling", 
  show_progress = FALSE)
forc$accuracy$individual

## nealmon(p = c(2, 0.5, -0.1), d = 17)

## Code generating oss_prec data set mentioned in Appendix

set.seed(1001)

gendata <- function(n) {
    trend <- c(1:n)
    z <- rnorm(12*n)
    fn.z <- nealmon(p = c(2, 0.5, -0.1), d = 17)
    y <- 2 + 0.1 * trend + mls(z, 0:16, 12) %*% fn.z + rnorm(n)
    list(y = as.numeric(y), z = z, trend = trend)
}

n <- c(50, 100, 200, 300, 500, 750, 1000)

mse <- function(x) {
    mean(residuals(x)^2)
}

bnorm <- function(x) {
    sqrt(sum((coef(x, midas = TRUE) - c(2, 0.1, nealmon(p = c(2, 0.5, -0.1), d = 17)))^2))
}

rep1 <- function(n) {
    dt <- gendata(round(1.25 * n))
    ni <- n
    ind <- 1:ni
    mind <- 1:(ni * 12)
    indt <- list(y = dt$y[ind], z = dt$z[mind], trend = dt$trend[ind])
    outdt <- list(y = dt$y[-ind], z = dt$z[-mind], trend = dt$trend[-ind])
    um <- midas_r(y ~ trend + mls(z, 0:16, 12), data = indt, start = NULL)
    nm <- midas_r(y ~ trend + mls(z, 0:16, 12, nealmon), data = indt, start = list(z = c(1, -1, 0)))
    am <- midas_r(y ~ trend + mls(z, 0:16, 12, almonp), data = indt, start = list(z = c(1, 0, 0, 0)))
    modl <- list(Unrestricted = um, Correct = nm, Incorrect = am)
    list(norms = sapply(modl, bnorm), 
         mse = sapply(modl, function(mod) mean((as.data.frame(forecast(mod, newdata = outdt)) - outdt$y)^2)))
}

repr <- function(n, R) {
    cc <- lapply(1:R, function(i) rep1(n))
    list(norms = t(sapply(cc, "[[", "norms")), mse = t(sapply(cc, "[[", "mse")))
}

res <- lapply(n, repr, R = 1000)

norms <- data.frame(n, t(sapply(lapply(res, "[[", "norms"), function(l) apply(l, 2, mean))))
mses <- data.frame(n, t(sapply(lapply(res, "[[", "mse"), function(l) apply(l, 2, mean))))

suppressMessages(library("reshape2"))

msd <- melt(mses, id = 1)
colnames(msd)[2] <- "Constraint"
nmd <- melt(norms, id = 1)
colnames(nmd)[2] <- "Constraint"

msd$Type <- "Mean squared error"
nmd$Type <- "Distance from true values"
oos_prec <- rbind(msd, nmd)
oos_prec$Type <- factor(oos_prec$Type, levels = c("Mean squared error", "Distance from true values"))

qplot(x = n, y = value, data = oos_prec, geom = "line", colour = Constraint, 
  ylab = "") + facet_wrap(~Type, scales = "free_y") + xlab("Sample size") + 
  theme_bw()
