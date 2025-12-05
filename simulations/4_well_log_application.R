# devtools::install_github("guillemr/robust-fpop")
library(robseg)

# install.packages(fastcpd, repos = c("https://doccstat.r-universe.dev", "https://cloud.r-project.org"))
library(fastcpd)

library(changepoint)


threshold_mood <- function(n, nbSeg, alpha = 0.01)
{
  m <- n/nbSeg - 1
  p_single <- 1 - (1 - alpha)^(1 / m)  # Dunn–Šidák
  qchisq(1 - p_single, df = nbSeg)
}

alpha <- 0.001
p_single <- alpha / (  m <- n/5 - 1)
qchisq(1 - p_single, df = 1)


y  <- well_log
n  <- length(well_log)
nbSeg <- 30  # <-- changed here


## --- Changepoint estimation ---


resPELT <- cpt.mean(
  y / est_sd,
  method    = "PELT",
  penalty   = "Manual",
  pen.value = 70
)


est_sd <- mad(y, constant = 1.4826)

resR <- Rob_seg.std(
  x          = y / est_sd,
  loss       = "Outlier",
  lambda     = 70,
  lthreshold = 2
)


resW <- SVP_costTEsts(
  y,
  gamma = 1.5 * sqrt((n / nbSeg)^3 / 12),
  test  = "WilcoxonCost"
)

resM <- SVP_costTEsts(
  y,
  gamma =  threshold_mood(n, 10),
  test  = "MedianMoodCost"
)



## --- Helper: piecewise constant (median) + vertical dashed lines ---

plot_piecewise_constant <- function(y, cps, main = "", col_line = 2) {
  n <- length(y)
  x <- seq_len(n)

  # clean and complete changepoint vector
  cps <- sort(unique(cps))
  cps <- cps[cps >= 1 & cps <= n]
  if (length(cps) == 0 || tail(cps, 1) != n) {
    cps <- c(cps, n)
  }

  starts <- c(1, cps[-length(cps)] + 1)

  # scatterplot of data
  plot(x, y, type = "p", pch = ".", cex = 3,
       main = main, xlab = "index", ylab = "y")

  # vertical dashed lines at changepoints
  abline(v = cps, col = "black", lty = 2)

  # median on each segment -> horizontal piecewise constant fit
  for (k in seq_along(cps)) {
    idx <- starts[k]:cps[k]
    med <- median(y[idx], na.rm = TRUE)
    lines(x[idx], rep(med, length(idx)), lwd = 2, col = col_line)
  }
}

## --- Three plots, one below the other ---

par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))


# 1) PELT segmentation
plot_piecewise_constant(y, resPELT@cpts,
                        main = "PELT",
                        col_line = 4)

# 2) Robust segmentation
plot_piecewise_constant(y, resR$t.est,
                        main = "Robust FPOP",
                        col_line = 3)

# 3) SVP_costTEsts segmentation
plot_piecewise_constant(y, resW$changepoints,
                        main = "SVP Wilcoson",
                        col_line = 2)

# 3) SVP_costTEsts segmentation
plot_piecewise_constant(y, resM$changepoints,
                        main = "SVP Median Mood",
                        col_line = 5)


length(resW$changepoints)
length(resM$changepoints)
