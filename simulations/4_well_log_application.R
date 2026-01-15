
# devtools::install_github("guillemr/robust-fpop")
library(robseg) # for Robust FPOP
library(DeCAFS) # for well log data
library(changepoint) # For PELT change un mean


#######

threshold_mood <- function(n, nbSeg, alpha = 0.01)
{
  m <- n/nbSeg - 1
  p_single <- 1 - (1 - alpha)^(1 / m)  # Dunn–Šidák
  qchisq(1 - p_single, df = nbSeg)
}


y <- (DeCAFS::oilWell)
n  <- length(y)
nbSeg <- 30  # <-- changed here


## variance estimation by MAD estimator

est_sd <- mad(y, constant = 1.4826)

## --- Changepoint estimation ---

system.time(
resPELT <- cpt.mean(
  y / est_sd,
  method    = "PELT",
  penalty   = "Manual",
  pen.value = 70
)
)

system.time(
resR <- Rob_seg.std(
  x          = y / est_sd,
  loss       = "Outlier",
  lambda     = 70,
  lthreshold = 2
)
)

system.time(
resW <- SVP_costTests(
  y,
  gamma = 1.5 * sqrt((n / nbSeg)^3 / 12),
  test  = "WilcoxonCost"
)
)
system.time(
resM <- SVP_costTests(
  y,
  gamma =  threshold_mood(n, 10),
  test  = "MedianMoodCost"
)
)


## --- Helper: piecewise constant (median) + vertical dashed lines ---

plot_piecewise_constant <- function(y,
                                    cps,
                                    main = "",
                                    col_line = 2)
{
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
  plot(x, y, type = "p", pch = ".", cex = 3,cex.main = 2,
       main = main, ylab = "", xlab = "")

  # vertical dashed lines at changepoints
  abline(v = cps, col = "black", lty = 2)

  # median on each segment -> horizontal piecewise constant fit
  for (k in seq_along(cps)) {
    idx <- starts[k]:cps[k]
    med <- median(y[idx], na.rm = TRUE)
    lines(x[idx], rep(med, length(idx)), lwd = 2, col = col_line)
  }
}

###################################################

out <- file.path("simulations", "plots", "4_segmentation_logdata.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

pdf(file = out, width = 10, height = 8)  # width/height in inches (like ggsave)
op <- par(mfrow = c(4, 1), mar = c(4, 4, 2, 0))
on.exit({ par(op); dev.off() }, add = TRUE)

# 1) PELT segmentation
plot_piecewise_constant(y, resPELT@cpts,
                        main = "PELT",
                        col_line = 4)

# 2) Robust segmentation
plot_piecewise_constant(y, resR$t.est,
                        main = "Robust FPOP",
                        col_line = 3)

# 3) SVP_costTests segmentation
plot_piecewise_constant(y, resW$changepoints,
                        main = "SVP Wilcoxon",
                        col_line = 2)

# 4) SVP_costTests segmentation
plot_piecewise_constant(y, resM$changepoints,
                        main = "SVP Median Mood",
                        col_line = 5)

dev.off()

length(resW$changepoints)
length(resM$changepoints)




