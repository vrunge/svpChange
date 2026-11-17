## Four independent Gaussian updown series, with change locations hidden.

source(file.path("simulations", "power_common.R"))

jumps <- c(0.6, 0.8, 1, 1.2)
n <- 1000L

series <- do.call(rbind, lapply(seq_along(jumps), function(i) {
  set_power_seed(20270101L + i)
  data.frame(
    time = seq_len(n),
    value = generate_power_signal(n, "updown", jumps[i], segments = 8L) +
      stats::rnorm(n),
    jump = paste0("Jump size = ", jumps[i])
  )
}))
series$jump <- factor(series$jump, paste0("Jump size = ", jumps))

plot <- ggplot2::ggplot(series, ggplot2::aes(time, value)) +
  ggplot2::geom_point(size = 0.5, alpha = 0.45, colour = "#245A81") +
  ggplot2::facet_wrap(~jump, nrow = 2L) +
  ggplot2::labs(x = "Time (index)", y = "Observed value") +
  paper_theme()

output <- file.path("simulations", "all_plots")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(output, "updown_four_jump_sizes.pdf"), plot,
                width = 10, height = 6)
