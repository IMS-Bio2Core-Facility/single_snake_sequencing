packages <- c("DropletUtils", "logger", "Matrix")
invisible(
  suppressPackageStartupMessages(
    lapply(packages, library, character.only = TRUE)
  )
)

set.seed(0)

# Configure logger
style <- layout_glue_generator("{time} :: {level} :: {namespace} :: {msg}")
logfile <- snakemake@log[[1]]
log_appender(appender_file(logfile))
log_layout(style)
log_threshold(INFO)

# Read in data
# readMM fails with floats, so we hack python
counts <- readMM(snakemake@input[["mtx"]])
log_success("Data read from {snakemake@input[['counts']]}")
log_info("Dimensions: {dim(counts)}")

totals <- colSums(counts)
totals <- totals[totals > 25]
totals <- totals[totals < 2000]

# The bound beneath which all droplets are "empty"
lower <- mean(totals) + (2 * sd(totals))
log_info("Lower threshold: {lower}")

# Tail end plot
png(snakemake@output[["tail"]], res = 300, height = 4, width = 4, units = "in")
invisible(hist(totals, breaks = 500, xlab = "nUMI", main = NULL, col = "black"))
invisible(abline(v = lower, col = "darkgreen", lty = 2))
invisible(dev.off())
log_success("Tail end plot saved at {snakemake@output[['tail']]}")

# Knee plot
bcrank <- barcodeRanks(counts)
uniq <- !duplicated(bcrank$rank)
png(snakemake@output[["knee"]], res = 300, height = 4, width = 4, units = "in")
invisible(plot(bcrank$rank[uniq], bcrank$total[uniq],
  log = "xy",
  xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2
))
invisible(abline(h = lower, col = "darkgreen", lty = 2))
invisible(dev.off())
log_success("Knee plot saved at {snakemake@output[['knee']]}")

# And run emptyDrops
results <- emptyDrops(
  counts,
  lower = lower, niters = snakemake@params[["niters"]]
)
write.csv(results, snakemake@output[["empty"]])
log_success("Data saved at {snakemake@output[['empty']]}")
