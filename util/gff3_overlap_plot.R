#! /usr/bin/Rscript

# Load required packages
library(methods)
library(RColorBrewer)

args <- commandArgs(T)
if (length(args) != 2) {
	cat("USAGE: ./summ_interval_class.R input_fn out_fn\n")
	q()
}

input_fn <- args[1]
out_fn <- args[2]

# Read in classification data
input_data <- read.table(input_fn, colClasses = c("character", "numeric"), sep = "\t", header = F)
colnames(input_data)= c("class", "count")

# plot
pdf(out_fn, title=out_fn)
cols <- brewer.pal(nrow(input_data), "Set1")
par(mar=c(11,4,4,3))
barplot(height = input_data$count, names.arg = input_data$class, xlab = "Modification class", beside = TRUE, las = 2, col = cols, ylim = c(0, max(input_data$count)*1.25), lwd = 2)
legend(legend = input_data$class, fill = cols, x = 'topright')
title(ylab = "Count", mgp = c(3, 0.8, 0))
par(lwd = 2)
pie(input_data$count, labels = "", radius=1, col=cols, lwd = 4)
legend(legend = input_data$class, fill = cols, x = 'topright')
dev.off()


