#! /usr/bin/Rscript

library(RColorBrewer)
library(colorspace)

args <- commandArgs(T)
if (length(args) < 2) {
  cat("USAGE: ./summ_mod_type_fromBed.R strict_input.bed out.pdf\n")
  q()
}

input_fn <- args[1]
out_fn <- args[2]
ignore_ref_eq_0 <- args[3]


# Read in classification data
mods <- read.table(input_fn, sep = "\t", header = F, as.is=T)

#assumes data are already filtered for true mods

# Summarize 
types = c(
"D", 
"i6A|t6A",
"m1A|m1I|ms2i6A",
"m1G",
"m2G|m22G",
"m3C",
"Y"
)
mod_type_list <- strsplit(mods[,5], ",")
names(mod_type_list) <- rownames(mods)
summ <- rbind(sapply(types, function (c) sum(sapply(mod_type_list, function(v) is.element(c, v)))))

rownames(summ) = "Type"

# plot
cols <- brewer.pal(7, "Set3")

pdf(out_fn, title=out_fn)
par(mar=c(11,4,4,3))
barplot(summ, xlab = "Modification Type", beside = TRUE, las = 2, col = "darkblue", ylim = c(0, max(summ)*1.25))
title(ylab = "Count", mgp = c(3, 0.8, 0))
legend("topright", legend = c("Modification Type"), fill = "darkblue")
pie(summ, labels=types, radius=1, col=cols)
dev.off()