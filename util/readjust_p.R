#!/usr/bin/env Rscript

argv = commandArgs(T)

if (length(argv) < 2) {
  cat("USAGE: readjust_p.R hamr_out_file adjusted_out_file")
  q()
}

infile = argv[1]
outfile = argv[2]

x = read.table(infile, sep="\t", as.is=T, header=T)
x$h1.padj = p.adjust(x$h1.p, method='BH')
x$h4.padj = p.adjust(x$h4.p, method='BH')
write.table(x, file=outfile, row.names=F, col.names=T, quote=F, sep="\t")

