#!/usr/bin/Rscript env

library('class')

if (length(commandArgs(T)) < 2) {
  cat("USAGE: classify_mods.R table model\n")
  q()
}

infn = commandArgs(T)[1]
modelfn = commandArgs(T)[2]

load(modelfn)

nucs = c('A','C','G','T')

x = read.table(infn, as.is=T, sep="\t", header=T)

x = x[x$sig == TRUE,]

for(i in 1:nrow(x)) {
 precursor = sub('T', 'U', x[i,'refnuc']) 
  x.p = pred.mod(modmodel2, precursor,
    (data.frame(x[i,nucs[nucs != x[i,'refnuc']]])))
  x[i,'pred.mod'] = x.p
}

write.table(x, file="", row.names=F, col.names=T, quote=F,sep="\t")



