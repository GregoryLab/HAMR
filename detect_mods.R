# Copyright (c) 2013 University of Pennsylvania

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


# get input arguments
args <- commandArgs(TRUE)
if (length(args) < 6)
{
  cat("USAGE: detect_mods.R in_table seq_err hypothesis_type max_p max_fdr\n")
  q()
}
infn = args[1] # input HAMR nucleotide freq table
seq.err = as.numeric(args[2]) # per-bp sequencing error rate
hypothesis = args[3] # H1 or H4
maxp = as.numeric(args[4]) # p-value threshold
maxq = as.numeric(args[5]) # q-value threshold
#refpercent = as.numeric(args[6]) # minimum percentage of ref nuc

nucs = c('A','C','G','T') #

# read HAMR nucleotide freq table
x = read.table(infn, as.is=T, sep="\t",
  col.names=c('chr', 'bp', 'strand',
    'refnuc', 'A', 'C', 'G', 'T', 'nonref'))

# select rows with A/C/G/T
u<-grepl('N',x$refnuc) 
x<-x[u==FALSE,] 

# add column with ref nuc counts
x$ref = rowSums(x[,nucs]) - x$nonref

hyps = c('AA', 'AC', 'AG', 'AT', 'CC', 'CG', 'CT', 'GG', 'GT', 'TT')
# num_positions x numGenotypes p-value table
hyp.ps = array(NA, dim=c(nrow(x), length(hyps)))
colnames(hyp.ps) = hyps

for(h in hyps) {
  correct.nucs = unique(unlist(strsplit(h,'')))
  err.nucs = nucs[!nucs %in% correct.nucs]
  correct.counts = rowSums(data.frame(x[,correct.nucs]))
  err.counts = rowSums(data.frame(x[,err.nucs]))
  # FIXME: H1: p=1-prob(seq.error) = prob of observing ref. nuc
  #	            (assuming homozygous reference)
  #        H4: p=prob(nuc1 | nuc2) (assumes biallelic locus)
  hyp.ps[,h] = pbinom(correct.counts,
                        rowSums(cbind(correct.counts, err.counts)),
                        1-seq.err,
                        lower.tail=T)
}

hyp.union.ps = array(NA, dim=c(nrow(x), 2))
colnames(hyp.union.ps) = paste("H", c(1,4), sep='')

# RR
hyp.union.ps[,'H1'] = sapply(1:nrow(x), function(i) {
  hyp.ps[i,sprintf("%s%s",x$refnuc[i],x$refnuc[i])]
})

# all
hyp.union.ps[,'H4'] = apply(hyp.ps, 1, max)

hyp.union.ps.adj = apply(hyp.union.ps, 2, p.adjust, method='BH')

x$h1.p = hyp.union.ps[,'H1']
x$h1.padj = hyp.union.ps.adj[,'H1']
x$h4.p = hyp.union.ps[,'H4']
x$h4.padj = hyp.union.ps.adj[,'H4']
if(hypothesis=='H1'){
  # FIXME: 1. 0.05/0.95 are hard-coded here; should be minq/1-minq?
  #        2. why H4 is involved in H1 test? this seems to be adhoc
  #        3. should H1 be dropped alltogether? It seem one needs to first classifify locus as homozygous or biallelic?
  u<-x[which(x$h4.padj>0.05 & x$h1.padj<.05),]
  h.edit<-sapply(1:nrow(u), function(i){
    n=match(u$bp[i],x$bp)
     
    ifelse(length(which(hyp.ps[n,]>.95))==1,1,0)
  })
  u$h.edit<-h.edit
  u$sig<-ifelse(u$h.edit==1,"TRUE", "FALSE")
  write.table(u,file="",row.names=F,col.names=T,quote=F,sep="\t")
}


if(hypothesis=='H4'){
  x$sig = hyp.union.ps[,hypothesis] < maxp &
    hyp.union.ps.adj[,hypothesis] < maxq

  write.table(x, file="", row.names=F, col.names=T, quote=F,sep="\t")
}
