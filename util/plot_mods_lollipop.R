###plot mods with splicegraphs

suppressPackageStartupMessages({
  library("RColorBrewer")
  library("gplots")
  library("ggplot2")
  library(trackViewer)
  library(rtracklayer)
  library(Gviz)
  library(VariantAnnotation)
})

args <- commandArgs(T)
if (length(args) < 2) {
  cat("USAGE: ./plot_mods_lollipop mods.bed transcripts.bed12\n")
  q()
}

mods_file <- args[1]
bed12_file <- args[2]
pdf_file <- args[3]

bed12_colnames = c("chrom", "start", "stop", "transcript", "gene", "strand", "CDS_start", "CDS_stop", "rgb", "exon_n", "exon_size", "exon_start")
mods_colnames = c("chrom", "start", "stop", "name", "type", "strand")

# Function to add splicegraph for each transcript isoform

feature_colors = rep(RColorBrewer::brewer.pal(n = 9, name = "Pastel1"),3)

# Function plot tracks for each transcript isoform and overlap modifications
bed12_splicegraph = function(mods, bed12) {
  #define total window to contain all transcript isoforms
  min_start = min(bed12[,"start"])
  max_stop = max(bed12[,"stop"])
  length = max_stop-min_start
  #define modification relative positions within this window
  total_types = factor(c("D","i6A|t6A", "m1A|m1I|ms2i6A","m1G", "m2G|m22G", "m3C","Y"))
  cols <- brewer.pal(7, "Set1")
  type2index = data.frame(index=1:7, row.names = total_types) 
  mods2index = type2index[mods$type,]
  pos = mods$start-min_start
  n_exons_all = c()
  # define GRanges object based on first transcript isoform
  chrom = bed12[1,"chrom"]
  offset = (bed12[1,"start"])-min_start
  starts = sapply(strsplit(as.character(bed12[1,"exon_start"]), ","), function(v) as.numeric(unlist(v)))+offset
  lengths = sapply(strsplit(as.character(bed12[1,"exon_size"]), ","), function(v) as.numeric(unlist(v)))
  n_exons = length(unlist(strsplit(as.character(bed12[1,"exon_size"]), ",")))
  n_exons_all = n_exons
  exon_numbers = (1:n_exons)
  exon_names=paste0("exon", 1:n_exons)
  features <- GRanges(chrom, IRanges(starts, width = lengths, names = exon_names), height = rep(0.05, n_exons))
  if (nrow(bed12)>1){
    for (i in 2:nrow(bed12)){
      offset = (bed12[1,"start"])-min_start
      starts = sapply(strsplit(as.character(bed12[i,"exon_start"]), ","), function(v) as.numeric(unlist(v)))+offset
      lengths = sapply(strsplit(as.character(bed12[i,"exon_size"]), ","), function(v) as.numeric(unlist(v)))
      n_exons = length(unlist(strsplit(as.character(bed12[i,"exon_size"]), ",")))
      n_exons_all[i] = n_exons
      exon_numbers = c(exon_numbers, (1:n_exons))
      exon_names=paste0("exon", 1:n_exons)
      more_features <- GRanges(chrom, IRanges(starts, width = lengths, names = exon_names), height = rep(0.05, n_exons))
      features = c(features, more_features)
    }
  }
  layer_names = rep(bed12$transcript, times=n_exons_all)
  exon_colors = layer_names
  levels(exon_colors) = feature_colors[1:nrow(bed12)]
  features$featureLayerID <- layer_names
  features$fill <- exon_colors
  #names(features) <- paste(layer_names, exon_numbers, sep="_")
  legend = cols
  names(legend) = total_types
  sample.gr <- GRanges(chrom, IRanges(pos, width=1, names=mods$type), color=cols[mods2index])
  levels(features$featureLayerID) = rev(levels(features$featureLayerID)) #reverses plot order to stay consistent with bed12 order
  lolliplot(sample.gr, features, legend = legend, ylab = bed12$transcript)
}

bed12_splicegraph_margin = function(mods, bed12, left, right) {
  #define total window to contain all transcript isoforms
  min_start = min(bed12[,"start"])
  max_stop = max(bed12[,"stop"])
  length = max_stop-min_start
  #define modification relative positions within this window
  total_types = factor(c("D","i6A|t6A", "m1A|m1I|ms2i6A","m1G", "m2G|m22G", "m3C","Y"))
  cols <- brewer.pal(7, "Set1")
  type2index = data.frame(index=1:7, row.names = total_types) 
  mods2index = type2index[mods$type,]
  pos = mods$start-min_start
  n_exons_all = c()
  # define GRanges object based on first transcript isoform
  chrom = bed12[1,"chrom"]
  n_exons = 1
  exon_numbers = c(1)
  exon_names="0"
  features <- GRanges(chrom, IRanges(start = -left, width = (length+right+left+1), names = exon_names), height = rep(0.05, 1))
  for (i in 1:nrow(bed12)){
    offset = (bed12[1,"start"])-min_start
    starts = sapply(strsplit(as.character(bed12[i,"exon_start"]), ","), function(v) as.numeric(unlist(v)))+offset
    lengths = sapply(strsplit(as.character(bed12[i,"exon_size"]), ","), function(v) as.numeric(unlist(v)))
    n_exons = length(unlist(strsplit(as.character(bed12[i,"exon_size"]), ",")))
    n_exons_all[i] = n_exons
    exon_numbers = c(exon_numbers, (1:n_exons))
    exon_names=paste0("exon", 1:n_exons)
    more_features <- GRanges(chrom, IRanges(starts, width = lengths, names = exon_names), height = rep(0.05, n_exons))
    features = c(features, more_features)
  }
  layer_names = c("0", rep(bed12$transcript, times=n_exons_all))
  exon_colors = layer_names
  levels(exon_colors) = feature_colors[1:nrow(bed12)+1]
  features$featureLayerID <- layer_names
  features$fill <- exon_colors
  #names(features) <- paste(layer_names, exon_numbers, sep="_")
  legend = cols
  names(legend) = total_types
  sample.gr <- GRanges(chrom, IRanges(pos, width=1, names=mods$type), color=cols[mods2index])
  lolliplot(sample.gr, features, legend = legend, ylab = bed12$transcript)
}


#########################

##plot mods for transcript of interest

bed12 = read.table(bed12_file)
colnames(bed12) = bed12_colnames
mods = read.table(mods_file)
colnames(mods) = mods_colnames

pdf(pdf_file, 8, 8)
bed12_splicegraph(mods, bed12)
# add transcipts labels. x and y parameters found through trial and error)
baseline = 0.225-0.051 #location for start of transcript id text
step = 0.051 #step between each transcript id 
x = 0.05
len = length(as.character(bed12$transcript))
grid.text(rev(as.character(bed12$transcript)), x = rep(x, len), y = seq(baseline,baseline+(step*(len-1)),step), gp = gpar(fontsize=10))
#bed12_splicegraph_margin(mods, bed12, 500, 500)
dev.off()
