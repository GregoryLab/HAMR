###plot mods with splicegraphs

##NOTE mods list will be plotted counterclockwise starting from 12 o`clock. Thus, sample 2 replicates should usually be in reverse order

suppressPackageStartupMessages({
  library("RColorBrewer")
  #library("gplots")
  #library("ggplot2")
  library(trackViewer)
  library(rtracklayer)
  library(Gviz)
  library(VariantAnnotation)
})

args <- commandArgs(T)
if (length(args) != 1) {
  cat("USAGE: ./plot_mods_lollipop filenames\n")
  q()
}

names <- args[1]

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
  #convert mods list to colors list
  modslist = strsplit(as.character(mods$mods), ",")
  for(i in seq_along(modslist)){
    for(j in seq_along(modslist[[i]])){
      modslist[[i]][[j]] = type2index[modslist[[i]][[j]],]
    }
  }
  for(i in seq_along(modslist)){
    modslist[[i]] = as.numeric(modslist[[i]])
    modslist[[i]] = cols[modslist[[i]]]
    modslist[[i]][which(is.na(modslist[[i]]))] = "#FFFFFF" #convert NAs to white space
  }
  # define GRanges object based on first transcript isoform
  pos = mods$start-min_start
  n_exons_all = c()
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
  # Define object for modifications
  sample.gr <- GRanges(chrom, IRanges(pos, width=1, names=mods$type))
  weights = as.numeric(unlist(strsplit(as.character(mods[1,"weights"]), ",")))
  sample.gr$value1 = weights[1]
  i = 1
  for (weight in weights[-1]) {
    i = i+1
    sample.gr$value1 = weights[i]
    value_name = paste("value", i, sep="")
    sample.gr[,value_name] = sample.gr[,"value1"]
  }
  sample.gr$color = modslist
  levels(features$featureLayerID) = rev(levels(features$featureLayerID)) #reverses plot order to stay consistent with bed12 order
  lolliplot(sample.gr, features, legend = legend, ylab = bed12$transcript, type='pie')
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
  #convert mods list to colors list
  modslist = strsplit(as.character(mods$mods), ",")
  for(i in seq_along(modslist)){
    for(j in seq_along(modslist[[i]])){
      modslist[[i]][[j]] = type2index[modslist[[i]][[j]],]
    }
  }
  for(i in seq_along(modslist)){
    modslist[[i]] = as.numeric(modslist[[i]])
    modslist[[i]] = cols[modslist[[i]]]
    modslist[[i]][which(is.na(modslist[[i]]))] = "#FFFFFF" #convert NAs to white space
  }  
  # define GRanges object based on first transcript isoform
  pos = mods$start-min_start
  n_exons_all = c()
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
  # Define object for modifications
  sample.gr <- GRanges(chrom, IRanges(pos, width=1, names=mods$type), color=cols[mods2index])
  weights = as.numeric(unlist(strsplit(as.character(mods[1,"weights"]), ",")))
  sample.gr$value1 = weights[1]
  i = 1
  for (weight in weights[-1]) {
    i = i+1
    sample.gr$value1 = weights[i]
    value_name = paste("value", i, sep="")
    sample.gr[,value_name] = sample.gr[,"value1"]
  }
  sample.gr$color = modslist
  levels(features$featureLayerID) = rev(levels(features$featureLayerID)) #reverses plot order to stay consistent with bed12 order
  lolliplot(sample.gr, features, legend = legend, ylab = bed12$transcript, type='pie')
}


#########################

##plot mods for transcript of interest

plot <- function(names) {
  mods_file = names[1]
  bed12_file = names[2]
  pdf_file = names[3]
  bed12_colnames = c("chrom", "start", "stop", "transcript", "gene", "strand", "CDS_start", "CDS_stop", "rgb", "exon_n", "exon_size", "exon_start")
  mods_colnames = c("chrom", "start", "stop", "weights", "mods", "strand")
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
}

names = read.table(names)
apply(names, 1, plot)

