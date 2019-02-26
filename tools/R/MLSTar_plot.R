#!/usr/bin/env Rscript

library("optparse")

## get options
option_list = list(
  make_option(c("-f_profile", "--file_profile"), type="character", help="folder path to download sequences", metavar="character"),
  make_option(c("-f_result", "--file_result"), type="character", help="folder path to download profile", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get help message
if (is.null(opt$file_profile)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## load library
library(MLSTar)

## get results
result_get <- read.delim(file=opt$file_result, header=TRUE, sep="\t")
# convert data frame, character, retype columns
resu <- as.data.frame(result_get, stringsAsFactors = FALSE)
resu[] <- lapply(resu, as.character)
resu <- resu[,c("ST",setdiff(names(resu),"ST"))]

## get profile
profile_get <- read.delim(file=opt$file_profile, header=TRUE, sep="\t")
# convert data frame, character
prof <- as.data.frame(profile_get, stringsAsFactors = FALSE)
prof[] <- lapply(profile_get, as.character)

#add attributes and convert to MLST object
out <- list(result = resu, profile = prof)
attr(out, 'infiles') <- "default"
attr(out, 'org') <- "default"
attr(out, 'scheme') <- "default"
attr(out, 'write') <- "default"
attr(out, 'pid') <- "default"
attr(out, 'scov') <- "default"
attr(out, 'class') <- 'mlst'

## generate plot and save to file
plot(out, 
     what='both',
     label=FALSE,
     alpha=0.4,
     col.axis='red',
     lty=3
)

##############################
## S3 method for class mlst
##  plot(x, type = "mst", what = "both", pt.size = 3,
##       label = FALSE, pf.col = "#E64B35FF", st.col = "#00A087FF",
##       nst.col = "white", alpha = 0.5, plot.igraph.args = list(),
##       plot.phylo.args = list(), tiplabels.args = list(), plot = TRUE,
##       ...)

##############################
##  Arguments
##  x:        An object of class mlst. 
##  type:     One of "mst" or "phylo", for plotting a Minimum Spanning Tree or a binary tree.
##  what:     One of "result", "profile", "both" (default). What should be plotted.
##  pt.size:  The size of the point. Default: 3.
##  label:    logical. Whether to plot node/tip labels or not. Default: FALSE.
##  pf.col:   The color of profile nodes/tips. Ignored if what="result".
##  st.col:   The color of result nodes/tips which have an ST assigned. Ignored if what="profile".
##  nst.col:  The color of result nodes/tips with no ST assigned. Default: "white". Ignored if what="profile".
##  alpha:    Color transparency, between [0,1]. See alpha.
##  plot.igraph.args: A list of arguments to be passed to plot.igraph. Used only if type="mst". Defaults try to keep aesthetics similar to phylo. Default: list(vertex.label = if (label) NULL else NA, vertex.size = pt.size).
##  plot.phylo.args: A list of arguments to be passed to plot.phylo. Used only if type="phylo". Defaults try to keep aesthetics similar to mst. Default: list(type=unrooted, show.tip.label = label).
##  tiplabels.args: A list of arguments to be passed to tiplabels. Used only if type="phylo". Defaults try to keep aesthetics similar to mst. Default: list(pch = 21, cex = pt.size/5).
##  plot: Default: TRUE.
##  ...: A list of arguments to be passed to par.

##############################
##  Details :
##  Distance is calculated using dist.gene function over the allele matrix. This distance metric 
##  only takes into account the number of differences between each pair of rows in the matrix. 
##  If type="mst", mst function is used to calculate a minimum spanning tree, and a graph is 
##  generated using graph.adjacency. If type="phylo", a "phylogeny" is inferred using the distance 
##  calculated above and nj function (neighbour-joining tree). Ethier a igraph or a phylo object is
##  returned invisibly so it can be further analysed using igraph/ape frameworks.
##  
##  It is worth noting that the result of this function is not strictly based on genetic information.
##  The distance used just counts the number of differences on the allele profiles, it doesnt use a 
##  genetic model. The clustering methods are also simple and are not based on complex evolution models. 
##  Despite the above, since alleles in mlst schemes usually differ in 1 or few nucleotides, the 
##  described methodology gives good enough results to have a general overview of your data.
##  
##############################
##  Value
##  A minimum spanning tree plot and an object of class igraph (invisible) if type="mst", or a binary 
##  tree plot and an object of class phylo (invisible) if type="phylo".
