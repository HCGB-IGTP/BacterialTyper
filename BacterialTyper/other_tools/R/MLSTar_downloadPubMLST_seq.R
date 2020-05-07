#!/usr/bin/env Rscript
library("optparse")
 
## get options
option_list = list(
	make_option(c("-s", "--species"), type="character", help="species name for PubMLST download", metavar="character"),
	make_option(c("-s", "--scheme"), type="integer", help="scheme id for PubMLST download", metavar="integer"),
	make_option(c("-d", "--dir_seq"), type="character", help="folder path to download sequences", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get help message
if (is.null(opt$species)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## get arguments

## load library
library(MLSTar, lib.loc=opt$lib.loc)
downloadPubmlst_seq(org=opt$species, scheme=opt$scheme, dir=opt$dir_seq)

