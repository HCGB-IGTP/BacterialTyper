#!/usr/bin/env Rscript

library("optparse")
 
## get options
option_list = list(
	make_option(c("-sp", "--species"), type="character", help="species name for PubMLST download", metavar="character"),
	make_option(c("-sc", "--scheme"), type="integer", help="scheme id for PubMLST download", metavar="integer"),
	make_option(c("-d_prf", "--dir_profile"), type="character", help="folder path to download profile", metavar="character")	
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
library(MLSTar)
downloadPubmlst_profile(org=opt$species, scheme=opt$scheme, dir=opt$dir_profile)

#############
## citation
#############
## Jolley KA, Bray JE, Maiden MCJ. A RESTful application programming interface for the PubMLST molecular typing and genome databases. Database: The Journal of Biological Databases and Curation. 2017;2017:bax060. doi:10.1093/database/bax060.
## MLSTar citation