#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-s", "--species"), type="character", help="organism name", metavar="character"),
  make_option(c("-o", "--output"), type="character", help="output file name", metavar="character"),
  make_option(c("-l", "--lib.loc"),type="character",help="Install path location", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get help message
if (is.null(opt$sp)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## load additional library
.libPaths(opt$lib.loc)
library(MLSTar)

## Check available PUBMLST schemes available

# Name the data frame
df <- data.frame(org=NA, scheme=NA, len=NA)

### loop through species
my_list_species <- as.list(strsplit(opt$sp, ",")[[1]])
for (j in 1:length(my_list_species))
{
  list_schemes <- listPubmlst_schemes(org=my_list_species[[j]])
  
  ### loop through scheme per species
  for (i in 1:length(list_schemes))
  {
    scheme_len <- lengths(list_schemes[i])
    df[nrow(df) + 1,] = list(my_list_species[[j]], i, scheme_len)
  }
}

## get results
write.table(df, file=opt$out, sep = ",", quote = FALSE) ## save to file



