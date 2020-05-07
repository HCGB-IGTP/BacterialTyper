#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-s", "--dir_seq"), type="character", help="folder path for sequences download", metavar="character"),
  make_option(c("-p", "--dir_profile"), type="character", help="folder path for profile download", metavar="character"),
  make_option(c("-f", "--file"), type="character", help="fasta file", metavar="character"),
  make_option(c("-d", "--dir"), type="character", help="folder path", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=2, help="threads", metavar="integer"),
  make_option(c("-n", "--name"), type="character", help="name", metavar="integer"),
  make_option(c("-l", "--lib.loc"),type="character",help="MLSTar library location", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## get message
if (is.null(opt$dir_seq)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

##
library(MLSTar, lib.loc=opt$lib.loc)

setwd(opt$dir)
nameDir = paste0(opt$name, '_alleles')

## get files
fastas_scheme <- list.files(path = opt$dir_seq, full.names = TRUE, pattern = "\\.fas$")
profile_scheme <- list.files(path = opt$dir_profile, full.names = TRUE)

## doMLST
x <- doMLST(
  infiles = opt$file,  # The fasta files

  schemeFastas = fastas_scheme,
  schemeProfile = profile_scheme,
  
  write = "new",    # write fasta files for alleles found
  fdir = nameDir,	## there is a bug when generating fdir
  
  n_threads = opt$threads
  
  ########### Additional #########
  # pid	: Percentage identity threshold to be consider as a novel allele. An integer <= 100. (Default: 90).
  # scov : Subject coverage threshold to be consider as a novel allele. A numeric between 0 and 1. Not recomended to set it below 0.7 . (Default 0.9) # @details  
)

res_file <- paste0(opt$dir, '/', opt$name, "_MLST_results.csv")

## get results
write.table(x$result, file=res_file,sep = ",", quote = FALSE) ## save to file
