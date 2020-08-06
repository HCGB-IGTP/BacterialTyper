#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-l", "--lib"),type="character",help="Library to install", metavar="character"),
  make_option(c("-p", "--path"),type="character",help="Install path location", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## load additional library
.libPaths(opt$path)
library(opt$lib, character.only = TRUE)