#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-l", "--lib.loc"),type="character",help="MLSTar library location", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(withr)
with_libpaths(new = opt$lib.loc, install_github('iferres/MLSTar'))
