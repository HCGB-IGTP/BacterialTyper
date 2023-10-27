#!/usr/bin/env Rscript
library("optparse")

## get options
option_list = list(
  make_option(c("-i", "--input_file"), type="character", help="Input csv file with BUSCO stats", metavar="character"),
  make_option(c("-o", "--out_folder"), type="character", help="Folder to store plot", metavar="character")
  
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## get message
if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("No arguments provided", call.=FALSE)
}

## load additional library
library(ggplot2)
library("grid")
library(reshape2)

my_input_file <- opt$input_file
my_output_folder <- opt$out_folder

######################################
# Original
# BUSCO summary figure
# @version 4.0.0
# @since BUSCO 2.0.0
# Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)
# Licensed under the MIT license. See LICENSE.md file.
######################################

## Original
my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442") # Colors
my_bar_height <- 0.75 # Bar height ratio
my_title <- "BUSCO Assessment Results" # Legend
my_family <- "sans" # Font
my_size_ratio <- 1 
labsize = 1 # Code to produce the graph
########################



# Read input csv
my_input <- read.csv(my_input_file, row.names = 1)

colnames(my_input) <- c('C', 'S', 'D', 'F', 'M', 'Total',  'c', 's', 'd', 'f', 'm', 't', 'database')

my_species <- as.factor(rep(rownames(my_input),each=4))
if (length(levels(my_species)) > 10){ labsize = 0.66 } ## Original

## create dataframe
my_df <- my_input[unique(my_species),c('s','d','f','m')]
my_df['names'] <- rownames(my_df) 
my_df_melt <- melt(my_df)
colnames(my_df_melt) <- c('names', 'category', 'percentage')

##
figure <- ggplot(data = my_df_melt) + 
  
  geom_bar(aes(x = names, y = percentage, fill = category), 
           position = position_stack(reverse = TRUE), data = my_df_melt, stat="identity", width=my_bar_height) + 
  
  ## Original
  coord_flip() + 
  theme_gray(base_size = 8) + 
  scale_y_continuous(labels = c("0","20","40","60","80","100"), 
                     breaks = c(0,20,40,60,80,100)) + 
  scale_fill_manual(values = my_colors,labels =c(" Complete (C) and single-copy (S)  ",
                                                 " Complete (C) and duplicated (D)",
                                                 " Fragmented (F)  ",
                                                 " Missing (M)")) +   
  ggtitle(my_title) + 
  xlab("") + 
  ylab("\n%BUSCOs") + 
  
  theme(plot.title = element_text(family=my_family, hjust=0.5, colour = "black", size = rel(2.2)*my_size_ratio, face = "bold")) + 
  theme(legend.position="top",legend.title = element_blank()) + 
  theme(legend.text = element_text(family=my_family, size = rel(1.2)*my_size_ratio)) + 
  theme(panel.background = element_rect(color="#FFFFFF", fill="white")) + 
  theme(panel.grid.minor = element_blank()) + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.text.x = element_text(family=my_family, colour = "black", size = rel(1.66)*my_size_ratio)) + 
  theme(axis.line = element_line(size=1*my_size_ratio, colour = "black")) + 
  theme(axis.ticks.length = unit(.85, "cm")) + 
  theme(axis.ticks.y = element_line(colour="white", size = 0)) + 
  theme(axis.ticks.x = element_line(colour="#222222")) + 
  theme(axis.ticks.length = unit(0.4, "cm")) + 
  theme(axis.title.x = element_text(family=my_family, size=rel(1.2)*my_size_ratio)) + 
  
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


for(i in rownames(my_input) ){
  figure <- figure + 
    annotate("text", label=paste("C:", my_input[i, 'C'], 
                                 " [S:", my_input[i, 'S'], 
                                 ", D:", my_input[i, 'D'], "], F:", my_input[i, 'F'], 
                                 ", M:", my_input[i, 'M'], ", n:", my_input[i, 't'], sep=""), 
             y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
}


my_input_file

pdf(file.path(my_output_folder, paste0(my_input[i,'database'], ".pdf")), paper = "A4r", width = 35, height = 12)
print(figure)
dev.off()

print("Done")
