#!/usr/bin/env Rscript

#' This R script plots a newick tree and adds the abundances as they are present in a Kallisto matrix. If no
#' Kallisto data is available, it only prints the tree. If Kallisto data is available, you can plot the heatmap
#' using a relative scale with the highest value set to 100 %, or a logarithmic scale. 
#' Special adjustments have been made for the analysis of the Glyphosate_gene_richness project.

library(argparse)
library(extrafont)
library(ggtree)
library(tidyverse)
library(reshape2)
library(data.table)

loadfonts(quiet = FALSE)
print(fonts())

parser <- ArgumentParser()
parser$add_argument("-t", "--tree", default = NULL, dest = "tree",
                    required = TRUE, help = "newick tree")
parser$add_argument("-k", "--kallisto", default = NULL, dest = "kallisto",
                    required = FALSE, help = "kallisto matrix (tab separated)")
parser$add_argument("-o", "--output", default = NULL, dest = "output",
                    required = TRUE, help = "output tree plot")
parser$add_argument("-s", "--scale_heatmap", default = "relative", dest = "scale",
					required = FALSE, help = "heatmap scale (logarithmic or relative)")
args <- parser$parse_args()


tree <- read.newick(args$tree)

tree_plot <- ggtree(tree) +
  geom_treescale(y = -2, offset = 0.75, fontsize = 3) +
  geom_tiplab(size = 2)


if(file.exists(args$kallisto)){
  
  stopifnot(args$scale == "logarithmic" | args$scale == "relative")
  kallisto <- read_delim(args$kallisto, delim = "\t")
  
  if(args$scale == "logarithmic"){
    
    kallisto_subset <- kallisto %>%
      select(-c(protein_name, tpm_B8, tpm_B9, tpm_B10)) %>%
      column_to_rownames("id") %>%
      rename(A1 = tpm_A1,
             A2 = tpm_A2,
             A3 = tpm_A3,
             A4 = tpm_A4,
             A5 = tpm_A5,
             A6 = tpm_A6,
             A7 = tpm_A7)

    log_kallisto <- apply(kallisto_subset, c(1, 2), function(x) log10(x+1))

    gheatmap(tree_plot, log_kallisto, offset = 2, low = "blue", high = "red", colnames_angle = 90)
    ggsave(args$output, device = "pdf", dpi = 300)
	
  } else if (args$scale == "relative"){
    
	kallisto_subset_for_melt <- kallisto[, c(1, 3:9)]
    kallisto_melt <- reshape2::melt(kallisto_subset_for_melt, id = "id")
    
	# turn data into long format to calculation of relative abundances
    kallisto_melt <- kallisto_melt[with(kallisto_melt, order(id)), ]
    kallisto_melt <- setDT(kallisto_melt)
    kallisto_melt[,relative_abundance := value / max(value) * 100, by = .(id)]
    kallisto_melt[is.na(relative_abundance), relative_abundance := 0]
    
	# remove absolute values before turning back into wide format
    kallisto_melt_relative <- kallisto_melt[, -"value"]
    kallisto_relative <- reshape2::dcast(kallisto_melt_relative, id ~ variable)
    kallisto_relative <- setDF(kallisto_relative)
    rel_kallisto <- kallisto_relative %>%
      column_to_rownames("id") %>%
      rename(A1 = tpm_A1,
             A2 = tpm_A2,
             A3 = tpm_A3,
             A4 = tpm_A4,
             A5 = tpm_A5,
             A6 = tpm_A6,
             A7 = tpm_A7)
    gheatmap(tree_plot, rel_kallisto, offset = 2, low = "white", high = "black", colnames_angle = 90)
    ggsave(args$output, device = "pdf", dpi = 300)
  }
} else {

  ggsave(tree_plot, file = args$output, 
					device = "pdf", 
					width = 10, 
					height = 10)
}
