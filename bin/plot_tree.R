#!/usr/bin/env Rscript

#' This R script plots a newick tree and adds the abundances as they are present in a Kallisto matrix. Special
#' adjustments have been made for the analysis of the Glyphosate_gene_richness project.

library(argparse)
library(extrafont)
library(ggtree)
library(tidyverse)

loadfonts(quiet = FALSE)
print(fonts())

parser <- ArgumentParser()
parser$add_argument("-t", "--tree", default = NULL, dest = "tree",
                    required = TRUE, help = "newick tree")
parser$add_argument("-k", "--kallisto", default = NULL, dest = "kallisto",
                    required = FALSE, help = "kallisto matrix (tab separated)")
parser$add_argument("-o", "--output", default = NULL, dest = "output",
                    required = TRUE, help = "output tree plot")
args <- parser$parse_args()


tree <- read.newick(args$tree)

if(file.exists(args$kallisto)){
  kallisto <- read_delim(args$kallisto, delim = "\t")

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

  tree_plot <- ggtree(tree) +
    geom_treescale(y = -2, offset = 0.75, fontsize = 3) +
    geom_tiplab(size = 2)

  gheatmap(tree_plot, log_kallisto, offset = 2, low = "blue", high = "red", colnames_angle = 90)
  ggsave(args$output, device = "pdf", dpi = 300)
} else {
  print("No kallisto data was provided, plotting tree without heatmap")
  
  tree_plot <- ggtree(tree) +
    geom_treescale(y = -2, offset = 0.75, fontsize = 3) +
    geom_tiplab(size = 2)
	
  ggsave(tree_plot, file = args$output, 
					device = "pdf", 
					width = 9, 
					height = 10, 
					dpi = 300)
}
