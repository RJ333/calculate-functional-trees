#!/usr/bin/env Rscript

#' This R script plots a newick tree and adds the abundances as they are present in a Kallisto matrix. Special
#' adjustments have been made for the analysis of the Glyphosate_gene_richness project.

library(argparse)
library(extrafont)
library(ggtree)
library(tidyverse)
library(reshape2)
library(data.table)  # add to conda.config

loadfonts(quiet = FALSE)
print(fonts())

parser <- ArgumentParser()
parser$add_argument("-t", "--tree", default = NULL, dest = "tree",
                    required = TRUE, help = "newick tree")
parser$add_argument("-k", "--kallisto", default = NULL, dest = "kallisto",
                    required = TRUE, help = "kallisto matrix (tab separated)")
parser$add_argument("-o", "--output", default = NULL, dest = "output",
                    required = TRUE, help = "output tree plot")
args <- parser$parse_args()


tree <- read.newick(args$tree)
#tree <- read.newick("/home/centos/scripts/calculate-functional-trees/results/RAxML_bestTree.phnJ")
kallisto <- read_delim(args$kallisto, delim = "\t")
#kallisto <- read_delim("/home/centos/scripts/calculate-functional-trees/input_files/kallisto_abundance_matrix.tsv.gz", delim = "\t")
kallisto_subset_for_melt <- kallisto[, c(1, 3:9)]

kallisto_melt <- reshape2::melt(kallisto_subset_for_melt, id = "id")

# library(dplyr)
# df1 <- df %>%
  # group_by(id,protein_name) %>%
  # mutate(relative_abundance = value/max(value)*100)

df1[is.na(df1)] <- 0
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



# kallisto_subset <- kallisto %>%
  # select(-c(protein_name, tpm_B8, tpm_B9, tpm_B10)) %>%
  # column_to_rownames("id") %>%
  # rename(A1 = tpm_A1,
         # A2 = tpm_A2,
         # A3 = tpm_A3,
         # A4 = tpm_A4,
         # A5 = tpm_A5,
         # A6 = tpm_A6,
         # A7 = tpm_A7)

# log_kallisto <- apply(kallisto_subset, c(1, 2), function(x) log10(x+1))

tree_plot <- ggtree(tree) +
  geom_treescale(y = -2, offset = 0.75, fontsize = 3) +
  geom_tiplab(size = 2)

#gheatmap(tree_plot, log_kallisto, offset = 2, low = "blue", high = "red", colnames_angle = 90)
gheatmap(tree_plot, rel_kallisto, offset = 2, low = "lightblue", high = "darkblue", colnames_angle = 90)
ggsave(args$output, device = "pdf", dpi = 300)
