#!/usr/bin/env nextflow

prokka_protein_fasta = file(params.prokka_protein_fasta)
prokka_gff = file(params.prokka_gff)
gene_name = params.gene_name

println """\
         workflow: build phylogenetic trees of proteins sequences
         ========================================================
         gene_name           : ${params.gene_name}
         prokka_protein_fasta: ${params.prokka_protein_fasta}
         prokka_gff          : ${params.prokka_gff}
         """
         .stripIndent()

process create_single_line_fasta {
    input:
      file prokka from prokka_protein_fasta

    output:
      file 'prokka.sl.faa' into single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka | tail -n +2 > prokka.sl.faa
    """
}

process subset_protein_fasta {
    publishDir 'results/'

    input:
      file prokka_single_line from single_line_fasta
      file gff_file from prokka_gff
    output:
      file 'prokka_subset.faa' into subset_prokka

    """
    cat $gff_file | \
      grep "gene=$gene_name" | \
      awk '\$3 == "CDS" {print \$0}' | \
      grep -A 1 "\$(awk -F '[\t;]' '\$3 == "CDS" {print \$9}' | grep -oP '^ID=\\K.*')" $prokka_single_line | \
      sed '/^--\$/d' > prokka_subset.faa
    """
}

