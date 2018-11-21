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
    /*
      This process joins all sequence lines and creates a single line fasta file.
    */

    input:
      file prokka from prokka_protein_fasta

    output:
      file 'prokka.sl.faa' into single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka | tail -n +2 > prokka.sl.faa
    """
}

process subset_gff {
    /*
      This process subsets all lines from the prokka gff file and only keeps those with the same gene name as set in
      the nextflow config file.
    */

    input:
      file gff_file from prokka_gff
    output:
      file 'prokka_subset.gff' into subset_gff

    """
    cat $gff_file | \
      grep "gene=$gene_name" | \
      awk '\$3 == "CDS" {print \$0}' > prokka_subset.gff
    """
}

process subset_protein_fasta {
    /*
      This process uses the subsetted gff file from the process subset_gff and extracts the gene IDs from those lines.
      A space is added at the end of the line (e.g. to parse ID1 but not ID10) and uses this file to obtain all
      protein sequences from the single line fasta file (from the process create_single_line_fasta) that belong to the
      gene of interest.
    */

    publishDir 'results/'

    input:
      file prokka_single_line from single_line_fasta
      file gff_subset_file from subset_gff
    output:
      file 'prokka_subset.faa' into subset_prokka

    """
    cat $gff_subset_file | \
      awk -F '[\t;]' '\$3 == "CDS" {print \$9}' | \
      grep -oP '^ID=\\K.*' | \
      sed 's/\$/ /' > gene_ids.txt
    grep -F -A 1 -f gene_ids.txt $prokka_single_line | \
      sed '/^--\$/d' > prokka_subset.faa
    """
}

