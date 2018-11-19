#!/usr/bin/env nextflow

prokka_input = file(params.prokka_input)

println """\
         workflow: build phylogenetic trees of proteins sequences
         ========================================================
         prokka_input: ${params.prokka_input}
         """
         .stripIndent()

process create_single_line_fasta {
    publishDir 'results/'

    input:
      file prokka from prokka_input

    output:
      file 'prokka.sl.faa' into single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka | tail -n +2 > prokka.sl.faa
    """
}

