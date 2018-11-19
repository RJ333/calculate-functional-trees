#!/usr/bin/env nextflow

prokka_input = file(params.prokka_input)

process create_single_line_fasta {
    publishDir 'results/'

    input:
      file prokka_input

    output:
      file single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka_input | tail -n +2 > single_line_fasta
    """
}

