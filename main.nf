#!/usr/bin/env nextflow

prokka_protein_fasta = file(params.prokka_protein_fasta)
prokka_gff = file(params.prokka_gff)
gene_name = params.gene_name
taxa_of_interest = file(params.tax_list)

println """\
         workflow: build phylogenetic trees of proteins sequences
         ========================================================
         gene_name           : ${params.gene_name}
         prokka_protein_fasta: ${params.prokka_protein_fasta}
         prokka_gff          : ${params.prokka_gff}
         taxa_of_interest    : ${params.tax_list}
         """
         .stripIndent()

process create_single_line_fasta {
    /*
      This process joins all sequence lines and creates a single line fasta file.
    */

    input:
      file prokka_protein_fasta

    output:
      file 'prokka.sl.faa' into single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka_protein_fasta | tail -n +2 > prokka.sl.faa
    """
}

process subset_gff {
    /*
      This process subsets all lines from the prokka gff file and only keeps those with the same gene name as set in
      the nextflow config file.
    */

    input:
      file prokka_gff
    output:
      file 'prokka_subset.gff' into subset_gff

    """
    cat $prokka_gff | \
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

    input:
      file single_line_fasta
      file subset_gff
    output:
      file 'prokka_subset.faa' into subset_prokka

    """
    cat $subset_gff | \
      awk -F '[\t;]' '\$3 == "CDS" {print \$9}' | \
      grep -oP '^ID=\\K.*' | \
      sed 's/\$/ /' > gene_ids.txt
    grep -F -A 1 -f gene_ids.txt $single_line_fasta | \
      sed '/^--\$/d' > prokka_subset.faa
    """
}

process download_sequences_from_uniprot {
    /*
      This process downloads the genes of interest of the taxa of interest (input file with TaxIDs) from UniProt-TrEMBL.
    */

    input:
      file taxa_of_interest
    output:
      file 'uniprot_seq.faa' into subset_target_genes

    """
    get_protein_sequences.py -i $taxa_of_interest -o uniprot_seq.faa -g $gene_name
    """
}

process merge_protein_sequences {
    /*
      This process merges the Prokka sequences and the UniProt sequences into one file.
    */

    input:
      file subset_prokka
      file subset_target_genes
    output:
      file 'merged_seqs.faa' into merged_proteins

    """
    cat $subset_prokka $subset_target_genes > merged_seqs.faa
    """
}

process remove_spaces_from_fasta_headers {
    /*
      This process removes spaces from the fasta headers and replaces them with "%".
    */

    input:
      file merged_proteins
    output:
      file 'merged_seqs.nospaces.faa' into merged_proteins_nospaces

    """
    sed '/^>/s/ /%/g' $merged_proteins > merged_seqs.nospaces.faa
    """
}

process remove_duplicates {
    /*
      This process removes the duplicated protein sequences from the merged protein fasta file.
    */

    input:
      file merged_proteins_nospaces
    output:
      file 'merged_seqs.nodup.faa' into merged_proteins_no_dup

    """
    cd-hit-dup -i $merged_proteins_nospaces -o merged_seqs.nodup.faa
    """
}

process run_msa {
    /*
      This process performes a multiple sequence alignment of the concatenated sequences with MAFFT.
    */

    input:
      file merged_proteins_no_dup
    output:
      file 'merged_seqs.msa.faa' into msa

    """
    mafft --auto $merged_proteins_no_dup > merged_seqs.msa.faa
    """
}

process fasta_to_phylip {
    /*
      This process converts fasta to phylip format.
    */

    input:
      file msa
    output:
      file 'merged_seqs.phylip' into phylip

    """
    fasta_to_phylip.py $msa merged_seqs.phylip
    """
}

process remove_special_chars_from_phylip {
    /*
      This process removes forbidden characters from the phylip format.
    */

    publishDir 'results/'

    input:
      file phylip
    output:
      file 'merged_seqs.clean.phylip' into phylip_clean

    """
    remove_special_characters_from_phylip.py $phylip > merged_seqs.clean.phylip
    """
}

