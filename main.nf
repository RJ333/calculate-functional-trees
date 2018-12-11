#!/usr/bin/env nextflow

prokka_protein_fasta = file(params.prokka_protein_fasta)
prokka_gff = file(params.prokka_gff)
gene_name = params.gene_name
heatmap_scale = params.heatmap_scale
taxa_of_interest = file(params.tax_list)
kallisto_matrix = file(params.kallisto_matrix)
raxml_randomseed = params.raxml_randomseed
raxml_model = params.raxml_model
raxml_algorithm = params.raxml_algorithm
raxml_runs = params.raxml_runs

println """\
         workflow: build phylogenetic trees of proteins sequences
         ========================================================
         gene_name           : ${params.gene_name}
         prokka_protein_fasta: ${params.prokka_protein_fasta}
         prokka_gff          : ${params.prokka_gff}
         taxa_of_interest    : ${params.tax_list}
         kallisto_matrix     : ${params.kallisto_matrix}
         heatmap_scale 		 : ${params.heatmap_scale}
		 raxml_randomseed    : ${params.raxml_randomseed}
         raxml_model         : ${params.raxml_model}
         raxml_algorithm     : ${params.raxml_algorithm}
         raxml_runs          : ${params.raxml_runs}
         """
         .stripIndent()

process create_single_line_fasta {
    /*
      This process joins all sequence lines and creates a single line fasta file.
    */

    input:
      file prokka_protein_fasta

    output:
      file "${prokka_protein_fasta.simpleName}.sl.faa" into single_line_fasta

    """
    perl -pe '/^>/ ? print "\n" : chomp' $prokka_protein_fasta | tail -n +2 > ${prokka_protein_fasta.simpleName}.sl.faa
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
      file "${prokka_gff.simpleName}_subset.gff" into subset_gff

    """
    cat $prokka_gff | \
      grep "gene=$gene_name" | \
      awk '\$3 == "CDS" {print \$0}' > ${prokka_gff.simpleName}_subset.gff
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
      file "${prokka_protein_fasta.simpleName}_subset.faa" into subset_prokka

    """
    cat $subset_gff | \
      awk -F '[\t;]' '\$3 == "CDS" {print \$9}' | \
      grep -oP '^ID=\\K.*' | \
      sed 's/\$/ /' > gene_ids.txt
    grep -F -A 1 -f gene_ids.txt $single_line_fasta | \
      sed '/^--\$/d' > ${prokka_protein_fasta.simpleName}_subset.faa
    """
}

process download_sequences_from_uniprot {
    /*
      This process downloads the genes of interest of the taxa of interest (input file with TaxIDs) from UniProt-TrEMBL.
    */

    publishDir 'results/'

    input:
      file taxa_of_interest
    output:
      file "uniprot_seq.faa" into subset_target_genes

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
      file "merged_seqs.faa" into merged_proteins

    """
    cat $subset_prokka $subset_target_genes > merged_seqs.faa
    """
}

process remove_duplicates {
    /*
      This process removes the duplicated protein sequences from the merged protein fasta file.
    */

    input:
      file merged_proteins 
    output:
      file "${merged_proteins.simpleName}.nodup.faa" into merged_proteins_no_dup

    """
    cd-hit-dup -i $merged_proteins -o ${merged_proteins.simpleName}.nodup.faa
    """
}

process run_msa {
    /*
      This process performes a multiple sequence alignment of the concatenated sequences with MAFFT.
    */

    input:
      file merged_proteins_no_dup
    output:
      file "${merged_proteins_no_dup.simpleName}.msa.faa" into msa

    """
    mafft --auto $merged_proteins_no_dup > ${merged_proteins_no_dup.simpleName}.msa.faa
    """
}

process fasta_to_phylip {
    /*
      This process converts fasta to phylip format.
    */

    input:
      file msa
    output:
      file "${msa.simpleName}.phylip" into phylip

    """
    fasta_to_phylip.py $msa ${msa.simpleName}.phylip
    """
}

process remove_special_chars_from_phylip {
    /*
      This process removes forbidden characters from the phylip format.
    */

    input:
      file phylip
    output:
      file "${phylip.simpleName}.clean.phylip" into phylip_clean

    """
    remove_special_characters_from_phylip.py $phylip > ${phylip.simpleName}.clean.phylip
    """
}

process run_raxml {
    /*
      This process runs raxml of the preprocessed phylip file.
    */

    cpus 28

    publishDir 'results/'

    input:
      file phylip_clean
    output:
      file "RAxML_bestTree.${gene_name}" into raxml_tree

    """
    raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -p $raxml_randomseed -s $phylip_clean -m $raxml_model -n $gene_name \
      -x $raxml_randomseed -f $raxml_algorithm -N $raxml_runs
    """
}

process plot_tree {
    /*
      This process uses ggtree to plot the tree together with the kallisto matrix.
    */

    publishDir 'results/'

    input:
      file raxml_tree
      file kallisto_matrix
    output:
      file "tree_${gene_name}.pdf" into tree_plot

    """
    plot_tree.R -t $raxml_tree -k $kallisto_matrix -s ${params.heatmap_scale} -o tree_${gene_name}.pdf
    """
}
