# Calculate trees with prokka and references

The idea of this program is threefold:

### First: 

Generate a workflow to perform a multiple sequence alignment of specific genes. The data is provided by prokka annotation and based on metagenomes. The tree based on this data helps to search for clusters within gene groups.

### Second: 

In a metagenomic dataset the taxonomic annotation is often lacking or unreliable. If we provide reference genes of the same group from isolated organisms, we can use sequence similarity as indication. Reference genes can be chosen based on the function you investigate and/or from the organisms you know which are present (e.g. by 16S data). It should be reminded that this is not a phylogenetic analysis: functional genes undergo evolutional selection and might be similar for this reason, not due to common ancestors.

### Third: 

If the experimental setup allows for it, we can later combine the tree information with abundance information. This allows us in the end to check whether specific gene clusters show a specific behaviour in abundance.

### Theoretical example: 

the phnJ gene is responsible for degradation of phosphonates in many bacteria. But not all bacteria containing it are capable of degrading the phosphonate glyphosate. To figure out if my glyphosate-inpacted metagenomes contain the "right" version of phnJ, we build trees of the phnJ gene and include reference genes from known glyphosate degraders. This tree is then linked to the phnJ abundance data to see if specific clusters become more abundant after glyphosate addition (arguebly those capable of degradation).

## How to run

```bash
nextflow run . --prokka_input <prokka fasta file>
```
