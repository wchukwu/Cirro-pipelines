# Variant Calling from FASTQ - Custom Genome (nf-core/sarek)

DNA variant calling based on GATK4 best practices (v3.5.0)
This configuration uses a custom FASTA sequence instead of an iGenomes reference,
with the goal of making it easy to align reads to small (microbial or plasmid) genomes.
The drawbacks of this approach are:

- The reference genome index must be created at runtime, which is not appropriate for very large genomes
- Variant annotation tools are not available

Workflow Engine: Nextflow


Category: DNA Sequencing


Documentation: [https://docs.cirro.bio/pipelines/catalog-dna-sequencing/#dna-variant-calling-nf-coresarek](https://docs.cirro.bio/pipelines/catalog-dna-sequencing/#dna-variant-calling-nf-coresarek)


Source: [nf-core/sarek](nf-core/sarek)


Version: `3.5.0`


Script: `main.nf`