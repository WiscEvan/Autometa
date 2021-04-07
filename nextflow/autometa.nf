#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LENGTH_FILTER; KMERS; KMER_COVERAGE; READ_COVERAGE; ORFS; MARKERS } from './common-tasks.nf'
include { TAXON_ASSIGNMENT } from './taxonomy-tasks.nf'
include { BINNING; UNCLUSTERED_RECRUITMENT } from './binning-tasks.nf'

workflow AUTOMETA {
  take:
    metagenome

  main:
    // Perform various annotations on provided metagenome
    LENGTH_FILTER(metagenome)

    // k-mer coverage vs. read coverage
    KMER_COVERAGE(LENGTH_FILTER.out.fasta)
    // READ_COVERAGE(LENGTH_FILTER.out.fasta, fwd_reads, rev_reads, se_reads)

    ORFS(LENGTH_FILTER.out.fasta)
    MARKERS(ORFS.out.prots)
    // Perform taxon assignment with filtered metagenome
    TAXON_ASSIGNMENT(LENGTH_FILTER.out.fasta, ORFS.out.prots)
    // Now perform binning with all of our annotations.
    KMERS(TAXON_ASSIGNMENT.out.bacteria)
    // KMERS(TAXON_ASSIGNMENT.out.archaea) ... for case of performing binning on archaea
    BINNING(KMERS.out.embedded, KMER_COVERAGE.out, LENGTH_FILTER.out.gc_content, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)
    // BINNING(KMERS.out.normalized, READ_COVERAGE.out.coverage, LENGTH_FILTER.out.gc_content, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)
    // Then unclustered recruitment of any unclustered contigs using binning assignments from above.
    UNCLUSTERED_RECRUITMENT(KMERS.out.normalized, KMER_COVERAGE.out, BINNING.out.binning, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)
    // UNCLUSTERED_RECRUITMENT(KMERS.out.normalized, READ_COVERAGE.out.coverage, BINNING.out.binning, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)


  emit:
    binning = BINNING.out.binning
    binning_main = BINNING.out.main
    recruitment = UNCLUSTERED_RECRUITMENT.out.binning
    recruitment_main = UNCLUSTERED_RECRUITMENT.out.main
    all_binning_results = BINNING.out.binning | mix(UNCLUSTERED_RECRUITMENT.out) | collect
}
