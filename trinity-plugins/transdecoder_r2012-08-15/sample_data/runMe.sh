#!/bin/bash -e

if [ -e test.genome.fasta.gz ] && [ ! -e test.genome.fasta ]; then
    gunzip -c test.genome.fasta.gz > test.genome.fasta
fi

if [ -e test.tophat.sam.gz ] && [ ! -e test.tophat.sam ]; then
    gunzip -c test.tophat.sam.gz > test.tophat.sam
fi

if [ -e transcripts.gtf.gz ] && [ ! -e transcripts.gtf ]; then
    gunzip -c transcripts.gtf.gz > transcripts.gtf
fi


## generate alignment gff3 formatted output
../util/cufflinks_gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3

## generate transcripts fasta file
../util/cufflinks_gtf_genome_to_cdna_fasta.pl transcripts.gtf test.genome.fasta > transcripts.fasta 

## find likely ORFs, only doing top strand here.
../transcripts_to_best_scoring_ORFs.pl -t transcripts.fasta -S

## convert to genome coordinates
../cdna_alignment_orf_to_genome_orf.pl best_candidates.eclipsed_orfs_removed.gff3 transcripts.gff3 transcripts.fasta > best_candidates.eclipsed_orfs_removed.genome.gff3


## make bed files for viewing with GenomeView

# covert cufflinks gtf to bed
../util/cufflinks_gtf_to_bed.pl transcripts.gtf > transcripts.bed

# convert the genome-based gene-gff3 file to bed
../util/gff3_file_to_bed.pl best_candidates.eclipsed_orfs_removed.genome.gff3 > best_candidates.eclipsed_orfs_removed.genome.bed



echo
echo
echo Done!  Coding region genome annotations provided as: best_candidates.eclipsed_orfs_removed.genome.gff3
echo
echo 

exit 0
