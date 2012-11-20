#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


my $BASEDIR = "$FindBin::Bin/../../";

my @conditions = qw(ds hs log plat);

## first, gunzip the inputs as needed.

my $reads_ALL_left_fq = "reads.ALL.left.fq";
my $reads_ALL_right_fq = "reads.ALL.right.fq";

my $REGENERATE_ALL_FQ = 1;
if (-s $reads_ALL_left_fq && -s $reads_ALL_right_fq) {
    $REGENERATE_ALL_FQ = 0;
}

my %condition_to_fastq;

## concatenate all entries before running Trinity

foreach my $condition (@conditions) {
    my $left_fq_file = "rnaseq_reads/Sp_${condition}.10k.left.fq";
    my $right_fq_file = "rnaseq_reads/Sp_${condition}.10k.right.fq";
 
    $condition_to_fastq{$condition}->{left} = $left_fq_file;
    $condition_to_fastq{$condition}->{right} = $right_fq_file;
    
    
    if (-s "$left_fq_file.gz" && ! -s $left_fq_file) {
        &process_cmd("gunzip -c $left_fq_file.gz > $left_fq_file");
    }

    if (-s "$right_fq_file.gz" && ! -s $right_fq_file) {
        &process_cmd("gunzip -c $right_fq_file.gz > $right_fq_file");
    }

    if ($REGENERATE_ALL_FQ) {

        &process_cmd("cat $left_fq_file >> $reads_ALL_left_fq");
        &process_cmd("cat $right_fq_file >> $reads_ALL_right_fq");
        
    }
 
}    
    
## Run Trinity:
my $cmd = "$BASEDIR/Trinity.pl --seqType fq --SS_lib_type RF --JM 1G "
    . " --left $reads_ALL_left_fq --right $reads_ALL_right_fq --CPU 4 "
    #. " --monitoring"
    ;


&process_cmd($cmd) unless (-s "trinity_out_dir/Trinity.fasta");;


## Align each of the read sets back against the Trinity assemblies
foreach my $condition (@conditions) {
    my $left_fq_file = $condition_to_fastq{$condition}->{left};
    my $right_fq_file = $condition_to_fastq{$condition}->{right};

    my $cmd = "$BASEDIR/util/alignReads.pl --target trinity_out_dir/Trinity.fasta "
        . " --left $left_fq_file --right $right_fq_file --SS_lib_type RF --seqType fq "
        . " --aligner bowtie --output ${condition}_bowtie -- -p 4 ";

    &process_cmd($cmd) unless (-s "${condition}_bowtie/${condition}_bowtie.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam");
    
}


## Run RSEM for each:
use Cwd;
my $workdir = cwd();

foreach my $condition (@conditions) {
    
    chdir($workdir) or die "Error, cannot CD to $workdir";
    
    my $bowtie_dir = "${condition}_bowtie";
    chdir $bowtie_dir or die "Error, cannot cd to $bowtie_dir";
    
    ## run RSEM
    my $cmd = "$BASEDIR/util/RSEM_util/run_RSEM.pl --transcripts target.fa "
        . " --name_sorted_bam ${condition}_bowtie.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam"
        . " --paired --SS_lib_type RF ";
    
    &process_cmd($cmd) unless (-s "RSEM.isoforms.results");
}

chdir($workdir) or die "Error, cannot CD to $workdir";


## Create the effective lengths file:
&process_cmd("cat ds_bowtie/RSEM.isoforms.results | cut -f1,3,4 > Trinity.trans_lengths.txt") unless (-s "Trinity.trans_lengths.txt");
&process_cmd("cat ds_bowtie/RSEM.genes.results | cut -f1,3,4 > Trinity.components_lengths.txt") unless (-s "Trinity.components_lengths.txt");


my @iso_links;
my @gene_links;
## Create the matrix of count data.
foreach my $condition (@conditions) {
    
    my $rsem_iso_file = "${condition}_bowtie/RSEM.isoforms.results";
    my $rsem_gene_file = "${condition}_bowtie/RSEM.genes.results";

    my $cmd = "ln -sf $rsem_iso_file ${condition}_trans";
    &process_cmd($cmd) unless (-e "${condition}_trans");
    
    push (@iso_links, "${condition}_trans");

    $cmd = "ln -sf $rsem_gene_file ${condition}_gene";
    &process_cmd($cmd) unless (-e "${condition}_gene");

    push (@gene_links, "${condition}_gene");
}

# make iso matrix
$cmd = "$BASEDIR/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl " . join (" " , @iso_links) . " > Trinity_trans.counts.matrix";
&process_cmd($cmd) unless (-s "Trinity_trans.counts.matrix");

# make the 'gene' matrix
$cmd = "$BASEDIR/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl " . join (" " , @gene_links) . " > Trinity_components.counts.matrix"; 
&process_cmd($cmd) unless (-s "Trinity_components.counts.matrix");


## run EdgeR
foreach my $target_type ("trans", "components") {
    
    chdir $workdir or die "Error, cannot cd to $BASEDIR";
    
    my $edgeR_dir = "edgeR_${target_type}";

    my $cmd = "$BASEDIR/Analysis/DifferentialExpression/run_EdgeR.pl --matrix Trinity_${target_type}.counts.matrix "
        . " --transcript_lengths Trinity.${target_type}_lengths.txt "
        . " --no_eff_length "
        . " --output $edgeR_dir ";
    
    &process_cmd($cmd) unless (-d $edgeR_dir);

    chdir $edgeR_dir or die "Error, cannot cd to $edgeR_dir";
    
    ## extract the diff. expressed transcripts.



}


print "Done.\n";
    



exit(0);


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}


