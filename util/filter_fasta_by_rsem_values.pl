#!/usr/bin/env perl

=head1 NAME

filter_fasta_by_rsem_values.pl - Use RSEM relative abundance values to filter a transcript assembly FASTA file

=head1 SYNOPSIS

USAGE: filter_fasta_by_rsem_values.pl 
            --rsem_output=/path/to/Trinity.RSEM.fpkm
            --fasta=/path/to/Trinity.fasta
            --output=/path/to/output.fasta
          [ --fpkm_cutoff=1200
            --perc_comp_fpkm_cutoff=1.0
          ]

=head1 OPTIONS

B<--rsem_output,-r>
    This file is the output file from util/RSEM_util/summarize_RSEM_fpkm.p

B<--fasta,-f>
    The FASTA file representing transcript assemblies from the primary Trinity run.

B<--output,-o>
    The output FASTA file to be created.

B<--fpkm_cutoff,-c>
    Optional.  Will filter transcripts, keeping those with FPKM values equal to or greater than this.

B<--perc_comp_fpkm_cutoff,-p>
    Optional.  Will filter transcripts, keeping those with % component FPKM values equal to or greater than this.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Trinity contains a pipeline for read alignment, visualization and relative abundance 
estimation here:

    http://trinityrnaseq.sourceforge.net/analysis/align_visualize_quantify.html

The product of this is the file 'Trinity.RSEM.fpkm', which contains calculated values
such as FPKM and %comp_FPKM.  See INPUT for more on this.

The user can use this file to filter the source transcript assembly FASTA file.

=head1  INPUT

The input RSEM file looks like this:

    #Total fragments mapped to transcriptome: 3233319.93999999
    transcript      length  eff_length      count   fraction        fpkm    %comp_fpkm
    comp3119_c0_seq1        370     71      14.59   5.02e-05        63.55   31.43
    comp3119_c0_seq2        349     50      22.41   8.69e-05        138.62  68.57

    comp2254_c0_seq1        2157    1858    9843.00 3.17e-03        1638.45 100.00

    comp186061_c0_seq1      765     466     19.00   2.08e-05        12.61   100.00

The corresponding input FASTA has headers like:

    >comp8_c0_seq1 len=523 path=[333:0-522]
    >comp9_c0_seq1 len=281 path=[259:0-280]
    >comp10_c0_seq1 len=635 path=[613:0-634]
    >comp10_c1_seq1 len=886 path=[1225:0-885]
    >comp11_c0_seq1 len=422 path=[1:0-421]

To do the actual filtering, you need to pass a --fpkm_cutoff or --perc_comp_fpkm_cutoff value.
If both are passed, both will be applied.


=head1  OUTPUT

The output is a FASTA file with the same headers and sequence, only filtered based
on the parameters passed.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'rsem_output|r=s',
                          'fasta|f=s',
                          'output|o=s',
                          'fpkm_cutoff|c=s',
                          'perc_comp_fpkm_cutoff|p=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

_log( "INFO: cutoffs: fpkm_cutoff=(" . ($options{fpkm_cutoff} || '') . "), perc_comp_fpkm_cutoff=(" . ($options{perc_comp_fpkm_cutoff} || '') . ")");

_log( "INFO: Opening RSEM file ($options{rsem_output})" );
my $rsem = load_rsem_output( $options{rsem_output} );

_log( "INFO: creating output file: ($options{output})" );
open(my $ofh, ">$options{output}") || logdie("ERROR: failed to create output file: $!");

_log( "INFO: reading input FASTA file: ($options{fasta})" );
open( my $ifh, $options{fasta} ) || logdie("ERROR: failed to read input FASTA file: $!");

my $keep = 0;

while ( my $line = <$ifh> ) {
    if ( $line =~ /^\>(\S+)/ ) {
        my $seq_id = $1;
        
        ## make sure we have cutoffs for this
        if ( ! exists $$rsem{$seq_id} ) {
            logdie("ERROR: found a transcript ID ($seq_id) in the FASTA that wasn't in the RSEM file");
        }
        
        ## optimism first
        $keep = 1;
        
        if ( defined $options{fpkm_cutoff} && $$rsem{$seq_id}{fpkm} < $options{fpkm_cutoff} ) {
            $keep = 0;
        }
        
        if ( defined $options{perc_comp_fpkm_cutoff} && $$rsem{$seq_id}{comp_fpkm} < $options{perc_comp_fpkm_cutoff} ) {
            $keep = 0;
        }
    }
    
    print $ofh $line if $keep;
}


exit(0);


sub load_rsem_output {
    my $file = shift;
    
    ## relative abundance data.  looks like:
    #   $rel{'comp3119_c0_seq1'} = { fpkm      => 63.55,
    #                                comp_fpkm => 31.43
    #                              }
    my %rel = ();
    
    open(my $ifh, $file) || logdie("ERROR: failed to read rsem_output file: $!");
    
    while (my $line = <$ifh>) {
        chomp $line;
        
        next if $line =~ /^\s*$/;
        my @cols = split("\t", $line);
        
        if ( scalar @cols == 7 && $cols[5] ne 'fpkm' ) {
            if ( exists $rel{$cols[0]} ) {
                logdie("ERROR: found more than one entry in the RSEM file for $cols[0]");
            } else {
                $rel{$cols[0]} = { fpkm => $cols[5], comp_fpkm => $cols[6]};
            }
        }
    }
    
    _log("INFO: Loaded RSEM values for " . scalar(keys %rel) . " transcripts");
    
    return \%rel;
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub logdie {
    my $msg = shift;
    
    print STDERR "$msg\n";
    print $logfh "$msg\n" if $logfh;
    
    exit(1);
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
+    my @required = qw( rsem_output fasta output );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            logdie("ERROR: --$option is a required option");
        }
    }
    
    ## either fpkm_cutoff or perc_comp_fpkm_cutoff need to be defined (or both)
    if ( ! defined $$options{fpkm_cutoff} && ! defined $$options{perc_comp_fpkm_cutoff} ) {
        logdie("ERROR: You must define either --fpkm_cutoff or --perc_comp_fpkm_cutoff for filtering");
    }
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
