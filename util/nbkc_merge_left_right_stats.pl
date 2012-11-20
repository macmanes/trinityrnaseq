#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 left.stats right.stats\n\n";

my $left_stats_file = $ARGV[0] or die $usage;
my $right_stats_file = $ARGV[1] or die $usage;

main: {

    open (my $left_fh, $left_stats_file) or die $!;
    open (my $right_fh, $right_stats_file) or die $!;

    while (1) {
        my $left_text = <$left_fh>;
        my $right_text = <$right_fh>;
        
        if ( (! $left_text) && (! $right_text)) {
            last;
        }
        elsif ($left_text xor $right_text) {
            die "Error, unequal number of lines read";
        }
        
        chomp $left_text;
        chomp $right_text;

        my ($left_median, $left_avg, $left_stdev, $left_pct, $left_acc, $left_cov) = split(/\t/, $left_text);
        my ($right_median, $right_avg, $right_stdev, $right_pct, $right_acc, $right_cov) = split(/\t/, $right_text);
        
        my $core_acc = $left_acc;
        if ($left_acc =~ /^(\S+)\/\d$/) {
            $core_acc = $1;
        }

        my $right_core = $right_acc;
        if ($right_acc =~ /^(\S+)\/\d$/) {
            $right_core = $1;
        }
        
        unless ($right_core eq $core_acc) {
            die "Error, core accs are not equivalent: [$core_acc] vs. [$right_core] reads.";
        }
        
        my $median_pair_cov = int( ($left_median+$right_median)/2 + 0.5);
        my $avg_pair_cov = int ( ($left_avg + $right_avg)/2 + 0.5);
        my $avg_stdev = int( ($left_stdev + $right_stdev)/2 + 0.5);
        
        if ($median_pair_cov == 0) { next; } ## ignore bad reads.
        

        my $avg_pair_pct = int( ($left_pct + $right_pct) / 2 + 0.5);

        print join("\t", $median_pair_cov, $avg_pair_cov, $avg_stdev, $avg_pair_pct, $core_acc) . "\n";
    }
    

    exit(0);
}
        
