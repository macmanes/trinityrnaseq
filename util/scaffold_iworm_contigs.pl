#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 nameSorted.sam iworm_fasta_file\n\n";

my $name_sorted_sam_file = $ARGV[0] or die $usage;
my $iworm_fasta = $ARGV[1] or die $usage;


my $SAM_OFH;


main: {

    
    open ($SAM_OFH, ">scaffolding_entries.sam") or die $!;
    
    my %iworm_acc_to_fasta_index;
    {
        
        my $counter = 0;
        open (my $fh, $iworm_fasta) or die "Error, cannot open file $iworm_fasta";
        while (<$fh>) {
            if (/^>(\S+)/) {
                my $acc = $1;
                $iworm_acc_to_fasta_index{$acc} = $counter; # starts at zero
                $counter++;
            }
        }
        close $fh;
    }
    

    my %paired_iworm_contigs;

    my $prev_core_acc = "";
    my %end_to_iworm;


    open (my $fh, $name_sorted_sam_file) or die "Error, cannot open file $name_sorted_sam_file";
    while (<$fh>) {
        my $line = $_;
        chomp;
        my @x = split(/\t/);
        my $read_acc = $x[0];
        my $iworm_acc = $x[2];

        if ($iworm_acc eq '*') { 
            # unmapped read
            next;
        }
        
        my $core_acc;
        my $frag_end;
        

        if ($read_acc =~ /^(\S+)\/([12])$/) { 
            
            $core_acc = $1;
            $frag_end = $2;
        }
        else {
            # must have mixed in a single read with the pairs...
            print STDERR "warning, ignoring read: $read_acc since cannot decipher if /1 or /2 of a pair.\n";
            next;
        }
        
        if ($core_acc ne $prev_core_acc) {

            &examine_frags(\%end_to_iworm, \%paired_iworm_contigs);
            %end_to_iworm = ();
            
        }

        $end_to_iworm{$frag_end}->{$iworm_acc} = $line;
        
        $prev_core_acc = $core_acc;

    }
    close $fh;

    ## get last one
    &examine_frags(\%end_to_iworm, \%paired_iworm_contigs);

    
    foreach my $pairing (reverse sort {$paired_iworm_contigs{$a}<=>$paired_iworm_contigs{$b}} keys %paired_iworm_contigs) {

        my ($iworm_acc_A, $iworm_acc_B) = split(/\t/, $pairing);
      
        my $index_A = $iworm_acc_to_fasta_index{$iworm_acc_A};
        unless (defined $index_A) {
            die "Error, no index for iworm acc: $iworm_acc_A";
        }
        my $index_B = $iworm_acc_to_fasta_index{$iworm_acc_B};
        unless (defined $index_B) {
            die "Error, no index for iworm acc: $iworm_acc_B";
        }
        
    
        print join("\t", $iworm_acc_A, $index_A, $iworm_acc_B, $index_B, $paired_iworm_contigs{$pairing}) . "\n";
    
    }
    

    close $SAM_OFH;
    
    exit(0);
}

####
sub examine_frags {
    my ($end_to_iworm_href, $paired_iworm_contigs_href) = @_;

    my @ends = keys %$end_to_iworm_href;

    unless (scalar @ends == 2) {
        ## no pairs
        return;
    }
    
    
    my @iworm_left = keys %{$end_to_iworm_href->{1}};
    my @iworm_right = keys %{$end_to_iworm_href->{2}};

    if (scalar(@iworm_left) == 1 && scalar(@iworm_right) == 1
        &&
        $iworm_left[0] ne $iworm_right[0]) {

        my $pair = join("\t", sort (@iworm_left, @iworm_right));
        

        #print STDERR "Got pairing: $pair\n";
        $paired_iworm_contigs_href->{$pair}++;

        foreach my $sam_line (values %{$end_to_iworm_href->{1}}, values %{$end_to_iworm_href->{2}}) {
            print $SAM_OFH $sam_line;
        }

        
    }
    
    
    return;
}

        
