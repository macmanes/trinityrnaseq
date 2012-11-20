#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/PerlLib");

use Sqlite_connect;

main: {


    unless (-s "TrinityFunctional.db") {
        die "Error, must first build Trinotate sqlite results database";
    }
    
    my $dbproc = &connect_to_db("TrinityFunctional.db");

    my @transcript_ids = &do_sql($dbproc, "select transcript_id from TrinityTranscript");

    # print header
    print join("\t", "#ID", "TopBlastHit", "Pfam", "SignalP", "TmHMM") . "\n";
    
    foreach my $id (@transcript_ids) {
        
        my $blast_info = &get_blast_results($dbproc, $id);
        my $pfam_info = &get_pfam_info($dbproc, $id);
        
        my $signalP_info = &get_signalP_info($dbproc, $id);
        
        my $TmHMM_info = &get_TmHMM_info($dbproc, $id);
        
        print join("\t", $id, $blast_info, $pfam_info, $signalP_info, $TmHMM_info) . "\n" if ($signalP_info ne ".");
    }
    
    exit(0);


}

####
sub get_TmHMM_info {
    my ($dbproc, $id) = @_;

    my $query = "select Score, PredHel, Topology from tmhmm where queryprotid = ? and PredHel != 'PredHel=0' ";
    my $result = &first_result_sql($dbproc, $query, $id);
    
    if ($result) {
        my $tmhmm_line = join("`", @$result);
        return($tmhmm_line);
    }
    else {
        return(".");
    }
}


####
sub get_signalP_info {
    my ($dbproc, $id) = @_;

    my $query = "select start, end, score, prediction from SignalP where query_prot_id = ?";
    my $result = &first_result_sql($dbproc, $query, $id);
    if ($result) {
        my $sigP_line = "sigP:" . join("`", @$result);
        return($sigP_line);
    }
    else {
        return(".");
    }
}


###
sub get_pfam_info {
    my ($dbproc, $id) = @_;

    my $query = "select pfam_id, HMMERDomain, HMMERTDomainDescription, QueryStartAlign, QueryEndAlign, ThisDomainEvalue "
        . " from HMMERDbase "
        . " where QueryProtID = ? and TCSeqeunceMet = 'YES' "
        . " order by QueryStartAlign ";
    

    my @results = &do_sql_2D($dbproc, $query, $id);

    if (@results) {
        my @encoded_hits;
        foreach my $result (@results) {
            my ($pfam_id, $domain, $domain_descr, $start, $end, $evalue) = @$result;
            
            my $hit = join("^", $pfam_id, $domain, $domain_descr, "$start-$end", "E:$evalue");
            push (@encoded_hits, $hit);
        }
        
        my $result_line = join("`", @encoded_hits);
        
        return($result_line);
    }
    else {
        return(".");
    }
}




####
sub get_blast_results {
    my ($dbproc, $id) = @_;

    my $query = "select FullAccesion, UniprotSearchString, QueryStart, QueryEnd, HitStart, HitEnd, PercentIdentity, Evalue "
        . "from BlastDbase where TrinityID = ?";

    my $result = &first_result_sql($dbproc, $query, $id);

    if ($result) {
        my ($FullAccesion, $UniprotSearchString, $QueryStart, $QueryEnd, $HitStart, $HitEnd, $PercentIdentity, $Evalue) = @$result;
        my $taxonomy_string = &get_taxonomy_string($dbproc, $UniprotSearchString);
        $taxonomy_string =~ s/[\`\s]+/ /g; 
        $taxonomy_string =~ s/^\s+|\s+$//g;
        
        my $description_line = &get_description_line($dbproc, $UniprotSearchString);
        
        $description_line =~ s/[\`\s]+/ /g; # using backtics as delimiters, and dont want tabs or newlines to disrupt formatting.
        $description_line =~ s/^\s+|\s+$//g;  
        
        ## encode the result
        my $ret_val = join("`", $FullAccesion, $UniprotSearchString, "Q:$QueryStart-$QueryEnd,H:$HitStart-$HitEnd", "$PercentIdentity%ID", "E:$Evalue", $description_line, $taxonomy_string);

        return($ret_val);
    }
    else {
        return(".");
    }

}

####
sub get_taxonomy_string {
    my ($dbproc, $uniprot_acc) = @_;

    my $query = "select t.TaxonomyValue from TaxonomyIndex t, UniprotIndex u "
        . " where t.NCBITaxonomyAccession = u.LinkId "
        . " and u.AttributeType = 'T' "
        . " and u.Accesion = ? ";

    my $result = &very_first_result_sql($dbproc, $query, $uniprot_acc);

    return($result);
}

####
sub get_description_line {
    my ($dbproc, $uniprot_acc) = @_;

    my $query = "select LinkId from UniprotIndex u where u.Accesion = ? and u.AttributeType = 'D'";
    my $description = &very_first_result_sql($dbproc, $query, $uniprot_acc);
    $description =~ s/^\s+|\s+$//g; # trim leading/trailing ws
    
    return($description);
}



        
