#!/usr/bin/env perl

use strict;
use warnings;
use threads;
no strict qw(subs refs);

use FindBin;
use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/PerlLibAdaptors");
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


open (STDERR, ">&STDOUT");  ## capturing stderr and stdout in a single stdout stream


# Site specific setup

my $CPU_MAX = 64;       # Set this to a high number if you really do want to use more processors.  I suspect it may only cause trouble, though.

my $IWORM_KMER_SIZE = 25;
my $MIN_IWORM_LEN = 25;

# option list:
my ($seqType, $left_file, $right_file, $single_file, $SS_lib_type, $min_contig_length,
    $group_pairs_distance, $jaccard_clip, $show_advanced_options,
    $output_directory, $prep_only
    );

# What is allowed for the options. Put string to be displayed in '%allowed'; this
#   will be showed to the user via  help and on error.   Keys are the variable names.
#   Actual hash to be used for checking is auto-generated. Fancy regex inside map 
#   is just to get rid of the syntaxical sugar 'or' in the display string.

my %allowed = 
    ( seqType       => 'cfa, cfq, fa, or fq'
    , kmer_method   => 'jellyfish, meryl, or inchworm'
    );

my %allowed_check;
foreach my $all (keys %allowed) {
    my %h = map { (my $s = $_) =~ s/^or //; $s => 1 } split ', ', $allowed{$all};
    $allowed_check{$all} = \%h;
}

# defaults:

$output_directory = "trinity_out_dir";

$min_contig_length = 200;
$group_pairs_distance = 500;
my $path_reinforcement_distance = 75;

my $NO_RUN_BUTTERFLY_FLAG = 0;
my $RERUN_BUTTERFLY_FLAG = 0;
my $bfly_opts = "";
my $bflyHeapSpaceMax = "20G";
my $bflyHeapSpaceInit = "1G";

my $min_kmer_cov = 1;
my $meryl_opts = "";
my $inchworm_cpu = 6;

my $min_percent_read_iworm_kmers = -1;

my $CPU = 2;
my $bflyCPU;
my $bflyCalculateCPU = 0;
my $bflyGCThreads;

my $long_reads = "";

my $max_number_of_paths_per_node = 10;
my $lenient_path_extension = 0;


## ADVANCED OPTIONS:
     # todo  add some.


my $no_meryl_flag = 0;

## Chrysalis opts
my $min_glue = 2;
my $min_iso_ratio = 0.05;
my $glue_factor = 0.05;
my $max_reads_per_graph = 20000000;
my $max_reads_per_loop = 1000000;
my $min_pct_read_mapping = 0;
my $NO_RUN_QUANTIFYGRAPH_FLAG = 0;
my $NO_RUN_CHRYSALIS_FLAG = 0;

my $help_flag;
my $SHOW_CITATION_FLAG = 0;

my $VERSION = "BLEEDING_EDGE"; 
my $show_version_flag = 0;

## Kmer methods
my $kmer_method = "";

## Jellyfish
my $max_memory;


## Grid computing options:
my $grid_computing_module;

## Performance monitoring options 
my $pm_logfile = "Trinity.timing";
my $pm_trinity_startstring;
my $pm_trinity_endstring;
my $pm_trinity_start=0;
my $pm_trinity_end=0;
my $pm_inchworm_start=0;
my $pm_inchworm_end=0;
my $pm_chrysalis_start=0;
my $pm_chrysalis_end=0;
my $pm_butterfly_start=0;
my $pm_butterfly_end=0;
my $pm_left_fa_size=0;
my $pm_right_fa_size=0;
my $pm_single_fa_size=0;
my $pm_trinity_fa_size=0;
my $pm_trinity_arguments="";
my $pm_inchworm_kmers=0;
my $pm_read_count=0;

my $run_with_collectl = 0;
my $collectl_output_directory = "collectl";
my $collectl_pid = 0;
my $collectl_out = "";
my $collectl_titlename = "";
my $start_dir = cwd();

## misc other opts, mostly for testing purposes
my $run_as_paired_flag = 0;  ## in case we have paired reads in single fasta file, already oriented.
my $weldmer_size = 48;
my $FORCE_INCHWORM_KMER_METHOD = 0; 

my $NO_TRIPLET_LOCK = 0;

# Note: For the Trinity logo below the backslashes are quoted in order to keep
#   them from quoting the character than follows them.  "\\" keeps "\ " from occuring.

my $usage = <<_EOUSAGE_;


###############################################################################
#
#     ______  ____   ____  ____   ____  ______  __ __
#    |      ||    \\ |    ||    \\ |    ||      ||  |  |
#    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
#    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
#      |  |  |    \\  |  | |  |  | |  |   |  |  |___, |
#      |  |  |  .  \\ |  | |  |  | |  |   |  |  |     |
#      |__|  |__|\\_||____||__|__||____|  |__|  |____/
#
###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( $allowed{seqType} )
#  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for 
#                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char 
#
#  If paired reads:
#      --left  <string>    :left reads
#      --right <string>    :right reads
#
#  Or, if unpaired reads:
#      --single <string>   :single reads   (note, if single file contains pairs, can use flag: --run_as_paired )
#
####################################
##  Misc:  #########################
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#  --long_reads <string>           :fasta file containing corrected pac bio reads
# 
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( "$output_directory" )
#  --CPU <int>                     :number of CPUs to use, default: $CPU
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=$min_contig_length)
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format
#                                   for reads).
#                                   (note: jaccard_clip is an expensive
#                                   operation, so avoid using it unless
#                                   necessary due to finding excessive fusion
#                                   transcripts w/o it.)
#  
#  --prep                          :Only prepare files (high I/O usage) and stop before kmer counting.
#
#  --no_cleanup                    :retain all intermediate input files.
#  --full_cleanup                  :only retain the Trinity fasta file, rename as \${output_dir}.Trinity.fasta
#
#  --cite                          :get the Trinity literature citation and those of tools leveraged within.
#  --monitoring                    :use collectl to monitor all steps of Trinity
#
#  --version                       :reports Trinity version ($VERSION) and exits.
#
####################################################
# Inchworm and K-mer counting-related options: #####
#
#  --min_kmer_cov <int>           :min count for K-mers to be assembled by
#                                  Inchworm (default: $min_kmer_cov)
#  --inchworm_cpu <int>           :number of CPUs to use for Inchworm, default is min(6, --CPU option)
#
###################################
# Chrysalis-related options: ######
#
#  --min_glue <int>               :min number of reads needed to glue two inchworm contigs
#                                  together. (default: $min_glue) 
#  --min_iso_ratio <float>        :min fraction of average kmer coverage between two iworm contigs
#                                  required for gluing.  (default: $min_iso_ratio)
#  --glue_factor <float>          :fraction of max (iworm pair coverage) for read glue support (default: $glue_factor)
#  --max_reads_per_graph <int>    :maximum number of reads to anchor within
#                                  a single graph (default: $max_reads_per_graph)
#  --max_reads_per_loop <int>     :maximum number of reads to read into
#                                  memory at once (default: $max_reads_per_loop)
#  --min_pct_read_mapping <int>   :minimum percent of a reads kmers that must map to an
#                                  inchworm bundle (aka. component)  default: 0
#
#  --no_run_chrysalis             :stop Trinity after Inchworm and before
#                                  running Chrysalis
#  --no_run_quantifygraph         :stop Trinity just before running the
#                                  parallel QuantifyGraph computes, to
#                                  leverage a compute farm and massively
#                                  parallel execution..
#
#####################################
###  Butterfly-related options:  ####
#
#  --bfly_opts <string>            :additional parameters to pass through to butterfly
#                                   (see butterfly documentation).
#  --max_number_of_paths_per_node <int>  :only most supported (N) paths are extended from node A->B,
#                                         mitigating combinatoric path explorations. (default: $max_number_of_paths_per_node)
#  --group_pairs_distance <int>    :maximum length expected between fragment pairs (default: $group_pairs_distance)
#                                   
#  --path_reinforcement_distance <int>   :minimum overlap of reads with growing transcript 
#                                        path (default: $path_reinforcement_distance)
#
#  --lenient_path_extension        :require minimal read overlap to allow for path extensions. 
#                                   (equivalent to --path_reinforcement_distance=1)
#
#  --no_triplet_lock               : do not lock triplet-supported nodes
#
#  --bflyHeapSpaceMax <string>     :java max heap space setting for butterfly
#                                   (default: $bflyHeapSpaceMax) => yields command
#                  'java -Xmx$bflyHeapSpaceMax -jar Butterfly.jar ... \$bfly_opts'
#  --bflyHeapSpaceInit <string>    :java initial hap space settings for
#                                   butterfly (default: $bflyHeapSpaceInit) => yields command
#                  'java -Xms$bflyHeapSpaceInit -jar Butterfly.jar ... \$bfly_opts'
#  --bflyGCThreads <int>           :threads for garbage collection
#                                   (default, not specified, so java decides)
#  --bflyCPU <int>                 :CPUs to use (default will be normal 
#                                   number of CPUs; e.g., $CPU)
#  --bflyCalculateCPU              :Calculate CPUs based on 80% of max_memory
#                                   divided by maxbflyHeapSpaceMax
#  --no_run_butterfly              :stops after the Chrysalis stage. You'll
#                                   need to run the Butterfly computes
#                                   separately, such as on a computing grid.
#                  Then, concatenate all the Butterfly assemblies by running:
#                  'find trinity_out_dir/ -name "\*allProbPaths.fasta" \
#                   -exec cat {} + > trinity_out_dir/Trinity.fasta'
#
#  --rerun_butterfly               :will reexecute butterfly commands if composed differently
#                                  from the previous execution.
#
#################################
# Grid-computing options: #######
#
#  --grid_computing_module <string>  : Perl module in $FindBin::RealBin/PerlLibAdaptors/ 
#                                      that implements 'run_on_grid()' 
#                                      for naively parallel cmds. (eg. 'BroadInstGridRunner')
#
#
###############################################################################
#
#  *Note, a typical Trinity command might be:
#        Trinity.pl --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6
#
#     see: $FindBin::RealBin/sample_data/test_Trinity_Assembly/
#          for sample data and 'runMe.sh' for example Trinity execution
#     For more details, visit: http://trinityrnaseq.sf.net
#
###############################################################################




_EOUSAGE_

    ;



=EXPERIMENTAL_OPTIONS

     ## DO NOT USE!  Not ready for prime-time yet, and doesn't seem to help.

#  Chyrsalis-related options:
#
#    --min_pcnt_read_iworm_kmers <int>      :min percentage of a read sequence that must be composed of inchworm kmers to be pursued 
#                                               by chrysalis (default: $min_percent_read_iworm_kmers)  note: off if < 0
#
#    --FORCE_INCHWORM_KMER_METHOD           :uses inchworm built-in kmer cataloger instead of jellyfish (not recommended)  

=cut



my $ROOTDIR = "$FindBin::RealBin";
my $UTILDIR = "$ROOTDIR/util";
my $INCHWORM_DIR = "$ROOTDIR/Inchworm";
my $CHRYSALIS_DIR = "$ROOTDIR/Chrysalis";
my $BUTTERFLY_DIR = "$ROOTDIR/Butterfly";
my $MERYL_DIR = "$ROOTDIR/trinity-plugins/kmer/meryl";
my $JELLYFISH_DIR = "$ROOTDIR/trinity-plugins/jellyfish";
my $FASTOOL_DIR = "$ROOTDIR/trinity-plugins/fastool";
my $COLLECTL_DIR = "$ROOTDIR/trinity-plugins/collectl/bin";
my $COREUTILS_DIR = "$ROOTDIR/trinity-plugins/coreutils/bin";


unless (@ARGV) {
    die "$usage\n";
}

# Log command line parameters for performance monitoring
foreach (@ARGV) {
    $pm_trinity_arguments = $pm_trinity_arguments . " " . $_;
};

my $NO_FASTOOL = 0;
my $NO_CLEANUP = 0;
my $FULL_CLEANUP = 0;

&GetOptions( 
    
    'h|help' => \$help_flag,
    
    ## general opts
    "seqType=s" => \$seqType,
    "left=s" => \$left_file,
    "right=s" => \$right_file,
    "single=s" => \$single_file,
    
    "SS_lib_type=s" => \$SS_lib_type,

    "long_reads=s" => \$long_reads,

    "output=s" => \$output_directory,
    
    "min_contig_length=i" => \$min_contig_length,

    "jaccard_clip" => \$jaccard_clip,
    
    "cite" => \$SHOW_CITATION_FLAG,
    
    'CPU=i' => \$CPU,

    'prep' => \$prep_only,    

    # Butterfly opts
    'no_run_butterfly'      => \$NO_RUN_BUTTERFLY_FLAG,
    'no_triplet_lock'       => \$NO_TRIPLET_LOCK,
    "group_pairs_distance=i" => \$group_pairs_distance,
    'bfly_opts=s'           => \$bfly_opts,
    'bflyHeapSpaceMax=s'    => \$bflyHeapSpaceMax,
    'bflyHeapSpaceInit=s'   => \$bflyHeapSpaceInit,
    'bflyGCThreads=i'       => \$bflyGCThreads,
    'bflyCPU=i'             => \$bflyCPU,
    'bflyCalculateCPU'      => \$bflyCalculateCPU,
    'max_number_of_paths_per_node=i' => \$max_number_of_paths_per_node,
    'lenient_path_extension' => \$lenient_path_extension,
    'path_reinforcement_distance=i' => \$path_reinforcement_distance,
    'rerun_butterfly' => \$RERUN_BUTTERFLY_FLAG,
    
    # Inchworm & kmer catalog opts

    'min_kmer_cov=i'        => \$min_kmer_cov,
    'inchworm_cpu=i'        => \$inchworm_cpu,
    'FORCE_INCHWORM_KMER_METHOD' => \$FORCE_INCHWORM_KMER_METHOD,              
    
    # Jellyfish
    'JM=s'          => \$max_memory, # in GB
    
    # Chrysalis -related opts
    'min_glue=i' => \$min_glue,
    'glue_factor=f' => \$glue_factor,
    'min_iso_ratio=f' => \$min_iso_ratio,
    'min_pcnt_read_iworm_kmers=i' => \$min_percent_read_iworm_kmers, 
    'no_run_quantifygraph' => \$NO_RUN_QUANTIFYGRAPH_FLAG,
    'max_reads_per_graph=i' => \$max_reads_per_graph,
    'max_reads_per_loop=i' => \$max_reads_per_loop,
    'no_run_chrysalis' => \$NO_RUN_CHRYSALIS_FLAG,         
    'min_pct_read_mapping=i' => \$min_pct_read_mapping,
    'weldmer_size=i' => \$weldmer_size,

             
    # Grid computing options
    'grid_computing_module=s' => \$grid_computing_module,         

    "show_advanced_options" => \$show_advanced_options,

             
    # misc
    'run_as_paired' => \$run_as_paired_flag,
    'no_fastool' => \$NO_FASTOOL,
    'no_cleanup' => \$NO_CLEANUP,
    'full_cleanup' => \$FULL_CLEANUP,
    'version' => \$show_version_flag,
    'monitoring' => \$run_with_collectl,

    );



if ($SHOW_CITATION_FLAG) {
    &show_lit_citation();
    exit(0);
}


if ($help_flag) {
    die "$usage\n";
}

if ($show_version_flag) {
    print "Trinity version: $VERSION\n";
    exit(1);
}

if ($NO_CLEANUP && $FULL_CLEANUP) {
    die "cannot set --no_cleanup and --full_cleanup as they contradict";
}


if (@ARGV) {
    die "Error, do not understand options: @ARGV\n";
}






## Check options set:

#     Subroutine takes variable *reference* plus name of variable. Lower-cases
#       variable value and checks to see if it one of the allowed ones.
#       'die' has new-line in order to keep line number from being shown to user.

sub check_option {
    my ($option, $name) = @_;
    $$option = lc $$option;
    if ($$option eq '') {
        die "Error, option '--$name' is required.\n";
    }
    if (!defined $allowed_check{$name}{$$option}) {
        die "Error, option '--$name' ($$option) not one of $allowed{$name}\n";
    }
}

check_option( \$seqType,     'seqType'     );


## load adaptors
my $grid_computing_method = "";
if ($grid_computing_module) {
    my $perl_lib_repo = "$FindBin::Bin/../PerlLibAdaptors";
    print STDERR "-importing module: $grid_computing_module\n";
    require "$grid_computing_module.pm" or die "Error, could not import perl module at run-time: $grid_computing_module";
    
    $grid_computing_method = $grid_computing_module . "::run_on_grid";
    
}

my $USE_FASTOOL = 1; # by default, using fastool for fastq to fasta conversion
if ($NO_FASTOOL) {
    $USE_FASTOOL = 0;
}


if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
        die "Error, unrecognized SS_lib_type value of $SS_lib_type. Should be: F, R, RF, or FR\n";
    }
}

unless ( ($left_file && $right_file) || $single_file ) {
    die "Error, need either options 'left' and 'right' or option 'single'\n";
}

if ($min_iso_ratio > 1) {
    die "Error, --min_iso_ratio should be <= 1 \n";
}

## keep the original 'xG' format string for the --JM option, then calculate the numerical value for max_memory
my $JM_string = $max_memory;    ## this one is used in the Chrysalis exec string
if ($max_memory) {
    $max_memory =~ /^([\d\.]+)G$/ or die "Error, cannot parse max_memory value of $max_memory.  Set it to 'xG' where x is a numerical value\n";
    
    $max_memory = $1;
    $max_memory *= 1024**3; # convert to from gig to bytes
}
else {
    die "Error, must specify max memory for jellyfish to use, eg.  --JM 10G \n";
}


## Try to remove stack limits
if ($^O eq "linux") {  # cannot set stacksize on newer macs for some reason...    
#    &try_unlimit_stacksize();
}

my $curr_limit_settings = `/bin/sh -c 'ulimit -a' `; 
unless ($curr_limit_settings && $curr_limit_settings =~ /\w/) {
    $curr_limit_settings = `/bin/csh -c limit`; # backup, probably not needed.
}

print "Current settings:\n$curr_limit_settings\n\n";


## Check Java version:
my $java_version = `java -version 2>&1 `;
unless ($java_version =~ /(java|openjdk) version \"1\.[7]\./) {
    die "Error, Trinity requires access to Java version 1.6.  Currently installed version is: $java_version";
}

# Give the variable with memory size and a user-oriented name

sub bfly_check {
    my ($mem, $name) = @_;
    my ($num, $type) = $mem =~ /^(\d+)([MG])$/;
    if (!defined $mem || !defined $type) {
        die "Error, $name must be set to a value of format: \\d+G or \\d+M  (eg. 1G or 1000M)\n  Currently: $mem\n";
    }
    return $type eq 'G' ? $num * 1024**3 : $num * 1024**2;
}

my $bflyHeapSpaceMaxBytes  = bfly_check($bflyHeapSpaceMax , 'bflyHeapSpaceMax' );
my $bflyHeapSpaceInitBytes = bfly_check($bflyHeapSpaceInit, 'bflyHeapSpaceInit');

if ($bflyHeapSpaceInitBytes > $bflyHeapSpaceMaxBytes) {
    die "Error, bflyHeapSpaceInit ($bflyHeapSpaceInit) must be less or equal to bflyHeapSpaceMax ($bflyHeapSpaceMax).\n";
}


if ($CPU > $CPU_MAX) {
    print STDERR "Warning, --CPU $CPU might be excessive.  Limiting it to $CPU_MAX for now.\n";
    $CPU = $CPU_MAX;
}

if ($inchworm_cpu > $CPU) {
    $inchworm_cpu = $CPU;
}

if ($bflyCalculateCPU && $max_memory) {
    $bflyCPU = int ($max_memory * 0.80 / $bflyHeapSpaceMaxBytes);
}

$bflyCPU = $CPU if !defined $bflyCPU;

if ($bflyCPU > $CPU_MAX) {
    print STDERR "Warning, --bflyCPU $bflyCPU might be excessive. Limiting it to $CPU_MAX for now.\n";
    $bflyCPU = $CPU_MAX;
}


if (defined($bflyGCThreads) && $bflyGCThreads > 32) {
    die "Error, you probably want fewer than $bflyGCThreads java garbage collection threads. Try a number less than 32.";
}


$ENV{OMP_NUM_THREADS} = $CPU; ## for Inchworm and Chrysalis

my $PAIRED_MODE = ( ($left_file && $right_file)  || $run_as_paired_flag) ? 1:0;
if ($PAIRED_MODE) {
    ## be sure we can find 'bowtie', since we use it as part of the iworm pair scaffolding step
    my $bowtie_path = `which bowtie`;
    if ($bowtie_path =~ /\w/) {
        print "Paired mode requires bowtie. Found bowtie at: $bowtie_path\n";
    }
    else {
        die "Error, cannot find path to bowtie, which is now needed as part of Chrysalis' read scaffolding step";
    }
}

sub perfmon_start {
    open (FILE, ">", "$output_directory/$pm_logfile") or die;
    print FILE "Statistics:\n";
    print FILE "===========\n";
    print FILE     "Trinity Version:      $VERSION\n";
    my $tempp="";
    $tempp=`ldd $INCHWORM_DIR/bin/inchworm | grep "libgomp"`;
    if  ($tempp eq "") {
	print FILE "Compiler:             Intel\n";
    } else {
	print FILE "Compiler:             GCC\n";
    }
    print FILE "Trinity Parameters:  $pm_trinity_arguments\n";
    $pm_trinity_startstring = `date`;
    $pm_trinity_start = `date +%s`;
    close (FILE);
}

sub perfmon_end {
    $pm_trinity_endstring = `date`;
    $pm_trinity_end = `date +%s`;
    my $timestamp = `date +%s`;
    if ( -e "$output_directory/$pm_logfile" ) {
	open (FILE, '>>', "$output_directory/$pm_logfile") or die;
	if ($PAIRED_MODE) {
		print FILE "Paired mode\n";
		print FILE " Input data\n";
	    if ($left_file && $right_file) {
		print FILE "  Left.fasta    $pm_left_fa_size MByte\n";
		print FILE "  Right.fasta   $pm_right_fa_size MByte\n";
	    } else {
		print FILE "  Single.fasta  $pm_single_fa_size MByte\n";
	    }
	} else {
	    print FILE "Unpaired read mode\n";
	    print FILE " Input data\n";
	    print FILE "  Single.fasta  $pm_single_fa_size MByte\n";
	}
    }
    $pm_inchworm_kmers = `cat $output_directory/inchworm.kmer_count`;
    print FILE "  Number of unique KMERs: $pm_inchworm_kmers";
    print FILE "  Number of reads:        $pm_read_count";
    print FILE " Output data\n";
    my $pm_temp = -s "$output_directory/Trinity.fasta" || 0;
    $pm_temp = $pm_temp / 1024 / 1024;
    my $pm_trinity_fa_size = sprintf('%.0f', $pm_temp);
    print FILE "  Trinity.fasta $pm_trinity_fa_size MByte\n\n";
    print FILE "Runtime\n";
    print FILE "=======\n";
    print FILE "Start:       $pm_trinity_startstring";
    print FILE "End:         $pm_trinity_endstring";
    my $pm_trinity_time = $pm_trinity_end - $pm_trinity_start;
    print FILE "Trinity.pl   $pm_trinity_time seconds\n";
    my $pm_inchworm_time = $pm_inchworm_end - $pm_inchworm_start;
    print FILE "  Inchworm   $pm_inchworm_time seconds\n";
    my $pm_chrysalis_time = $pm_chrysalis_end - $pm_chrysalis_start;
    print FILE "  Chrysalis  $pm_chrysalis_time seconds\n";
    my $pm_butterfly_time = $pm_butterfly_end - $pm_butterfly_start;
    print FILE "  Butterfly  $pm_butterfly_time seconds\n";
    my $pm_rest_time = $pm_trinity_time - $pm_butterfly_time - $pm_chrysalis_time - $pm_inchworm_time;
    print FILE "  Rest       $pm_rest_time seconds\n";
    close (FILE);
}

main: {
    $ENV{OMP_NUM_THREADS} = $inchworm_cpu;
    unless ($NO_RUN_BUTTERFLY_FLAG) {
        print STDERR "-since butterfly will eventually be run, lets test for proper execution of java\n";
        &test_java_failure_capture();
    }
    
    
    ## create complete paths for input files:
    $left_file = &create_full_path($left_file) if $left_file;
    $right_file = &create_full_path($right_file) if $right_file;
    $single_file = &create_full_path($single_file) if $single_file;
    $output_directory = &create_full_path($output_directory);
    $long_reads = &create_full_path($long_reads) if $long_reads;
    
    
    unless (-d $output_directory) {
        
        mkdir $output_directory or die "Error, cannot mkdir $output_directory";
    }
    
    chdir ($output_directory) or die "Error, cannot cd to $output_directory";
    
    if ($run_with_collectl){
        $collectl_output_directory = "$start_dir/collectl";
        `rm -rf $collectl_output_directory `;
        $collectl_output_directory = &create_full_path($collectl_output_directory);
        unless (-d $collectl_output_directory) {
            mkdir $collectl_output_directory or die "Error, cannot mkdir $collectl_output_directory";
        }
        my $user = $ENV{USER};
        my $cmd = "cd $collectl_output_directory && ${COLLECTL_DIR}/collectl  -i5:5 -sZ --procfilt U$user -F0 -f $collectl_output_directory/y &";
        &process_cmd($cmd);
        #&process_cmd("sync; sleep 1"); #try to make sure that collectl has started
        
        print STDERR "CHECKING FOR COLLECTL PID\n";

        while (! $collectl_pid) {
            $collectl_pid = `ps h -C collectl -o pid | head -n1`;
            chomp $collectl_pid;
            print STDERR "COLLECTL PID = $collectl_pid\n";
            sleep(1);
        }
    }
    
    &perfmon_start();
    ## create inchworm file name
    my $inchworm_file = "inchworm.K$IWORM_KMER_SIZE.L$MIN_IWORM_LEN";
    unless ($SS_lib_type) {
        $inchworm_file .= ".DS";
    }
    $inchworm_file .= ".fa";
    $inchworm_file = &create_full_path($inchworm_file);
    
    my $trinity_target_fa = ($single_file) ? "single.fa" : "both.fa"; 
    my $inchworm_target_fa = $trinity_target_fa; # change this later if we have long_reads
    

    ## Don't prep the inputs if Inchworm already exists.... Resuming earlier operations.
    my $inchworm_finished_checkpoint_file = "$inchworm_file.finished";
    if (-s $inchworm_file && -e $inchworm_finished_checkpoint_file) {
        print "\n\n#######################################################################\n"
            . "Inchworm file: $inchworm_file detected.\n"
            . "Skipping Inchworm Step, Using Previous Inchworm Assembly\n"
            . "#######################################################################\n\n";
        sleep(2);
    }
    else {
        
        ## Prep data for Inchworm
        
        if ($left_file && $right_file) {

            unless (-s $trinity_target_fa && !-e "left.fa" && !-e "right.fa") {
                
                my ($left_SS_type, $right_SS_type);
                if ($SS_lib_type) {
                    ($left_SS_type, $right_SS_type) = split(//, $SS_lib_type);
                }
                print("Converting input files. (in parallel)");
                my $thr1;
                my $thr2;
                if (!(-s "left.fa")) {
                    $thr1 = threads->create('prep_seqs', $left_file, $seqType, "left", $left_SS_type);
                } else {
                    $thr1 = threads->create(sub { print ("left file exists, nothing to do");});
                }
                if (!(-s "right.fa")) {
                    $thr2 = threads->create('prep_seqs', $right_file, $seqType, "right", $right_SS_type);
                } else {
                    $thr2 = threads->create(sub { print ("right file exists, nothing to do");});
                }
                $thr1->join();
                $thr2->join();
                
                print("Done converting input files.");
                ## Calculate input file sizes for performance monitoring
                my $pm_temp = -s "$left_file";
                $pm_temp = $pm_temp / 1024 / 1024;
                $pm_left_fa_size = sprintf('%.0f', $pm_temp);
                $pm_temp = -s "$right_file";
                $pm_temp = $pm_temp / 1024 / 1024;
                $pm_right_fa_size = sprintf('%.0f', $pm_temp);
                
                &process_cmd("cat left.fa right.fa > $trinity_target_fa") unless (-s $trinity_target_fa && (-s $trinity_target_fa == ((-s "left.fa") + (-s "right.fa"))));
                unless (-s $trinity_target_fa == ((-s "left.fa") + (-s "right.fa"))){
                    die "$trinity_target_fa is smaller (".(-s $trinity_target_fa)." bytes) than the combined size of left.fa and right.fa (".((-s "left.fa") + (-s "right.fa"))." bytes)\n";
                }
                
                unlink ("left.fa", "right.fa"); # no longer needed now that we have 'both.fa', which is needed by chryaslis
            }
        }
        elsif ($single_file) {
            
            &prep_seqs($single_file, $seqType, "single", $SS_lib_type) unless (-s "single.fa");
            ## Calculate input file sizes for performance monitoring
            my $pm_temp = -s "$single_file";
            $pm_temp = $pm_temp / 1024 / 1024;
            $pm_single_fa_size = sprintf('%.0f', $pm_temp);
        }
        
        else {
            die "not sure what to do. "; # should never get here.
        }
    
        if ($long_reads) {
            $inchworm_target_fa .= ".wLongReads.fa";
            &process_cmd("cat $long_reads $trinity_target_fa > $inchworm_target_fa");
        }
            
    }
    
    if ($prep_only){
	print "Data has been prepared. Exiting now as per user request\n";
	exit();
    }
    
    #################
    ## Inchworm step:
    $pm_inchworm_start = `date +%s`;
    unless (-s $inchworm_file && -e $inchworm_finished_checkpoint_file) {
                    

        &run_inchworm($inchworm_file, $inchworm_target_fa, $SS_lib_type, $kmer_method);
        &process_cmd("touch $inchworm_finished_checkpoint_file");
    }
    $pm_inchworm_end = `date +%s`;

    
    unless (-s $inchworm_file) {

        if ($FULL_CLEANUP && -e $inchworm_file && -e $inchworm_finished_checkpoint_file) {
            ## GG-trinity mode, clean-up gracefully
            &process_cmd("rm -rf $output_directory");
            exit(0);
        }
        else {
            die "Error, no Inchworm output is detected at: $inchworm_file";
        }
    }

    
    if ($jaccard_clip) {

        eval {
            
            if ($jaccard_clip && $left_file && $right_file) {
                $inchworm_file = &run_jaccard_clip_left_right($inchworm_file, $left_file, $right_file, $seqType, $SS_lib_type);
            }
            elsif ($jaccard_clip && $single_file) {
                $inchworm_file = &run_jaccard_clip_single_but_really_paired($inchworm_file, $single_file, $seqType, $SS_lib_type);
            }
        };

        if ($@) {
            if ($FULL_CLEANUP) {
                ## GG-trinity mode, clean up gracefully
                &process_cmd("rm -rf $output_directory");
                exit(0);
            }
            else {
                die "Error, jaccard-clip failed: $@";
            }
        }
    }
    
    
    if ($NO_RUN_CHRYSALIS_FLAG) {
        print "\n\n\n";
        print "#########################################################################\n";
        print "Inchworm is complete.  --no_run_chrysalis was specified, so stopping here.\n";
        print "#########################################################################\n\n\n";
    
        exit(0);
    }
    $ENV{OMP_NUM_THREADS} = $CPU;
    ##################
    ## Chrysalis step:
    
    if ($min_percent_read_iworm_kmers > 0) {
        
        ###  EXPERIMENTAL:  DO NOT USE!
        
        $trinity_target_fa = &extract_reads_with_iworm_kmers($trinity_target_fa, $inchworm_file, $min_percent_read_iworm_kmers, $SS_lib_type);
        
    }
    
    ## butterfly commands can be reparameterized for exploring different assembly requirements
    ## chrysalis will just run or resume depending on what's already been processed.
    $pm_chrysalis_start = `date +%s`;
    my $butterfly_cmds = &run_chrysalis($inchworm_file, $inchworm_target_fa,
                                        $min_contig_length, $group_pairs_distance, $SS_lib_type, $trinity_target_fa);
    $pm_chrysalis_end = `date +%s`;

    print "Butterfly_cmds: $butterfly_cmds\n";
    
    if ($butterfly_cmds && -s $butterfly_cmds) {

        if ($NO_RUN_BUTTERFLY_FLAG) {
            
            print "\n\nYou've opted to run butterfly commands independently from this script, such as on a computing grid.\n\n";
            print "Butterfly commands to execute are available here:\n"
                . "\t$butterfly_cmds\n\n";
            print "After executing Butterfly commands, concatenate all Butterfly outputs by running:\n"
                . "\t\tfind $output_directory/ -name \"\*allProbPaths.fasta\" -exec cat {} + > $output_directory/Trinity.fasta\n\n\n";
            
        }
        else {
            
            ## Run Butterfly
            
            print "Inchworm and Chrysalis complete.  Butterfly commands to execute are provided here:\n"
                . $butterfly_cmds . "\n\n";
            
            
            print STDERR "---------------------------------------------------------------\n"
                . "-------------------- Butterfly --------------------------------\n"
                . "-- (Reconstruct transcripts from reads and de Bruijn graphs) --\n"
                . "---------------------------------------------------------------\n\n";
            
            $pm_butterfly_start = `date +%s`;
            if ($grid_computing_module) {
                my @bfly_cmds = `cat $butterfly_cmds`;
                chomp @bfly_cmds;
                
                &$grid_computing_method(@bfly_cmds);
                
            }
            else {
                my $cmd = "$INCHWORM_DIR/bin/ParaFly -c $butterfly_cmds -shuffle -CPU $bflyCPU -failed_cmds failed_butterfly_commands.$$.txt -v ";  # shuffle them since the first ones are usually the longest-running ones.
                &process_cmd($cmd);
            }
            $pm_butterfly_end = `date +%s`;

            ## capture results:
            # my $cmd = 'find ./chrysalis -name "*allProbPaths.fasta" -exec cat {} + > Trinity.fasta.tmp';
            # no longer scan the file system... we know which files should exist
            my $cmd = "$UTILDIR/print_butterfly_assemblies.pl ./chrysalis/component_file_listing.txt > Trinity.fasta.tmp";
            &process_cmd($cmd);
            
        }

    }
     
    if ($FULL_CLEANUP) {
        print "Fully cleaning up.\n";
        $output_directory =~ s|/+$||g; # remove any trailing directory slash
    
        if (-s "Trinity.fasta.tmp") {
            rename("Trinity.fasta.tmp", "$output_directory.Trinity.fasta") or die "Error, cannot rename Trinity.fasta.tmp to $output_directory.Trinity.fasta";
        }
        &process_cmd("rm -rf $output_directory");
        
        
        print "\n\n";
        print "###################################################################\n";
        print "Butterfly assemblies are written to $output_directory.Trinity.fasta\n";
        print "###################################################################\n\n\n";
        
        
    }
    else {
        
        if (-s "Trinity.fasta.tmp") {
            rename("Trinity.fasta.tmp", "Trinity.fasta") or die "Error, cannot rename Trinity.fasta.tmp to Trinity.fasta"; # now that process has finished.
        }
        
        if (-s "Trinity.fasta") {
        
            print "\n\n";
            print "###################################################################\n";
            print "Butterfly assemblies are written to $output_directory/Trinity.fasta\n";
            print "###################################################################\n\n\n";
        }
        else {
            die "ERROR, no butterfly assemblies written to $output_directory/Trinity.fasta";
        }
        
    }
    
    &perfmon_end();
    # finish monitoring and create collectl statistics
    if ($run_with_collectl){
        &process_cmd("kill -15 $collectl_pid");
	#&process_cmd("sync");
	#&process_cmd("sleep 3");
        chdir ($collectl_output_directory) or die "Error, cannot cd to $collectl_output_directory";
        &process_cmd("$COLLECTL_DIR/make_data_files.sh");
        &process_cmd("$COLLECTL_DIR/timetable.sh");
        $collectl_titlename = "${VERSION} ${CPU} ${left_file}${single_file}";
        &process_cmd("$COLLECTL_DIR/plot.sh \"$collectl_titlename\" ${CPU}");
    }
    exit(0);
}


####
sub run_chrysalis {
    my ($inchworm_file, $reads_file,
        $min_contig_length, $group_pairs_distance, $SS_lib_type, $pairs_fa) = @_;
    
    
    my $butterfly_cmds = &create_full_path("chrysalis/butterfly_commands");
    
    my $quantify_graph_cmds = &create_full_path("chrysalis/quantifyGraph_commands");
    
    my $adjusted_butterfly_cmds = "$butterfly_cmds.adj";
    
    my $chrysalis_finished_checkpoint = "chrysalis/chrysalis.finished";
    
    if (-s $butterfly_cmds && -e $chrysalis_finished_checkpoint) {
        
        print "###################################################################\n";
        print "#### Chrysalis results already exist. Not rerunning Chrysalis. ####\n";
        print "###################################################################\n\n\n"; 
        
        sleep(2);
    
    }
    else {
        ## run Chrysalis
        
        my $cmd = "$CHRYSALIS_DIR/Chrysalis -i $reads_file -iworm $inchworm_file -o chrysalis -cpu $CPU "
            . " -min_glue $min_glue -min_iso_ratio $min_iso_ratio -glue_factor $glue_factor -weldmer_size $weldmer_size "
            . " -min $min_contig_length -dist $group_pairs_distance -max_reads $max_reads_per_graph "
            . " -sort_buffer_size $JM_string -max_mem_reads $max_reads_per_loop ";
        
        if ($SS_lib_type) {
            $cmd .= " -strand 1 ";
        }
        
        if ($PAIRED_MODE) {
            $cmd .= " -paired ";
            $cmd .= " -reads_for_pairs $pairs_fa ";
        }
        
    
        if ($min_pct_read_mapping) {
            $cmd .= " -min_pct_read_mapping $min_pct_read_mapping ";
        }
        

        $cmd .= " -butterfly $BUTTERFLY_DIR/Butterfly.jar ";
        
        if ($NO_CLEANUP) {
            $cmd .= " -no_cleanup ";
        }
        
        $cmd .= " 2>&1 ";
        
        eval {
            
            &process_cmd($cmd);
            
        };
        
        
        if ($@) {
            
            my $errmsg = "$curr_limit_settings\n";
            $errmsg .= "Error, the Chrysalis process failed:\n$@\n";
            croak $errmsg;
        }
        
           
        print "Chrysalis initial stage completed successfully.\n";
        &process_cmd("touch $chrysalis_finished_checkpoint");
    }
    
    unless (-s $butterfly_cmds) {
        
        if ($FULL_CLEANUP) {
            ## Trinity GG mode
            return("");
        }
        
        
        croak "\n\n*************\nError, chrysalis did not report butterfly commands file: $butterfly_cmds\n"
            . "Could it be that your read data is too sparse and no components were capable of generating minimal length contigs?\n"
            . "*****************\n\n";
    }
    
        
    ## Rewrite the Butterfly commands

    ## add additional butterfly opts.
    open (my $fh, $butterfly_cmds) or die "Error, cannot read file $butterfly_cmds";
    open (my $ofh, ">$adjusted_butterfly_cmds") or die "Error, cannot write to $adjusted_butterfly_cmds";
    my @bfly_cmds;
    while (<$fh>) {
        my $line = $_;
        chomp $line;
        $line =~ s/^java /java -Xmx$bflyHeapSpaceMax -Xms$bflyHeapSpaceInit / or die "Error, cannot modify command";
        
        if (defined($bflyGCThreads)) {
            $line =~ s/^java /java -XX:ParallelGCThreads=$bflyGCThreads / or die "Error, cannot modify command";
        }
        
        if ($bfly_opts) {
            $line .= " $bfly_opts ";
        }
        $line .= " --max_number_of_paths_per_node=$max_number_of_paths_per_node ";
        if ($lenient_path_extension) {
            $line .= " --lenient_path_extension ";
        }
        else {
            $line .= " --path_reinforcement_distance=$path_reinforcement_distance ";
            $line .= " --triplet-lock " unless ($NO_TRIPLET_LOCK);
        }
        
        ## check to see if this butterfly job already finished successfully (don't rerun it)
        my @line_pts = split(/\s+/, $line);
        my ($comp_entry) = grep { /comp\d+$/ } @line_pts;
        unless ($comp_entry) {
            die "Error, couldn't decipher the component entry from $line";
        }
        my $bfly_finished_file = "$comp_entry.bfly.finished";
        if (-e $bfly_finished_file && ! $RERUN_BUTTERFLY_FLAG) {
            print STDERR "Butterfly already processed for component: $line, not executing again.\n";
        }
        else {
            push (@bfly_cmds, $line);
        }
    }
    
        
    foreach my $bfly_cmd (@bfly_cmds) {
        print $ofh $bfly_cmd . "\n";
        #print STDERR $line . "\n";
    }
    
    close $ofh;
    close $fh;
    

     
    # see if we need to run the quantifyGraph commands:
    if ($NO_RUN_QUANTIFYGRAPH_FLAG) {

        print "#############################################################################\n";
        print "## Ceasing Trinity prior to execution of massively parallel operations.\n";
        print "##\n";
        print "## To complete Trinity, execute the following sets of commands:\n";
        print "##\n";
        print "## First, run the Chrysalis QuantifyGraph commands in parallel:\n";
        print "##    $quantify_graph_cmds\n";
        print "##\n";
        print "## Then, execute all the Butterfly commands:\n";
        print "##    $adjusted_butterfly_cmds\n";
        print "##\n";
        print "## And, finally, concatenate all Butterfly assemblies into a single output file:\n";
        print "##\n";
        print "##     find $output_directory/ -name \"\*allProbPaths.fasta\" -exec cat {} + > $output_directory/Trinity.fasta\n";
        print "##\n";
        print "##############################################################################\n";
        print "\n\n";
        
        exit(0);
    }
    else {

        
        my $quantify_graph_cmds_finished = &create_full_path("chrysalis/quantifyGraph_commands.run.finished");
        if (! -e $quantify_graph_cmds_finished) {
            ## run it
            
            print STDERR "---------------------------------------------------\n"
                       . "----------- Chrysalis: QuantifyGraph --------------\n"
                       . "-- (Integrate mapped reads into de Bruijn graph) --\n"
                       . "---------------------------------------------------\n\n";
            
            
            if ($grid_computing_module) {
                my @quantify_graph_cmds = `cat $quantify_graph_cmds`;
                chomp @quantify_graph_cmds;
                
                

                &$grid_computing_method(@quantify_graph_cmds);
            }
            else {

                my $cmd = "$INCHWORM_DIR/bin/ParaFly -c $quantify_graph_cmds -CPU $CPU -failed_cmds failed_quantify_graph_commands.$$.txt -v -shuffle ";
                &process_cmd($cmd);
            }
            
            # write checkpoint
            &process_cmd("touch $quantify_graph_cmds_finished");
        }
                
        return($adjusted_butterfly_cmds);
    
    }
}


####
sub run_inchworm {
    my ($inchworm_outfile, $reads, $strand_specific_flag, $kmer_method) = @_;
    
    
    ## get count of number of reads to be assembled.
    my $read_count_file = "$reads.read_count";
    if (! -s $read_count_file) {
        my $count_of_reads = `grep '>' $reads | wc -l`;
        $pm_read_count = $count_of_reads;
        open (my $ofh, ">$read_count_file") or die $!;
        print $ofh $count_of_reads;
        close $ofh;
    }
    

    my $inchworm_cmd;
    
    my @tmp_files; # to be deleted after successful inchworm run.

    
    #####################################################
    ## Using Jellyfish kmer method
    #####################################################

    if (! $FORCE_INCHWORM_KMER_METHOD) {

        my $jelly_kmer_fa_file = "jellyfish.kmers.fa";
        my $jelly_finished_checkpoint_file = "jellyfish.$min_kmer_cov.finished";
        unless (-e $jelly_finished_checkpoint_file) {
            

            print STDERR "-------------------------------------------\n"
                       . "----------- Jellyfish  --------------------\n"
                       . "-- (building a k-mer catalog from reads) --\n"
                       . "-------------------------------------------\n\n";

            
            my $read_file_size = -s $reads;
            
            my $jelly_hash_size = int( ($max_memory - $read_file_size)/7); # decided upon by Rick Westerman
            
            
            if ($jelly_hash_size < 100e6) {
                $jelly_hash_size = 100e6; # seems reasonable for a min hash size as 100M
            }
            
            my $cmd = "$JELLYFISH_DIR/bin/jellyfish count -t $CPU -m $IWORM_KMER_SIZE -s $jelly_hash_size ";
            
            unless ($SS_lib_type) {
                ## count both strands
                $cmd .= " --both-strands ";
            }
            
            $cmd .= " $reads";
            
            &process_cmd($cmd);
            
            my @kmer_db_files;
            
            if (-s $jelly_kmer_fa_file) {
                unlink($jelly_kmer_fa_file) or die "Error, cannot unlink $jelly_kmer_fa_file";
            }
            
            foreach my $file (<mer_counts_*>) {
                my $cmd = "$JELLYFISH_DIR/bin/jellyfish dump -L $min_kmer_cov $file >> $jelly_kmer_fa_file";

                &process_cmd($cmd);
                push (@tmp_files, $file); # don't retain the individual jelly kmer files.
            }
            
            ## if got this far, consider jellyfish done.
            &process_cmd("touch $jelly_finished_checkpoint_file");
        }

        
        $inchworm_cmd = "$INCHWORM_DIR/bin/inchworm --kmers $jelly_kmer_fa_file --run_inchworm -K $IWORM_KMER_SIZE -L $MIN_IWORM_LEN --monitor 1 ";

        # hold on to the jellyfish file - we might use it for other applications.
        #push (@tmp_files, $jelly_finished_checkpoint_file, $jelly_kmer_fa_file) unless $NO_CLEANUP;
        
    }
    else {
        
        ######################################################
        ## Using Inchworm kmer method (original, slow method)
        ######################################################
                
        $inchworm_cmd = "$INCHWORM_DIR/bin/inchworm --reads $reads --run_inchworm -K $IWORM_KMER_SIZE -L $MIN_IWORM_LEN --monitor 1 ";
        if ($min_kmer_cov > 1) {
            $inchworm_cmd .= " --minKmerCount $min_kmer_cov ";
        }
    }
    

    ## finish constructing the inchworm command to execute
    
    unless ($strand_specific_flag) {
        $inchworm_cmd .= " --DS ";
    }
    
    
    #$inchworm_cmd .= " 2>inchworm.log > $inchworm_outfile.tmp";
    $inchworm_cmd .= " > $inchworm_outfile.tmp"; 
    
    print STDERR "----------------------------------------------\n"
               . "--------------- Inchworm ---------------------\n"
               . "-- (Linear contig construction from k-mers) --\n"
               . "----------------------------------------------\n\n";

    
    eval {
                
        &process_cmd($inchworm_cmd);;
    };
    
    if ($@) {
        
        print STDERR "$@\n";
        print "** The inchworm process failed.";
        print STDERR "\n\nIf it indicates bad_alloc(), then Inchworm ran out of memory.  You'll need to either reduce the size of your data set or run Trinity on a server with more memory available.\n\n";
        exit(1);
    }
    
    rename("$inchworm_outfile.tmp", $inchworm_outfile) or die "Error, cannot rename $inchworm_outfile.tmp to $inchworm_outfile"; # now we know for sure it's done.
    

    ## remove the tmp files
    foreach my $tmp_file (@tmp_files) {
        unlink($tmp_file);
    }
        
    return;
    
}

####
sub prep_seqs {
    my ($initial_file, $seqType, $file_prefix, $SS_lib_type) = @_;

    if ($seqType eq "fq") {
        # make fasta
        
        my $perlcmd = "$UTILDIR/fastQ_to_fastA.pl -I $initial_file ";
        my $fastool_cmd = "$FASTOOL_DIR/fastool";
        if ($SS_lib_type && $SS_lib_type eq "R") {
            $perlcmd .= " --rev ";
            $fastool_cmd .= " --rev ";
        }
        $fastool_cmd .= " --illumina-trinity --to-fasta $initial_file > $file_prefix.fa";
        $perlcmd .= " > $file_prefix.fa";  
        
        my $cmd = ($USE_FASTOOL) ? $fastool_cmd : $perlcmd;
        
        &process_cmd($cmd) unless (-e "$file_prefix.fa");
    }
    elsif ($seqType eq "fa") {
        if ($SS_lib_type && $SS_lib_type eq "R") {
            my $cmd = "$UTILDIR/revcomp_fasta.pl $initial_file > $file_prefix.fa";
            &process_cmd($cmd) unless (-e "$file_prefix.fa");
        }
        else {
            ## just symlink it here:
            my $cmd = "ln -s $initial_file $file_prefix.fa";
            &process_cmd($cmd) unless (-e "$file_prefix.fa");
        }
    }
    elsif (($seqType eq "cfa") | ($seqType eq "cfq")) {
        # make double-encoded fasta
        my $cmd = "$UTILDIR/csfastX_to_defastA.pl -I $initial_file ";
        if ($SS_lib_type && $SS_lib_type eq "R") {
            $cmd .= " --rev ";
        }
        $cmd .= "> $file_prefix.fa";
        &process_cmd($cmd) unless (-e "$file_prefix.fa");
  }
    return;
}



###
sub create_full_path {
    my ($file) = @_;

    my $cwd = cwd();
    if ($file !~ m|^/|) { # must be a relative path
        $file = $cwd . "/$file";
    }
    
    return($file);
}



####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system($cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

    return;
}


####
sub run_jaccard_clip_left_right {
    my ($inchworm_file, $left_file, $right_file, $seqType, $SS_lib_type) = @_;

    my $output_file = "$inchworm_file.clipped.fa";

    if (-s $output_file) {
        print STDERR "###### WARNING: $output_file already exists, skipping the jaccard-clip step, using already existing output: $output_file\n";
        return($output_file);
    }
    
    my $cmd = "$UTILDIR/inchworm_transcript_splitter.pl --iworm $inchworm_file "
        . " --left $left_file --right $right_file --seqType $seqType --CPU $CPU ";

    if ($SS_lib_type) {
        $cmd .= " --SS_lib_type $SS_lib_type ";
    }
    
    &process_cmd($cmd);



    unless (-s $output_file) {
        croak "Error, jaccard clipping didn't produce the expected output file: $output_file";
    }

    return($output_file);
}



####
sub run_jaccard_clip_single_but_really_paired {
    my ($inchworm_file, $single_file, $seqType, $SS_lib_type) = @_;

    my $output_file = "$inchworm_file.clipped.fa";

    if (-s $output_file) {
        print STDERR "###### WARNING: $output_file already exists, skipping the jaccard-clip step, using already existing output: $output_file\n";
        return($output_file);
    }
    
    my $cmd = "$UTILDIR/inchworm_transcript_splitter.pl --iworm $inchworm_file "
        . " --single_but_really_paired $single_file --seqType $seqType --CPU $CPU ";
    
    if ($SS_lib_type) {
        $cmd .= " --SS_lib_type $SS_lib_type ";
    }
    
    &process_cmd($cmd);



    unless (-s $output_file) {
        croak "Error, jaccard clipping didn't produce the expected output file: $output_file";
    }

    return($output_file);
}

####
sub test_java_failure_capture {
    
    print "#######################################\n";
    print "Running Java Tests\n";
    
    my $java_prog = `which java`;
    unless ($java_prog) {
        die "Error, cannot find 'java'.  Please be sure it is available within your \${PATH} setting and then try again.";
    }
    

    my $cmd = "java -jar $UTILDIR/ExitTester.jar 0";
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "Error encountered in testing for running of a simple java application. ";
        print "$@\n\n";
        print STDERR "Please check your java configuration.\n";
        exit(1);
        
    }
    
    $cmd = "java -jar $UTILDIR/ExitTester.jar 1";
    eval {
        &process_cmd($cmd);
    };

    if ($@) {
        print "-we properly captured the java failure status, as needed.  Looking good.\n";
    }
    else {
        print STDERR "-we are unable to properly capture java failure status.  Please be sure that java (or any wrapper around java that's being used) can properly capture and propagate failure status before proceeding.\n";
        exit(1);
    }

    print "Java tests succeeded.\n";
    print "###################################\n\n";
    
    return;
}


####
sub extract_reads_with_iworm_kmers {
    my ($trinity_target_fa, $inchworm_file, $min_percent_read_containing_kmers, $SS_lib_type) = @_;

    my $extracted_reads_file = "$trinity_target_fa." . $min_percent_read_containing_kmers . "pcnt.iworm_extracted";

    my $cmd = "$INCHWORM_DIR/bin/pull_reads_with_kmers "
        . "--target $inchworm_file "
        . "--reads $trinity_target_fa "
        . "--min_percent_read_containing_kmers $min_percent_read_containing_kmers ";
    
    unless ($SS_lib_type) {
        $cmd .= " --DS ";
    }
    
    $cmd .= " > $extracted_reads_file ";
    
    if (-s $extracted_reads_file) {
        print STDERR "-warning, iworm kmer-extracted reads file already exists: $extracted_reads_file.  Re-using it.\n";
    }
    else {

        &process_cmd($cmd);
    }

    return($extracted_reads_file);
}


sub try_unlimit_stacksize {

    # from Ryan Thompson
    eval "use BSD::Resource; setrlimit(RLIMIT_STACK, RLIM_INFINITY, RLIM_INFINITY); ";
    
    if( $@ ) {
        warn <<"EOF";
        
            $@

            Unable to set unlimited stack size. Please install the BSD::Resource
            Perl module to allow this script to set the stack size, or set it
            yourself in your shell before running Trinity (ignore this warning if
            you have set the stack limit in your shell). See the following URL for
            more information:

            http://trinityrnaseq.sourceforge.net/trinity_faq.html#ques_E

EOF
;
    }
    else {
        print "Successfully set unlimited stack size.\n";
        print "###################################\n\n";
    }
    return;
}


####
sub show_lit_citation {
    
    print "\n\n";
    print "###########################################################################################\n\n";

    print "------------------------------------------------------------------------------------------\n";
    print "----- Trinity-based Transcript Reconstruction --------------------------------------------\n";
    print "------------------------------------------------------------------------------------------\n\n"; 
    
    
    print "* Trinity:\n"
        . "Full-length transcriptome assembly from RNA-Seq data without a reference genome.\n"
        . "Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,\n"
        . "Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F,\n"
        . "Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.\n"
        . "Nature Biotechnology 29, 644–652 (2011)\n"
        . "Paper: http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html\n"
        . "Code:  http://trinityrnaseq.sf.net\n\n";
    
    print "------------------------------------------------------------------------------------------\n";
    print "----- Tools Below are Used Within Trinity Accordingly ------------------------------------\n";
    print "------------------------------------------------------------------------------------------\n\n"; 
    
    print "* Fastool (for fast fastQ-to-fastA conversion):\n"
        . "Francesco Strozzi\n"
        . "Code: https://github.com/fstrozzi/Fastool\n\n";
    
    print "* Meryl (for fast K-mer counting):\n"
        . "Brian Walenz\n"
        . "Code: http://kmer.sourceforge.net/meryl\n\n";
    
    print "* Jellyfish (for fast K-mer counting):\n"
        . "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers.\n"
        . "Guillaume Marcais and Carl Kingsford.\n"
        . "Bioinformatics (2011) 27(6): 764-770\n"
        . "Paper: http://bioinformatics.oxfordjournals.org/content/27/6/764.long\n"
        . "Code: http://www.cbcb.umd.edu/software/jellyfish\n\n";
    



    print "* RSEM (for abundance estimation):\n"
        . "RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome.\n"
        . "Bo Li and Colin N Dewey.\n"
        . "BMC Bioinformatics. 2011; 12: 323.\n"
        . "Paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3163565\n"
        . "Code: http://deweylab.biostat.wisc.edu/rsem\n\n";
    

    print "###########################################################################################\n\n";
    print "\n\n";

    
    
    return;
}

