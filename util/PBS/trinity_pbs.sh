#!/bin/bash 
set -e  # turn on exit on error
##################################################################################################################################
##########################                                                                ########################################
##########################     Trinity PBS job submission with multi part dependencies    ########################################
##########################                                                                ########################################
##################################################################################################################################
### Author: Josh Bowden, CSIRO
### Version 0.1
### Based on script by Alexie Papanicolaou for submission to PBS pro.
###
###		
###		Script to split the Trinity workflow into multiple stages so as to efficiently request 
###		and use appropriate resources (walltime and number of cores) on a computer cluster / supercomputer.
###		Currently creates scripts for PBS Torque or PBSpro
###	    
###		trinity_pbs script install instructions:
###			1. Copy all trinity_pbs.* files into a directory (we will call it "TRINITY_PBS_DIR").
###			2. Add TRINITY_PBS_DIR to the PATH i.e. export or set PATH=TRINITY_PBS_DIR:$PATH  (perhaps export PATH in .bashrc file)
###			3. Change the "RUNDIR" variable found below to point to the directory also. i.e. RUNDIR=TRINITY_PBS_DIR
###			4. Set MEMDIRIN to name of a node-local filesystem so a network drive is not needed unecesarily for Scripts 4b and 5b
###			5. Set MODTRINITY to the path to trinity executables and/or modules to be loaded so Trinity.pl can be run.
###			6. Set PBSTYPE to --pbspro or --pbs, dependent on the system present.
###		That should be all that is needed from an admin perspective (besides making scripts accessible and exectable for users)			
###
###		Users need make a copy of TRINITY.CONFIG and then modify variables in it. See TRINITY.CONFIG for further details.
###		
###		The current script does the following.
###		Part 1. Reads data from TRINITY.CONF and creates the input directory, data file names, output data directory and Trinity.pl command line 
###			User inputs from config file :
###			JOBPREFIX				A string of less than 11 characters long. PBS will use this as a jobname prefix.
###			DATADIRECTORY  			Where input data exists
###			OUTPUTDIR  				Where user wants output data to go - requires a lot of space even for small datatsets 
###			STANDARD_JOB_DETAILS 	the Trinity.pl command line
###			ACCOUNT  				Account details of user (if required by PBS system being used)
###			STANDARD_JOB_DETAILS	Other Trinity runtime input options
###	   		
###		Part 2. Writes scripts to run Trinity.pl in 6 stages: 
###			3 intial (Inchworm, and 2 x Chrysalis stages: Chrysalis::GraphFromFasta and Chrysalis::ReadsFromTranscripts) 
###			2 parallel stages (Chrysalis::QuantifyGraph and Butterfly) which use PBS array job for execution across a cluster.
###			1 collection of results as Trinity.Fasta.
###	      			
###			Information input from command line filename for stage 'x' : 
###			WALLTIME_Px  			Amount of time stage requires
###			MEM_Px					The amount of memory the stage requires
###			NCPU_Px     			The number of CPUs the stage may use
###			PBSQUEUE _Px				The PBS (for --pbspro only) queue name
###			NUMPERARRAYITEM_Px		The number of massively parallel jobs per "pbs array" job
###	          		                                     
###		Part 3. Runs scripts dependant upon what stage has been detected as completed, using PBS job dependencies  
###                                                                                                      
###			Command line usage:
###				To start (or re-start) an analysis:
###					>trinity_pbs.sh TRINITY.CONFIG
###				To stop previously started PBS jobs on the queue:
###					>trinity_pbs.sh --rm OUTPUTDIR 
###				Where:
###					TRINITY.CONFIG = user specific job details
###					OUTPUTDIR   = is path to output data directory
### 
### Output job script submission files. These are saved in the output directory (OUTPUTDIR) and can be modified/re-run if any job fails.
### *_run.sh       Runs all the following scripts - with job dependencies and only the jobs that still need to be run.
### *_p1.sh        Runs Inchworm stage. Does not scale well past a single socket. Only request at most the number of cores on a single CPU.
### *_p2.sh        Runs Chrysalis::GrapghFromFasta clustering of Inchworm output. Should scale to number of cores on node
### *_p3.sh        Runs Chrysalis::ReadsToTranscripts. I/O limited. Try to use local filesystem (not implemented)
### *_p4a.sh       Creates a PBS array job to run Chrysalis::QuantifyGraph parallel tasks
### *_p4b.sh       QuantifyGraph array job. "NUMPERARRAYITEM" tasks from the file /chrysalis/quantifyGraph_commands are run for each job in the array
### *_p5a.sh       Creates a PBS array job to run Butterfly parallel tasks
### *_p5b.sh       Butterfly array job. "NUMPERARRAYITEM" tasks from the file /chrysalis/butterfly_commands.adj are run for each job in the array
### *_p6           Collects together into Trinity.fasta file. If this completes then full computation is finished. Otherwise re-run *_run.sh 
###                to start off at last completed stage. At present leaves all data on temporary area of shared network drive and 
###                copies Trinity.fatsa to home directory (with specific job prefix in filename). 
### * = $JOBPREFIX. $JOBPREFIX should not be > 10 characters long

##################################################################################################################################
########    Set RUNDIR to the directory where trinity_pbs.sh and ancillary scripts are installed
##################################################################################################################################
RUNDIR="/home2/bow355/bin/trinity_pbs"
#RUNDIR="/home/bow355/ACP_butterfly/trinity_pbs"

##################################################################################################################################
########  Set cluster specific name for compute node local filesystem
##################################################################################################################################
MEMDIRIN="\$TMPDIR"  # available on Barrine  
# MEMDIRIN="\$MEMDIR"  # available on Bragg-l

##################################################################################################################################
########     Select the PBS type availble on your system "--pbspro" or "--pbs" (for PBS Torque)
##################################################################################################################################
## PBS Pro
PBSTYPE="--pbspro"
## PBS torque:
##PBSTYPE="--pbs"

##################################################################################################################################
########     Set system specific PATHS and load system specific modules (if available) 
########     Please note, (as of 15/10/2012) Trinity requires java/1.6 so ensure no other version is loaded/available
##################################################################################################################################
MODTRINITY="module load mpt/2.00 perl/5.15.8 bowtie/12.7 jellyfish/1.1.5 samtools/1.18 trinity/2012-06-08; module unload java/1.7.0_02; module load java/1.6.0_22-sun;

## Another example that sets PATH directly (for systems without 'modules')
#MODTRINITY=" trinity/2012-06-08; export PATH=/home/bow355/trinityrnaseq_r2012-06-08:/home/bow355/trinityrnaseq_r2012-06-08/Inchworm/bin:/home/bow355/trinityrnaseq_r2012-06-08/Chrysalis:/home/bow355/trinityrnaseq_r2012-06-08/trinity-plugins/bwa-0.5.7-patched_multi_map:/home/bow355/trinityrnaseq_r2012-06-08/trinity-plugins/fastool:\$PATH ;"

## That should be all the admin modifications needed (unless the PBS system has a different syntax).
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################


##################################################################################################################################
#########    Load files that contains functions 
##################################################################################################################################
# Modify function F_GETNODESTRING in file trinity_pbs.header so that a correct PBS header is returned to suit your PBS cluster
if [ -e "$RUNDIR"/trinity_pbs.header ] ; then
  	source "$RUNDIR"/trinity_pbs.header 
else
	echo "$1 requires file \"trinity_pbs.header\" to be present in: "
	echo "$RUNDIR"
	exit 1
fi

##################################################################################################################################
#########    Load input config file
##################################################################################################################################
if [ -e "$1" ] ; then
  	source "$1"
else
	echo "Error: Input file does not exist: "$1"  "
	exit 1
fi


##################################################################################################################################
## Common variables to PBS and PBSpro
## and other needed variables that a user should not need to modify
##################################################################################################################################
PBSUSER="#PBS -M "$UEMAIL""
HASHBANG="#!/bin/bash"

##################################################################################################################################
## PBS torque and PBSpro have some differences. 
## Organise these here and also check further on (line 184) and change NODETYPE to match cluster system
## MODTRINITY will also be different on different clusters - it sets up the paths to the required executables
##################################################################################################################################
if [[ "$PBSTYPE" = "--pbspro" ]] ; then	
	JOBARRAY="-J"
	JOBARRAY_ID="\$PBS_ARRAY_INDEX"
	AFTEROKARRAY="afterok"	
elif [[ "$PBSTYPE" = "--pbs" ]] ; then
	## PBS torque:
	JOBARRAY="-t"
	JOBARRAY_ID="\$PBS_ARRAYID"
	AFTEROKARRAY="afterokarray"	
else # no paramaters present
	F_USAGE	
	exit 0
fi

##############################################################################################################################################
##########     Part 1: Set up file names for input directory and for output data dir and Trinity.pl command line          ####################
##############################################################################################################################################
	echo ""
	## Ensure JOBPREFIX is not greatr than 11 characters as PBS-pro can not handle > 15 characters for total job name length
	echo "submitting trinity jobs with prefix: "
	echo "		$JOBPREFIX"	

	###### Set input data directory - $DATADIR is CSIRO specific	
	echo "Input directory: "
	echo "		$DATADIRECTORY"

	###### Set output data directory (OUTPUTDIR)
	echo "Output directory: Scripts and output data will be written to:"
	echo "		$OUTPUTDIR" 
	mkdir -p "$OUTPUTDIR"
	cd "$OUTPUTDIR"

	
	### Modify STANDARD_JOB_DETAILS for analysis specific input to Trinity.pl 
	echo "The following trinity command line will be run:"
	echo "$STANDARD_JOB_DETAILS" 
	echo "" 
	if [[ "$1" = "--pbs" ]] ; then
		echo "		Use: \"qpeek -t PBS_JOBID\" "
		echo "		To view stdout and stderr from each separate job while they are running" 
	fi
	
###########################################################################################################################################################		
###		Do some checking that files exist etc. (User should not modify)  
###		This sets the $DS variable
F_CHECKFILES "$FILENAMEINPUT"    "$DATADIRECTORY"   "$FILENAMESINGLE"  "$FILENAMELEFT"   "$FILENAMERIGHT"

###########################################################################################################################################################	
#################################                                                                                       ###################################
#################################       Part 2: Create the shell scripts to be run via the PBS batch system             ###################################
#################################       Users should modify WALLTIME and MEM dependent upon dataset size                ###################################
#################################       and NCPU to appropriate value for compute node cpu resources                    ###################################
#################################       Check: "Trinity RNA-seq Assembler Performance Optimisation" (Henschel 2012)     ###################################
#################################       for current best practice.                                                      ###################################
#################################       N.B. On busy clusters it may be best not to try to request                      ###################################
#################################       a full nodes resources. i.e. If 8 cores per node are present,                   ###################################
#################################       only request half of these.                                                     ###################################
###########################################################################################################################################################

###########################################################################################################################################################
##############################  Script 1: Write script to run Inchworm                                         ############################################
#############################   MEM should equal JFMEM, which is the amount of memory requested for Jellyfish ############################################

JOBNAME1="$JOBPREFIX"_p1
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P1" "$NCPU_P1" "$PBSQUEUE_P1" "$WALLTIME_P1" "$JOBNAME1" "$ACCOUNT" "$PBSUSER" "$MODTRINITY" )

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p1"

#############################################################################################################################################################
##############################  Script 2: Chrysalis::GraphFromFasta                                             #############################################
##############################  This script has a dependency on part 1 completion without error.                #############################################
# It would be good to force an exit(0) before ReadsToTranscripts after checkpoint file /chrysalis/GraphFromIwormFasta.finished is 
# written (i.e. add --no_run_readstotrans to Trinity and pass through to Chrysalis) as the script 3  can be started directly after.
# This may be less of an issue with the new (fast) version of Trinity::GraphFromFasta (since version 2012-06-08).


JOBNAME2="$JOBPREFIX"_p2 
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P2" "$NCPU_P2" "$PBSQUEUE_P2" "$WALLTIME_P2" "$JOBNAME2" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p2"

###########################################################################################################################################################
##############################  Script 3: Script to run Chrysalis::ReadsToTranscripts                                                   ###################
##############################  ReadsToTranscripts can be slow due to reads from disk,                                                  ###################
##############################  This script has a dependency on part 2 completion with error.                                            ###################
##############################  Section script is skipped if enough time was given in Part 2.                                           ###################


JOBNAME3="$JOBPREFIX"_p3 
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P3" "$NCPU_P3" "$PBSQUEUE_P3" "$WALLTIME_P3" "$JOBNAME3" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p3"

##########################################################################################################################################################
############################## Script 4a: Write script to call the Chrysalis QuantifyGraph array job                                    ##################
############################## N.B. SLOTLIMIT="%x" indicates 'slot limit' i.e. the number of concurrent jobs to execute in an array     ##################
############################## (Not available in PBSpro )                                                                               ##################


#SLOTLIMIT="%64"   # available for PBS Torque
JOBNAME4="$JOBPREFIX"_p4a
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" 1gb 1 "$PBSQUEUE_P4" 00:02:00 "$JOBNAME4" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p4a"

###########################################################################################################################################################
##############################   Script Array part 4b: Write script to be run as an array Job . Runs Chrysalis::QuantifyGraph            ##################
##############################   This scipt has a dependency on part 4a being run. If an array component fails it will email user.       ##################
##############################   Files named quantifyGraph_commands_X are written with subset of total commands (X is the array ID from the PBS system)  ##
##############################   quantifyGraph_commands_X.fin will be present in chrysalis output directory if all commands have been executed  ###########

NCPU="1"
PBSUSER=""
JOBNAMEARRAY1="$JOBPREFIX"_p4b
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P4" "$NCPU" "$PBSQUEUE_P4" "$WALLTIME_P4" "$JOBNAMEARRAY1" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p4b"

###########################################################################################################################################################
############################  Script 5a: Writes script to call the butterfly array job and then collate the results (part 6)             ##################
############################  This script has a dependency on part 4b array jobs completion without error.                               ##################

JOBNAME5=""$JOBPREFIX"_p5a"
JOBNAME6=""$JOBPREFIX"_p6"
PBSUSER="#PBS -M "$UEMAIL""
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" 2gb 1 "$PBSQUEUE_P5" 00:05:00 "$JOBNAME5" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p5a"

###########################################################################################################################################################
##############################   Script Array 5b: Write script to run Butterfly Array Job                                                ##################
##############################   This script has a dependency on part 5a and 4b being run.                                               ##################
##############################   Files named butterfly_commands.adj_X are written with subset of total commands (X is the array ID from the PBS system) ###
##############################   butterfly_commands.adj_X.fin will be present in chrysalis output directory if all commands have been executed ############

NCPU="1"
JOBNAMEARRAY2=""$JOBPREFIX"_p5b"
PBSUSER=""
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" "$MEM_P5" "$NCPU" "$PBSQUEUE_P5" "$WALLTIME_P5" "$JOBNAMEARRAY2" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p5b"

############################################################################################################################################################
###############################  Script 6: Script to gather results from parallel sections                                               ###################
###############################  This scipt has a dependency on part 5b array jobs completion without error.                             ###################
###############################  Collate results script - final results end in Trinity.fasta                                             ###################
###############################  This script is started by Script 5a and is dependant on 5b completing                                   ###################

PBSQUEUE="any"
PBSUSER="#PBS -M "$UEMAIL""
NODESCPUS=$(F_GETNODESTRING "$PBSTYPE" 4gb 1 "$PBSQUEUE" "$WALLTIME_P6" "$JOBNAME6" "$ACCOUNT" "$PBSUSER" "$MODTRINITY")

F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.p6"

############################################################################################################################################################
###############                                                                                                                          ###################
###############  Part 3:  Write main control script that executes scripts that were created above.                                       ###################
###############                                                                                                                          ###################
############################################################################################################################################################
F_WRITESCRIPT "$0" ""$RUNDIR"/trinity_pbs.cont"

###########################################################################################################################################################
######### Exit here if we only want to test script creation of the above code                                                            ##################
# echo "Exiting without submitting scripts" ; exit 0 

###############  Run the script written in Part 3
###############  Checks to see what is current stage of calculation and executes scripts created in above code                           ###################
./"$JOBPREFIX"_run.sh



exit 0

