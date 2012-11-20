################################################################################################################
########################     README file                                                ########################
########################     Trinity PBS job submission with multi part dependencies    ########################
########################     Author: Josh Bowden, CSIRO                                 ########################
########################     Version 0.1                                                ########################
########################     Based on script by Alexie Papanicolaou, CSIRO              ########################
################################################################################################################

DESCRIPTION: trinity_pbs.sh is a BASH shell script used for submission of Trinity jobs to clusters that use PBS Torque or PBS Pro.
The set of scripts stages the parts of the Trinity workflow into 6 stages:
	1/ Inchworm
	2,3/ Chrysalis::GraphFromFasta and Chrysalis::ReadsToTranscripts
	4/ Chrysalis::QuantifyGraph
	5/ Butterfly
	6/ Gather together resulting transcripts into "Trinity.fasta" file
	
Each stage is submitted as a PBS job, with dependencies i.e. the following stage will only execute after successfull completion of the stage before (an exception to this is stage 2 and 3). Due to their parallel nature, stages 4 and 5 are submitted as 'array jobs', with each job made up of a user defined number of subtasks.

ADMINISTRATION setup instructions:
	trinity_pbs script install instructions:
		1. Copy all trinity_pbs.* files into a directory (we will call it "TRINITY_PBS_DIR") in a user accessible, read only, area.
		2. Add TRINITY_PBS_DIR to the PATH i.e. export or set PATH=TRINITY_PBS_DIR:$PATH  (perhaps export PATH in .bashrc file)
		In the file "trinity_pbs.sh", do the following:
		3. Change the "RUNDIR" variable to point to the directory installation directory. i.e. RUNDIR=TRINITY_PBS_DIR
		4. Set MEMDIRIN to the name of a node-local filesystem so a network drive is not needed unecesarily for Scripts 4b and 5b
		5. Set MODTRINITY to the path to the trinity executables and/or modules to be loaded so Trinity.pl can be run.
		6. Set PBSTYPE to --pbspro or --pbs, dependent on the system present.
	That should be all that is needed from an admin perspective. These scripts have been tested on PBS Torque 3.0.6 and PBS Pro 11.0.2. There may be PBS system specific changes due to PBS version incompatabilities.
	
USER setup instructions:
		Users should make a copy of TRINITY.CONFIG (possibly a copy for each job they want to run) and then modify variables in it to suit the system the job is being run on (number of CPUs, amount of memory) and the expected runtime of the Trinity process (which is a function of the datset size). Further instructions are provided within the TRINITY.CONFIG file.
		


