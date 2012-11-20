##################################################################################################################################
##########################                                                                ########################################
##########################     Trinity PBS job submission with multi part dependencies    ########################################
##########################                                                                ########################################
##################################################################################################################################
### Author: Josh Bowden, CSIRO
### Version 0.1
### Create Trinity.fasta file p6 Script
##################################################################################################################################


JOBSTRING6=""$HASHBANG"
"$NODESCPUS"

 cd "$OUTPUTDIR"
 # this collects output data and saves to Trinity.fasta
 find chrysalis/ -name "*allProbPaths.fasta" -exec cat {} \; > Trinity.fasta
 
 # Some other logical steps we could do in this stage:
 # Move output to backed up area (home directory):
 # cp Trinity.fasta ~/"$JOBPREFIX"_Trinity.fasta
 # Remove the working directory:
 # 
 "
		
	# Write the JOBSTRING6 to a file for later execution
	echo "${JOBSTRING6}" | cat -> ""$JOBNAME6".sh"

