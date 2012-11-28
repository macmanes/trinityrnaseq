#!/usr/bin/env python

#    preprocess.py
#    November 2012  Matthew MacManes (macmanes@gmail.com)
#
#   This wrapper is free software: you can redistribute it and/or modify
#
#  v.0.0.7
#
#ISSUES
#
#clean version!
#1st version on git
#

import sys
import subprocess
import optparse
import shutil
import os
from datetime import datetime, date, time
import os.path
import glob

print ""
print ""
print ""
print "***********************************************************************"
print "***   preproc.py                                 "
print "***   To run this program, you must have AllPathsLG (R43869 or newer) "
print "***   installed and in your $PATH                           "
print "***********************************************************************"
print ""
print "Execution of this script with the '--full True' option selected will produce 5 files: corr* files have been Error Corrected but not normalized"
print "norm.corr* files have been Error Corrected and normalized. Unpaired reads have been concatenated to the /1 file"
print "These norm.corr* files are apprproate for assembly using trinity.pl"
print ""


##########################################
## date function
##########################################

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

##########################################
## Options
##########################################
def getOptions():
    parser = optparse.OptionParser(usage="usage: python %prog -m [int]  -p [33 or 64] -p [num threads] -c [max_cov] -l left.fq -r right.fq -u unpaired.fq]",
                          version="%prog 0.0.7")
    parser.add_option("-m", "--memory",
                      dest="memory",
                      default="4",
                      metavar='memory',
                      help="How much memory can you allocate to processes")
    parser.add_option("-p", "--phred",
                      dest="phred",
                      metavar='phred',
                      default="33",
                      help="Casava 1.8+ should be -p33",)
    parser.add_option("-t", "--threads",
                      dest="threads",
                      metavar='[INT]',
                      default="2",
                      help="Number of threads to use",)
    parser.add_option("-l", "--left",
                      dest="left",
                      metavar='file.fq',
                      default="",
                      help="left reads",)
    parser.add_option("-u", "--unpaired",
                      dest="unpaired",
                      metavar='unpaired reads',
                      default="",
                      help="full path to unpaired reads")
    parser.add_option("-r", "--right",
                      dest="right",
                      metavar='file.fq',
                      default="",
                      help="right reads",)
    parser.add_option("-o", "--out",
                      dest="out",
                      metavar='how will reads be named',
                      default="read",
                      help="how will reads be named",)
    parser.add_option("-c", "--max_cov",
                      dest="maxcov",
                      metavar='maxcov',
                      default="40",
                      help="for normalization: what is the max coverage: Recommended 30-100X")
    parser.add_option("-k", "--max_kmer",
                      dest="maxkmer",
                      metavar='maxkmer',
                      default="0",
                      help="max kmer freq to mark for removal")
    parser.add_option("-e", "--EC_K",
                      dest="eck",
                      metavar='eck',
                      default="25",
                      help="kmer for error correction, whould always be 25")
    parser.add_option("-F", "--fill_frags",
                      dest="FF",
                      metavar='FF',
                      default="False",
                      help="use Fill Fragment Module")
    parser.add_option("-H", "--haploidify",
                      dest="hap",
                      metavar='hap',
                      default="False",
                      help="Attempt to haploidify")
    parser.add_option("--error_correct",
                      dest="dec",
                      metavar='dec',
                      default="False",
                      help="Do you ONLY want to error correct reads?")                                      
    parser.add_option("--full",
                      dest="full",
                      metavar='full',
                      default="False",
                      help="Do you want to error correct and normalize?")                                        
    (options, args) = parser.parse_args()
    return options

##########################################
## Error Correct
##########################################
BASE_DIR = os.path.join( os.path.dirname(sys.argv[0]), '../../' )
def error_corr(options):
        ec = subprocess.Popen(['ErrorCorrectReads.pl', 
		'MAX_MEMORY_GB=', options.memory, 
		'THREADS=', options.threads, 
		'PHRED_ENCODING=', options.phred, 
		'PAIRED_READS_A_IN=', options.left,  
		'PAIRED_READS_B_IN=', options.right, 
		'UNPAIRED_READS_IN=', options.unpaired,
		'FE_MAX_KMER_FREQ_TO_MARK=', options.maxkmer, 
		'EC_K=', options.eck, 
		'HAPLOIDIFY=', options.hap, 
		'FILL_FRAGMENTS=', options.FF, 
		'FF_K=25', 
		'READS_OUT=', options.out])
        output = ec.communicate()
        assert ec.returncode == 0, output[0] + "Error Correction failed\n"
##########################################
## Fastool
##########################################
def fastool1(options):
	log1 = open("corr.left.fa", "w")
        fastool1 = subprocess.Popen([BASE_DIR + 'trinity-plugins/fastool/fastool', 
	'--append', '/1', 
	'--to-fasta', options.out + '.paired.A.fastq'], stdout=log1)
        output = fastool1.communicate()
        assert fastool1.returncode == 0, output[0] + "FasTool1 Failed\n"
def fastool2(options):
	log2 = open("corr.right.fa", "w")
        fastool2 = subprocess.Popen([BASE_DIR + 'trinity-plugins/fastool/fastool', 
	'--append', '/2', 
	'--to-fasta', options.out + '.paired.B.fastq'], stdout=log2)
        output = fastool2.communicate()
        assert fastool2.returncode == 0, output[0] + "FasTool2 Failed\n"
def fastool3(options):
	log3 = open("corr.unpaired.fa", "w")
        fastool3 = subprocess.Popen([BASE_DIR + 'trinity-plugins/fastool/fastool', 
	'--append', '/1', 
	'--to-fasta', options.out + '.unpaired.fastq'], stdout=log3)
        output = fastool3.communicate()
        assert fastool3.returncode == 0, output[0] + "FasTool3 Failed\n"
def fastool4(options):
	log4 = open("corr.single.fa", "w")
        fastool4 = subprocess.Popen([BASE_DIR + 'trinity-plugins/fastool/fastool', 
	'--append', '/1', 
	'--to-fasta', options.out + '.fastq'], stdout=log4)
        output = fastool4.communicate()
        assert fastool4.returncode == 0, output[0] + "FasTool4 Failed\n"
##########################################
## Normalize Pairs
##########################################
def normalize_pairs(options):
        normalize_pairs = subprocess.Popen([BASE_DIR + 'util/normalize_by_kmer_coverage.pl', 
	'--seqType', 'fa', 
	'--max_cov', options.maxcov, 
	'--left', 'corr.left.fa', 
	'--right', 'corr.right.fa', 
	'--JM', options.memory + 'G', 
	'--JELLY_CPU', options.threads])
        output = normalize_pairs.communicate()
        assert normalize_pairs.returncode == 0, output[0] + "Pair Normalization failed\n"
##########################################
## Normalize UnPaired Reads
##########################################
def normalize_unpaired(options):
        normalize_unpaired = subprocess.Popen([BASE_DIR + 'util/normalize_by_kmer_coverage.pl', 
	'--seqType', 'fa', 
	'--max_cov', options.maxcov, 
	'--single', 'corr.unpaired.fa', 
	'--JM', options.memory + 'G', 
	'--JELLY_CPU', options.threads])
        output = normalize_unpaired.communicate()
        assert normalize_unpaired.returncode == 0, output[0] + "UnPaired Normalization failed\n"
##########################################
## Normalize Single Reads
##########################################
def normalize_single(options):
        normalize_single = subprocess.Popen([BASE_DIR + 'util/normalize_by_kmer_coverage.pl', 
	'--seqType', 'fa', 
	'--max_cov', options.maxcov, 
	'--single', 'corr.single.fa', 
	'--JM', options.memory + 'G', 
	'--JELLY_CPU', options.threads])
        output = normalize_single.communicate()
        assert normalize_single.returncode == 0, output[0] + "Single-end Normalization failed\n"
##########################################
## Concatenate
##########################################
def concatenate(options):
	destination = open('norm.corr.left.fa','wb')
	shutil.copyfileobj(open('corr.left.fa','rb'), destination)
	shutil.copyfileobj(open('corr.unpaired.fa','rb'), destination)
	destination.close()

#########################################
## alignment depend
##########################################
def checkAllPaths():
        try:
            p = subprocess.Popen(['RunAllPathsLG'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError:
            print "Could not find AllPaths"
            print "Make sure that it is properly installed on your $PATH"
            sys.exit(1)
#########################################
## Do this if --full True
##########################################
def main():
        options = getOptions()
        checkAllPaths()    
        if options.full=='True':
            print >> sys.stderr,"\nRunning Error Correction Module, If Necessary: [%s] \n" % (right_now())
            if os.path.exists(options.out + '.paired.B.fastq'):
                    print >> sys.stderr,"\Error Correction apprears to be complete! I'm going straight to the Normalization Steps. \n"
            else: 
                    error_corr(options)
            print ""
            print "*********************************************************************************"	
            print "*********************************************************************************"
            print "Running FasTool conversion: [%s] \n" % (right_now())
            print "*********************************************************************************"
            print "*********************************************************************************"	
            print ""
            if os.path.exists(options.out + '.paired.A.fastq'):
                    fastool1(options)
                    os.remove(options.out + '.paired.A.fastq')
            else: 
                    print >> sys.stderr,"working.."		
            if os.path.exists(options.out + '.paired.B.fastq'):
                    fastool2(options)
                    os.remove(options.out + '.paired.B.fastq')
            else: 
                    print >> sys.stderr,"working.."		
            if os.path.exists(options.out + '.unpaired.fastq'):
                    fastool3(options)
                    os.remove(options.out + '.unpaired.fastq')
            else: 
                    print >> sys.stderr,"working.."			
            if os.path.exists(options.out + '.fastq'):
                    os.remove(options.out + '.fastq')
            else: 
                    print >> sys.stderr,"working.."			
            print ""
            print "*********************************************************************************"	
            print "*********************************************************************************"
            print "Running Normalization Steps for Paired Reads: [%s] \n" % (right_now())
            print "*********************************************************************************"
            print "*********************************************************************************"	
            print ""
            if os.path.exists('corr.left.fa'):
                    normalize_pairs(options)
                    for filename in glob.glob('corr.right.fa.normalized_K25_C*_pctSD*.fa'):
                            os.rename(filename,'norm.corr.right.fa')
            else: 
                    print >> sys.stderr,"working.."		
            print ""
            print "*********************************************************************************"	
            print "*********************************************************************************"
            print "Running Normalization Steps for UnPaired Reads: [%s] \n" % (right_now())
            print "*********************************************************************************"
            print "*********************************************************************************"	
            print ""
            if os.path.exists(options.out + '.unpaired.fa'):
                    normalize_unpaired(options)
            else: 
                    print >> sys.stderr,"working.."		
            print ""
            print "*********************************************************************************"	
            print "*********************************************************************************"
            print "Concatenate and Clean up: [%s] \n" % (right_now())
            print "*********************************************************************************"
            print "*********************************************************************************"	
            print ""
            if os.path.exists('corr.unpaired.fa'):
                    concatenate(options)	
                    for filename in glob.glob('corr.left.fa.normalized_K25_C*_pctSD*.fa'):
                            os.remove(filename) # Need to make sure if ok is something other than C40
            else: 
                    print >> sys.stderr,"should have concat.."	
                    for filename in glob.glob('corr.left.fa.normalized_K25_C*_pctSD*.fa'):
                            os.rename(filename,'norm.corr.left.fa')
            os.remove(options.out + '.fastq.ids')
            shutil.rmtree('normalized_reads')
            print ""
            print "*********************************************************************************"	
            print "*********************************************************************************"
            print "\nDone.. Happy Assembling! [%s] \n" % (right_now())
            print "*********************************************************************************"
            print "*********************************************************************************"	
            print ""
#########################################
## Do this if --error_corr True
##########################################
        else:
            if options.dec=='True':
                print >> sys.stderr,"\nRunning Error Correction Module: [%s] \n" % (right_now())
                if os.path.exists(options.out + '.paired.B.fastq'):
                        print >> sys.stderr,"\Error Correction apprears to be complete! \n"
                else: 
                        error_corr(options)
                print ""
                print "*********************************************************************************"	
                print "*********************************************************************************"
                print "Running FasTool conversion: [%s] \n" % (right_now())
                print "*********************************************************************************"
                print "*********************************************************************************"	
                print ""
                if os.path.exists(options.out + '.paired.A.fastq'):
                        fastool1(options)
                        os.remove(options.out + '.paired.A.fastq')
                else: 
                        print >> sys.stderr,"working.."		
                if os.path.exists(options.out + '.paired.B.fastq'):
                        fastool2(options)
                        os.remove(options.out + '.paired.B.fastq')
                else: 
                        print >> sys.stderr,"working.."		
                if os.path.exists(options.out + '.unpaired.fastq'):
                        fastool3(options)
                        os.remove(options.out + '.unpaired.fastq')
                else: 
                        print >> sys.stderr,"working.."			
                if os.path.exists(options.out + '.fastq'):
                        os.remove(options.out + '.fastq')
                else: 
                        print >> sys.stderr,"working.."
            else:
                print >> sys.stderr,"\n\n\n*** You must specify either '--full True' if you want to error correct AND normalize, or '--error_corr True' if you want to do just error correction. \n\n\n"
if __name__ == "__main__":
        main()
