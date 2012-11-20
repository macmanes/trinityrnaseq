#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "irke_common.hpp"
#include "argProcessor.hpp"

unsigned int IRKE_COMMON::MONITOR = 0;


int KMER_SIZE = 25; 

int MAX_THREADS = 6;

void populate_kmer_counter (KmerCounter& kcounter, string& kmers_fasta_file);
vector<int> compute_kmer_coverage (string& sequence, KmerCounter& kcounter);
int median_coverage (vector<int>);
string usage (ArgProcessor args);
int sum ( vector<int>& vals);
float average( vector<int>& vals);
float stDev( vector<int>& vals);

int main (int argc, char* argv[]) {

    ArgProcessor args(argc, argv);
    
    if (args.isArgSet("--help")
        ||
        (! (args.isArgSet("--reads") && args.isArgSet("--kmers")) )
               
        ) {
        
        cerr << usage(args) << endl << endl;
        exit(1);
    }
    

    string reads_fasta_file = args.getStringVal("--reads");
    
    string kmers_fasta_file = args.getStringVal("--kmers");
    bool is_DS = (! args.isArgSet("--SS"));

    if (args.isArgSet("--kmer_size")) {
        KMER_SIZE = args.getIntVal("--kmer_size");
        if (KMER_SIZE < 20) {
            cerr << "Error, min kmer size is 20";
            exit(2);
        }
    }
    
    if (args.isArgSet("--monitor")) {
        IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
    }
    
    if (omp_get_max_threads() > MAX_THREADS) {
        omp_set_num_threads(MAX_THREADS);
    }
    
    KmerCounter kcounter (KMER_SIZE, is_DS);

    populate_kmer_counter(kcounter, kmers_fasta_file);
    
    
    Fasta_reader fasta_reader(reads_fasta_file);
    
    ofstream* filewriter = NULL;
    ofstream* covwriter = NULL;    
    bool write_coverage_info = args.isArgSet("--capture_coverage_info");
    
    
    while (true) {
        Fasta_entry fe = fasta_reader.getNext();
        string sequence = fe.get_sequence();
        if (sequence == "") break;
        
        
        string header = fe.get_header();

        vector<int> kmer_coverage = compute_kmer_coverage(sequence, kcounter);
        
        int median_cov = median_coverage(kmer_coverage);
        float average_cov = average(kmer_coverage);
        float stdev = stDev(kmer_coverage);
        float pct_stdev_of_avg = stdev/average_cov*100;



        stringstream stats_text;
        stats_text << median_cov << "\t" 
                   << average_cov << "\t"
                   << stdev << "\t" 
                   << pct_stdev_of_avg << "\t"
                   << fe.get_accession();


        
        if (write_coverage_info) {
            // add the coverage info
            stats_text << "\t";
            for (int i = 0; i < kmer_coverage.size(); i++) {
                stats_text<< kmer_coverage[i];
                if (i != kmer_coverage.size() - 1) {
                    stats_text<< ",";
                }
            }
        }
        
        stats_text << endl;
        
        
        cout << stats_text.str();
        
    }
    
            
    return(0);
}

void populate_kmer_counter (KmerCounter& kcounter, string& kmers_fasta_file) {

    // code largely copied from IRKE.cpp

    int i, myTid;
    unsigned long sum, 
        *record_counter = new unsigned long[omp_get_max_threads()];
    unsigned long start, end;

    // init record counter
    for (int i = 0; i < omp_get_max_threads(); i++) {
        record_counter[i] = 0;
    }

    cerr << "-reading Kmer occurences..." << endl;
    start = time(NULL);

    Fasta_reader fasta_reader(kmers_fasta_file);

    #pragma omp parallel private (myTid)
    {
        myTid = omp_get_thread_num();
        record_counter[myTid] = 0;
        
        while (true) {
            Fasta_entry fe = fasta_reader.getNext();
            if (fe.get_sequence() == "") break;

            record_counter[myTid]++;

            if (IRKE_COMMON::MONITOR) {
                if (myTid == 0 && record_counter[myTid] % 100000 == 0)
                    {
                        sum = record_counter[0];
                        for (i=1; i<omp_get_num_threads(); i++)
                            sum+= record_counter[i];
                        cerr << "\r [" << sum/1000000 << "M] Kmers parsed.     ";
                    }
            }


            string seq = fe.get_sequence();
            if (seq.length() != KMER_SIZE) {
                cerr << "ERROR: kmer " << seq << " is not of length: " << KMER_SIZE << endl;
                continue;
            }
            
            kmer_int_type_t kmer = kcounter.get_kmer_intval(seq);
            int count = atoi(fe.get_header().c_str());
            kcounter.add_kmer(kmer, count);
        }
    }
    end = time(NULL);

    sum = record_counter[0];
    for (i=1; i<omp_get_max_threads(); i++)
        sum+= record_counter[i];
    delete [] record_counter;

    cerr << endl << " done parsing " << sum << " Kmers, " << kcounter.size() << " added, taking " << (end-start) << " seconds." << endl;

    return;
    

}

vector<int> compute_kmer_coverage(string& sequence, KmerCounter& kcounter) {
    
    vector<int> coverage;

    if (IRKE_COMMON::MONITOR) {
        cerr << "processing sequence: " << sequence << endl;
    }

    for (int i = 0; i <= (int) sequence.length() - KMER_SIZE; i++) {
        
        // cerr << "i: " << i << ", <= " << sequence.length() - KMER_SIZE << endl;

        string kmer = sequence.substr(i, KMER_SIZE);
     
        if (IRKE_COMMON::MONITOR >= 2) {
            for (int j = 0; j <= i; j++) {
                cerr << " ";
            }
            cerr << kmer << endl;
        }
        
   
        int kmer_count = 0;
        if (! contains_non_gatc(kmer)) {

            kmer_count = kcounter.get_kmer_count(kmer);
            
        }
                
        coverage.push_back(kmer_count);
        
    }

    return(coverage);
    
}
 

int median_coverage (vector<int> coverage) {
    

    int num_pts = coverage.size();
    if (num_pts == 0) {
        return(0);
    }
    else if (num_pts == 1) {
        return(coverage[0]);
    }

    sort(coverage.begin(), coverage.end());
    
    if (num_pts % 2 > 0) {
        // odd value
        int midpt = num_pts / 2;
        return(coverage[midpt]);
    }
    else {
        // even 
        int midpt = num_pts / 2;
        int median = (coverage[midpt-1] + coverage[midpt])/(float)2;
        return(median);
    }
}


string usage(ArgProcessor args) {

    stringstream usage_info;

    usage_info << endl << endl;

    usage_info << "Usage: " << endl
               << "  --reads  <str>             " << ":fasta file containing reads" << endl
               << "  --kmers  <str>             " << ":fasta file containing kmers" << endl
               << "* optional:" << endl
               << "  --kmer_size <int>          " << ":default = 25" << endl
               << "  --DS                             " << ":double-stranded RNA-Seq mode (not strand-specific)" << endl
               << "  --capture_coverage_info    " << ":writes coverage info file." << endl
               << "  --monitor <int>            " << ":verbose level for debugging" << endl
               << endl;
    

    return(usage_info.str());
}

int sum ( vector<int>& vals) {

    int sum_vals = 0;
    for (int i = 0; i < (int) vals.size(); i++) {
        sum_vals += vals[i];
    }

    return(sum_vals);
}

float average ( vector<int>& vals) {

    int sum_vals = sum(vals);
    float avg = (float)sum_vals / vals.size();

    return(avg);
}

float stDev( vector<int>& vals) {

    float avg = average(vals);

    int num_vals = vals.size();

    float sum_avg_diffs_sqr = 0;
    for (int i = 0; i < (int) vals.size(); i++) {
        int val = vals[i];
        float delta = val - avg;
        float sqr =  pow(delta,2);
        sum_avg_diffs_sqr += sqr;
    }

    float stdev = sqrt( sum_avg_diffs_sqr/(num_vals-1) );

    return(stdev);
}
