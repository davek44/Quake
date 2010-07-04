#ifndef EDIT_H
#define EDIT_H

#include <string>
#include <vector>
#include <iostream>

using namespace::std;

////////////////////////////////////////////////////////////////////////////////
// options
////////////////////////////////////////////////////////////////////////////////
extern char* fastqf;
extern char* file_of_fastqf;
extern int threads;
extern unsigned int chunks_per_thread;

////////////////////////////////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////////////////////////////////
void combine_output(string fqf, string mid_ext);
void combine_output_paired(string fqf1, string fqf2, string mid_ext);
void chunkify_fastq(string fqf, vector<streampos> & starts, vector<unsigned long long> & counts);
void guess_quality_scale(string fqf);
vector<string> parse_fastq(vector<string> & fastqfs, vector<int> & pairedend_codes);
void unzip_fastq(string & fqf);
void zip_fastq(string fqf);
vector<string> split(string s, char c);
vector<string> split(string);
#endif
