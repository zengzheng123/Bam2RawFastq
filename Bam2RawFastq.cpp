//
//  main.cpp
//  Bam2RawFastq
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 10/8/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include "api/BamReader.h"
#include <gzstream.h>


using namespace std;
using namespace BamTools;

const string VERSION = "Bam2RawFastq 1.0.0";


string bam_input_file;
string fastq_input_list;
string fastq_output_prefix;
bool pair_end_read = true;

vector< pair<string, string> > fastq_input_pairs;

void printUsage(string msg = "")
{
    cout << endl;
    cout << VERSION << endl;
    cout << "Usage: " << endl;
    cout << "\t--bam                   <string>                        Input bam file" << endl;
    cout << "\t--fastq_list            <string>                        Comma seperated list of raw fastq gzipped files. For pair end read, the list need to be in correct order: R1_file1,R2_file1,R1_file2, R2_file2,..." << endl;
    cout << "\t--output                <string>                        Output file prefix, for PE read, the output files will be <prefix>_R1.fastq.gz & <prefix>_R2.fastq.gz. For single end read, the output file will be <prefix>.fastq.gz" << endl;
    cout << "\t--single_end                                            The script assumes the the data is pair end read by default, use this option if the data is single end read" << endl;
    cout << "\t--help                                                  Print command line usage" << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"bam",                     required_argument,      0,     'b'},
    {"fastq_list",              required_argument,      0,     'l'},
    {"output",                  required_argument,      0,     'o'},
    {"single_end",              no_argument,            0,     's'},
    {"help",                    no_argument,            0,     'h'},
    {0, 0, 0, 0}
};


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "b:l:o:sh", long_options, &option_index);
        switch(next_option)
        {
            case 'b':
                bam_input_file = optarg;
                break;
            case 'l':
                fastq_input_list = optarg;
                break;
            case 'o':
                fastq_output_prefix = optarg;
                break;
            case 's':
                pair_end_read = false;
                break;
            case 'h':
                printUsage();
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR] Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    
    if(bam_input_file.empty())
        printUsage("[ERROR] Please specify input bam file");
    if(fastq_input_list.empty())
        printUsage("[ERROR] Please specify raw fastq files with --fastq_list");
    if(fastq_output_prefix.empty())
        printUsage("[ERROR] please specify output file prefix");
#ifdef _DEBUG
    cerr << "[DEBUG] Parsing options complete." << endl;
#endif
}


void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cerr << "[DEBUG] Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // use to skip empty item
        parsed_item.push_back(item);
    }
}


void prepare_input_fastq_list()
{
    vector<string> fastq_vec;
    split(fastq_input_list, ',', fastq_vec);
    if(pair_end_read)
    {
        if(fastq_vec.size() % 2 != 0)
        {
            cerr << "[ERROR]: Data is specified as pair end read, but number raw fastq files provided with --fastq_list is not an even number" << endl;
            exit(1);
        }
        for(size_t i = 0; i < fastq_vec.size(); i += 2)
            fastq_input_pairs.push_back(make_pair(fastq_vec[i], fastq_vec[i+1]) );
    }
    else
    {
        for(size_t i = 0; i < fastq_vec.size(); i++)
            fastq_input_pairs.push_back(make_pair(fastq_vec[i], "") );
    }
}



void print_read(map<string, int> &readname_hash, string fastq_input_file1, string fastq_input_file2, ogzstream *output_fastq1_fs, ogzstream *output_fastq2_fs)
{
    cout << "[INFO]: processing raw fastq files: " << fastq_input_file1 << ", " << fastq_input_file2 << endl;
    igzstream input_fs1;
    input_fs1.open(fastq_input_file1.c_str());
    if(!input_fs1)
    {
        cerr << "[ERROR]: failed to open input file: " << fastq_input_file1 << endl;
        exit(1);
    }
    
    igzstream input_fs2;
    if(pair_end_read)
    {
        input_fs2.open(fastq_input_file2.c_str());
        if(!input_fs2)
        {
            cerr << "[ERROR]: failed to open input file: " << fastq_input_file2 << endl;
            exit(1);
        }
    }
    
    string end1_line1, end1_line2, end1_line3, end1_line4, end2_line1, end2_line2, end2_line3, end2_line4;
    string read_name1, read_name2;
    while(true)
    {
        getline(input_fs1, end1_line1);
        getline(input_fs1, end1_line2);
        getline(input_fs1, end1_line3);
        getline(input_fs1, end1_line4);
        if(input_fs1.eof())
            break;
        
        char name1_tmp[10000];
        sscanf(end1_line1.c_str(), "@%[^/\t ]", name1_tmp);
        read_name1 = name1_tmp;

        if(pair_end_read)
        {
            getline(input_fs2, end2_line1);
            getline(input_fs2, end2_line2);
            getline(input_fs2, end2_line3);
            getline(input_fs2, end2_line4);
            if(input_fs2.eof())
                break;
            
            char name2_tmp[10000];
            sscanf(end2_line1.c_str(), "@%[^/\t ]", name2_tmp);
            read_name2 = name2_tmp;
        }
        
        if(pair_end_read && read_name1 != read_name2)
        {
            cerr << "[ERROR]: Input pair end fastq files are not synchronized: " << read_name1 << ", " << read_name2 << endl;
            exit(1);
        }
        
        if(readname_hash.find(read_name1) != readname_hash.end())
        {
            (*output_fastq1_fs) << end1_line1 << endl;
            (*output_fastq1_fs) << end1_line2 << endl;
            (*output_fastq1_fs) << end1_line3 << endl;
            (*output_fastq1_fs) << end1_line4 << endl;
            if(pair_end_read)
            {
                (*output_fastq2_fs) << end2_line1 << endl;
                (*output_fastq2_fs) << end2_line2 << endl;
                (*output_fastq2_fs) << end2_line3 << endl;
                (*output_fastq2_fs) << end2_line4 << endl;
            }
            readname_hash[read_name1] = 1;
        }
    }
    input_fs1.close();
    input_fs2.close();
}


void convert_read()
{
    BamReader my_bam_reader;
    if(!my_bam_reader.Open(bam_input_file))
    {
        cerr << "[ERROR] fail to open input bam file: " << bam_input_file << endl;
        exit(1);
    }
    
    string fastq_output_file1;
    string fastq_output_file2;
    ogzstream output_PE1_fs;
    ogzstream output_PE2_fs;
    ogzstream *output_PE1_fs_ptr = NULL;
    ogzstream *output_PE2_fs_ptr = NULL;
    
    if(pair_end_read)
    {
        fastq_output_file1 = fastq_output_prefix + "_R1.fastq.gz";
        fastq_output_file2 = fastq_output_prefix + "_R2.fastq.gz";
    }
    else
    {
        fastq_output_file1 = fastq_output_prefix + ".fastq.gz";
    }


    output_PE1_fs.open(fastq_output_file1.c_str());
    if(!output_PE1_fs)
    {
        cerr << "[ERROR] fail to open output file: " << fastq_output_file1 << endl;
        exit(1);
    }
    output_PE1_fs_ptr = &output_PE1_fs;

    if(pair_end_read)
    {
        output_PE2_fs.open(fastq_output_file2.c_str());
        if(!output_PE2_fs)
        {
            cerr << "[ERROR] fail to open output file: " << fastq_output_file2 << endl;
            exit(1);
        }
        output_PE2_fs_ptr = &output_PE2_fs;
    }
    cout << "[INFO]: Processing BAM file: " << bam_input_file << endl;
    map<string, int> readname_hash;
    BamAlignment my_bam_alignment;
    while(my_bam_reader.GetNextAlignment(my_bam_alignment))
    {
        if(my_bam_alignment.IsMapped())
            readname_hash[my_bam_alignment.Name] = 0;
    }

    for(size_t i = 0; i < fastq_input_pairs.size(); i++)
    {
        print_read(readname_hash, fastq_input_pairs[i].first, fastq_input_pairs[i].second, output_PE1_fs_ptr, output_PE2_fs_ptr);
    }
    
    for(map<string, int>::iterator it = readname_hash.begin(); it != readname_hash.end(); it++)
    {
        if(it->second == 0)
            cerr << "[Warning]: could not find original raw fastq data for BAM entry: " << it->first << endl;
    }
    my_bam_reader.Close();
    output_PE1_fs.close();
    output_PE2_fs.close();
}


int main(int argc, const char * argv[]) {
    parseOption(argc, argv);
    prepare_input_fastq_list();
    convert_read();
    return 0;
}

























