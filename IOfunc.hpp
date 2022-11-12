#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// struct containing command line parameters and other globals
struct Args {
  std::string filename = "";
  std::string outname = "";
  int format = 0; // input file format: 0 fasta | 1 fastq
  int variant = 0; // bwt variant: 0 optBWT | 1 inputBWT | 2 inputBWT + SAP
  bool verbose = false; // verbosity level
};

// function to load a fasta file
template<typename uint_s>
void load_fasta(const char *filename, std::vector<char>& Text, std::vector<uint_s>& onset,
                uint_s& sum, uint_s& ns, bool concat){
    // open the input stream
    std::ifstream input(filename);
    if(!input.good()) { // exit if we cannot open the stream
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    // resize the input vector
    input.seekg(0, std::ios::end);
    Text.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    // check input file size
    #if M64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        if(Text.size() > pow(2,32) - 1){ 
            // the input file is too big
            std::cerr << "Error, the file size is > 4.29 GB, please use ./optsais64. exiting..." << std::endl;
            exit(-1);
        }
    #endif 
    // start reading line by line
    std::string line, DNA_sequence;
    sum = 0, ns = 0;
    while(std::getline(input, line)) {
        // skip empty lines
        if(line.empty()){ continue; }
        // header of a new sequence
        if (line[0] == '>') {
            ns++; // increase sequence count
            size_t seqlen = DNA_sequence.size();
            // insert previous DNA sequence
            if(seqlen > 0){
                // copy new sequence in Text
                memcpy(&Text[sum],&DNA_sequence[0],seqlen);
                // increase current text size
                sum += seqlen;
                // if asked compute the concatenation
                if( concat ){
                    // add separator
                    Text[sum++] = 1;
                }
                // insert new offset
                onset.push_back(sum);
            }
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else {
            // add new line
            DNA_sequence += line;
        }
    }
    // insert last sequence
    memcpy(&Text[sum],&DNA_sequence[0],DNA_sequence.size());
    sum += DNA_sequence.size();
    if( concat ){ Text[sum++] = 1; }
    onset.push_back(sum);
    // resize Text to the correct size
    Text.resize(sum);
    Text.shrink_to_fit();
    // close stream
    input.close();
}

// function to load a fastq file
template<typename uint_s>
void load_fastq(const char *filename, std::vector<char>& Text, std::vector<uint_s>& onset,
                uint_s& sum, uint_s& ns, bool concat){
    // open the input stream
    std::ifstream input(filename);
    if(!input.good()) { // exit if we cannot open the stream
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    // resize the input vector
    input.seekg(0, std::ios::end);
    Text.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    // check input file size
    #if P64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        if(Text.size() > pow(2,32) - 1){ 
            // the input file is too big
            std::cerr << "Error, the file size is > 4.29 GB, please use ./optsais64. exiting..." << std::endl;
            exit(-1);
        }
    #endif 
    // start reading line by line
    std::string line, DNA_sequence;
    sum = 0, ns = 0;
    char last = ' ';
    while(std::getline(input, line)) {
        // skip empty lines
        if(line.empty()){ continue; }
        // header of a new sequence
        if (line[0] == '@') {
            last = '@'; // last identifier seen
            ns++; // increase sequence count
            size_t seqlen = DNA_sequence.size();
            // insert previous DNA sequence
            if(seqlen > 0){
                // copy new sequence in Text
                memcpy(&Text[sum],&DNA_sequence[0],seqlen);
                // increase current text size
                sum += seqlen;
                // if asked compute the concatenation
                if( concat ){
                    // add separator
                    Text[sum++] = 1;
                }
                // insert new offset
                onset.push_back(sum);
            }
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else if (line[0] == '+'){
            // qualities line beginning
            last = '+';
        }
        else {
            // add new line
            if(last == '@') DNA_sequence += line;
        }
    }
    // insert last sequence
    memcpy(&Text[sum],&DNA_sequence[0],DNA_sequence.size());
    sum += DNA_sequence.size();
    if( concat ){ Text[sum++] = 1; }
    onset.push_back(sum);
    // resize Text to the correct size
    Text.resize(sum);
    Text.shrink_to_fit();
    // close stream
    input.close();
}
// function to load a fasta file
template<typename uint_s>
void load_fasta_conc(const char *filename, std::vector<char>& Text, uint_s& sum, uint_s& ns){
    // open the input stream
    std::ifstream input(filename);
    if(!input.good()) { // exit if we cannot open the stream
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    // resize the input vector
    input.seekg(0, std::ios::end);
    Text.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    // check input file size
    #if M64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        if(Text.size() > pow(2,32) - 1){ 
            // the input file is too big
            std::cerr << "Error, the file size is > 4.29 GB, please use ./optsais64. exiting..." << std::endl;
            exit(-1);
        }
    #endif 
    // start reading line by line
    std::string line, DNA_sequence;
    sum = 0, ns = 0;
    while(std::getline(input, line)) {
        // skip empty lines
        if(line.empty()){ continue; }
        // header of a new sequence
        if (line[0] == '>') {
            ns++; // increase sequence count
            size_t seqlen = DNA_sequence.size();
            // insert previous DNA sequence
            if(seqlen > 0){
                // copy new sequence in Text
                memcpy(&Text[sum],&DNA_sequence[0],seqlen);
                // increase current text size
                sum += seqlen;
                // add separator
                Text[sum++] = 1;
            }
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else {
            // add new line
            DNA_sequence += line;
        }
    }
    // insert last sequence
    memcpy(&Text[sum],&DNA_sequence[0],DNA_sequence.size());
    sum += DNA_sequence.size();
    Text[sum++] = 1;
    // resize Text to the correct size
    Text.resize(sum);
    Text.shrink_to_fit();
    // close stream
    input.close();
}

// function to load a fastq file
template<typename uint_s>
void load_fastq_conc(const char *filename, std::vector<char>& Text, uint_s& sum, uint_s& ns){
    // open the input stream
    std::ifstream input(filename);
    if(!input.good()) { // exit if we cannot open the stream
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    // resize the input vector
    input.seekg(0, std::ios::end);
    Text.resize(input.tellg());
    input.seekg(0, std::ios::beg);
    // check input file size
    #if M64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        if(Text.size() > pow(2,32) - 1){ 
            // the input file is too big
            std::cerr << "Error, the file size is > 4.29 GB, please use ./optsais64. exiting..." << std::endl;
            exit(-1);
        }
    #endif 
    // start reading line by line
    std::string line, DNA_sequence;
    sum = 0, ns = 0;
    char last = ' ';
    while(std::getline(input, line)) {
        // skip empty lines
        if(line.empty()){ continue; }
        // header of a new sequence
        if (line[0] == '@') {
            last = '@'; // last identifier seen
            ns++; // increase sequence count
            size_t seqlen = DNA_sequence.size();
            // insert previous DNA sequence
            if(seqlen > 0){
                // copy new sequence in Text
                memcpy(&Text[sum],&DNA_sequence[0],seqlen);
                // increase current text size
                sum += seqlen;
                // add separator
                Text[sum++] = 1;
            }
            // the current DNA sequence
            DNA_sequence.clear();
        }
        else if (line[0] == '+'){
            // qualities line beginning
            last = '+';
        }
        else {
            // add new line
            if(last == '@') DNA_sequence += line;
        }
    }
    // insert last sequence
    memcpy(&Text[sum],&DNA_sequence[0],DNA_sequence.size());
    sum += DNA_sequence.size();
    Text[sum++] = 1;
    // resize Text to the correct size
    Text.resize(sum);
    Text.shrink_to_fit();
    // close stream
    input.close();
}

// function to load file as a single text
template<typename uint_s>
void load_text(const char *filename, std::vector<char>& Text, uint_s& size){
    // open the input stream
    std::ifstream input(filename);
    if(!input.good()) { // exit if we cannot open the stream
        std::cerr << "Error opening " << filename << ". exiting..." << std::endl;
        exit(-1);
    }
    // resize the input vector
    input.seekg(0, std::ios::end);
    size = input.tellg();
    std::cout << size << std::endl;
    Text.resize(size);
    std::cout << size << std::endl;
    input.seekg(0, std::ios::beg);
    // check input file size
    #if M64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        if(Text.size() > pow(2,32) - 1){ 
            // the input file is too big
            std::cerr << "Error, the file size is > 4.29 GB, please use ./optsais64. exiting..." << std::endl;
            exit(-1);
        }
    #endif 
    // read all the file
    input.read(reinterpret_cast<char*>(&Text[0]), size); 
    // close stream
    input.close();
}