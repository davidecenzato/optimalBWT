#include <string>
#include <iostream>
#include <getopt.h>

/*
Description:
Tool for computing the optimal BWT (optBWT), the SAP-array, and the input order BWT
using the optimal SAIS algorithm (optSAIS).
   "Computing the optimal BWT of very large string collections"
   by Davide Cenzato, Veronica Guerrini, Zsuzsanna Lipták, and Giovanna Rosone
   submitted paper 
   
The input file cannot contain the characters < 2 which are used internally by the algorithm.
Input file smaller than 4.29 GB will take 5n + n/8 bytes in RAM, while input files larger than
4.29 GB will take 9n + n/8 bytes.
*/

#include "computeTransform.h"
#include "external/malloc_count/malloc_count.h" // memory counter

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-p \tconstruct the optimalBWT, def. True " << std::endl
        << "\t-d \tconstruct the BWT (input order) and SAP-array, def. False " << std::endl
        << "\t-e \tconstruct the BWT (input order), def. True " << std::endl
        << "\t-f \ttake in input a fasta file, def. True " << std::endl
        << "\t-q \ttake in input a fastq file, def. False " << std::endl 
        << "\t-v \tset verbose mode, def. False " << std::endl
        << "\t-o O\tbasename for the output files, def. <input filename>" << std::endl;

  exit(-1);
}
// function for parsing the input arguments
void parseArgs( int argc, char** argv, Args& arg ) {
  int c;
  extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");
 
  std::string sarg;
  while ((c = getopt( argc, argv, "edpfqho:v") ) != -1) {
    switch(c) {
      case 'p':
        arg.variant = 0; break;
        // compute the ebwt
      case 'e':
        arg.variant = 1; break;
        // compute the dollar ebwt
      case 'd':
        arg.variant = 2; break;
        // compute the optimal bwt
      case 'f':
        arg.format = 0; break;
        // take in input a fasta file
      case 'q':
        arg.format = 1; break;
        // take in input a fastq file
      case 'v':
        arg.verbose = true; break;
        // verbose mode
      case 'o':
        sarg.assign( optarg );
        arg.outname.assign( sarg ); break;
        // store the output files path
      case 'h':
        print_help(argv); exit(-1);
        // fall through
      default:
        std::cout << "Unknown option. Use -h for help." << std::endl;
        exit(-1);
    }
  }
  // the only input parameter is the file name
  if (argc == optind+1) {
    arg.filename.assign( argv[optind] );
  }
  else {
    std::cout << "Invalid number of arguments" << std::endl;
    print_help(argv);
  }
  // set output files basename
  if(arg.outname == "") arg.outname = arg.filename;
}

int main(int argc, char** argv)
{
  // translate command line arguments
  Args arg;
  parseArgs(argc, argv, arg);
  
  // select the BWT variant to compute
  switch(arg.variant) {
    case 0:
      if(arg.verbose) std::cout << "Computing the optimal BWT of : " << arg.filename << std::endl; 
      compute_optbwt(arg);
      break;
    case 1:
      if(arg.verbose) std::cout << "Computing the BWT (input order) : " << arg.filename << std::endl; 
      compute_bwt(arg);
      break;
    case 2:
      if(arg.verbose) std::cout << "Computing the BWT (input order) and SAP-array of : " << arg.filename << std::endl;
      compute_bwt_sap(arg);
      break;
    default:
      std::cout << "Error! select a valid BWT variant. exiting..." << std::endl;
      exit(1);
  }

  return 0;
}