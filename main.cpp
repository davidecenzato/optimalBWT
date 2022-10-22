#include <string>
#include <iostream>
#include <getopt.h>

#include "computeTransform.h"
#include "external/malloc_count/malloc_count.h" // memory counter

// function that prints the instructions for using the tool
void print_help(char** argv) {
  std::cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << std::endl;
  std::cout << "  Options: " << std::endl
        << "\t-e \tconstruct the mdollarBWT and SAP-array, def. True " << std::endl
        << "\t-d \tconstruct the optimalBWT (internal memory), def. False " << std::endl
        << "\t-f \ttake in input a fasta file, def. True " << std::endl
        << "\t-q \ttake in input a fastq file, def. False " << std::endl 
        //<< "\t-n \ttake in input a text file, def. False " << std::endl 
        //<< "\t-s \twrite the whole (circular) suffix array, def. False " << std::endl
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
  while ((c = getopt( argc, argv, "edfqsho:n") ) != -1) {
    switch(c) {
      case 'e':
        arg.variant = 0; break;
        // compute the ebwt
      case 'd':
        arg.variant = 1; break;
        // compute the dollar ebwt
      case 'p':
        arg.variant = 2; break;
        // compute the optimal bwt
      case 'f':
        arg.format = 0; break;
        // take in input a fasta file
      case 'q':
        arg.format = 1; break;
        // take in input a fastq file
      case 'n':
        arg.format = 2; break;
        // take in input a nls file
      case 's':
        arg.format = true; break;
        // write the suffix array
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
  
  // select what to compute
  switch(arg.variant) {
    case 0:
      std::cout << "Computing the mdollarBWT and SAP-array of : " << arg.filename << std::endl;
      compute_bwt_sap_unsigned(arg);
      break;
    case 1:
      std::cout << "Computing the optimalBWT (internal memory) of : " << arg.filename << std::endl;
      compute_optbwt_unsigned(arg);
      break;
    default:
      std::cout << "error" << std::endl;
  }

  return 0;
}