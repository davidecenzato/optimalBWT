#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path

Description = """
optimal BWT is a tool that computes the optimal BWT of string collections.
"""

dirname = os.path.dirname(os.path.abspath(__file__))
optsais_exe    =  os.path.join(dirname, "optsais")
optsais64_exe  =  os.path.join(dirname, "optsais64")
bcr_exe        =  os.path.join(dirname, "external/BCR_LCP_GSA/BCR_LCP_GSA")
permute_exe    =  os.path.join(dirname, "permute")

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',             help='input file name', type=str)
    parser.add_argument('output',            help='output file name', type=str)
    parser.add_argument('-a', '--algorithm', help='select which algorithm to use ((sais | bcr) def. sais)', default="sais", type=str)
    parser.add_argument('-f', '--fasta',     help='take in input a fasta file (sais only, def. True)', action='store_true')
    parser.add_argument('-q', '--fastq',     help='take in input a fastq file (sais only, def. False)', action='store_true')
    parser.add_argument('-v', '--verbose',   help='verbose (def. False)',action='store_true')
    parser.add_argument('-b', '--buffer',    help='set memory buffer size in MB (BCR only, def. 10)', default=10, type=str)
    args = parser.parse_args()

    logfile_name = args.input + ".log"
    # get main optbwt directory
    args.optbwt_dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)
    with open(logfile_name,"a") as logfile:

        # ---------- start computation #
        start = time.time()
        # ---------- SAIS based algorithm #
        if args.algorithm == "sais":
            file_size = os.path.getsize(args.input)
            command = ""
            if( file_size > 2**(32)-1 ):
                if(args.verbose): print("Running optSAIS in 64 bit mode.")
                command = "{exe} {input}".format(exe = os.path.join(args.optbwt_dir,optsais64_exe), input=args.input)
            else:
                if(args.verbose): print("Running optSAIS in 32 bit mode.")
                command = "{exe} {input}".format(exe = os.path.join(args.optbwt_dir,optsais_exe), input=args.input)
            # check input format
            if( args.fastq ):
                command += " -q"
            else:
                command += " -f"
            # check verbosity
            if( args.verbose ):
                command += " -v"
            # add output file path
            command += " -o " + args.output
            # run the command
            if(args.verbose): print(command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            # print time
            if(args.verbose): print("Elapsed time: {0:.4f}".format(time.time()-start))
        # ---------- BCR based algorithm #
        else:
            if(args.verbose): print("Running BCR based algorithm.")
            command = "{exe} {input} {output}".format(exe = bcr_exe, input = args.input, output = args.output)
            # run the command
            if(args.verbose): print(command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            # permute characters
            if(args.verbose): print("Permuting the BWT characters.")
            command = "{exe} {bwt} {sap} {buffer}".format(exe = permute_exe, bwt = args.output+".ebwt", sap = args.output+".bwt.red_sap",
                                                          buffer = args.buffer)
            if(args.verbose): print(command)
            if(execute_command(command,logfile,logfile_name)!=True):
                return
            if(args.verbose): print("Elapsed time: {0:.4f}".format(time.time()-start))
            if(args.verbose): print("Remove temporary files.")
            os.rename(args.output+".ebwt.optbwt",args.output+".optbwt")
            os.remove(args.output+".ebwt")
            os.remove(args.output+".bwt.red_sap")

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,env=None):
  try:
    #subprocess.run(command.split(),stdout=logfile,stderr=logfile,check=True,env=env)
    subprocess.check_call(command.split(),stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True
  

if __name__ == '__main__':
    main()
