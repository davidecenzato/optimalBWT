# optimalBWT
optimalBWT is a tool that computes the optimal BWT using either a SAIS- or BCR-based algorithm.

# Usage

```
usage: optimalBWT.py [-h] [-a ALGORITHM] [-f] [-q] [-v] [-b BUFFER] input output

optimal BWT is a tool that computes the optimal BWT of string collections.

positional arguments:
  input                 input file name
  output                output file name

options:
  -h, --help            show this help message and exit
  -a ALGORITHM, --algorithm ALGORITHM
                        select which algorithm to use ((sais | bcr) def. sais)
  -f, --fasta           take in input a fasta file (sais only, def. True)
  -q, --fastq           take in input a fastq file (sais only, def. False)
  -v, --verbose         verbose (def. False)
  -b BUFFER, --buffer BUFFER
                        set memory buffer size in MB (BCR only, def. 10)
```
When using the SAIS-based algorithm you can choose your input format between fasta and fastq format. When using the BCR-based algorithm you can choose the size of the input buffer for the algorithm permuting the characters in the SAP-intervals.

### Requirements

The optimalBWT tool requires
* A modern Python 3 release version 3.7 or higher.
* A modern C++11 compiler such as `g++` version 4.9 or higher.

# Example

### Download and Compile

```console
git clone https://github.com/davidecenzato/optimalBWT.git
cd optimalBWT
git submodule update --init --recursive

make
make install_bcr
```

### Run on Example Data

```console
// Construct the optimal BWT of fasta file using the SAIS-based algorithm
python3 input.fasta output --algorithm sais --fasta --verbose 

// Construct the optimal BWT of fasta file using the BCR-based algorithm
python3 input.fasta output --algorithm bcr --verbose -b 10
```

# External resources

* [malloc_count](https://github.com/bingmann/malloc_count)
* [BCR_LCP_GSA](https://github.com/giovannarosone/BCR_LCP_GSA.git)

# Citation 

If you use this tool in an academic setting, please cite this work as follows:

### optimalBWT
    @unpublished  {optBWT,
      author    = {Davide Cenzato and
                   Veronica Guerrini and
                   {\relax Zs}uzsanna Lipt{\'{a}}k and
                   Giovanna Rosone},
      title     = {Computing the optimal BWT of very large string collections},
      booktitle = {Submitted for revision},
      year      = {2022}
    }

# Authors

### Theoretical results:

* Davide Cenzato
* Veronica Guerrini
* Zsuzsanna Lipták
* Giovanna Rosone

### Experimental results:

* Davide Cenzato
* Veronica Guerrini
