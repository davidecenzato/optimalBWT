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
```

### Run on Example Data

```console
// Construct the optimal BWT of fasta file using the SAIS-based algorithm
python3 optimalBWT.py input.fasta output --algorithm sais --fasta --verbose 

// Construct the optimal BWT of fasta file using the BCR-based algorithm
python3 optimalBWT.py input.fasta output --algorithm bcr --verbose -b 10
```

# External resources

* [malloc_count](https://github.com/bingmann/malloc_count)
* [BCR_LCP_GSA](https://github.com/giovannarosone/BCR_LCP_GSA.git)

# Authors

### Theoretical results:

* Davide Cenzato
* Veronica Guerrini
* Zsuzsanna Lipták
* Giovanna Rosone

### Experimental results:

* Davide Cenzato
* Veronica Guerrini

# Reference and citation 

[1] Davide Cenzato, Veronica Guerrini, Zsuzsanna Lipták, Giovanna Rosone: Computing the optimal BWT of very large string collections. DCC 2023: 71-80 ([go to the paper](https://doi.org/10.1109/DCC55655.2023.00015))

If you use optimalBWT in an academic setting, please cite this work as follows:

### optimalBWT
    @inproceedings{CenzatoGLR23,
      author       = {Davide Cenzato and
                      Veronica Guerrini and
                      Zsuzsanna Lipt{\'{a}}k and
                      Giovanna Rosone},
      title        = {Computing the optimal {BWT} of very large string collections},
      booktitle    = {In Proc. of the 33rd Data Compression Conference, {DCC} 2023, 2023},
      pages        = {71--80},
      year         = {2023},
      doi          = {10.1109/DCC55655.2023.00015}
    }
