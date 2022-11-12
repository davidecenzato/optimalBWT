# optimalBWT
optimalBWT is a tool that computes the optimal BWT using the optSAIS and BCR algorithms.

# Usage

```
Usage: ./optsais <input filename> [options]
  Options:
        -p      construct the optimalBWT, def. True
        -d      construct the BWT (input order) and SAP-array, def. False
        -e      construct the BWT (input order), def. True
        -f      take in input a fasta file, def. True
        -q      take in input a fastq file, def. False
        -v      set verbose mode, def. False
        -o O    basename for the output files, def. <input filename>
```
When computing the optimalBWT you can choose your input format between fasta and fastq format, and set a path for the output.

### Requirements

The optimalBWT tool requires:
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
// Construct the optimal BWT of fasta file using the optSAIS algorithm
./optsais file.fasta -p -f -v

// Compute the input order BWT and the SAP-array, then run Bentley et al. algorithm
./ optsais file.fasta -d -f -v
./ permute file.fasta.inputbwt file.fasta.sap 10
```

# External resources

* [malloc_count](https://github.com/bingmann/malloc_count)

# Citation 

If you use this tool in an academic setting, please cite this work as follows:

### optimalBWT
    @unpublished  {optBWT,
      author    = {Davide Cenzato and
                   Veronica Guerrini and
                   {\relax Zs}uzsanna Lipt{\'{a}}k and
                   Giovanna Rosone},
      title     = {Computing the optimal BWT of very large string collections},
      booktitle = {Submitted for revision at DCC 2023},
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

### Implementation:

* [Davide Cenzato](https://github.com/davidecenzato) 
