# Introduction
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/equirep/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/equirep/badges/downloads.svg)](https://anaconda.org/bioconda/equirep)

# Introduction

EquiRep is a tool for finding tandem repeats in sequences. The data and scripts for reproducing the results of the paper are available [here](https://github.com/Shao-Group/EquiRep-test).

# Dependency

EquiRep requires C++11 or later to compile and run successfully. Make sure your compiler supports C++11 by checking its version.

# Installation

EquiRep can be easily installed via 
[conda](http://bioconda.github.io/recipes/equirep/README.html).

If you wish to install it from source code, download the source code of latest EquiRep from [here](https://github.com/Shao-Group/EquiRep/releases/download/v1.0.0/EquiRep-1.0.0.tar.gz).

Use the following commands to build EquiRep:
```
cd EquiRep
./configure
make
```
The executable file `EquiRep` will appear at `src/EquiRep`.

# Usage

EquiRep processes an input FASTA file and generates an output FASTA file with the repeat units.

The usage of EquiRep is:
```
EquiRep <input_file> <output_file_prefix>
```
Arguments:

`<input_file>` - Path to the input FASTA file.

`<output_file_prefix>` - Prefix for the output FASTA file.

# Running EquiRep on a small example

A small example of input data `input.fasta` is available in the `example` directory.

Commands to enter `example` directory and run EquiRep using `input.fasta` as input:
```
cd ./example
../src/EquiRep input.fasta output
```

An output file named `output.fasta` will appear in the `example` directory.

