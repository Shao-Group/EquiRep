# Introduction

EquiRep is a tool for finding tandem repeats in sequences.

# Installation

Clone the git repository of EquiRep using the command:

```
git clone https://github.com/Shao-Group/EquiRep.git
```

Or download the source code of latest EquiRep from [here](https://github.com/Shao-Group/EquiRep/releases/download/v1.0.0/EquiRep-1.0.0.tar.gz).

Use the following commands to build EquiRep:
```
cd EquiRep/src
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

# Data

All simulated and real data used in the experiments are available in the `data` directory.

`simulated_data`: `dat3` contains the randomly simulated data of different unit lengths and copy numbers. `dat_aax2` contains simulated sequences with 2 recurring kmers. `dat_aax3` contains simulated sequences with 3 recurring kmers. `error_10` and `error_20` inside each of these folders represent error rates of 10% and 20% respectively.

`HOR_data`: `hor_repeats.fasta` contains the 13 Higher Order Repeat sequences from human chromosome 5. `hor_combined.fasta` contains the concatenated sequences (x), concatenated sequences with flanking regions (axa), and error rate of 1%, 5%, 10% applied to them (x_err1, x_err5, x_err10, axa_err1, axa_err5, axa_err10).

`RCA_data`: `RCA_101.fasta` contains selected 101 Nanopore long read sequences from human prostate tissue (GEO, accession number: GSE141693).