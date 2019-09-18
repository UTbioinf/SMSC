## dependency

* Csh: required by MUMmer
* Blasr: please install this and make sure you can execute the command `blasr -h`
* Python 2.7

The following dependencies are included in this repo for easier installation:

* MUMmer (with bug fixed)
* Lemon
* mhap
* muscle

## Install

Run `python install.py`, this will install `smsc` to the `build/local/bin` directory inside this repository. If you want to install it to `<local>/bin/`, run `python install.py --prefix <local>`

### Known Issue

* `gcc-9.1`, or the clang linked to this version of gcc has a problem compiling `Lemon` (see Issue #2). If this is the case, please either use a lower version gcc (gcc-4.8, gcc-6, gcc-7, gcc-8.3 have been tested), or install Lemon by yourself. If you choose to install Lemon by yourself, please use the command `ls $(dirname $(which dimacs-solver))/../share/lemon/cmake/LEMONConfig.cmake` to make sure the `LEMONConfig.cmake` file exists.


# run

Go to the directory of `smsc` (default: `build/local/bin`)

Type `./smsc -h` to see the running options

# output

* `<prefix>.final.bin`: the scaffolding in binary format, including the counting for each nucleotide **A simplified parser will be added to this repository**
* `<prefix>.final.fasta`: the scaffolding in fasta format
* `<prefix>.dis`: disassembling locations
    * Each line is a triple delimited by a space *draft_id start_loc end_loc*
    * *draft_id*: the index of the draft assembly (start from 0)
    * *start_loc*: the start location of the corresponding draft assembly (start from 1)
    * *end_loc*: the end location of the corresponding draft assembly (right included)
* `<prefix>.path`: reassembling path
    * There are multiple blocks in this file. Each block is a path
    * In each block, the first line is tuple delimited by a space *path_id number_of_verticess*
    * In the following *number_of_vertices* lines, each line is a quadruple *draft_id start_loc end_loc direction*
    * The first three entry are identical to that in `<prefix>.dis`
    * *direction*: if 0, then it's we use the forward direction. If 1, we use the reverse complement.
    * Notice: *start_loc* and *end_loc* always denote the original locations in the draft assembly. So if *direction=1*, then we should extract the corresponding subsequence in the draft assembly and then compute the reverse complement.

# sample data

The following is a sample of a simulated data.

* The *E. Coli* K12 (ground truth) can be downloaded [here](http://www.ncbi.nlm.nih.gov/nuccore/NC_000913)
* The manually mutated genome is in `sample/E_coli_k12.mut.fasta`
* The error-corrected Nanopore long reads can be downloaded [here](http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/2015.07.07/ecoli_ONT_Nanocorr_Corrected_reads.fa.gz)

A simple way to run of the program would be
```
smsc -p <prefix> <draft_assembly> <long_reads>
```

The program can also run in multithreading mode my setting `-j <# of threads>`.

The default alignment tool is Nucmer. If `-E` is set, the alignment tool would be `Blasr`.

For example, suppose that the sample data `E_coli_k12.mut.fasta` and the Nanopore file `Nanopore.fasta` are in the current directory, the command 
```
smsc -j 8 -p Nanopore E_coli_k12.mut.fasta Nanopore.fasta
``` 
will use 8 threads to produce a corrected-assembly. The result is in the file `Nanopore.final.fasta`. There is also a binary file `Nanopore.final.bin` containing the counting of each nucleotide in each scaffolds of final output.


---

# Update 2.2.1

* **N.B.:** We restrict the running instances of MHAP to at most 8 on purpose. The mhap.jar keeps crashing if there are 16 instances.

# Update 2.1.0
* deal with N's in the genome sequence

# Update 2.0.1
* new path cover algorithm
* Usage: as before (change path in `run_nucmer.sh`)

---
