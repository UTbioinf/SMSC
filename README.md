# Update 2.2.1

* **WARNING:** restrict the running instances of MHAP to at most 8. The mhap.jar keeps crashing if there are 16 instances.

Since 2.2.1, we use cmake to compile SMSC.

## dependency

* MUMmer: There is a bug in MUMmer and you need to fix it by your self. `cd <MUMmer directory>/src/tigr` and open `mgaps.cc`, in the `main()` function, change `char line [MAX_LINE];` and `char save [MAX_LINE];` to `char line [MAX_LINE] = {0};` and `char save [MAX_LINE] = {0};`, respectively. This bug occurs if the input file or input pipe is empty.
* Blasr
* Lemon
* mhap

## Install

Before installing the program, you need to modify the `config.sh` by following steps:

* `cd smsc/src.2.2.1`
* open `config.sh`, and modify the configurations to the absolute directory where you install MUMmer 

Then go to the `smsc` folder, do the following things:

```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PWD/local -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

The `smsc` will then be installed in `local/bin`

**Notice:** the `make install` is required so that the script `run_nucmer.sh` will be accessed correctly by `smsc`. To change the install location, simply change `$PWD/local` to the location you prefer.

# run

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

# Update 2.1.0
* deal with N's in the genome sequence

# Update 2.0.1
* new path cover algorithm
* Usage: as before (change path in `run_nucmer.sh`)

---
