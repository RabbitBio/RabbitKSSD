

RabbitKSSD is a high-speed genome distance estimation tool with a series of function mudules for genome alalysis, includeing:
* `sketch`, generate sketch files.
* `dist`, compute pairwise distances between reference and query genomes.
* `set operation`, operations on sketches, including subtraction and union.


## Installation
RabbitKSSD version 1.0.0 can only support 64-bit Linux Systems.

### Dependancy
* cmake v.3.0 or later
* c++14
* [zlib](https://zlib.net/)
* [RabbitFX](https://github.com/RabbitBio/RabbitFX)

### Install
* Pay attention to the `--recursive` option for submodule `RabbitFX`.
* The `install.sh` is used for compiling and installing RabbitKSSD.
* The `init_shuffle.sh` is used for generating the shuffle files and save in the `shuf_file/` directory.
```bash
git clone --recursive https://github.com/RabbitBio/RabbitKSSD.git
cd RabbitKSSD 
# Compile and install RabbitKSSD
./install.sh
# Init generating the shuffle files
./init_shuffle.sh
```

### Compile from the source code
```bash
git clone --recursive https://github.com/RabbitBio/RabbitKSSD.git
cd RabbitKSSD

# Compile rabbitFX library
cd RabbitFX &&
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install && 
cd ../../ &&

# Compile RabbitKSSD
mkdir -p build
cd build
cmake ..
make -j8 && make install
cd ..

# Init generating the shuffle files
./init_shuffle.sh
```

## Usage
RabbitKSSD has several sub-commands, including `shuffle`, `sketch`, `alldist`, `dist`, `union`, `sub`, `convert`, `merge`, and `info`.\
All sub-commands support the `-h` or `--help` option for help information.\
The sub-command `convert` is used for converting the sketch format between the Kssd and RabbitKSSD formats.\
The sub-command `merge` is used for merging multiple sketch files into one sketch file.\
The sub-command `info` is used for formatting the sketch files in a human-readable format.\
Other sub-commands are introduced in detail in the following sub-sections.
```bash
# For help informations:
./rabbit_kssd -h
# It will print:
rabbit_kssd: accelerating Kssd-based genome distance estimation on modern multi-core architectures
Usage: ./rabbit_kssd [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit

Subcommands:
  shuffle                     generate the shuffle file for sketching usage
  sketch                      compute sketches for the input genome list
  alldist                     compute all-vs-all distances for one input dataset
  dist                        compute the all-vs-all distances between reference genomes and query datasets
  union                       compute the set union from multiple sketches
  sub                         subtract the reference sketch from the query sketches
  convert                     convert the sketches between Kssd format and RabbitKSSD format
  merge                       merge multiple sketch files into one single sketch file
  info                        get the information of the sketch file
```

### shuffle
**The sub-command `shuffle` is used for generating the shuffle files to create sketches.**

`-h`: For help information.\
`-k`: Half-length of k-mer, `-k x` means using k-mer of length `2x`.\
`-s`: Half-length -f k-mer substring, `-s x` means the whole space is the collection of all `2x-mer`.\
`-l`: Dimension reduction level. `-l x` means the expected rate of dimensionality reduction is $16^x$.\
Users can use `init_shuffle.sh` to generate the shuffle files automatically.

Recommend: 
* `-k 8` for bacteria; `-k 10` for mammals or metagenomics; `-k 9` for other genome size in-between.
* `-s 6` is the default setting, usually no need to change this setting; `s < k`.
* `-l 3` for bacteria; `-l 4` or `-l 5` for mammals; `l < s`.

```bash
# For help information
./rabbit_kssd shuffle -h
# Example:
./rabbit_kssd shuffle -k 10 -s 6 -l 3 -o shuf_file/L3K10.shuf
```


### sketch 
**The sub-command `sketch` is used for generating sketches from the input genome list.**

`-h`: For help information.\
`-L L3K10.shuf`: Load the shuffle file `L3K10.shuf` (k-mer length is 20, dimension reduction level is 3, which means the expected sketch size is $1/16^3$ of the genome size)\
`-i bacteria.list`: Input genome list file `bacteria.list`, one genome path per line.\
`-o bacteria.sketch`: Output sketch file `bacteria.sketch`. The hash dictionary file `bacteria.sketch.dict` and hash offset file `bacteria.sketch.index` are also generated for distance computation (`alldist` and `dist` sub-commands).\
`-t threads`: Set the thread number, default all CPUs of the platform.

```bash
# For help information
./rabbit_kssd sketch -h
# Example: 
./rabbit_kssd sketch -L shuf_file/L3K10.shuf -i bacteria.list -o bacteria.sketch
```

### alldist
**The sub-command `alldist` is used for computing the all-vs-all distances for a single genome dataset.**

`-h`: For help information.\
`-i bacteria.list`: Input genome list file `bacteria.list`, one genome path per line.\
`-i bacteria.sketch`: Input sketch file `bacteria.sketch`, which is generated by `sketch` sub-command.\
`-L L3K10.shuf`: Load the shffle file `L3K10.shuf` for generating sketches.\
`-D 0.05`: Set the maximum output distance threshold, distances over `0.05` are omitted.\
`-o bacteria.alldist`: Set the output distance file `bacteria.alldist`. When the file is over `4GB`, the results are stored to sub-files in `bacteria.alldist.dir/` directory, and an appending index file `bacteria.alldist.index` is generated for retrieving the result.\
`-t threads`: Set the thread number, default all CPUs of the platform.
```bash
# For help information:
./rabbit_kssd alldist -h
```
```bash
# Example: computing the all-vs-all distances of bacteria, input as a genome list
./rabbit_kssd alldist -i bacteria.list -L shuf_file/L3K10.shuf -D 0.05 -o bacteria.alldist
```
is equal to:

```bash
# Generating the sketch file:
./rabbit_kssd sketch -L shuf_file/L3K10.shuf -i bacteria.list -o bacteria.sketch
# Computing the all-vs-all distances:
./rabbit_kssd alldist -i bacteria.sketch -D 0.05 -o bacteria.alldist
```

### dist
**The sub-command `dist` is used for computing the all-vs-all distances between the reference and query genome datasets.**

`-h`: For help information.\
`-r ref.list`: Input reference genome list file `ref.list`, one genome per line.\
`-r ref.sketch`: Input reference sketch file `ref.sketch`, which is generated by `sketch` sub-command.\
`-q query.list`: Input query genome list file `query.list`, one genome per line.\
`-q query.sketch`: Input query sketch file `query.sketch`, which is generated by `sketch` sub-command.\
`-L L3K10.shuf`: Load the shuffle file `L3K10.shuf` for generating sketches.\
`-D 0.05`: Set the maximum output distance threshold, distances over `0.05` are omitted.\
`-o ref_query.dist`: Set the output distance file `ref_query.dist`. When the file is over `4GB`, the results are stored to sub-files in `ref_query.dist.dir/` directory, and an appending index file `ref_query.dist.index` is generated for retrieving the result.\
`-N 3`: Set the maximum number of nearest neighbor reference output for each query genome.\
`-t threads`: Set the thread number, default all CPUs of the platform.

```bash
# For help information:
./rabbit_kssd dist -h
```
Compute the all-vs-all distances between reference and query genome datasets.
```bash
# Example: computing the all-vs-all distances between reference and query genome datasets, input as genome list
./rabbit_kssd dist -r ref.list -q query.list -L shuf_file/L3K10.shuf -D 0.05 -o ref_query.dist
```
is equal to:
```bash
# Generating the sketch file:
./rabbit_kssd sketch -L shuf_file/L3K10.shuf -i ref.list -o ref.sketch
./rabbit_kssd sketch -L shuf_file/L3K10.shuf -i query.list -o query.sketch
# Computing the all-vs-all distances:
./rabbit_kssd dist -r ref.sketch -q query.sketch -D 0.05 -o ref_query.dist
```
Find the nearest neighbor reference for query dataset.
```bash
# Example: computing the nearest 3 reference genomes for each query genome
# Input as genome list:
./rabbit_kssd dist -r ref.list -q query.list -L shuf_file/L3K10.shuf -N 3 -o ref_query_nearest_3.dist

# Input as sketch file:
./rabbit_kssd dist -r ref.sktch -q query.sketch -N 3 -o ref_query_nearest_3.dist
```

### union
**The sub-command `union` is used for computing the set union from multiple sketches.** 

`-h`: For help information.\
`-i ref.sketch`: Input sketch file `ref.sketch`, which includes repeated hash values from multiple sketches.\
`-o ref_union.sketch`: Output sketch file `ref_union.sketch`, which includes a union of hash value set.

```bash
# For help information:
./rabbit_kssd union -h
# Example: computing the union of sketch set for ref.sketch
./rabbit_kssd union -i ref.sketch -o ref_union.sketch
```


### sub 
**The sub-command `sub` is used for subtracting the reference sketch hash set from each query sketch.**

`-h`: For help information.\
`--rs human.sketch`: Set the reference sketch file `human.sketch` to be subtracted from the query sketches.\
`--qs 1000Genomes.sketch`: Set the query sketch file `1000Genomes.sketch`.\
`-o`: Output sketch file `1000Genomes_sub.sketch`, which consists of query sketches after subtraction.\
`-t threads`: Set the thread number, default all CPUs of the platform.

```bash
# For help information:
./rabbit_kssd sub -h
# Example:
./rabbit_kssd sub --rs human.sketch --qs 1000Genomes.sketch -o 1000Genomes_sub.sketch
```


## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
Xiaoming Xu, Zekun Yin, Lifeng Yan, Huiguang Yi, Hua Wang, Bertil Schmidt, Weiguo Liu, RabbitKSSD: accelerating genome distance estimation on modern multi-core architectures, Bioinformatics, 2023;, btad695, https://doi.org/10.1093/bioinformatics/btad695
