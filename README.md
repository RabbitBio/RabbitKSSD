
RabbitKSSD is a fast and memory-effcient genome distance computation tool based on the k-mer substring space sampling (kssd).
It enables processing of large-scale datasets in practical time.



## Installation
RabbitKSSD version 1.0 can only support 64-bit Linux Systems.

### Dependancy
* cmake v.3.0 or later
* c++14
* [zlib](https://zlib.net/)

### Compile from the source code automatically
```bash
git clone --recursive https://github.com/RabbitBio/RabbitKSSD.git
cd RabbitKSSD 
./install.sh
```

### Compile from the source code
```bash
git clone --recursive https://github.com/RabbitBio/RabbitKSSD.git
cd RabbitKSSD

#make rabbitFX library
cd RabbitFX &&
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install && 
cd ../../ &&

#compile RabbitKSSD
mkdir -p build
cd build
cmake ..
make -j8 && make install
cd ..
```

## Usage
RabbitKSSD has several sub-command including sketch, alldist, dist, merge, sub, and info.

### sketch
The subcommand ***sketch*** is used for generation sketches for the input genome list.
```bash
# example: (dafault use the total number of threads, set manually by -t options)
# input as genome list, one genome per line
./kssd sketch -L shuf_file/L3K10.shuf -i bacteria.list -o bacteria.sketch
```

### alldist
The subcommand ***alldist*** is used for computing the all-to-all distance for an input dataset.
```bash
# examples:(-m set the maximum distance threshold, distance over the threshold will not be outputed)
# input as a sketch (suggested)
./kssd alldist -i bacteria.sketch -m 0.05 -o bacteria.alldist

# input as the genome list
./kssd alldist -i bacteria.list -m 0.05 -o bacteria.alldist
```

### dist
The subcommand ***dist*** is used for computing the distance between all pairwise reference and query genomes.
```bash
# examples: (-m set the maximum distance threshold, distance over the threshold will not be outputed)
# input as sketches (suggested)
./kssd dist -r reference.sketch -q query.sketch -m 0.05 -o reference_query.dist

# input as genome lists
./kssd dist -r reference.list -q reference.list -m 0.05 -o reference_query.dist
```

### info
The subcommand ***info*** is used for parsing the sketch contents into human-readable format.
```bash
# example: 
./kssd info -i bacteria.sketch -o bacteria.info
```
### merge
The subcommand ***merge*** is used for getting the union of the hash sets in multi-sketches.
```bash
# example:
./kssd merge -i bacteria.sketch -o bacteria.merge.sketch
```
### sub 
The subcommand ***sub*** is used for the set substraction.
```bash
# example:
./kssd sub --rs bacteria.merge.sketch --qs million.sketch -o million_sub_bact.sketch
```


## Example:

## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
The paper of RabbitKSSD is under drafting.

