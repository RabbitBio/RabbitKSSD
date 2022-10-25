
RabbitKSSD is a fast and memory-effcient genome distance computation tool based on the k-mer substring space sampling (kssd).
It enables processing of large-scale datasets in practical time.



## Installation
RabbitKSSD version 1.0 can only support 64-bit Linux Systems.

### Dependancy
* cmake v.3.0 or later
* c++14
* [zlib](https://zlib.net/)

### Compile from the source code
```bash
git clone --recursive https://github.com/RabbitBio/RabbitKSSD.git
cd RabbitKSSD && mkdir build && cd build
cmake ..
make && make install
cd ..
```

## Usage
**need to update**
```bash
usage: kssd [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-i] <string> [-o] <string>
-h         : this help message
-k <int>   : set kmer size, default 21, for both clust-mst and clust-greedy
-s <int>   : set sketch size, default 1000, for both clust-mst and clust-greedy
-c <int>   : set sampling ratio to compute variable sketchSize, sketchSize = genomeSize/samplingRatio, only support with MinHash sketch function, for clust-greedy
-d <double>: set the distance threshold, default 0.05 for both clust-mst and clust-greedy
-t <int>   : set the thread number, default take full usage of platform cores number, for both clust-mst and clust-greedy
-l         : input is a file list, not a single genome file. Lines in the input file list specify paths to genome files, one per line, for both clust-mst and clust-greedy
-i <string>: path of input file. One file list or single genome file. Two input file with -f and -E option
-f         : two input files, genomeInfo and MSTInfo files for clust-mst; genomeInfo and sketchInfo files for clust-greedy 
-E         : two input files, genomeInfo and sketchInfo for clust-mst
-F <string>: set the sketch function, including MinHash and KSSD, default MinHash, for both clust-mst and clust-greedy
-o <string>: path of output file, for both clust-mst and clust-greedy
-e         : not save the intermediate file generated from the origin genome file, such as the GenomeInfo, MSTInfo, and SketchInfo files, for both clust-mst and clust-greedy

```

## Example:
```bash
#input is a file list, one genome path per line:
./kssd bacteria.list


```

## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
The paper of RabbitKSSD is under drafting.
