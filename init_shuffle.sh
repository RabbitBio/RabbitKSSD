#!/usr/bin/bash

mkdir -p shuf_file

./rabbit_kssd shuffle -k 11 -s 6 -l 4 -o shuf_file/L4K11.shuf
./rabbit_kssd shuffle -k 10 -s 6 -l 4 -o shuf_file/L4K10.shuf
./rabbit_kssd shuffle -k 9 -s 6 -l 4 -o shuf_file/L4K9.shuf
./rabbit_kssd shuffle -k 8 -s 6 -l 4 -o shuf_file/L4K8.shuf

./rabbit_kssd shuffle -k 11 -s 6 -l 3 -o shuf_file/L3K11.shuf
./rabbit_kssd shuffle -k 10 -s 6 -l 3 -o shuf_file/L3K10.shuf
./rabbit_kssd shuffle -k 9 -s 6 -l 3 -o shuf_file/L3K9.shuf
./rabbit_kssd shuffle -k 8 -s 6 -l 3 -o shuf_file/L3K8.shuf

./rabbit_kssd shuffle -k 10 -s 6 -l 2 -o shuf_file/L2K10.shuf
./rabbit_kssd shuffle -k 9 -s 6 -l 2 -o shuf_file/L2K9.shuf
./rabbit_kssd shuffle -k 8 -s 6 -l 2 -o shuf_file/L2K8.shuf
./rabbit_kssd shuffle -k 7 -s 6 -l 2 -o shuf_file/L2K7.shuf


