#!/usr/bin/bash


./kssd shuffle -k 11 -s 6 -l 4 -o shuf_file/L4K11.shuf
./kssd shuffle -k 10 -s 6 -l 4 -o shuf_file/L4K10.shuf
./kssd shuffle -k 9 -s 6 -l 4 -o shuf_file/L4K9.shuf
./kssd shuffle -k 8 -s 6 -l 4 -o shuf_file/L4K8.shuf

./kssd shuffle -k 11 -s 6 -l 3 -o shuf_file/L3K11.shuf
./kssd shuffle -k 10 -s 6 -l 3 -o shuf_file/L3K10.shuf
./kssd shuffle -k 9 -s 6 -l 3 -o shuf_file/L3K9.shuf
./kssd shuffle -k 8 -s 6 -l 3 -o shuf_file/L3K8.shuf

./kssd shuffle -k 10 -s 6 -l 2 -o shuf_file/L2K10.shuf
./kssd shuffle -k 9 -s 6 -l 2 -o shuf_file/L2K9.shuf
./kssd shuffle -k 8 -s 6 -l 2 -o shuf_file/L2K8.shuf
./kssd shuffle -k 7 -s 6 -l 2 -o shuf_file/L2K7.shuf

./kssd shuffle -k 9 -s 6 -l 1 -o shuf_file/L1K9.shuf
./kssd shuffle -k 8 -s 6 -l 1 -o shuf_file/L1K8.shuf
./kssd shuffle -k 7 -s 6 -l 1 -o shuf_file/L1K7.shuf
./kssd shuffle -k 6 -s 5 -l 1 -o shuf_file/L1K6.shuf

