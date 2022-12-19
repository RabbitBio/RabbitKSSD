#!/usr/bin/bash 


if [ $# -lt 1 ]; then
  echo "Usage: $0 SRA_accession.list"
	echo "The SRA_accession.list is the list files of SRA accession id, one id per line"
	echo ""
  exit
fi

accessionList=$1
mkdir -p tmpSketchDir
rm tmpSketchDir/*

cat $accessionList | while read line
do
	fastq-dump $line
	ls ${line}.fastq > ${line}.list
	./kssd sketch -L shuf_file/L3K10.shuf -i ${line}.list -o tmpSketchDir/${line}.sketch -q
	rm ${line}.fastq ${line}.list
done

ls tmpSketchDir/*.sketch >tmpSketch.list
./kssd merge -i tmpSketch.list -o ${accessionList}.sketch

rm tmpSketch.list
rm -rf tmpSketchDir



