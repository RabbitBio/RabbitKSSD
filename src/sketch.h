#ifndef _SKETCH_H_
#define _SKETCH_H_
#include <iostream>
#include <vector>
#include <stdint.h>
#include <string>
#include "common.h"

using namespace std;

typedef struct fileInfo
{
	string fileName;
	uint64_t fileSize;
} fileInfo_t;

typedef struct sketch
{
	string fileName;
	string seqName;
	string comment;
	int id;
	vector<uint32_t> hashSet;

} sketch_t;

typedef struct sketchInfo
{
	int half_k;
	int half_subk;
	int drlevel;
	int genomeNumber;
} sketchInfo_t;

bool cmpSketch(sketch_t s1, sketch_t s2);

bool isSketchFile(string inputFile);
bool sketchFile(string inputFile, bool isReference, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, string outputFile);
void saveSketches(vector<sketch_t> sketches, sketchInfo_t info, string outputFile);
void readSketches(vector<sketch_t>& sketches, sketchInfo_t& info, string inputFile);
void transSketches(vector<sketch_t> sketches, sketchInfo_t info, string dictFile, string indexFile, int numThreads);
void printSketches(vector<sketch_t> sketches, string outputFile);
void printInfos(vector<sketch_t> sketches, string outputFile);
#endif
