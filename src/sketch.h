#ifndef _SKETCH_H_
#define _SKETCH_H_
#include <iostream>
#include <vector>
#include <stdint.h>
#include <string>
#include "common.h"

using namespace std;
typedef struct sketch
{
	string fileName;
	string seqName;
	string comment;
	int id;
	vector<uint32_t> hashSet;

} sketch_t;

bool isSketchFile(string inputFile);
bool sketchFile(string inputFile, bool isReference, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, string outputFile);
void saveSketches(vector<sketch_t> sketches, string outputFile);
void readSketches(vector<sketch_t>& sketches, string inputFile);
void printSketches(vector<sketch_t> sketches, string outputFile);
void printInfos(vector<sketch_t> sketches, string outputFile);
#endif
