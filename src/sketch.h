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
	vector<uint64_t> hashSet;

} sketch_t;

bool sketchFile(string inputFile, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches);
void saveSketches(vector<sketch_t> sketches, string outputFile);
#endif
