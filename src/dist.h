#ifndef _DIST_H_
#define _DIST_H_

#include <iostream>
#include <vector>
#include "sketch.h"
#include <stdint.h>

using namespace std;

typedef struct setResult
{
	int common;
	int size0;
	int size1;
	double jaccard;
} setResult_t;

setResult_t getJaccard(vector<uint32_t> list1, vector<uint32_t> list2);


void tri_dist(vector<sketch_t> sketches, string outputFile, int kmer_size, double maxDist, int numThreads);
void dist(vector<sketch_t> ref_sketches, vector<sketch_t> query_sketches, string outputFile, int kmer_size, double maxDist, int numThreads);




#endif
