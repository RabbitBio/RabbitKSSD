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

struct DistInfo
{
	string refName;
	int common;
	int refSize;
	double jorc;
	double dist;
};
struct cmpDistInfo
{
	bool operator()(DistInfo d1, DistInfo d2){
		return d1.dist < d2.dist;
	}
};


void index_tridist(vector<sketch_t>& sketches, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads);

void tri_dist(vector<sketch_t>& sketches, string outputFile, int kmer_size, double maxDist, int numThreads);
void index_dist(vector<sketch_t>& ref_sketches, string refSketchOut, vector<sketch_t>& query_sketches, string outputFile, int kmer_size, double maxDist, int maxNeighbor, bool isNeighbor, int isContainment, int numThreads);
void dist(vector<sketch_t>& ref_sketches, vector<sketch_t>& query_sketches, string outputFile, int kmer_size, double maxDist, int numThreads);


#endif
