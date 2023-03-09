#ifndef _SKETCH_H_
#define _SKETCH_H_
#include <iostream>
#include <vector>
#include <stdint.h>
#include <string>
#include "common.h"
#include <fstream>

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
	int id;
	int half_k;
	int half_subk;
	int drlevel;
	int genomeNumber;
} sketchInfo_t;

//for converting Kssd sketch into RabbitKSSD sketch format
typedef struct co_dirstat
{
  unsigned int shuf_id;
	bool koc;
	int kmerlen;
	int dim_rd_len;
	int comp_num;
	int infile_num;
	uint64_t all_ctx_ct;
} co_dstat_t;


bool existFile(string fileName);
bool isFastaList(string inputList);
bool isFastqList(string inputList);
bool isFastaGZList(string inputList);
bool isFastqGZList(string inputList);

bool cmpSketch(sketch_t s1, sketch_t s2);

//for result accuracy testing
bool cmpSketchName(sketch_t s1, sketch_t s2);

bool isSketchFile(string inputFile);
bool sketchFastaFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, sketchInfo_t& info, string outputFile);
bool sketchFastqFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, int leastQual, int leastNumKmer, vector<sketch_t>& sketches, sketchInfo_t& info, string outputFile);
void saveSketches(vector<sketch_t>& sketches, sketchInfo_t& info, string outputFile);
void readSketches(vector<sketch_t>& sketches, sketchInfo_t& info, string inputFile);
void transSketches(vector<sketch_t>& sketches, sketchInfo_t& info, string dictFile, string indexFile, int numThreads);
void printSketches(vector<sketch_t>& sketches, string outputFile);
void printInfos(vector<sketch_t>& sketches, string outputFile);
void convertSketch(vector<sketch_t>& sketches, sketchInfo_t& info, string inputDir, int numThreads);
void convert_from_RabbitKSSDSketch_to_KssdSketch(vector<sketch_t>& sketches, sketchInfo_t& info, string outputDir, int numThreads);
#endif
