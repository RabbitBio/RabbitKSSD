#ifndef _SUB_COMMAND_H_
#define _SUB_COMMAND_H_

#include <string>
#include "common.h"

using namespace std;

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads);
void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads);
void command_merge(string sketchFile, string outputFile, int threadNumber);
void command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads);
void command_sketch(string refList, string outputFile, kssd_parameter_t kssd_parameter, int threads);
void command_info(string sketchFile, string outputFile);

#endif
