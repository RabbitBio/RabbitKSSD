#ifndef _SUB_COMMAND_H_
#define _SUB_COMMAND_H_

#include <string>
#include "common.h"

using namespace std;

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, double maxDist, int isContainment, int threads);
void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, double maxDist, int maxNeighbor, bool isNeighbor, int isContainment, int threads);
void command_union(string sketchFile, string outputFile, int threadNumber);
void command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads);
void command_sketch(string refList, bool isReference, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, int threads);
void command_info(string sketchFile, bool isDetail, string outputFile);
void command_convert(string inputDir, bool to_Kssd_sketch, bool isQuery, string outputFile, int threads);
void command_merge(string inputList, string outputFile, int threads);

void new_command_union(string sketchFile, string outputFile, int threadNumber);
void new_command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads);
#endif
