#include "subCommand.h"
#include "sketch.h"
#include "dist.h"
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_set>

bool cmp(uint64_t a, uint64_t b){
	return a < b;
}

bool cmpSketch(sketch_t s1, sketch_t s2){
	return s1.id < s2.id;
}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	vector<sketch_t> sketches;
	bool success = sketchFile(refList, threads, kssd_parameter, sketches);
	cerr << "finish the sketch generation " << success << endl;
	cerr << "the size of sketches is: " << sketches.size() << endl;
	double t3 = get_sec();
	cerr << "===================time of generator sketches is: " << t3 - t2 << endl;

	std::sort(sketches.begin(), sketches.end(), cmpSketch);

	double t4 = get_sec();
	cerr << "===================time of sort sketches order is: " << t4 - t3 << endl;

	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < sketches.size(); i++)
	{
		std::sort(sketches[i].hashSet.begin(), sketches[i].hashSet.end(), cmp);
	}
	double t5 = get_sec();
	cerr << "===================time of sort each sketches is: " << t5 - t4 << endl;


	string refHashOut = refList + ".hash";
	saveSketches(sketches, refHashOut);

	double t6 = get_sec();
	cerr << "===================time of save sketches into hash.out is: " << t6 - t5 << endl;

	//compute the pair distance
	
	cerr << "start the distance computing" << endl;
	//string dist_output = "result.out";

	tri_dist(sketches, outputFile, kmer_size, maxDist, threads);

	//FILE * fp0 = fopen(dist_output.c_str(), "w");
	double t7 = get_sec();
	cerr << "===================time of get total distance matrix file is: " << t7 - t6 << endl;
}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	vector<sketch_t> ref_sketches;
	bool success0 = sketchFile(refList, threads, kssd_parameter, ref_sketches);
	cerr << "finish the ref_sketches generation " << success0 << endl;
	cerr << "the size of ref_sketches is: " << ref_sketches.size() << endl;
	double t3 = get_sec();
	cerr << "===================time of generator reference sketches is: " << t3 - t2 << endl;

	vector<sketch_t> query_sketches;
	bool success1 = sketchFile(queryList, threads, kssd_parameter, query_sketches);
	double t3_1 = get_sec();
	cerr << "===================time of generator query sketches is: " << t3_1 - t3 << endl;

	std::sort(ref_sketches.begin(), ref_sketches.end(), cmpSketch);
	std::sort(query_sketches.begin(), query_sketches.end(), cmpSketch);

	double t4 = get_sec();
	cerr << "===================time of sort sketches order is: " << t4 - t3_1 << endl;

	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < ref_sketches.size(); i++)
	{
		std::sort(ref_sketches[i].hashSet.begin(), ref_sketches[i].hashSet.end(), cmp);
	}
	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i = 0; i < query_sketches.size(); i++)
	{
		std::sort(query_sketches[i].hashSet.begin(), query_sketches[i].hashSet.end(), cmp);
	}

	double t5 = get_sec();
	cerr << "===================time of sort each sketches is: " << t5 - t4 << endl;

	string refHashOut = refList + ".hash";
	string queryHashOut = queryList + ".hash";
	saveSketches(ref_sketches, refHashOut);
	saveSketches(query_sketches, queryHashOut);

	double t6 = get_sec();
	cerr << "===================time of save sketches into hash.out is: " << t6 - t5 << endl;

	//compute the pair distance
	
	cerr << "start the distance computing" << endl;
	//string dist_output = "result.out";

	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);

	//FILE * fp0 = fopen(dist_output.c_str(), "w");
	double t7 = get_sec();
	cerr << "===================time of get total distance matrix file is: " << t7 - t6 << endl;
}

void command_merge(string sketchFile, string outputFile, int threadNumber){
	
	vector<sketch_t> sketches;
	readSketches(sketches, sketchFile);
	//for(int i = 0; i < sketches.size(); i++){
	//	cout << sketches[i].fileName << endl;
	//	//for(int j = 0; j < sketches[i].hashSet.size(); j++){
	//	//	cout << sketches[i].hashSet[j] << '\t';
	//	//}
	//	//cout << endl;
	//}
	//exit(0);
	string totalName("");
	unordered_set<uint64_t> mergedSet;
	for(int i = 0; i < sketches.size(); i++){
		totalName += sketches[i].fileName + '\n';
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			uint64_t curHash = sketches[i].hashSet[j];
			mergedSet.insert(curHash);
		}
	}
	totalName = totalName.substr(0, totalName.length()-1);
	//cout << totalName << endl;
	cerr << "the size of merged hash set is: " << mergedSet.size() << endl;
	vector<uint64_t> mergedArr;
	for(auto x : mergedSet){
		mergedArr.push_back(x);
	}

	vector<sketch_t> mergedSketches;
	sketch_t s;
	s.id = 0;
	s.fileName = totalName;
	s.hashSet = mergedArr;
	mergedSketches.push_back(s);
	saveSketches(mergedSketches, outputFile);
}
	









