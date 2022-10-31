#include "subCommand.h"
#include "sketch.h"
#include "dist.h"
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_set>


void command_sketch(string refList, string outputFile, kssd_parameter_t kssd_parameter, int threads){
	vector<sketch_t> sketches;
	bool success = sketchFile(refList, threads, kssd_parameter, sketches);
	if(!isSketchFile(outputFile)){
		outputFile = outputFile + ".sketch";
	}
	saveSketches(sketches, outputFile);
	cerr << "save the sketches into: " << outputFile << endl;
}

void command_info(string sketchFile, string outputFile){
	vector<sketch_t> sketches;
	readSketches(sketches, sketchFile);
	cerr << "the number of genome is: " << sketches.size() << endl;
	printSketches(sketches, outputFile);


}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	vector<sketch_t> sketches;
	if(isSketchFile(refList)){
		readSketches(sketches, refList);
	}
	else{
		bool success = sketchFile(refList, threads, kssd_parameter, sketches);
		string refSketchOut = refList + ".sketch";
		saveSketches(sketches, refSketchOut);
		cerr << "finish the sketch generation " << success << endl;
	}
	cerr << "the size of sketches is: " << sketches.size() << endl;

	double t3 = get_sec();
	if(isSketchFile(refList))
		cerr << "===================time of read sketches from file is " << t3 - t2 << endl;
	else
		cerr << "===================time of computing sketches and save sketches into file is " << t3 - t2 << endl;

	//compute the pair distance
	
	tri_dist(sketches, outputFile, kmer_size, maxDist, threads);

	//FILE * fp0 = fopen(dist_output.c_str(), "w");
	double t4 = get_sec();
	cerr << "===================time of get total distance matrix file is: " << t4 - t3 << endl;
}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	vector<sketch_t> ref_sketches;
	if(isSketchFile(refList)){
		readSketches(ref_sketches, refList);
	}
	else{
		bool success0 = sketchFile(refList, threads, kssd_parameter, ref_sketches);
		string refHashOut = refList + ".sketch";
		saveSketches(ref_sketches, refHashOut);
		cerr << "finish the ref_sketches generation " << success0 << endl;
	}
	cerr << "the size of ref_sketches is: " << ref_sketches.size() << endl;
	double t3 = get_sec();
	if(isSketchFile(refList))
		cerr << "===================time of read reference sketches from file is: " << t3 - t2 << endl;
	else 
		cerr << "===================time of computing reference sketches and save sketches into file is: " << t3 - t2 << endl;

	vector<sketch_t> query_sketches;
	if(isSketchFile(queryList)){
		readSketches(query_sketches, queryList);
	}
	else{
		bool success1 = sketchFile(queryList, threads, kssd_parameter, query_sketches);
		string queryHashOut = queryList + ".sketch";
		saveSketches(query_sketches, queryHashOut);
		cerr << "finish the query_sketches generation " << success1 << endl;
	}
	cerr << "the size of query_sketches is: " << query_sketches.size() << endl;

	double t4 = get_sec();
	if(isSketchFile(queryList))
		cerr << "===================time of read query sketches from file is: " << t4 - t3 << endl;
	else
		cerr << "===================time of computing query sketches and  save sketches into file is: " << t4 - t3 << endl;


	//compute the pair distance
	
	//string dist_output = "result.out";

	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);

	//FILE * fp0 = fopen(dist_output.c_str(), "w");
	double t5 = get_sec();
	cerr << "===================time of get total distance matrix file is: " << t5 - t4 << endl;
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
	unordered_set<uint32_t> mergedSet;
	for(int i = 0; i < sketches.size(); i++){
		totalName += sketches[i].fileName + '\n';
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			uint32_t curHash = sketches[i].hashSet[j];
			mergedSet.insert(curHash);
		}
	}
	totalName = totalName.substr(0, totalName.length()-1);
	//cout << totalName << endl;
	cerr << "the size of merged hash set is: " << mergedSet.size() << endl;
	vector<uint32_t> mergedArr;
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

void command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads){
	vector<sketch_t> refSketches;
	readSketches(refSketches, refSketchFile);
	vector<sketch_t> querySketches;
	readSketches(querySketches, querySketchFile);
	cerr << "the refSketches size is: " << refSketches.size() << endl;
	cerr << "the querySketches size is: " << querySketches.size() << endl;

	unordered_set<uint32_t> refHashSet;
	for(int i = 0; i < refSketches.size(); i++){
		for(auto x : refSketches[i].hashSet){
			refHashSet.insert(x);
		}
	}
	cerr << "the size of refHashSet is: " << refHashSet.size() << endl;

	vector<sketch_t> subSketches;
	for(int i = 0; i < querySketches.size(); i++){
		sketch_t s;
		s.fileName = querySketches[i].fileName;
		s.id = querySketches[i].id;
		vector<uint32_t> newHashSet;
		for(auto x : querySketches[i].hashSet){
			if(refHashSet.count(x) == 0){
				newHashSet.push_back(x);
			}
		}
		s.hashSet = newHashSet;
		subSketches.push_back(s);
	}
	printSketches(subSketches, "hello.out");
	
	saveSketches(subSketches, outputFile);
}
	









