#include "subCommand.h"
#include "sketch.h"
#include "dist.h"
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_set>
#include "robin_hood.h"
#include <err.h>


void command_sketch(string refList, string outputFile, kssd_parameter_t kssd_parameter, int threads){
	vector<sketch_t> sketches;
	bool isReference = true;
	bool success = sketchFile(refList, isReference, threads, kssd_parameter, sketches, outputFile);
}

void command_info(string sketchFile, string outputFile){
	vector<sketch_t> sketches;
	readSketches(sketches, sketchFile);
	cerr << "the number of genome is: " << sketches.size() << endl;
	//printInfos(sketches, outputFile);
	printSketches(sketches, outputFile);
}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	double t3;
	string refSketchOut;
	vector<sketch_t> sketches;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(sketches, refList);
		t3 = get_sec();
		cerr << "===================time of read sketches from file is " << t3 - t2 << endl;
	}
	else{
		refSketchOut = refList + ".sketch";
		bool success = sketchFile(refList, true, threads, kssd_parameter, sketches, refSketchOut);
		t3 = get_sec();
		cerr << "===================time of computing sketches and save sketches into file is " << t3 - t2 << endl;
		//cerr << "finish the sketch generation " << success << endl;
	}
	//cerr << "the size of sketches is: " << sketches.size() << endl;

	//compute the pair distance
	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	index_dist(sketches, refSketchOut, outputFile, kmer_size, maxDist, threads);

}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int threads){
	double t2 = get_sec();
	double t3, t4;
	vector<sketch_t> ref_sketches;
	if(isSketchFile(refList)){
		readSketches(ref_sketches, refList);
		t3 = get_sec();
		cerr << "===================time of read reference sketches from " << refList << " is: " << t3 - t2 << endl;
	}
	else{
		string refHashOut = refList + ".sketch";
		bool success0 = sketchFile(refList, true, threads, kssd_parameter, ref_sketches, refHashOut);
		t3 = get_sec();
		cerr << "===================time of computing reference sketches and save sketches into " << refHashOut << " is: " << t3 - t2 << endl;
		//cerr << "finish the ref_sketches generation " << success0 << endl;
	}
	//cerr << "the size of ref_sketches is: " << ref_sketches.size() << endl;

	vector<sketch_t> query_sketches;
	if(isSketchFile(queryList)){
		readSketches(query_sketches, queryList);
		t4 = get_sec();
		cerr << "===================time of read query sketches from " << queryList << " is: " << t4 - t3 << endl;
	}
	else{
		string queryHashOut = queryList + ".sketch";
		bool success1 = sketchFile(queryList, false, threads, kssd_parameter, query_sketches, queryHashOut);
		t4 = get_sec();
		cerr << "===================time of computing query sketches and save sketches into " << queryHashOut << " is: " << t4 - t3 << endl;
		//cerr << "finish the query_sketches generation " << success1 << endl;
	}
	//cerr << "the size of query_sketches is: " << query_sketches.size() << endl;


	//compute the pair distance
	
	dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);

	//double t5 = get_sec();
	//cerr << "===================time of get total distance matrix file is: " << t5 - t4 << endl;
}

void command_merge(string sketchFile, string outputFile, int threadNumber){
	
	vector<sketch_t> sketches;
	if(!isSketchFile(sketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", sketchFile.c_str());
	}
	readSketches(sketches, sketchFile);
	//exit(0);
	string totalName("");
	robin_hood::unordered_set<uint32_t> mergedSet;
	//unordered_set<uint32_t> mergedSet;
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

	std::sort(mergedArr.begin(), mergedArr.end());

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
	if(!isSketchFile(refSketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", refSketchFile.c_str());
	}
	readSketches(refSketches, refSketchFile);
	vector<sketch_t> querySketches;
	if(!isSketchFile(querySketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", querySketchFile.c_str());
	}
	readSketches(querySketches, querySketchFile);
	cerr << "the refSketches size is: " << refSketches.size() << endl;
	cerr << "the querySketches size is: " << querySketches.size() << endl;

	robin_hood::unordered_set<uint32_t> refHashSet;
	//unordered_set<uint32_t> refHashSet;
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
		vector<uint32_t> newHashArr;
		for(auto x : querySketches[i].hashSet){
			if(refHashSet.count(x) == 0){
				newHashArr.push_back(x);
			}
		}
		std::sort(newHashArr.begin(), newHashArr.end());
		s.hashSet = newHashArr;
		subSketches.push_back(s);
	}
	//printSketches(subSketches, "hello.out");
	
	saveSketches(subSketches, outputFile);
}
	








