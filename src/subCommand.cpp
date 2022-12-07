#include "subCommand.h"
#include "sketch.h"
#include "dist.h"
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_set>
#include "robin_hood.h"
#include <err.h>

void command_convert(string inputDir, bool isQuery, string outputFile, int threads){
	double t0 = get_sec();
	vector<sketch_t> sketches;
	sketchInfo_t info;
	convertSketch(sketches, info, inputDir, threads);
	double t1 = get_sec();
	cerr << "===================time of converting sketches from Kssd to RabbitKSSD format is " << t1 - t0 << endl;
	saveSketches(sketches, info, outputFile);
	double t2 = get_sec();
	cerr << "===================time of saving sketches to file is " << t2 - t1 << endl;

	if(!isQuery){
		double tstart = get_sec();
		string dictFile = outputFile + ".dict";
		string indexFile = outputFile + ".index";
		transSketches(sketches, info, dictFile, indexFile, threads);
		double tend = get_sec();
		cerr << "===============the time of transSketches is: " << tend - tstart << endl;
	}
}


void command_sketch(string refList, bool isQuery, string outputFile, kssd_parameter_t kssd_parameter, int threads){
	vector<sketch_t> sketches;
	sketchInfo_t info;
	//bool isReference = true;
	bool success = sketchFile(refList, isQuery, threads, kssd_parameter, sketches, outputFile);
}

void command_info(string sketchFile, bool isDetail, string outputFile){
	vector<sketch_t> sketches;
	sketchInfo_t info;
	readSketches(sketches, info, sketchFile);
	cerr << "the number of genome is: " << sketches.size() << endl;
	if(isDetail)
		printSketches(sketches, outputFile);
	else
		printInfos(sketches, outputFile);
}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int isContainment, int threads){
	double t2 = get_sec();
	double t3;
	string refSketchOut;
	vector<sketch_t> sketches;
	sketchInfo_t info;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(sketches, info, refList);
		t3 = get_sec();
		cerr << "===================time of read sketches from file is " << t3 - t2 << endl;
	}
	else{
		refSketchOut = refList + ".sketch";
		bool success = sketchFile(refList, true, threads, kssd_parameter, sketches, refSketchOut);
		t3 = get_sec();
		cerr << "===================time of computing sketches and save sketches into file is " << t3 - t2 << endl;
	}
	//cerr << "the size of sketches is: " << sketches.size() << endl;

	//compute the pair distance
	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	index_tridist(sketches, refSketchOut, outputFile, kmer_size, maxDist, isContainment, threads);

}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, double maxDist, int isContainment, int threads){
	double t2 = get_sec();
	double t3, t4;
	string refSketchOut, querySketchOut;
	vector<sketch_t> ref_sketches;
	sketchInfo_t info;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(ref_sketches, info, refList);
		t3 = get_sec();
		cerr << "===================time of read reference sketches from " << refList << " is: " << t3 - t2 << endl;
	}
	else{
		refSketchOut = refList + ".sketch";
		bool success0 = sketchFile(refList, true, threads, kssd_parameter, ref_sketches, refSketchOut);
		t3 = get_sec();
		cerr << "===================time of computing reference sketches and save sketches into " << refSketchOut << " is: " << t3 - t2 << endl;
	}
	//cerr << "the size of ref_sketches is: " << ref_sketches.size() << endl;

	vector<sketch_t> query_sketches;
	if(isSketchFile(queryList)){
		querySketchOut = queryList;
		readSketches(query_sketches, info, queryList);
		t4 = get_sec();
		cerr << "===================time of read query sketches from " << queryList << " is: " << t4 - t3 << endl;
	}
	else{
		querySketchOut = queryList + ".sketch";
		bool success1 = sketchFile(queryList, false, threads, kssd_parameter, query_sketches, querySketchOut);
		t4 = get_sec();
		cerr << "===================time of computing query sketches and save sketches into " << querySketchOut << " is: " << t4 - t3 << endl;
	}
	//cerr << "the size of query_sketches is: " << query_sketches.size() << endl;

	//dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);
	index_dist(ref_sketches, refSketchOut, query_sketches, outputFile, kmer_size, maxDist, isContainment, threads);

	//double t5 = get_sec();
	//cerr << "===================time of get total distance matrix file is: " << t5 - t4 << endl;
}

void command_merge(string sketchFile, string outputFile, int threads){
	
	vector<sketch_t> sketches;
	if(!isSketchFile(sketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", sketchFile.c_str());
	}
	#ifdef Timer_inner
	double t0 = get_sec();
	#endif
	sketchInfo_t info;
	readSketches(sketches, info, sketchFile);
	
	string totalName = outputFile + " merge sketches";
	vector<uint32_t> mergedArr;

	#ifdef Timer_inner
	double t1 = get_sec();
	cerr << "-----time of read sketches is: " << t1 - t0 << endl;
	#endif

	size_t dictSize = (1LLU << 32) / 64;
	uint64_t * dict = (uint64_t*)malloc(dictSize * sizeof(uint64_t));
	memset(dict, 0, dictSize * sizeof(uint64_t));
	for(int i = 0; i < sketches.size(); i++){
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			uint32_t curHash = sketches[i].hashSet[j];
			dict[curHash/64] |= (0x8000000000000000LLU >> (curHash % 64));
		}
	}
	for(int i = 0; i < dictSize; i++){
		if(dict[i]){
			for(int j = 0; j < 64; j++){
				if((0x8000000000000000LLU >> j) & dict[i]){
					uint32_t hash = 64 * i + j;
					mergedArr.push_back(hash);
				}
			}
		}
	}

	#ifdef Timer_inner
	double t2 = get_sec();
		cerr << "-----time of set union by dictionary is: " << t2 - t1 << endl;
	#endif

	std::sort(mergedArr.begin(), mergedArr.end());
	#ifdef Timer_inner
	double t3 = get_sec();
		cerr << "-----time of sort sketch hash values is: " << t3 - t2 << endl;
	#endif

	vector<sketch_t> mergedSketches;
	sketch_t s;
	s.id = 0;
	s.fileName = totalName;
	s.hashSet = mergedArr;
	mergedSketches.push_back(s);
	//info.genomeNumber = 1;
	saveSketches(mergedSketches, info, outputFile);
	#ifdef Timer_inner
	double t4 = get_sec();
	cerr << "-----time of save sketches is: " << t4 - t3 << endl;
	#endif

	//string dictFile = outputFile + ".dict";
	//string indexFile = outputFile + ".index";
	//transSketches(mergedSketches, info, dictFile, indexFile, threads);
}

void command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads){
	#ifdef Timer_inner
	double t0 = get_sec();
	#endif
	vector<sketch_t> refSketches;
	if(!isSketchFile(refSketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", refSketchFile.c_str());
	}
	sketchInfo_t info;
	readSketches(refSketches, info, refSketchFile);

	#ifdef Timer_inner
	double t1 = get_sec();
	cerr << "-----time of read reference sketch: " << refSketchFile << " is: " << t1 - t0 << endl;
	#endif
	vector<sketch_t> querySketches;
	if(!isSketchFile(querySketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", querySketchFile.c_str());
	}
	readSketches(querySketches, info, querySketchFile);

	#ifdef Timer_inner
	double t2 = get_sec();
	cerr << "-----time of read query sketch: " << querySketchFile << " is: " << t2 - t1 << endl;
	#endif
	//cerr << "the refSketches size is: " << refSketches.size() << endl;
	//cerr << "the querySketches size is: " << querySketches.size() << endl;
	
	double t3;
	vector<sketch_t> subSketches;

if(1){
	size_t dictSize = (1LLU << 32) / 64;
	uint64_t * dict = (uint64_t*)malloc(dictSize * sizeof(uint64_t));
	memset(dict, 0, dictSize * sizeof(uint64_t));
	for(int i = 0; i < refSketches.size(); i++){
		for(int j = 0; j < refSketches[i].hashSet.size(); j++){
			uint32_t curHash = refSketches[i].hashSet[j];
			dict[curHash/64] |= (0x8000000000000000LLU >> (curHash % 64));
		}
	}

	#ifdef Timer_inner
	t3 = get_sec();
	cerr << "-----time of generate reference dictionary is: " << t3 - t2 << endl;
	#endif

	#pragma omp parallel for num_threads(threads)
	for(int i = 0; i < querySketches.size(); i++){
		sketch_t s;
		s.fileName = querySketches[i].fileName;
		s.id = querySketches[i].id;
		vector<uint32_t> newHashArr;
		for(int j = 0; j < querySketches[i].hashSet.size(); j++){
			uint32_t curHash = querySketches[i].hashSet[j];
			uint64_t isIn = dict[curHash/64] & (0x8000000000000000 >> (curHash % 64));
			if(!isIn){
				newHashArr.emplace_back(curHash);
			}
		}
		//std::sort(newHashArr.begin(), newHashArr.end());
		s.hashSet = newHashArr;
		#pragma omp critical
		{
			subSketches.push_back(s);
		}
	}
	//std::sort(subSketches.begin(), subSketches.end(), cmpSketch);
}

else{
	robin_hood::unordered_set<uint32_t> refHashSet;
	//unordered_set<uint32_t> refHashSet;
	for(int i = 0; i < refSketches.size(); i++){
		for(auto x : refSketches[i].hashSet){
			refHashSet.insert(x);
		}
	}
	cerr << "the size of refHashSet is: " << refHashSet.size() << endl;

	#ifdef Timer_inner
	t3 = get_sec();
	cerr << "-----time of generate reference dictionary is: " << t3 - t2 << endl;
	#endif

	#pragma omp parallel for num_threads(threads)
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
		//std::sort(newHashArr.begin(), newHashArr.end());
		s.hashSet = newHashArr;
		#pragma omp critical
		{
		subSketches.push_back(s);
		}
	}
	//std::sort(subSketches.begin(), subSketches.end(), cmpSketch);
}

	#ifdef Timer_inner
	double t4 = get_sec();
	cerr << "-----time of multithreading substraction sketches is: " << t4 - t3 << endl;
	#endif
	
	saveSketches(subSketches, info, outputFile);

	#ifdef Timer_inner
	double t5 = get_sec();
	cerr << "-----time of save sketches is: " << t5 - t4 << endl;
	#endif
	//string dictFile = outputFile + ".dict";
	//string indexFile = outputFile + ".index";
	//transSketches(subSketches, info, dictFile, indexFile, threads);
}
	








