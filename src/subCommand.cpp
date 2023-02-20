#include "subCommand.h"
#include "sketch.h"
#include "dist.h"
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <unordered_set>
#include "robin_hood.h"
#include <err.h>
#include <queue>
#include <omp.h>

void command_convert(string inputDir, bool to_Kssd_sketch, bool isQuery, string outputFile, int threads){
	if(!to_Kssd_sketch){
		double t0 = get_sec();
		vector<sketch_t> sketches;
		sketchInfo_t info;
		convertSketch(sketches, info, inputDir, threads);
		double t1 = get_sec();
		cerr << "===================time of converting sketches from Kssd to RabbitKSSD format is " << t1 - t0 << endl;
		if(!isSketchFile(outputFile)){
			outputFile = outputFile + ".sketch";
		}
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
	else{
		if(!isSketchFile(inputDir)){
			cerr << "need input RabbitKSSD sketch file: " << inputDir << endl;
			exit(1);
		}
		vector<sketch_t> sketches;
		sketchInfo_t info;
		readSketches(sketches, info, inputDir);
		convert_from_RabbitKSSDSketch_to_KssdSketch(sketches, info, outputFile, threads);
	}
}


void command_sketch(string refList, bool isQuery, string outputFile, kssd_parameter_t kssd_parameter, int leastNumKmer, int threads){
	vector<sketch_t> sketches;
	sketchInfo_t info;
	//bool isReference = true;
	//bool success = sketchFile(refList, isQuery, threads, kssd_parameter, sketches, outputFile);
	//bool success = sketchFastaFile(refList, isQuery, threads, kssd_parameter, sketches, outputFile);
	bool success;
	if(isFastaList(refList)){
		success = sketchFastaFile(refList, isQuery, threads, kssd_parameter, sketches, outputFile);
	}
	else if(isFastqList(refList)){
		success = sketchFastqFile(refList, isQuery, threads, kssd_parameter, leastNumKmer, sketches, outputFile);
	}
	else{
		cerr << "the input file list for sketching must be list of fasta and fastq file" << endl;
		exit(1);
	}
}

void command_info(string sketchFile, bool isDetail, string outputFile){
	//vector<sketch_t> sketches;
	sketchInfo_t info;
	//readSketches(sketches, info, sketchFile);
	FILE * fp = fopen(sketchFile.c_str(), "rb+");
	fread(&info, sizeof(sketchInfo_t), 1, fp);
	int sketchNumber = info.genomeNumber;
	cerr << "the number of genome is: " << sketchNumber << endl;
	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	fread(genomeNameSize, sizeof(int), sketchNumber, fp);
	fread(hashSetSize, sizeof(int), sketchNumber, fp);
	FILE * fp1 = fopen(outputFile.c_str(), "w+");
	fprintf(fp1, "the number of sketches are: %d\n", sketchNumber);
	int maxNameLength=1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	for(int i = 0; i < sketchNumber; i++){
		int curLength = genomeNameSize[i];
		if(curLength > maxNameLength){
			maxNameLength = curLength;
			curName = new char[maxNameLength+1];
		}
		//char * curName = new char[curLength+1];
		int nameLength = fread(curName, sizeof(char), curLength, fp);
		if(nameLength != curLength){
			cerr << "error: the read nameLength is not equal to the saved nameLength, exit!" << endl;
			exit(0);
		}
		fprintf(fp1, "%s\t%d\n", curName, hashSetSize[i]);
		int curSize = hashSetSize[i];
		if(curSize > maxHashSize){
			maxHashSize = curSize;
			curPoint = new uint32_t[maxHashSize];
		}
		//uint32_t * curPoint = new uint32_t[curSize];
		int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp);
		if(hashSize != curSize){
			cerr << "error: the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
			exit(0);
		}
		if(isDetail){
			for(int j = 0; j < curSize; j++){
				fprintf(fp1, "%u\t", curPoint[j]);
				if(j % 10 == 9) fprintf(fp1, "\n");
			}
			fprintf(fp1, "\n");
		}
	}
	fclose(fp);
	fclose(fp1);
	//if(isDetail)
	//	printSketches(sketches, outputFile);
	//else
	//	printInfos(sketches, outputFile);
}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, int leastNumKmer, double maxDist, int isContainment, int threads){
	double t2 = get_sec();
	double t3;
	string refSketchOut;
	vector<sketch_t> sketches;
	sketchInfo_t info;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(sketches, info, refList);
		string indexFile = refSketchOut + ".index";
		string dictFile = refSketchOut + ".dict";
		if(!existFile(indexFile) || !existFile(dictFile)){
			transSketches(sketches, info, dictFile, indexFile, threads);
		}
		t3 = get_sec();
		cerr << "===================time of read sketches from file is " << t3 - t2 << endl;
	}
	else{
		refSketchOut = refList + ".sketch";
		bool success;
		if(isFastaList(refList)){
			success = sketchFastaFile(refList, false, threads, kssd_parameter, sketches, refSketchOut);
		}
		else if(isFastqList(refList)){
			success = sketchFastqFile(refList, false, threads, kssd_parameter, leastNumKmer, sketches, refSketchOut);
		}
		else{
			cerr << "the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}

		t3 = get_sec();
		cerr << "===================time of computing sketches and save sketches into file is " << t3 - t2 << endl;
	}
	//cerr << "the size of sketches is: " << sketches.size() << endl;

	//compute the pair distance
	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	index_tridist(sketches, refSketchOut, outputFile, kmer_size, maxDist, isContainment, threads);

}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int kmer_size, int leastNumKmer, double maxDist, int maxNeighbor, bool isNeighbor, int isContainment, int threads){
	double t2 = get_sec();
	double t3, t4;
	string refSketchOut, querySketchOut;
	vector<sketch_t> ref_sketches;
	sketchInfo_t info;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(ref_sketches, info, refList);
		string indexFile = refSketchOut + ".index";
		string dictFile = refSketchOut + ".dict";
		if(!existFile(indexFile) || !existFile(dictFile)){
			transSketches(ref_sketches, info, dictFile, indexFile, threads);
		}
		t3 = get_sec();
		cerr << "===================time of read reference sketches from " << refList << " is: " << t3 - t2 << endl;
	}
	else{
		refSketchOut = refList + ".sketch";
		bool success0;
		if(isFastaList(refList)){
			success0 = sketchFastaFile(refList, false, threads, kssd_parameter, ref_sketches, refSketchOut);
		}
		else if(isFastqList(refList)){
			success0 = sketchFastqFile(refList, false, threads, kssd_parameter, leastNumKmer, ref_sketches, refSketchOut);
		}
		else{
			cerr << "the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}


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
		bool success1;
		if(isFastaList(queryList)){
			success1 = sketchFastaFile(queryList, true, threads, kssd_parameter, query_sketches, querySketchOut);
		}
		else if(isFastqList(queryList)){
			success1 = sketchFastqFile(queryList, true, threads, kssd_parameter, leastNumKmer, query_sketches, querySketchOut);
		}
		else{
			cerr << "the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}
		
		t4 = get_sec();
		cerr << "===================time of computing query sketches and save sketches into " << querySketchOut << " is: " << t4 - t3 << endl;
	}
	//cerr << "the size of query_sketches is: " << query_sketches.size() << endl;

	//dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);
	index_dist(ref_sketches, refSketchOut, query_sketches, outputFile, kmer_size, maxDist, maxNeighbor, isNeighbor, isContainment, threads);

	//double t5 = get_sec();
	//cerr << "===================time of get total distance matrix file is: " << t5 - t4 << endl;
}

void command_union(string sketchFile, string outputFile, int threads){
	
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

	//std::sort(mergedArr.begin(), mergedArr.end());
	#ifdef Timer_inner
	double t3 = get_sec();
	//cerr << "-----time of sort sketch hash values is: " << t3 - t2 << endl;
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

void new_command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads){
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
	double t2 = get_sec();
	cerr << "-----time of generate reference dictionary is: " << t2 - t1 << endl;
	#endif

	vector<sketch_t> querySketches;
	if(!isSketchFile(querySketchFile)){
		err(errno, "error: %s is not sketch file, need input sketch file\n", querySketchFile.c_str());
	}
	//readSketches(querySketches, info, querySketchFile);
	sketchInfo_t query_info;
	FILE* fp_query = fopen(querySketchFile.c_str(), "rb");
	if(!fp_query){
		cerr << "error open the file: " << querySketchFile << endl;
		exit(1);
	}
	fread(&query_info, sizeof(sketchInfo_t), 1, fp_query);
	int query_num = query_info.genomeNumber;
	int * genomeNameSize = new int[query_num];
	int * hashSetSize = new int[query_num];
	fread(genomeNameSize, sizeof(int), query_num, fp_query);
	fread(hashSetSize, sizeof(int), query_num, fp_query);

	int maxNameLength = 1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];

	FILE* fp_result = fopen(outputFile.c_str(), "w+");
	if(!fp_result){
		cerr << "error open the output file: " << outputFile << endl;
		exit(1);
	}
	int * out_genomeNameSize = new int [query_num];
	int * out_hashSetSize = new int[query_num];
	fwrite(&query_info, sizeof(sketchInfo_t), 1, fp_result);
	fwrite(out_genomeNameSize, sizeof(int), query_num, fp_result);
	fwrite(out_hashSetSize, sizeof(int), query_num, fp_result);

	vector<sketch_t> subSketches;
	queue<sketch_t> sketch_queue;
	int index = 0;
	int produced_num = 0;
	int solved_num = 0;
	omp_lock_t queue_lock;
	omp_lock_t vector_lock;
	omp_init_lock(&queue_lock);
	omp_init_lock(&vector_lock);
	sketch_t s_arr[threads];
	#pragma omp parallel num_threads(threads)
	{
		int tid = omp_get_thread_num();
		//cerr << "DEBUG: the tid is: " << tid << endl;
		if(tid == 0)//producer
		{
			//cerr << "DEBUG: enter the producer " << endl;
			for(index = 0; index < query_num; index++){
				int curLength = genomeNameSize[index];
				if(curLength > maxNameLength){
					maxNameLength = curLength;
					curName = new char[maxNameLength+1];
				}
				int nameLength = fread(curName, sizeof(char), curLength, fp_query);
				string genomeName;
				genomeName.assign(curName, curName + curLength);

				int curSize = hashSetSize[index];
				if(curSize > maxHashSize){
					maxHashSize = curSize;
					curPoint = new uint32_t[maxHashSize];
				}
				int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp_query);

				vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
				sketch_t s;
				s.fileName = genomeName;
				s.id = index;
				s.hashSet=curHashSet;
				//#pragma omp critical
				omp_set_lock(&queue_lock);
				{
					sketch_queue.push(s);
					produced_num++;
					//cerr << "DEBUG: index is: " << index << endl;
				}
				omp_unset_lock(&queue_lock);
			}//end for
		}
		else//consumer
		{
			//while(solved_num < query_num)//may conflict when using multiple threads.
			while(1){
				omp_set_lock(&queue_lock);
				{
					if(sketch_queue.empty()){
						//cerr << "DEBUG: sketch_queue is empty " << tid << endl;
						omp_unset_lock(&queue_lock);
						if(solved_num >= query_num)	break;
						continue;
					}
					s_arr[tid] = sketch_queue.front();
					sketch_queue.pop();
				}
				omp_unset_lock(&queue_lock);
				
				vector<uint32_t> newHashArr;
				//cerr << "DEBUG: s.hashSet.size() is: " << s.hashSet.size() << endl;
				for(int j = 0; j < s_arr[tid].hashSet.size(); j++){
					uint32_t curHash = s_arr[tid].hashSet[j];
					uint64_t isIn = dict[curHash/64] & (0x8000000000000000 >> (curHash % 64));
					if(!isIn){
						newHashArr.emplace_back(curHash);
					}
				}
				sketch_t new_s;
				new_s.fileName = s_arr[tid].fileName;
				new_s.id = s_arr[tid].id;
				new_s.hashSet = newHashArr;
				omp_set_lock(&vector_lock);
				{
					//subSketches.push_back(new_s);
					//solved_num++;
					////cerr << "DEBUG: solved_num is: " << solved_num << endl;
					const char * namePoint = new_s.fileName.c_str();
					int curNameLength = new_s.fileName.length();
					uint32_t * curPoint = new_s.hashSet.data();
					int curHashSetSize = new_s.hashSet.size();
					fwrite(namePoint, sizeof(char), curNameLength, fp_result);
					fwrite(curPoint, sizeof(uint32_t), curHashSetSize, fp_result);
					out_genomeNameSize[solved_num] = curNameLength;
					out_hashSetSize[solved_num] = curHashSetSize;
					solved_num++;
				}
				omp_unset_lock(&vector_lock);
			}
		}
	}

	rewind(fp_result);
	fwrite(&query_info, sizeof(sketchInfo_t), 1, fp_result);
	fwrite(out_genomeNameSize, sizeof(int), query_num, fp_result);
	fwrite(out_hashSetSize, sizeof(int), query_num, fp_result);

	#ifdef Timer_inner
	double t3 = get_sec();
	cerr << "-----time of substraction is: " << t3 - t2 << endl;
	#endif

	//saveSketches(subSketches, info, outputFile);

	#ifdef Timer_inner
	double t4 = get_sec();
	cerr << "-----time of save sketches is: " << t4 - t3 << endl;
	#endif

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


void command_merge(string inputList, string outputFile, int threads){
	ifstream ifs(inputList);
	string line;
	vector<string> fileList;
	while(getline(ifs, line)){
		if(!isSketchFile(line)){
			cerr << "the file: " << line << " is not a sketch file in the list file: " << inputList << endl;
			exit(1);
		}
		fileList.push_back(line);
	}
	ifs.close();

	int totalSketchNumber = 0;
	vector<int> totalGenomeNameSize;
	vector<int> totalHashSetSize;
	sketchInfo_t resInfo;

	vector<sketch_t> resSketches;
	//readSketches(resSketches, resInfo, fileList[0]);
	FILE * fp0 = fopen(fileList[0].c_str(), "rb+");
	fread(&resInfo, sizeof(sketchInfo_t), 1, fp0);
	int sketchNumber0 = resInfo.genomeNumber;
	int * genomeNameSize0 = new int[sketchNumber0];
	int * hashSetSize0 = new int[sketchNumber0];
	fread(genomeNameSize0, sizeof(int), sketchNumber0, fp0);
	fread(hashSetSize0, sizeof(int), sketchNumber0, fp0);
	totalSketchNumber += sketchNumber0;
	for(int j = 0; j < sketchNumber0; j++){
		totalGenomeNameSize.emplace_back(genomeNameSize0[j]);
		totalHashSetSize.emplace_back(hashSetSize0[j]);
	}

	for(int i = 1; i < fileList.size(); i++){
		FILE * fpTmp = fopen(fileList[i].c_str(), "rb+");
		sketchInfo_t curInfo;
		fread(&curInfo, sizeof(sketchInfo_t), 1, fpTmp);
		if(curInfo.half_k != resInfo.half_k || curInfo.half_subk != resInfo.half_subk || curInfo.drlevel != resInfo.drlevel){
			cerr << "the sketch info in the list file: " << inputList << " is not same, cannot merge" << endl;
			exit(1);
		}
		int curSketchNumber = curInfo.genomeNumber;
		int * genomeNameSize = new int[curSketchNumber];
		int * hashSetSize = new int[curSketchNumber];
		fread(genomeNameSize, sizeof(int), curSketchNumber, fpTmp);
		fread(hashSetSize, sizeof(int), curSketchNumber, fpTmp);
		totalSketchNumber += curSketchNumber;
		for(int j = 0; j < curSketchNumber; j++){
			totalGenomeNameSize.emplace_back(genomeNameSize[j]);
			totalHashSetSize.emplace_back(hashSetSize[j]);
		}
		fclose(fpTmp);
	}
	//cerr << "finish the totalSketchNumber calculation" << endl;

	FILE* fp1 = fopen(outputFile.c_str(), "w+");
	resInfo.genomeNumber = totalSketchNumber;
	fwrite(&resInfo, sizeof(sketchInfo_t), 1, fp1);
	fwrite(totalGenomeNameSize.data(), sizeof(int), totalSketchNumber, fp1);
	fwrite(totalHashSetSize.data(), sizeof(int), totalSketchNumber, fp1);
	fclose(fp1);
	
	FILE* fp = fopen(outputFile.c_str(), "a+");
	for(int i = 0; i < fileList.size(); i++){
		vector<sketch_t> curSketches;
		sketchInfo_t curInfo;
		readSketches(curSketches, curInfo, fileList[i]);

		int sketchNumber = curSketches.size();
		for(int i = 0; i < sketchNumber; i++){
			const char * namePoint = curSketches[i].fileName.c_str();
			fwrite(namePoint, sizeof(char),curSketches[i].fileName.length(), fp);
			uint32_t * curPoint = curSketches[i].hashSet.data();
			fwrite(curPoint, sizeof(uint32_t),curSketches[i].hashSet.size(), fp);
		}
	}
	fclose(fp);

}
	








