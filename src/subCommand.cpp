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
			cerr << "ERROR: command_convert(), need input RabbitKSSD sketch file: " << inputDir << endl;
			exit(1);
		}
		vector<sketch_t> sketches;
		sketchInfo_t info;
		readSketches(sketches, info, inputDir);
		convert_from_RabbitKSSDSketch_to_KssdSketch(sketches, info, outputFile, threads);
	}
}


void command_sketch(string refList, bool isQuery, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, int threads){
	vector<sketch_t> sketches;
	sketchInfo_t info;
	bool success;
	if(isFastaList(refList) || isFastaGZList(refList)){
		success = sketchFastaFile(refList, isQuery, threads, kssd_parameter, sketches, info, outputFile);
	}
	else if(isFastqList(refList) || isFastqGZList(refList)){
		success = sketchFastqFile(refList, isQuery, threads, kssd_parameter, leastQual, leastNumKmer, sketches, info, outputFile);
	}
	else{
		cerr << "ERROR: command_sketch(), the input file list for sketching must be list of fasta and fastq file in normal format or gz format" << endl;
		exit(1);
	}
	if(!success){
		cerr << "ERROR: command_sketch(), generate sketches failed" << endl;
		exit(1);
	}
}

void command_info(string sketchFile, bool isDetail, string outputFile){
	//vector<sketch_t> sketches;
	sketchInfo_t info;
	//readSketches(sketches, info, sketchFile);
	FILE * fp = fopen(sketchFile.c_str(), "rb");
	if(!fp){
		fprintf(stderr, "ERROR: command_info(), cannot open file: %s\n", sketchFile.c_str());
		exit(1);
	}
	int read_sketch_info = fread(&info, sizeof(sketchInfo_t), 1, fp);
	int sketchNumber = info.genomeNumber;
	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	cerr << "the number of genome is: " << sketchNumber << endl;
	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	int read_genome_name_size = fread(genomeNameSize, sizeof(int), sketchNumber, fp);
	int read_hash_set_size = fread(hashSetSize, sizeof(int), sketchNumber, fp);
	if(read_sketch_info != 1 || read_genome_name_size != sketchNumber || read_hash_set_size != sketchNumber){
		cerr << "ERROR: command_info(), mismatched read sketch_info, genome_name_size, hash_set_size" << endl;
		exit(1);
	}

	FILE * fp1 = fopen(outputFile.c_str(), "w+");
	fprintf(fp1, "the number of sketches are: %d\n", sketchNumber);
	int maxNameLength=1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	uint64_t * curPoint64 = new uint64_t[maxHashSize];
	for(int i = 0; i < sketchNumber; i++){
		int curLength = genomeNameSize[i];
		if(curLength > maxNameLength){
			maxNameLength = curLength;
			curName = new char[maxNameLength+1];
		}
		int nameLength = fread(curName, sizeof(char), curLength, fp);
		if(nameLength != curLength){
			cerr << "ERROR: command_info(), the read nameLength is not equal to the saved nameLength, exit!" << endl;
			exit(1);
		}
		string genomeName;
		genomeName.assign(curName, curName + curLength);
		fprintf(fp1, "%s\t%d\n", genomeName.c_str(), hashSetSize[i]);
		int curSize = hashSetSize[i];
		if(curSize > maxHashSize){
			maxHashSize = curSize;
			curPoint = new uint32_t[maxHashSize];
			curPoint64 = new uint64_t[maxHashSize];
		}
		int hashSize;
		if(use64){
			hashSize = fread(curPoint64, sizeof(uint64_t), curSize, fp);
		}
		else
			hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp);
		if(hashSize != curSize){
			cerr << "ERROR: command_info(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
			exit(1);
		}
		if(isDetail){
			for(int j = 0; j < curSize; j++){
				if(use64)
					fprintf(fp1, "%lu\t", curPoint64[j]);
				else
					fprintf(fp1, "%u\t", curPoint[j]);
				if(j % 10 == 9) fprintf(fp1, "\n");
			}
			fprintf(fp1, "\n");
		}
	}
	fclose(fp);
	fclose(fp1);
	delete [] genomeNameSize;
	delete [] hashSetSize;
	delete [] curName;
	delete [] curPoint;
	delete [] curPoint64;
}

void command_alldist(string refList, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, double maxDist, int isContainment, int threads){
	double t2 = get_sec();
	double t3;
	if(maxDist < 0.0){
		cerr << "ERROR: command_alldist(), maxDist must be > 0" << endl;
		cerr << "Use -D to set the maxDist" << endl;
		exit(1);
	}
	int kmer_size = kssd_parameter.half_k * 2;
	string refSketchOut;
	vector<sketch_t> sketches;
	sketchInfo_t info;
	if(isSketchFile(refList)){
		refSketchOut = refList;
		readSketches(sketches, info, refList);
		kmer_size = info.half_k * 2;
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
			success = sketchFastaFile(refList, false, threads, kssd_parameter, sketches, info, refSketchOut);
		}
		else if(isFastqList(refList)){
			success = sketchFastqFile(refList, false, threads, kssd_parameter, leastQual, leastNumKmer, sketches, info, refSketchOut);
		}
		else{
			cerr << "ERROR: command_alldist(), the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}
		if(!success){
			cerr << "ERROR: command_alldist(), filed sketching genomes" << endl;
			exit(1);
		}

		t3 = get_sec();
		cerr << "===================time of computing sketches and save sketches into file is " << t3 - t2 << endl;
	}
	//cerr << "the size of sketches is: " << sketches.size() << endl;

	//compute the pair distance
	//tri_dist(sketches, outputFile, kmer_size, maxDist, threads);
	index_tridist(sketches, info, refSketchOut, outputFile, kmer_size, maxDist, isContainment, threads);

}


void command_dist(string refList, string queryList, string outputFile, kssd_parameter_t kssd_parameter, int leastQual, int leastNumKmer, double maxDist, int maxNeighbor, bool isNeighbor, int isContainment, int threads){
	double t2 = get_sec();
	double t3, t4;
	if(maxDist < 0.0){
		cerr << "ERROR: command_dist(), maxDist must be > 0" << endl;
		cerr << "Use -D to set the maxDist" << endl;
		exit(1);
	}
	if(maxNeighbor < 0){
		cerr << "ERROR: command_dist(), maxNeighbor must be > 0" << endl;
		cerr << "Use -N to set the maxNeighbor " << endl;
		exit(1);
	}
	int kmer_size = kssd_parameter.half_k * 2;
	string refSketchOut, querySketchOut;

	vector<sketch_t> ref_sketches;
	sketchInfo_t ref_info;
	if(isSketchFile(refList))
	{
		refSketchOut = refList;
		readSketches(ref_sketches, ref_info, refList);

		kmer_size = ref_info.half_k * 2;
		string indexFile = refSketchOut + ".index";
		string dictFile = refSketchOut + ".dict";
		if(!existFile(indexFile) || !existFile(dictFile)){
			transSketches(ref_sketches, ref_info, dictFile, indexFile, threads);
		}
		cerr << "the ref_sketch size is: " << ref_sketches.size() << endl;
		t3 = get_sec();
		#ifdef Timer
		cerr << "===================time of read reference sketches from " << refList << " is: " << t3 - t2 << endl;
		#endif
	}
	else
	{
		refSketchOut = refList + ".sketch";
		bool success0;
		if(isFastaList(refList)){
			success0 = sketchFastaFile(refList, false, threads, kssd_parameter, ref_sketches, ref_info, refSketchOut);
		}
		else if(isFastqList(refList)){
			success0 = sketchFastqFile(refList, false, threads, kssd_parameter, leastQual, leastNumKmer, ref_sketches, ref_info, refSketchOut);
		}
		else{
			cerr << "ERROR: command_dist(), the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}
		if(!success0){
			cerr << "ERROR: command_dist(), failed sketching reference genomes" << endl;
			exit(1);
		}

		t3 = get_sec();
		#ifdef Timer
		cerr << "===================time of computing reference sketches and save sketches into " << refSketchOut << " is: " << t3 - t2 << endl;
		#endif
	}

	vector<sketch_t> query_sketches;
	sketchInfo_t query_info;
	if(isSketchFile(queryList))
	{
		querySketchOut = queryList;
		readSketches(query_sketches, query_info, queryList);
		t4 = get_sec();
		#ifdef Timer
		cerr << "===================time of read query sketches from " << queryList << " is: " << t4 - t3 << endl;
		#endif
	}
	else{
		querySketchOut = queryList + ".sketch";
		bool success1;
		if(isFastaList(queryList)){
			success1 = sketchFastaFile(queryList, true, threads, kssd_parameter, query_sketches, query_info, querySketchOut);
		}
		else if(isFastqList(queryList)){
			success1 = sketchFastqFile(queryList, true, threads, kssd_parameter, leastQual, leastNumKmer, query_sketches, query_info, querySketchOut);
		}
		else{
			cerr << "ERROR: command_dist(), the input file list for sketching must be list of fasta and fastq file" << endl;
			exit(1);
		}
		if(!success1){
			cerr << "ERROR: command_dist(), failed sketching query genomes" << endl;
			exit(1);
		}
		
		t4 = get_sec();
		#ifdef Timer
		cerr << "===================time of computing query sketches and save sketches into " << querySketchOut << " is: " << t4 - t3 << endl;
		#endif
	}
	if(query_info.id != ref_info.id){
		fprintf(stderr, "ERROR: command_dist(), the sketch infos between reference and query files are not match\n");
		fprintf(stderr, "try to use the same shuffle file to generate sketches of the reference and query datasets\n");
		exit(1);
	}

	//dist(ref_sketches, query_sketches, outputFile, kmer_size, maxDist, threads);
	index_dist(ref_sketches, ref_info, refSketchOut, query_sketches, outputFile, kmer_size, maxDist, maxNeighbor, isNeighbor, isContainment, threads);
}

void command_union(string sketchFile, string outputFile, int threads){
	vector<sketch_t> sketches;
	if(!isSketchFile(sketchFile)){
		err(errno, "ERROR: command_union, %s is not sketch file, need input sketch file\n", sketchFile.c_str());
	}
	#ifdef Timer_inner
	double t0 = get_sec();
	#endif

	sketchInfo_t query_info;
	FILE* fp_query = fopen(sketchFile.c_str(), "rb");
	if(!fp_query){
		fprintf(stderr, "ERROR: command_union(): cannot open the file: %s\n", sketchFile.c_str());
		exit(1);
	}
	int read_query_info = fread(&query_info, sizeof(sketchInfo_t), 1, fp_query);
	bool use64 = query_info.half_k - query_info.drlevel > 8 ? true : false;
	int query_num = query_info.genomeNumber;
	int progress_bar_size = get_progress_bar_size(query_num);
	int * genomeNameSize = new int[query_num];
	int * hashSetSize = new int[query_num];
	int read_genome_name_size = fread(genomeNameSize, sizeof(int), query_num, fp_query);
	int read_hash_set_size = fread(hashSetSize, sizeof(int), query_num, fp_query);
	if(read_query_info != 1 || read_genome_name_size != query_num || read_hash_set_size != query_num){
		cerr << "ERROR: command_union(), mismatched read query_info, genome_name_size, hash_set_size" << endl;
		exit(1);
	}

	uint64_t remain_memory = get_total_system_memory();
	uint64_t available_memory = remain_memory * 0.7;

	size_t dictSize = use64 ? (1LLU << 4*(query_info.half_k-query_info.drlevel)) / 64 : (1LLU << 32) / 64;
	size_t dict_memory = dictSize * sizeof(uint64_t);
	int available_consumer = available_memory / dict_memory;
	int consumer_num = std::min(available_consumer, threads-1);
	if(use64)
		cerr << "-----use hash64, and the dictSize is: " << dictSize << endl;
	else
		cerr << "-----not use hash64, and the dictSize is: " << dictSize << endl;
	cerr << "the consumer_num is: " << consumer_num << endl;
	uint64_t * dictArr[consumer_num];
	#pragma omp parallel for num_threads(consumer_num)
	for(int i = 0; i < consumer_num; i++){
		//dictArr[i] = (uint64_t*)malloc(dictSize * sizeof(uint64_t));
		dictArr[i] = new uint64_t[dictSize];
		memset(dictArr[i], 0, dictSize * sizeof(uint64_t));
	}

	int maxNameLength = 1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	uint64_t * curPoint64 = new uint64_t[maxHashSize];
	string totalName = sketchFile + " merged sketches";
	vector<uint32_t> mergedArr;
	vector<uint64_t> mergedArr64;

	queue<sketch_t> sketch_queue;
	int index = 0;
	int produced_num = 0;
	int solved_num = 0;
	omp_lock_t queue_lock;
	omp_lock_t vector_lock;
	omp_init_lock(&queue_lock);
	omp_init_lock(&vector_lock);
	sketch_t s_arr[threads];

	cerr << "the total genome number in sketch file is: " << query_num << endl;
	#pragma omp parallel num_threads(consumer_num+1)
	{
		int tid = omp_get_thread_num();
		if(tid == 0)//producer
		{
			for(index = 0; index < query_num; index++){
				int curLength = genomeNameSize[index];
				if(curLength > maxNameLength){
					maxNameLength = curLength;
					curName = new char[maxNameLength+1];
				}
				int nameLength = fread(curName, sizeof(char), curLength, fp_query);
				if(nameLength != curLength){
					cerr << "ERROR: command_union(), mismatched read name length" << endl;
					exit(1);
				}
				string genomeName;
				genomeName.assign(curName, curName + curLength);

				sketch_t s;
				s.fileName = genomeName;
				s.id = index;
				if(use64){
					int curSize = hashSetSize[index];
					if(curSize > maxHashSize){
						maxHashSize = curSize;
						curPoint64 = new uint64_t[maxHashSize];
					}
					int hashSize = fread(curPoint64, sizeof(uint64_t), curSize, fp_query);
					if(hashSize != curSize){
						cerr << "ERROR: command_union(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
						exit(1);
					}
					vector<uint64_t> curHashSet64(curPoint64, curPoint64 + curSize);
					s.hashSet64 = curHashSet64;
				}
				else{
					int curSize = hashSetSize[index];
					if(curSize > maxHashSize){
						maxHashSize = curSize;
						curPoint = new uint32_t[maxHashSize];
					}
					int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp_query);
					if(hashSize != curSize){
						cerr << "ERROR: command_union(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
						exit(1);
					}
					vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
					s.hashSet=curHashSet;
				}
				//#pragma omp critical
				omp_set_lock(&queue_lock);
				{
					sketch_queue.push(s);
					produced_num++;
				}
				omp_unset_lock(&queue_lock);
			}//end for
		}
		else//consumer
		{
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

				sketch_t new_s;
				new_s.fileName = s_arr[tid].fileName;
				new_s.id = s_arr[tid].id;
				if(use64){
					for(size_t j = 0; j < s_arr[tid].hashSet64.size(); j++){
						uint64_t curHash64 = s_arr[tid].hashSet64[j];
						uint64_t dict_id = curHash64 / 64;
						uint64_t dict_offset = 0x8000000000000000LLU >> (curHash64 % 64);
						dictArr[tid-1][dict_id] |= dict_offset;
					}

				}
				else{
					for(size_t j = 0; j < s_arr[tid].hashSet.size(); j++){
						uint32_t curHash = s_arr[tid].hashSet[j];
						uint64_t dict_id = curHash / 64;
						uint64_t dict_offset = 0x8000000000000000LLU >> (curHash % 64);
						dictArr[tid-1][dict_id] |= dict_offset;
					}
				}
				omp_set_lock(&vector_lock);
				solved_num++;
				if(solved_num % progress_bar_size == 0)	cerr << "finished: " << solved_num << endl;
				omp_unset_lock(&vector_lock);
			}
		}
	}

	double tx = get_sec();
	#pragma omp parallel for num_threads(threads) schedule(static)
	for(size_t j = 0; j < dictSize; j++){
		for(int i = 1; i < consumer_num; i++){
			dictArr[0][j] |= dictArr[i][j];
		}
	}
	double ty = get_sec();
	cerr << "the test time is: " << ty - tx << endl;

	#ifdef Timer_inner
	double t1 = get_sec();
		cerr << "-----time of multithreads generate the dictionary is: " << t1 - t0 << endl;
	#endif

	vector<sketch_t> mergedSketches;
	sketch_t s;
	s.id = 0;
	s.fileName = totalName;
	if(use64)
	{
		for(size_t i = 0; i < dictSize; i++){
			if(dictArr[0][i]){
				for(size_t j = 0; j < 64; j++){
					if((0x8000000000000000LLU >> j) & dictArr[0][i]){
						uint64_t hash = 64 * i + j;
						mergedArr64.push_back(hash);
					}
				}
			}
		}
		s.hashSet64 = mergedArr64;
	}
	else
	{
		for(size_t i = 0; i < dictSize; i++){
			if(dictArr[0][i]){
				for(size_t j = 0; j < 64; j++){
					if((0x8000000000000000LLU >> j) & dictArr[0][i]){
						uint32_t hash = 64 * i + j;
						mergedArr.push_back(hash);
					}
				}
			}
		}
		s.hashSet = mergedArr;
	}
	mergedSketches.push_back(s);

	#ifdef Timer_inner
	double t2 = get_sec();
		cerr << "-----time of get hash array by dictionary is: " << t2 - t1 << endl;
	#endif

	query_info.genomeNumber = 1;
	saveSketches(mergedSketches, query_info, outputFile);
	#ifdef Timer_inner
	double t3 = get_sec();
	cerr << "-----time of save sketches is: " << t3 - t2 << endl;
	#endif
	for(int i = 0; i < consumer_num; i++)	delete [] dictArr[i];
	delete [] genomeNameSize;
	delete [] hashSetSize;
	delete [] curPoint;
	delete [] curPoint64;
}

void command_sub(string refSketchFile, string querySketchFile, string outputFile, int threads){
	#ifdef Timer_inner
	double t0 = get_sec();
	#endif
	vector<sketch_t> refSketches;
	if(!isSketchFile(refSketchFile)){
		err(errno, "ERROR: command_sub(), %s is not sketch file, need input sketch file\n", refSketchFile.c_str());
	}
	sketchInfo_t info;
	readSketches(refSketches, info, refSketchFile);

	#ifdef Timer_inner
	double t1 = get_sec();
	cerr << "-----time of read reference sketch: " << refSketchFile << " is: " << t1 - t0 << endl;
	#endif

	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	uint64_t* dict;
	if(use64){
		size_t hash_dimention = 1LLU << (4 * (info.half_k - info.drlevel));
		size_t dictSize = hash_dimention / 64;
		dict = new uint64_t[dictSize];
		memset(dict, 0, dictSize * sizeof(uint64_t));
		for(size_t i = 0; i < refSketches.size(); i++){
			for(size_t j = 0; j < refSketches[i].hashSet64.size(); j++){
				uint64_t curHash = refSketches[i].hashSet64[j];
				dict[curHash/64] |= (0x8000000000000000LLU >> (curHash % 64));
			}
		}
	}
	else{
		size_t dictSize = (1LLU << 32) / 64;
		dict = new uint64_t[dictSize];
		memset(dict, 0, dictSize * sizeof(uint64_t));
		for(size_t i = 0; i < refSketches.size(); i++){
			for(size_t j = 0; j < refSketches[i].hashSet.size(); j++){
				uint32_t curHash = refSketches[i].hashSet[j];
				dict[curHash/64] |= (0x8000000000000000LLU >> (curHash % 64));
			}
		}
	}

	#ifdef Timer_inner
	double t2 = get_sec();
	cerr << "-----time of generate reference dictionary is: " << t2 - t1 << endl;
	#endif

	vector<sketch_t> querySketches;
	if(!isSketchFile(querySketchFile)){
		err(errno, "ERROR: command_sub(): %s is not sketch file, need input sketch file\n", querySketchFile.c_str());
	}
	//readSketches(querySketches, info, querySketchFile);
	sketchInfo_t query_info;
	FILE* fp_query = fopen(querySketchFile.c_str(), "rb");
	if(!fp_query){
		fprintf(stderr, "ERROR: command_sub(): cannot open the file: %s\n", querySketchFile.c_str());
		exit(1);
	}
	int read_query_info = fread(&query_info, sizeof(sketchInfo_t), 1, fp_query);
	if(query_info.id != info.id){
		fprintf(stderr, "ERROR: command_sub(): the sketch infos between subtraction reference and query sketches are not same\n");
		exit(1);
	}
	int query_num = query_info.genomeNumber;
	int * genomeNameSize = new int[query_num];
	int * hashSetSize = new int[query_num];
	int read_genome_name_size = fread(genomeNameSize, sizeof(int), query_num, fp_query);
	int read_hash_set_size = fread(hashSetSize, sizeof(int), query_num, fp_query);
	if(read_query_info != 1 || read_genome_name_size != query_num || read_hash_set_size != query_num){
		cerr << "ERROR: command_sub(), mismatched read query_info, genome_name_size, hash_set_size" << endl;
		exit(1);
	}

	int maxNameLength = 1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	uint64_t * curPoint64 = new uint64_t[maxHashSize];

	FILE* fp_result = fopen(outputFile.c_str(), "w+");
	if(!fp_result){
		fprintf(stderr, "ERROR: command_sub(): cannot open the file: %s\n", outputFile.c_str());
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
	int progress_bar_size = get_progress_bar_size(query_num);
	cerr << "=====total: " << query_num << endl;
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
				if(nameLength != curLength){
					cerr << "ERROR: command_sub(), mismatch read name length" << endl;
					exit(1);
				}
				string genomeName;
				genomeName.assign(curName, curName + curLength);

				sketch_t s;
				s.fileName = genomeName;
				s.id = index;
				if(use64){
					int curSize = hashSetSize[index];
					if(curSize > maxHashSize){
						maxHashSize = curSize;
						curPoint64 = new uint64_t[maxHashSize];
					}
					int hashSize = fread(curPoint64, sizeof(uint64_t), curSize, fp_query);
					if(hashSize != curSize){
						cerr << "ERROR: command_sub(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
						exit(1);
					}
					vector<uint64_t> curHashSet64(curPoint64, curPoint64 + curSize);
					s.hashSet64 = curHashSet64;
				}
				else{
					int curSize = hashSetSize[index];
					if(curSize > maxHashSize){
						maxHashSize = curSize;
						curPoint = new uint32_t[maxHashSize];
					}
					int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp_query);
					if(hashSize != curSize){
						cerr << "ERROR: command_sub(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
						exit(1);
					}
					vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
					s.hashSet=curHashSet;
				}
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
						omp_unset_lock(&queue_lock);
						if(solved_num >= query_num)	break;
						continue;
					}
					s_arr[tid] = sketch_queue.front();
					sketch_queue.pop();
				}
				omp_unset_lock(&queue_lock);
				sketch_t new_s;
				new_s.fileName = s_arr[tid].fileName;
				new_s.id = s_arr[tid].id;

				if(use64){
					vector<uint64_t> newHashArr64;
					for(size_t j = 0; j < s_arr[tid].hashSet64.size(); j++){
						uint64_t curHash = s_arr[tid].hashSet64[j];
						uint64_t isIn = dict[curHash/64] & (0x8000000000000000 >> (curHash % 64));
						if(!isIn){
							newHashArr64.emplace_back(curHash);
						}
					}
					new_s.hashSet64 = newHashArr64;
				}
				else{
					vector<uint32_t> newHashArr;
					for(size_t j = 0; j < s_arr[tid].hashSet.size(); j++){
						uint32_t curHash = s_arr[tid].hashSet[j];
						uint64_t isIn = dict[curHash/64] & (0x8000000000000000 >> (curHash % 64));
						if(!isIn){
							newHashArr.emplace_back(curHash);
						}
					}
					new_s.hashSet = newHashArr;
				}
				omp_set_lock(&vector_lock);
				{
					const char * namePoint = new_s.fileName.c_str();
					int curNameLength = new_s.fileName.length();
					out_genomeNameSize[solved_num] = curNameLength;
					fwrite(namePoint, sizeof(char), curNameLength, fp_result);
					int curHashSetSize;
					if(use64){
						uint64_t * curPoint64 = new_s.hashSet64.data();
						curHashSetSize = new_s.hashSet64.size();
						fwrite(curPoint64, sizeof(uint64_t), curHashSetSize, fp_result);
					}
					else{
						uint32_t * curPoint = new_s.hashSet.data();
						curHashSetSize = new_s.hashSet.size();
						fwrite(curPoint, sizeof(uint32_t), curHashSetSize, fp_result);
					}
					out_hashSetSize[solved_num] = curHashSetSize;
					solved_num++;
					if(solved_num % progress_bar_size == 0) cerr << "=====finished: " << solved_num << endl;
				}
				omp_unset_lock(&vector_lock);
			}
		}
	}

	rewind(fp_result);
	fwrite(&query_info, sizeof(sketchInfo_t), 1, fp_result);
	fwrite(out_genomeNameSize, sizeof(int), query_num, fp_result);
	fwrite(out_hashSetSize, sizeof(int), query_num, fp_result);

	delete [] curName;
	delete [] curPoint;
	delete [] curPoint64;
	delete [] dict;
	delete [] genomeNameSize;
	delete [] hashSetSize;
	delete [] out_genomeNameSize;
	delete [] out_hashSetSize;

	#ifdef Timer_inner
	double t3 = get_sec();
	cerr << "-----time of substraction is: " << t3 - t2 << endl;
	#endif
}

void command_merge(string inputList, string outputFile, int threads){
	ifstream ifs(inputList);
	string line;
	vector<string> fileList;
	while(getline(ifs, line)){
		if(!isSketchFile(line)){
			cerr << "ERROR: command_merge(), the file: " << line << " is not a sketch file in the list file: " << inputList << endl;
			exit(1);
		}
		fileList.push_back(line);
	}
	ifs.close();

	int totalSketchNumber = 0;
	vector<int> totalGenomeNameSize;
	vector<int> totalHashSetSize;
	sketchInfo_t resInfo;
	bool use64;

	vector<sketch_t> resSketches;
	FILE * fp0 = fopen(fileList[0].c_str(), "rb");
	if(!fp0){
		fprintf(stderr, "ERROR: command_merge(), cannot open file: %s\n", fileList[0].c_str());
		exit(1);
	}
	int read_res_info0 = fread(&resInfo, sizeof(sketchInfo_t), 1, fp0);
	use64 = resInfo.half_k - resInfo.drlevel > 8 ? true : false;
	int sketchNumber0 = resInfo.genomeNumber;
	int * genomeNameSize0 = new int[sketchNumber0];
	int * hashSetSize0 = new int[sketchNumber0];
	int read_genome_name_size0 = fread(genomeNameSize0, sizeof(int), sketchNumber0, fp0);
	int read_hash_set_size0 = fread(hashSetSize0, sizeof(int), sketchNumber0, fp0);
	totalSketchNumber += sketchNumber0;
	for(int j = 0; j < sketchNumber0; j++){
		totalGenomeNameSize.emplace_back(genomeNameSize0[j]);
		totalHashSetSize.emplace_back(hashSetSize0[j]);
	}
	fclose(fp0);
	if(read_res_info0 != 1 || read_genome_name_size0 != sketchNumber0 || read_hash_set_size0 != sketchNumber0){
		cerr << "ERROR: command_merge(), mismatched read res_info0, genome_name_size0, hash_set_size0" << endl;
		exit(1);
	}

	for(size_t i = 1; i < fileList.size(); i++){
		FILE * fpTmp = fopen(fileList[i].c_str(), "rb");
		sketchInfo_t curInfo;
		int read_cur_info = fread(&curInfo, sizeof(sketchInfo_t), 1, fpTmp);
		if(curInfo.id != resInfo.id){
			fprintf(stderr, "ERROR: command_merge(), mismatched sketch info\n");
			exit(1);
		}
		int curSketchNumber = curInfo.genomeNumber;
		int * genomeNameSize = new int[curSketchNumber];
		int * hashSetSize = new int[curSketchNumber];
		int read_genome_name_size = fread(genomeNameSize, sizeof(int), curSketchNumber, fpTmp);
		int read_hash_set_size = fread(hashSetSize, sizeof(int), curSketchNumber, fpTmp);
		totalSketchNumber += curSketchNumber;
		for(int j = 0; j < curSketchNumber; j++){
			totalGenomeNameSize.emplace_back(genomeNameSize[j]);
			totalHashSetSize.emplace_back(hashSetSize[j]);
		}
		fclose(fpTmp);
		if(read_cur_info != 1 || read_genome_name_size != curSketchNumber || read_hash_set_size != curSketchNumber){
			cerr << "ERROR: command_merge(), mismatched read cur_info, genome_name_size, hash_set_size" << endl;
			exit(1);
		}
	}

	FILE* fp1 = fopen(outputFile.c_str(), "w+");
	resInfo.genomeNumber = totalSketchNumber;
	fwrite(&resInfo, sizeof(sketchInfo_t), 1, fp1);
	fwrite(totalGenomeNameSize.data(), sizeof(int), totalSketchNumber, fp1);
	fwrite(totalHashSetSize.data(), sizeof(int), totalSketchNumber, fp1);
	fclose(fp1);
	
	FILE* fp = fopen(outputFile.c_str(), "a+");
	for(size_t i = 0; i < fileList.size(); i++){
		vector<sketch_t> curSketches;
		sketchInfo_t curInfo;
		readSketches(curSketches, curInfo, fileList[i]);

		size_t sketchNumber = curSketches.size();
		for(size_t i = 0; i < sketchNumber; i++){
			const char * namePoint = curSketches[i].fileName.c_str();
			fwrite(namePoint, sizeof(char),curSketches[i].fileName.length(), fp);
			if(use64){
				uint64_t * curPoint64 = curSketches[i].hashSet64.data();
				fwrite(curPoint64, sizeof(uint64_t),curSketches[i].hashSet64.size(), fp);
			}
			else{
				uint32_t * curPoint = curSketches[i].hashSet.data();
				fwrite(curPoint, sizeof(uint32_t),curSketches[i].hashSet.size(), fp);
			}
		}
	}
	fclose(fp);
}
	








