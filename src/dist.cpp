#include "dist.h"
#include <omp.h>
#include <cstdio>
#include <math.h>
#include <sys/stat.h>
#include <immintrin.h>
#include <cstring>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <atomic>
#include <stdint.h>
#include <queue>
#include <algorithm>
#include "robin_hood.h"

int u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3);
int u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t &a, uint64_t &b);
setResult_t getJaccard(vector<uint32_t> list1, vector<uint32_t> list2);
setResult_t vgetJaccard(vector<uint32_t> list1, vector<uint32_t> list2);

void index_tridist(vector<sketch_t>& sketches, sketchInfo_t& info, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads){

	#ifdef Timer
	double t0 = get_sec();
	#endif
	string indexFile = refSketchOut + ".index";
	string dictFile = refSketchOut + ".dict";
	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
	uint32_t* sketchSizeArr = NULL;
	size_t* offset = NULL;
	uint32_t* indexArr = NULL;
	//uint64_t* dict;
	if(use64)
	{
		cerr << "-----use hash64 in index_tridist() " << endl;
		size_t hash_number;
		FILE* fp_index = fopen(indexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: index_tridist(), cannot open index file: " << indexFile << endl;
			exit(1);
		}
		int read_hash_num = fread(&hash_number, sizeof(size_t), 1, fp_index);
		uint64_t * hash_arr = new uint64_t[hash_number];
		uint32_t * hash_size_arr = new uint32_t[hash_number];
		size_t read_hash_arr = fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		size_t read_hash_size_arr = fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		if(read_hash_num != 1 || read_hash_arr != hash_number || read_hash_size_arr != hash_number){
			cerr << "ERROR: index_tridist(), error read hash_number, hash_arr, and hash_size_arr" << endl;
			exit(1);
		}
		
		fclose(fp_index);

		FILE* fp_dict = fopen(dictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: index_tridist(), cannot open dict file: " << dictFile << endl;
			exit(1);
		}
		uint32_t max_hash_size = 1LLU << 24;
		uint32_t* cur_point = new uint32_t[max_hash_size];
		for(size_t i = 0; i < hash_number; i++){
			uint32_t cur_hash_size = hash_size_arr[i];
			if(cur_hash_size > max_hash_size){
				max_hash_size = cur_hash_size;
				cur_point = new uint32_t[max_hash_size];
			}
			uint32_t hash_size = fread(cur_point, sizeof(uint32_t), cur_hash_size, fp_dict);
			if(hash_size != cur_hash_size){
				cerr << "ERROR: index_tridist(), the read hash number is not equal to the saved hash number information" << endl;
				exit(1);
			}
			vector<uint32_t> cur_genome_arr(cur_point, cur_point + cur_hash_size);
			uint64_t cur_hash = hash_arr[i];
			hash_map_arr.insert({cur_hash, cur_genome_arr});
		}
		delete [] cur_point;
		delete [] hash_arr;
		delete [] hash_size_arr;
		fclose(fp_dict);
	}
	else
	{
		cerr << "-----not use hash64 in index_tridist() " << endl;
		size_t hashSize;
		uint64_t totalIndex;
		FILE * fp_index = fopen(indexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: index_tridist(), cannot open the index sketch file: " << indexFile << endl;
			exit(1);
		}
		int read_hash_size = fread(&hashSize, sizeof(size_t), 1, fp_index);
		int read_total_index = fread(&totalIndex, sizeof(uint64_t), 1, fp_index);
		//sketchSizeArr = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
		sketchSizeArr = new uint32_t[hashSize];
		size_t read_sketch_size_arr = fread(sketchSizeArr, sizeof(uint32_t), hashSize, fp_index);

		//offset = (size_t*)malloc(hashSize * sizeof(size_t));
		offset = new size_t[hashSize];
		uint64_t totalHashNumber = 0;
		for(size_t i = 0; i < hashSize; i++){
			totalHashNumber += sketchSizeArr[i];
			offset[i] = sketchSizeArr[i];
			if(i > 0) offset[i] += offset[i-1];
		}
		if(totalHashNumber != totalIndex){
			cerr << "ERROR: index_tridist(), mismatched total hash number" << endl;
			exit(1);
		}
		fclose(fp_index);

		//cerr << "the hashSize is: " << hashSize << endl;
		//cerr << "totalIndex is: " << totalIndex << endl;
		//cerr << "totalHashNumber is: " << totalHashNumber << endl;
		//cerr << "offset[n-1] is: " << offset[hashSize-1] << endl;;

		//indexArr = (uint32_t*)malloc(totalHashNumber * sizeof(uint32_t));
		indexArr = new uint32_t[totalHashNumber];
		FILE * fp_dict = fopen(dictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: index_tridist(), cannot open the dictionary sketch file: " << dictFile << endl;
			exit(1);
		}
		size_t read_index_arr = fread(indexArr, sizeof(uint32_t), totalHashNumber, fp_dict);
		if(read_hash_size != 1 || read_total_index != 1 || read_sketch_size_arr != hashSize || read_index_arr != totalHashNumber){
			cerr << "ERROR: index_tridist(), error read hash_size, total_index, sketch_size_arr, index_arr" << endl;
			exit(1);
		}
	}

	#ifdef Timer
	double t1 = get_sec();
	cerr << "===================time of read index and offset sketch file is: " << t1 - t0 << endl;
	#endif
	size_t numRef = sketches.size();

	vector<FILE*> fpArr;
	vector<FILE*> fpIndexArr;
	vector<string> dist_file_list;
	vector<string> dist_index_list;
	//vector<int*> intersectionArr;
	int** intersectionArr = new int*[numThreads];
	
	string folderPath = outputFile + ".dir";
	string command0 = "mkdir -p " + folderPath;
	int status = system(command0.c_str());
	if(!status){
		cerr << "success create: " << folderPath << endl;
	}

	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = folderPath + '/' + outputFile + '.' + to_string(i);
		dist_file_list.push_back(tmpName);
		FILE * fp0 = fopen(tmpName.c_str(), "w+");
		fpArr.push_back(fp0);

		string tmpIndexName = outputFile + ".index." + to_string(i);
		dist_index_list.push_back(tmpIndexName);
		FILE * fp1 = fopen(tmpIndexName.c_str(), "w+");
		fpIndexArr.push_back(fp1);

		//int * arr = (int*)malloc(numRef * sizeof(int));
		//int* arr = new int[numRef];
		//intersectionArr.push_back(arr);
		intersectionArr[i] = new int[numRef];
	}
	
	//cerr << "before generate the intersection " << endl;
	
	int progress_bar_size = get_progress_bar_size(numRef);
	cerr << "=====total: " << numRef << endl;
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < numRef; i++){
		if(i % progress_bar_size == 0) cerr << "=====finish: " << i << endl;
		int tid = omp_get_thread_num();
		fprintf(fpIndexArr[tid], "%s\t%s\n", sketches[i].fileName.c_str(), dist_file_list[tid].c_str());
		memset(intersectionArr[tid], 0, numRef * sizeof(int));
		if(use64){
			for(size_t j = 0; j < sketches[i].hashSet64.size(); j++){
				uint64_t hash64 = sketches[i].hashSet64[j];
				//if(!(dict[hash64/64] & (0x8000000000000000LLU >> (hash64 % 64))))	continue;
				if(hash_map_arr.count(hash64) == 0) continue;
				//for(auto x : hash_map_arr[hash64])
				for(size_t k = 0; k < hash_map_arr[hash64].size(); k++){
					size_t cur_index = hash_map_arr[hash64][k];
					intersectionArr[tid][cur_index]++;
					//cerr << hash64 << '\t' << cur_index << endl;
				}
			}
		}
		else{
			for(size_t j = 0; j < sketches[i].hashSet.size(); j++){
				uint32_t hash = sketches[i].hashSet[j];
				if(sketchSizeArr[hash] == 0) continue;
				size_t start = hash > 0 ? offset[hash-1] : 0;
				size_t end = offset[hash];
				for(size_t k = start; k < end; k++){
					size_t curIndex = indexArr[k];
					intersectionArr[tid][curIndex]++;
				}
			}
		}

		string strBuf("");
		for(size_t j = i+1; j < numRef; j++){
			int common = intersectionArr[tid][j];
			int size0, size1;
			if(use64){
				size0 = sketches[i].hashSet64.size();
				size1 = sketches[j].hashSet64.size();
			}
			else{
				size0 = sketches[i].hashSet.size();
				size1 = sketches[j].hashSet.size();
			}
			if(!isContainment){
				int denom = size0 + size1 - common;
				double jaccard;
				if(size0 == 0 || size1 == 0)
					jaccard = 0.0;
				else
					jaccard = (double)common / denom;
				double mashD;
				if(jaccard == 1.0)
					mashD = 0.0;
				else if(jaccard == 0.0)
					mashD = 1.0;
				else
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				if(mashD < maxDist){
					strBuf += sketches[j].fileName + '\t' + sketches[i].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(jaccard) + '\t' + to_string(mashD) + '\n';
				}
					//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[j].fileName.c_str(), sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
			else{
				int denom = std::min(size0, size1);
				double containment;
				if(size0 == 0 || size1 == 0)
					containment = 0.0;
				else
					containment = (double)common / denom;
				double AafD;
				if(containment == 1.0)
					AafD = 0.0;
				else if(containment == 0.0)
					AafD = 1.0;
				else
					AafD = (double)-1.0 / kmer_size * log(containment);
				if(AafD < maxDist)
					strBuf += sketches[j].fileName + '\t' + sketches[i].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(containment) + '\t' + to_string(AafD) + '\n';
					//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[j].fileName.c_str(), sketches[i].fileName.c_str(), common, size0, size1, containment, AafD);
			}
		}
		fprintf(fpArr[tid], "%s", strBuf.c_str());
		strBuf = "";
	}


	//cerr << "finished multithread computing" << endl;

	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
		fclose(fpIndexArr[i]);
		delete [] intersectionArr[i];
	}
	//cerr << "finished fpArr fclose" << endl;

	#ifdef Timer
	double t2 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t2 - t1 << endl;
	#endif

	uint64_t totalSize = 0;
	uint64_t maxSize = 1LLU << 32; //max total distance file size 4GB
	bool isMerge = false;
	for(int i = 0; i < numThreads; i++){
		struct stat cur_stat;
		stat(dist_file_list[i].c_str(), &cur_stat);
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
	}
	if(totalSize <= maxSize)	isMerge = true;

	if(isMerge){
		FILE * cofp;
		FILE * com_cofp = fopen(outputFile.c_str(), "w");
		cerr << "-----save the output distance file: " << outputFile << endl;
		fprintf(com_cofp, " genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++)
		{
			cofp = fopen(dist_file_list[i].c_str(), "rb+");
			while(1)
			{
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
				//cerr << "the lengthRead is: " << lengthRead << endl;
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp);
			remove(dist_file_list[i].c_str());
			remove(dist_index_list[i].c_str());
		}
		remove(folderPath.c_str());

		free(bufRead);
		fclose(com_cofp);
	}
	else{
		cerr << "-----the output distance file is too big to merge into one single file, saving the result into directory: " << folderPath << endl;
		FILE * cofp1;
		string outputIndexFile = outputFile + ".index";
		cerr << "-----save the index between genomes and distance sub-files into: " << outputIndexFile << endl;
		FILE * com_cofp1 = fopen(outputIndexFile.c_str(), "w+");
		fprintf(com_cofp1, "genomeName\tdistFileName\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++){
			cofp1 = fopen(dist_index_list[i].c_str(), "rb+");
			while(1){
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp1);
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp1);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp1);
			remove(dist_index_list[i].c_str());
		}
		free(bufRead);
		fclose(com_cofp1);
	}

	#ifdef Timer
	double t3 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t3 - t2 << endl;
	#endif

}

void tri_dist(vector<sketch_t>& sketches, string outputFile, int kmer_size, double maxDist, int numThreads){
	
	double t00 = get_sec();
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < sketches.size(); i++){
		std::sort(sketches[i].hashSet.begin(), sketches[i].hashSet.end());
	}

	double t0 = get_sec();
	cerr << "time of multithreading sort sketches hash values is: " << t0 - t00 << endl;
	vector<string> dist_file_list;
	vector<FILE*> fpArr;
	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = outputFile + to_string(i);
		dist_file_list.push_back(tmpName);

		FILE * fp0 = fopen(tmpName.c_str(), "w");
		fpArr.push_back(fp0);
	}
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < sketches.size(); i++)
	{
		//if(i % 300 == 0) cerr << "finish: " << i << endl;
		int tid = omp_get_thread_num();
		for(size_t j = i+1; j < sketches.size(); j++)
		{
			//setResult_t tmpResult = getJaccard(sketches[i].hashSet, sketches[j].hashSet);
			setResult_t tmpResult = vgetJaccard(sketches[i].hashSet, sketches[j].hashSet);
			double jaccard = tmpResult.jaccard;
			int common = tmpResult.common;
			int size0 = tmpResult.size0;
			int size1 = tmpResult.size1;
			double mashD;
			if(jaccard == 1.0) 
				mashD = 0.0;
			else if(jaccard == 0.0) 
				mashD = 1.0;
			else 
				mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
			//fprintf(fp0, "jaccard between[%d] and [%d] is: %lf\n", i, j, jaccard);
			if(mashD < maxDist)
				fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[j].fileName.c_str(), sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
		}
	}
	//close the subResultFile
	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
	}
	double t1 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

	//struct stat cof_stat;
	FILE * cofp;
	FILE * com_cofp = fopen(outputFile.c_str(), "w");
	//fprintf(fpArr[tid], " %s\t%s\t%d\t%d\t%lf\t%lf\n", sketches[i].fileName.c_str(), sketches[j].fileName.c_str(), common, denom, jaccard, mashD);
	fprintf(com_cofp, " genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");
	int bufSize = 1 << 24;
	int lengthRead = 0;
	char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
	for(int i = 0; i < numThreads; i++)
	{
		cofp = fopen(dist_file_list[i].c_str(), "rb+");

		while(1)
		{
			lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
			//cerr << "the lengthRead is: " << lengthRead << endl;
			fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
			if(lengthRead < bufSize) break;
		}
		fclose(cofp);

		remove(dist_file_list[i].c_str());
	}

	free(bufRead);
	fclose(com_cofp);
	double t2 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t2 - t1 << endl;

}

void index_dist(vector<sketch_t>& ref_sketches, sketchInfo_t& ref_info, string refSketchOut, vector<sketch_t>& query_sketches, string outputFile, int kmer_size, double maxDist, uint64_t maxNeighbor, bool isNeighbor, int isContainment, int numThreads){

	#ifdef Timer
	double t0 = get_sec();
	#endif
	string refIndexFile = refSketchOut + ".index";
	string refDictFile = refSketchOut + ".dict";
	bool use64 = ref_info.half_k - ref_info.drlevel > 8 ? true : false;
	robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
	//uint64_t * dict;
	uint32_t* refSketchSizeArr = NULL;
	size_t* refOffset = NULL;
	int* refIndexArr = NULL;
	if(use64){
		cerr << "use 64 in index_dist() " << endl;
		size_t hash_number;
		FILE* fp_index = fopen(refIndexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: index_dist(), cannot open index file: " << refIndexFile << endl;
			exit(1);
		}
		int read_hash_number = fread(&hash_number, sizeof(size_t), 1, fp_index);
		uint64_t* hash_arr = new uint64_t[hash_number];
		uint32_t* hash_size_arr = new uint32_t[hash_number];
		size_t read_hash_arr = fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		size_t read_hash_size_arr = fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		if(read_hash_number != 1 || read_hash_arr != hash_number || read_hash_size_arr != hash_number){
			cerr << "ERROR: index_dist(), error read hash_number, hash_arr, hash_size_arr" << endl;
			exit(1);
		}

		FILE* fp_dict = fopen(refDictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: index_dist(), cannot open dict file: " << refDictFile << endl;
			exit(1);
		}
		uint32_t max_hash_size = 1LLU << 24;
		uint32_t* cur_point = new uint32_t[max_hash_size];
		for(size_t i = 0; i < hash_number; i++){
			uint32_t cur_hash_size = hash_size_arr[i];
			if(cur_hash_size > max_hash_size){
				max_hash_size = cur_hash_size;
				cur_point = new uint32_t[max_hash_size];
			}
			uint32_t hash_size = fread(cur_point, sizeof(uint32_t), cur_hash_size, fp_dict);
			if(hash_size != cur_hash_size){
				cerr << "ERROR: index_dist(), the read hash number is not equal to the saved hash number information" << endl;
				exit(1);
			}
			vector<uint32_t> cur_genome_arr(cur_point, cur_point + cur_hash_size);
			uint64_t cur_hash = hash_arr[i];
			hash_map_arr.insert({cur_hash, cur_genome_arr});
		}
		delete [] cur_point;
		delete [] hash_arr;
		delete [] hash_size_arr;
		fclose(fp_dict);
	}
	else
	{
		size_t refHashSize;
		uint64_t refTotalIndex;
		FILE* fp_index = fopen(refIndexFile.c_str(), "rb");
		int read_ref_hash_size = fread(&refHashSize, sizeof(size_t), 1, fp_index);
		int read_ref_total_index = fread(&refTotalIndex, sizeof(uint64_t), 1, fp_index);
		refSketchSizeArr = new uint32_t[refHashSize];
		size_t read_ref_sketch_size_arr = fread(refSketchSizeArr, sizeof(uint32_t), refHashSize, fp_index);
		fclose(fp_index);

		refOffset = new size_t[refHashSize];
		uint64_t refTotalHashNumber = 0;
		for(size_t i = 0; i < refHashSize; i++){
			refTotalHashNumber += refSketchSizeArr[i];
			refOffset[i] = refSketchSizeArr[i];
			if(i > 0) refOffset[i] += refOffset[i-1];
		}
		if(refTotalHashNumber != refTotalIndex){
			cerr << "ERROR: index_dist(), mismatched the total hash number of refSketch" << endl;
			exit(1);
		}
		//cerr << "the refHashSize is: " << refHashSize << endl;
		//cerr << "the refTotalIndex is: " << refTotalIndex << endl;
		//cerr << "refTotalHashNumber is: " << refTotalHashNumber << endl;
		//cerr << "the refOffset[n-1] is: " << refOffset[refHashSize-1] << endl;
		
		refIndexArr = new int[refTotalHashNumber];
		FILE* fp_dict = fopen(refDictFile.c_str(), "rb");
		size_t read_ref_index_arr = fread(refIndexArr, sizeof(int), refTotalHashNumber, fp_dict);
		fclose(fp_dict);
		if(read_ref_hash_size != 1 || read_ref_total_index != 1 || read_ref_sketch_size_arr != refHashSize || read_ref_index_arr != refTotalHashNumber){
			cerr << "ERROR: index_dist(), error read ref_hash_size, ref_total_index, ref_sketch_size_arr, ref_index_arr" << endl;
			exit(1);
		}
	}
	#ifdef Timer
	double t1 = get_sec();
	cerr << "===================time of read the index and dict sketch file is: " << t1 - t0 << endl;
	#endif
	size_t numRef = ref_sketches.size();
	size_t numQuery = query_sketches.size();

	string folderPath = outputFile + ".dir";
	string command0 = "mkdir -p " + folderPath;
	int status = system(command0.c_str());
	if(!status){
		cerr << "success create: " << folderPath << endl;
	}

	vector<FILE*> fpArr;
	vector<FILE*> fpIndexArr;
	vector<string> dist_file_list;
	vector<string> dist_index_list;
	int** intersectionArr = new int*[numThreads];
	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = folderPath + '/' + outputFile + '.' + to_string(i);
		dist_file_list.push_back(tmpName);
		FILE * fp0 = fopen(tmpName.c_str(), "w+");
		fpArr.push_back(fp0);
		string tmpIndexName = outputFile + ".index." + to_string(i);
		dist_index_list.push_back(tmpIndexName);
		FILE * fp1 = fopen(tmpIndexName.c_str(), "w+");
		fpIndexArr.push_back(fp1);

		intersectionArr[i] = new int[numRef];
	}
	//cerr << "before generate the intersection " << endl;

	int progress_bar_size = get_progress_bar_size(numQuery);
	cerr << "=====total: " << numQuery << endl;
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < numQuery; i++){
		if(i % progress_bar_size == 0) cerr << "=====finish: " << i << endl;
		int tid = omp_get_thread_num();
		fprintf(fpIndexArr[tid], "%s\t%s\n", query_sketches[i].fileName.c_str(), dist_file_list[tid].c_str());
		memset(intersectionArr[tid], 0, numRef *sizeof(int));
		if(use64){
			for(size_t j = 0; j < query_sketches[i].hashSet64.size(); j++){
				uint64_t hash64 = query_sketches[i].hashSet64[j];
				//if(!(dict[hash64/64] & (0x8000000000000000LLU >> (hash64 % 64))))	continue;
				if(hash_map_arr.count(hash64) == 0) continue;
				for(size_t k = 0; k < hash_map_arr[hash64].size(); k++){
					size_t cur_index = hash_map_arr[hash64][k];
					intersectionArr[tid][cur_index]++;
				}
			}
		}
		else{
			for(size_t j = 0; j < query_sketches[i].hashSet.size(); j++){
				uint32_t hash = query_sketches[i].hashSet[j];
				if(refSketchSizeArr[hash] == 0) continue;
				size_t start = hash > 0 ? refOffset[hash-1] : 0;
				size_t end = refOffset[hash];
				for(size_t k = start; k < end; k++){
					size_t curIndex = refIndexArr[k];
					intersectionArr[tid][curIndex]++;
				}
			}
		}
		
		//double nearDist = 1.0;
		//int nearCommon = 0;
		//int nearRSize = 0;
		//int nearRid = 0;
		//double nearJaccard = 0.0;
		//double nearContainment = 0.0;

		string strBuf("");
		int size0, size1;
		priority_queue<DistInfo, vector<DistInfo>, cmpDistInfo> curQueue;
		for(size_t j = 0; j < numRef; j++){
			int common = intersectionArr[tid][j];
			if(use64){
				size0 = ref_sketches[j].hashSet64.size();
				size1 = query_sketches[i].hashSet64.size();
			}
			else{
				size0 = ref_sketches[j].hashSet.size();
				size1 = query_sketches[i].hashSet.size();
			}
			if(!isContainment){
				int denom = size0 + size1 - common;
				double jaccard;
				if(size0 == 0 || size1 == 0)
					jaccard = 0.0;
				else
					jaccard = (double)common / denom;
				double mashD;
				if(jaccard == 1.0)
					mashD = 0.0;
				else if(jaccard == 0.0)
					mashD = 1.0;
				else
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				if(mashD <= maxDist){
					if(isNeighbor){
						DistInfo t0;
						t0.refName = ref_sketches[j].fileName;
						t0.common = common;
						t0.refSize = size0;
						t0.jorc = jaccard;
						t0.dist = mashD;

						if(curQueue.size() < maxNeighbor){
							curQueue.push(t0);
						}
						else if(mashD < curQueue.top().dist){
							curQueue.push(t0);
							curQueue.pop();
						}
					}
					else{
						strBuf += query_sketches[i].fileName + '\t' + ref_sketches[j].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(jaccard) + '\t' + to_string(mashD) + '\n';

					}
				}
			}
			else{
				int denom = std::min(size0, size1);
				double containment;
				if(size0 == 0 || size1 == 0)
					containment = 0.0;
				else
					containment = (double)common / denom;
				double AafD;
				if(containment == 1.0)
					AafD = 0.0;
				else if(containment == 0.0)
					AafD = 1.0;
				else
					AafD = (double)-1.0 / kmer_size * log(containment);
				if(AafD <= maxDist){
					if(isNeighbor){
						DistInfo t0;
						t0.refName = ref_sketches[j].fileName;
						t0.common = common;
						t0.refSize = size0;
						t0.jorc = containment;
						t0.dist = AafD;
						if(curQueue.size() < maxNeighbor){
							curQueue.push(t0);
						}
						else if(AafD < curQueue.top().dist){
							curQueue.push(t0);
							curQueue.pop();
						}
					}
					else{
						strBuf += query_sketches[i].fileName + '\t' + ref_sketches[j].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(containment) + '\t' + to_string(AafD) + '\n';
					}
				}
			}
		}
		if(isNeighbor){
			while(!curQueue.empty()){
				DistInfo t0 = curQueue.top();
				strBuf += query_sketches[i].fileName + '\t' + t0.refName + '\t' + to_string(t0.common) + '|' + to_string(t0.refSize) + '|' + to_string(size1) + '\t' + to_string(t0.jorc) + '\t' + to_string(t0.dist) + '\n';
				curQueue.pop();
			}
		}
		fprintf(fpArr[tid], "%s", strBuf.c_str());
		strBuf = "";
	}//end this query genome

	//cerr << "finished multithread computing" << endl;

	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
		fclose(fpIndexArr[i]);
		delete [] intersectionArr[i];
	}

	//cerr << "finished fpArr fclose" << endl;

	#ifdef Timer
	double t2 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t2 - t1 << endl;
	#endif

	uint64_t totalSize = 0;
	uint64_t maxSize = 1LLU << 32; //max total distance file size 4GB
	bool isMerge = false;
	for(int i = 0; i < numThreads; i++){
		struct stat cur_stat;
		stat(dist_file_list[i].c_str(), &cur_stat);
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
	}
	if(totalSize <= maxSize)	isMerge = true;

	if(isMerge){
		FILE * cofp;
		FILE * com_cofp = fopen(outputFile.c_str(), "w");
		cerr << "-----save the output distance file: " << outputFile << endl;
		fprintf(com_cofp, " genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++)
		{
			cofp = fopen(dist_file_list[i].c_str(), "rb+");
			while(1)
			{
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
				//cerr << "the lengthRead is: " << lengthRead << endl;
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp);
			remove(dist_file_list[i].c_str());
			remove(dist_index_list[i].c_str());
		}
		remove(folderPath.c_str());

		free(bufRead);
		fclose(com_cofp);
	}
	else{
		cerr << "-----the output distance file is too big to merge into one single file, saving the result into directory: " << folderPath << endl;
		FILE * cofp1;
		string outputIndexFile = outputFile + ".index";
		cerr << "-----save the index between genomes and distance sub-files into: " << outputIndexFile << endl;
		FILE * com_cofp1 = fopen(outputIndexFile.c_str(), "w+");
		fprintf(com_cofp1, "genomeName\tdistFileName\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++){
			cofp1 = fopen(dist_index_list[i].c_str(), "rb+");
			while(1){
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp1);
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp1);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp1);
			remove(dist_index_list[i].c_str());
		}
		free(bufRead);
		fclose(com_cofp1);
	}
	#ifdef Timer
	double t3 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t3 - t2 << endl;
	#endif

}

void dist(vector<sketch_t>& ref_sketches, vector<sketch_t>& query_sketches, string outputFile, int kmer_size, double maxDist, int numThreads){

	double t00 = get_sec();
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < ref_sketches.size(); i++){
		std::sort(ref_sketches[i].hashSet.begin(), ref_sketches[i].hashSet.end());
	}
	
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < query_sketches.size(); i++){
		std::sort(query_sketches[i].hashSet.begin(), query_sketches[i].hashSet.end());
	}
	double t0 = get_sec();
	cerr << "the time of sort reference and query sketch is: " << t0 - t00 << endl;
	vector<string> dist_file_list;
	vector<FILE*> fpArr;
	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = outputFile + to_string(i);
		dist_file_list.push_back(tmpName);

		FILE * fp0 = fopen(tmpName.c_str(), "w");
		fpArr.push_back(fp0);
	}

	int refSize = ref_sketches.size();
	int querySize = query_sketches.size();
	if(refSize >= querySize){
		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(int i = 0; i < refSize; i++){
			int tid = omp_get_thread_num();
			for(int j = 0; j < querySize; j++)
			{
				//setResult_t tmpResult = getJaccard(ref_sketches[i].hashSet, query_sketches[j].hashSet);
				setResult_t tmpResult = vgetJaccard(ref_sketches[i].hashSet, query_sketches[j].hashSet);
				double jaccard = tmpResult.jaccard;
				int common = tmpResult.common;
				int size0 = tmpResult.size0;
				int size1 = tmpResult.size1;
				double mashD;
				if(jaccard == 1.0) 
					mashD = 0.0;
				else if(jaccard == 0.0) 
					mashD = 1.0;
				else 
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				//fprintf(fp0, "jaccard between[%d] and [%d] is: %lf\n", i, j, jaccard);
				if(mashD < maxDist)
					//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", ref_sketches[i].fileName.c_str(), query_sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", query_sketches[j].fileName.c_str(), ref_sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
		}
	}
	else{
		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(int i = 0; i < querySize; i++)
		{
			int tid = omp_get_thread_num();
			for(int j = 0; j < refSize; j++)
			{
				setResult_t tmpResult = getJaccard(query_sketches[i].hashSet, ref_sketches[j].hashSet);
				double jaccard = tmpResult.jaccard;
				int common = tmpResult.common;
				int size0 = tmpResult.size0;
				int size1 = tmpResult.size1;
				double mashD;
				if(jaccard == 1.0) 
					mashD = 0.0;
				else if(jaccard == 0.0) 
					mashD = 1.0;
				else 
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				//fprintf(fp0, "jaccard between[%d] and [%d] is: %lf\n", i, j, jaccard);
				if(mashD < maxDist)
					//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", ref_sketches[j].fileName.c_str(), query_sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", query_sketches[i].fileName.c_str(), ref_sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
		}
	}

	//close the subResultFile
	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
	}
	double t1 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

	//struct stat cof_stat;
	FILE * cofp;
	FILE * com_cofp = fopen(outputFile.c_str(), "w");
	//fprintf(fpArr[tid], " %s\t%s\t%d\t%d\t%lf\t%lf\n", sketches[i].fileName.c_str(), sketches[j].fileName.c_str(), common, denom, jaccard, mashD);
	fprintf(com_cofp, " referenceGenome\tqueryGenome\tcommon|size0|size1\tjaccard\tmashD\n");
	int bufSize = 1 << 24;
	int lengthRead = 0;
	char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
	for(int i = 0; i < numThreads; i++)
	{
		cofp = fopen(dist_file_list[i].c_str(), "rb+");
		while(1)
		{
			lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
			//cerr << "the lengthRead is: " << lengthRead << endl;
			fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
			if(lengthRead < bufSize) break;
		}
		fclose(cofp);

		remove(dist_file_list[i].c_str());
	}

	free(bufRead);
	fclose(com_cofp);
	double t2 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t2 - t1 << endl;
}



setResult_t getJaccard(vector<uint32_t> list1, vector<uint32_t> list2)
{
	double jaccard = 1.0;
	int size0 = list1.size();
	int size1 = list2.size();
	int i = 0, j = 0;
	int common = 0;
	while(i < size0 && j < size1)
	{
		if(list1[i] < list2[j]){
			i++;
		}
		else if(list1[i] > list2[j]){
			j++;
		}
		else{
			i++;
			j++;
			common++;
		}

	}
	//cout << "the common is: " << common << endl;
	int denom = size0 + size1 - common;
	jaccard = (double)common / denom;

	setResult_t result{common, size0, size1, jaccard};

	return result;
}

setResult_t vgetJaccard(vector<uint32_t> list1, vector<uint32_t> list2)
{
	int size1 = list1.size();
	int size2 = list2.size();
	int common = u32_intersect_vector_avx2(list1.data(), size1, list2.data(), size2, size1+size2);
	int denom = size1 + size2 - common;
	double jaccard = (double)common / denom;
	setResult_t result{common, size1, size2, jaccard};

	return result;
}


int u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t &a, uint64_t &b){
    uint64_t counter=0;
    const uint32_t *end1 = list1+size1, *end2 = list2+size2;
    a = 0;
    b = 0;
    while(list1 != end1 && list2 != end2 ){
        if(*list1 < *list2){
            list1++;
            a++;
            size3--;
        }else if(*list1 > *list2){
            list2++; 
            b++;
            size3--;
        }else{
            //result[counter++] = *list1;
            counter++;
            list1++; list2++; 
            a++;
            b++;
            size3--;
        }
        if(size3 == 0) break;
    }
    return counter;
}

int u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3){
    //assert(size3 <= size1 + size2);
    int count=0;
		//int a, b;
		uint64_t a = 0; 
		uint64_t b = 0;
		//int *i_a = &a;
		//int *i_b = &b;
		////int * i_a, * i_b;
    //*i_a = 0;
    //*i_b = 0;
    uint64_t st_a = (size1 / 8) * 8;
    uint64_t st_b = (size2 / 8) * 8;

    uint64_t i_a_s, i_b_s;

    if(size3 <= 16){
        count += u32_intersect_scalar_stop(list1, size1, list2, size2, size3, a, b);
        return count;
    }

    //uint64_t stop = size3 - 16;
    while(a < st_a && b < st_b){

        uint32_t a_max = list1[a+7];
        uint32_t b_max = list2[b+7];

        __m256i v_a = _mm256_loadu_si256((__m256i*)&(list1[a]));
        __m256i v_b = _mm256_loadu_si256((__m256i*)&(list2[b]));

        a += (a_max <= b_max) * 8;
        b += (a_max >= b_max) * 8;


        /*constexpr*/ const int32_t cyclic_shift = _MM_SHUFFLE(0,3,2,1); //rotating right
        /*constexpr*/ const int32_t cyclic_shift2= _MM_SHUFFLE(2,1,0,3); //rotating left
        /*constexpr*/ const int32_t cyclic_shift3= _MM_SHUFFLE(1,0,3,2); //between
        __m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a, v_b);
        __m256 rot1 = _mm256_permute_ps((__m256)v_b, cyclic_shift);
        __m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a, (__m256i)rot1);
        __m256 rot2 = _mm256_permute_ps((__m256)v_b, cyclic_shift3);
        __m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a, (__m256i)rot2);
        __m256 rot3 = _mm256_permute_ps((__m256)v_b, cyclic_shift2);
        __m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a, (__m256i)rot3);

        __m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);

        __m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
        __m256 rot5 = _mm256_permute_ps(rot4, cyclic_shift);
        __m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
        __m256 rot6 = _mm256_permute_ps(rot4, cyclic_shift3);
        __m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
        __m256 rot7 = _mm256_permute_ps(rot4, cyclic_shift2);
        __m256i cmp_mask8 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);

        __m256i cmp_mask = _mm256_or_si256(
                _mm256_or_si256(
                    _mm256_or_si256(cmp_mask1, cmp_mask2),
                    _mm256_or_si256(cmp_mask3, cmp_mask4)
                    ),
                _mm256_or_si256(
                    _mm256_or_si256(cmp_mask5, cmp_mask6),
                    _mm256_or_si256(cmp_mask7, cmp_mask8)
                    )
                );
        int32_t mask = _mm256_movemask_ps((__m256)cmp_mask);

        count += _mm_popcnt_u32(mask);

        //if(*i_a + *i_b - count >= stop){
        //    //count -= _mm_popcnt_u32(cmp0);
        //    //*i_a -= (a_max <= b_max) * 16;
        //    //*i_b -= (a_max >= b_max) * 16;
        //    break;
        //}

    }
    count += u32_intersect_scalar_stop(list1+a, size1-a, list2+b, size2-b, size3 - (a+b - count), i_a_s, i_b_s);

    a += i_a_s;
    b += i_b_s;
    return count;
}












