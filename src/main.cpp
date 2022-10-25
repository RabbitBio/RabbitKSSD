#include <iostream>
#include <errno.h>
#include <err.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <stdint.h>

#include "shuffle.h"
#include "common.h"
#include "parameter.h"
#include "sketch.h"
#include "dist.h"


#include <algorithm>
#include <omp.h>
#include <cmath>


#include <sys/stat.h>


int numThreads = 48;

using namespace std;

bool cmp(uint64_t a, uint64_t b){
	return a < b;
}



bool cmpSketch(sketch_t s1, sketch_t s2){
	return s1.id < s2.id;
}

#ifdef DIST_INDEX
void get_dist_by_index(vector< vector<uint64_t> > final_hash_index, vector<sketch_t> sketches, int kmer_size);
#endif

int main(int argc, char * argv[]){
	double t0 = get_sec();
	if(argc < 2)
		err(errno, "argument number %d must be greater than 1", argc);

	//cout << "start main: " << endl;
	//global arguments
	string inputFile = argv[1];
	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	int kmer_size = 2 * half_k;

	//arguments for sketch query
	bool isList = true;
	double t1 = get_sec();
	cerr << "===================time of init parameters is: " << t1 - t0 << endl;

	//step1: get the shuffle array
	bool readShuffleFile = true;

	int * shuffled_dim;
	if(readShuffleFile){
		string shufFile = "shuf_file/bact.shuf";
		shuffled_dim = read_shuffle_dim(shufFile);
	}
	else{
		shuffled_dim = generate_shuffle_dim(half_subk);
	}

	kssd_parameter_t kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);

	double t2 = get_sec();
	cerr << "===================time of read shuffle file is: " << t2 - t1 << endl;


	//step2: generate the query sketches for each genome or sequence.
	//uint64_t ** coList = (uint64_t **)malloc(fileList.size() * sizeof(uint64_t*));

	vector<sketch_t> sketches;
	//vector< vector<uint64_t> > hashList;
	
	bool success = sketchFile(inputFile, numThreads, kssd_parameter, sketches);
	cerr << "finish the sketch generation " << success << endl;
	cerr << "the size of sketches is: " << sketches.size() << endl;



	double t3 = get_sec();
	cerr << "===================time of generator sketches is: " << t3 - t2 << endl;

	std::sort(sketches.begin(), sketches.end(), cmpSketch);
	//for(int i = 0; i < fileList.size(); i++){
	//	vector<uint64_t> hashArr;
	//	for(int j = 0; j < hashSize; j++){
	//		if(coList[i][j] != 0){
	//			hashArr.push_back(coList[i][j]);
	//		}
	//	}
	//	hashList.push_back(hashArr);
	//}
	double t4 = get_sec();
	cerr << "===================time of sort sketches order is: " << t4 - t3 << endl;

	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int i = 0; i < sketches.size(); i++)
	{
		std::sort(sketches[i].hashSet.begin(), sketches[i].hashSet.end(), cmp);
	}
	double t5 = get_sec();
	cerr << "===================time of sort each sketches is: " << t5 - t4 << endl;

	string hash_output = "hash.out";
	FILE * fp = fopen(hash_output.c_str(), "w+");
	for(int i = 0; i < sketches.size(); i++){
		//cout << "the sketch size of file " << i << " is: " << hashList[i].size() << endl;
		//cerr << "the sketch size of " << sketches[i].fileName << " is: " << sketches[i].hashSet.size() << endl;
		//fprintf(stderr, "the sketch size of file %s is: %d\n", i, sketches[i].fileName.c_str(), sketches[i].hashSet.size());
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			fprintf(fp, "%llu\t", sketches[i].hashSet[j]);
			if(j % 10 == 9) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fprintf(fp, "\n");
		//fprintf(fp, "\n");
	}
	fclose(fp);

	double t6 = get_sec();
	cerr << "===================time of save sketches into hash.out is: " << t6 - t5 << endl;

	//compute the pair distance
	
	cerr << "start the distance computing" << endl;
	string dist_output = "result.out";

	tri_dist(sketches, dist_output, kmer_size, numThreads);

	//FILE * fp0 = fopen(dist_output.c_str(), "w");
	double t7 = get_sec();
	cerr << "===================time of get total distance matrix file is: " << t7 - t6 << endl;


	return 0;
}

//========================================================================================
//========================================================================================
//========================================================================================

#ifdef DIST_INDEX
	uint64_t hash_range = 1 << (4 * (half_k - drlevel));
	vector< vector<uint64_t> > hash_index_arr[numThreads];
	//uint64_t ** hash_index_arr = (uint64_t **)malloc(numThreads * sizeof(uint64_t *));
	for(int i = 0; i < numThreads; i++)
	{
		hash_index_arr[i].resize(hash_range);
	}
	
	//double t22 = get_sec();
	//cerr << "===================time of multiThread generate sketches and index  is: " << t22 - t2 << endl;

	vector< vector<uint64_t> > final_hash_index;
	final_hash_index.resize(hash_range);
	for(int i = 0; i < hash_range; i++)
	{
		for(int j = 0; j < numThreads; j++)
		{
			final_hash_index[i].insert(final_hash_index[i].end(), hash_index_arr[j][i].begin(), hash_index_arr[j][i].end());
		}
	}

	double t33 = get_sec();
	cerr << "===================time of merge multi index is: " << t33 - t22 << endl;

	get_dist_by_index(final_hash_index, sketches, kmer_size);
	double t44 = get_sec();
	cerr << "===================time of get distance matrix by index is: " << t44 - t33 << endl;
	exit(0);
	//return 0;//end main
#endif
	




#ifdef DIST_INDEX
void get_dist_by_index(vector< vector<uint64_t> > final_hash_index, vector<sketch_t> sketches, int kmer_size)
{
	int genomeNumber = sketches.size();
	int dist_size = genomeNumber * genomeNumber;
	uint32_t * common_matrix = (uint32_t*)malloc(dist_size * sizeof(uint32_t));
	for(int i = 0; i < dist_size; i++)
	{
		common_matrix[i] = 0;
	}
	//memset(common_matrix, 0, dist_size);
	cerr << "start to generate the distance matrix " << endl;

	for(int i = 0; i < sketches.size(); i++)
	{
		//int dist_offset = i * genomeNumber;
		for(int j = 0; j < sketches[i].hashSet.size(); j++)
		{
			uint64_t ind = sketches[i].hashSet[j];
			int common_genome_number = final_hash_index[ind].size();
			for(int k = 0; k < common_genome_number; k++)
			{
				int common_gid = final_hash_index[ind][k];
				common_matrix[i * genomeNumber + common_gid]++;
			}
		}
	}//end for loop to compute the common_matrix
	cerr << "end to generate the distance matrix " << endl;

	//double time_1 = get_sec();
	string output = "result.out";
	vector<string> dist_file_list;
	vector<FILE*> fpArr;
	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = output + to_string(i);
		dist_file_list.push_back(tmpName);

		FILE * fp0 = fopen(tmpName.c_str(), "w");
		fpArr.push_back(fp0);
	}
	
	//FILE * fp = fopen(output.c_str(), "w");
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int i = 0; i < sketches.size(); i++)
	{
		int tid = omp_get_thread_num();
		for(int j = 0; j < sketches.size(); j++)
		{
			int common = common_matrix[i * genomeNumber + j];
			int union_number = sketches[i].hashSet.size() + sketches[j].hashSet.size() - common;
			double jaccard = (double)common / union_number;
			double mashDist;
			if(jaccard == 0.0)
				mashDist = 1.0;
			else if(jaccard == 1.0)
				mashDist = 0.0;
			else
				mashDist = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
			fprintf(fpArr[tid], "%s\t%s\t%d\t%lf\t%lf\n", sketches[i].fileName.c_str(), sketches[j].fileName.c_str(), common, jaccard, mashDist);
		}
	}
	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
	}
	//fclose(fp);
	double time_2 = get_sec();
	//cerr << "===================time of output distance.out in index is: " << time_2 - time_1 << endl;

	cerr << "close the result.out" << endl;

	//free(common_matrix);
	//vector< vector<uint63_t> >().swap(final_hash_index);
	//vector<sketch_t>().swap(sketches);

}
#endif





