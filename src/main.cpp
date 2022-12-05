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
#include "sketch.h"
#include "dist.h"
#include "subCommand.h"


#include <algorithm>
#include <omp.h>
#include <cmath>
#include "CLI11.hpp"


#include <sys/stat.h>
#include <sys/sysinfo.h>


using namespace std;

int main(int argc, char * argv[]){
	
	double t0 = get_sec();

	CLI::App app{"kssd: k-mer substring space sampling for analyzing large-scale genome databases"};
	app.require_subcommand(1);
	CLI::App * sketch = app.add_subcommand("sketch", "computing sketches for the input genome list");
	CLI::App * info = app.add_subcommand("info", "get the information of the sketches");
	CLI::App * alldist = app.add_subcommand("alldist", "computing all to all distances for the input genome list");
	CLI::App * dist = app.add_subcommand("dist", "computing the distances between reference genomes and query datasets");
	CLI::App * merge = app.add_subcommand("merge", "merging the sketches into one sketch containing a set of hash values");
	CLI::App * sub = app.add_subcommand("sub", "subtracting the specific sketch from the query sketches");
	CLI::App * convert = app.add_subcommand("convert", "converting the sketches from Kssd format to RabbitKSSD format");

	string refList = "default";
	string queryList = "default";
	double maxDist = 1.0;
	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	int threads = get_nprocs_conf();
	string shuf_file = "shuf_file/bact.shuf";
	string outputFile = "result.out";
	int isContainment = 0;
	bool isQuery = false;
	string inputDir = "default";

	auto sketch_option_i = sketch->add_option("-i, --input", refList, "list of input genome path, one genome per line");
	auto sketch_option_k = sketch->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto sketch_option_s = sketch->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto sketch_option_l = sketch->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto sketch_option_L = sketch->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto sketch_option_o = sketch->add_option("-o, --output", outputFile, "set the output file");
	auto sketch_option_t = sketch->add_option("-t, --threads", threads, "set the thread number, default all cpus of the platform");
	auto sketch_flag_q = sketch->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	sketch_option_k->excludes(sketch_option_L);
	sketch_option_s->excludes(sketch_option_L);
	sketch_option_l->excludes(sketch_option_L);
	sketch_option_i->required();
	sketch_option_o->required();

	auto convert_option_i = convert->add_option("-i, --input", inputDir, "input Kssd sketches directory, including cofiles.stat, combco.index.0, combco.0");
	auto convert_option_o = convert->add_option("-o, --output", outputFile, "output sketches file in RabbitKSSD format");
	auto convert_option_t = convert->add_option("-t, --threads", threads, "set the thread number, default all cpus of the platform");
	auto convert_flag_q = convert->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	convert_option_i->required();
	convert_option_o->required();


	
	auto alldist_option_i = alldist->add_option("-i, --input", refList, "list of input genome path, one genome per line");
	auto alldist_option_m = alldist->add_option("-m, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	auto alldist_option_k = alldist->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto alldist_option_s = alldist->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto alldist_option_l = alldist->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto alldist_option_L = alldist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto alldist_option_o = alldist->add_option("-o, --output", outputFile, "set the output file");
	auto alldist_option_t = alldist->add_option("-t, --threads", threads, "set the thread number");
	auto alldist_option_M = alldist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	alldist_option_k->excludes(alldist_option_L);
	alldist_option_s->excludes(alldist_option_L);
	alldist_option_l->excludes(alldist_option_L);
	alldist_option_i->required();
	alldist_option_o->required();

	auto dist_option_r = dist->add_option("-r, --reference", refList, "list of reference genome path, one genome per line");
	auto dist_option_q = dist->add_option("-q, --query", queryList, "list of query genome path, one genome per line");
	auto dist_option_m = dist->add_option("-m, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	auto dist_option_k = dist->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto dist_option_s = dist->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto dist_option_l = dist->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto dist_option_L = dist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto dist_option_o = dist->add_option("-o, --output", outputFile, "set the output file");
	auto dist_option_t = dist->add_option("-t, --threads", threads, "set the thread number");
	auto dist_option_M = dist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	dist_option_k->excludes(dist_option_L);
	dist_option_s->excludes(dist_option_L);
	dist_option_l->excludes(dist_option_L);
	dist_option_r->required();
	dist_option_q->required();
	dist_option_o->required();

	string sketchFile = "default";
	auto info_option_i = info->add_option("-i, --input", sketchFile, "input sketch file to get the infomation");
	auto info_option_o = info->add_option("-o, --output", outputFile, "the genome name and hash values in each sketch");
	info_option_i->required();
	info_option_o->required();

	auto merge_option_i = merge->add_option("-i, --input", sketchFile, "sketch file including hash values from multi-sketches, which is saved from the distance computing");
	auto merge_option_o = merge->add_option("-o, --output", outputFile, "result file for merged hash values");
	auto merge_option_t = merge->add_option("-t, --threads", threads, "set the thread number");
	merge_option_i->required();
	merge_option_o->required();

	string refSketchFile = "default";	
	string querySketchFile = "default";	
	auto sub_option_r = sub->add_option("--rs", refSketchFile, "the union sketches for reference sketches, which need to be subtracted from the query sketches");
	auto sub_option_q = sub->add_option("--qs", querySketchFile, "the query sketches file");
	auto sub_option_o = sub->add_option("-o, --output", outputFile, "result file for merged hash values");
	auto sub_option_t = sub->add_option("-t, --threads", threads, "set the thread number");
	sub_option_r->required();
	sub_option_q->required();
	sub_option_o->required();
	

	CLI11_PARSE(app, argc, argv);

	if(app.got_subcommand(info)){
		cerr << "-----run the subcommand: info" << endl;
		command_info(sketchFile, outputFile);
		return 0;
	}
	else if(app.got_subcommand(merge)){
		cerr << "-----run the subcommand: merge" << endl;
		command_merge(sketchFile, outputFile, threads);
		return 0;
	}
	else if(app.got_subcommand(sub)){
		cerr << "-----run the subcommand: sub" << endl;
		command_sub(refSketchFile, querySketchFile, outputFile, threads);
		return 0;
	}
	else if(app.got_subcommand(convert)){
		cerr << "-----run the subcommand: convert" << endl;
		command_convert(inputDir, isQuery, outputFile, threads);
		return 0;
	}

	double t1 = get_sec();
	//cerr << "===================time of init parameters is: " << t1 - t0 << endl;

	int kmer_size = 2 * half_k;
	int * shuffled_dim;
	
	//if(*sketch_option_L || *alldist_option_L || *dist_option_L){
	//	cerr << "---read the shuffle file: " << shuf_file << endl;
	//	shuffled_dim = read_shuffle_dim(shuf_file);
	//}
	//else{
	//	cerr << "---generate the shuffle file: " << endl;
	//	shuffled_dim = generate_shuffle_dim(half_subk);
	//}

	//kssd_parameter_t kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);

	//double t2 = get_sec();
	//cerr << "===================time of read shuffle file is: " << t2 - t1 << endl;

	if(app.got_subcommand(sketch)){
		cerr << "-----run the subcommand: sketch" << endl;
		if(*sketch_option_L){
			cerr << "---read the shuffle file: " << shuf_file << endl;
			shuffled_dim = read_shuffle_dim(shuf_file);
		}
		else{
			cerr << "---generate the shuffle file: " << endl;
			shuffled_dim = generate_shuffle_dim(half_subk);
		}
		kssd_parameter_t kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);

		command_sketch(refList, isQuery, outputFile, kssd_parameter, threads);
	}
	else if(app.got_subcommand(alldist)){
		cerr << "-----run the subcommand: alldist" << endl;
		kssd_parameter_t kssd_parameter;
		if(!isSketchFile(refList)){
			if(*alldist_option_L){
				cerr << "---read the shuffle file: " << shuf_file << endl;
				shuffled_dim = read_shuffle_dim(shuf_file);
			}
			else{
				cerr << "---generate the shuffle file: " << endl;
				shuffled_dim = generate_shuffle_dim(half_subk);
			}
			kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);
		}
		command_alldist(refList, outputFile, kssd_parameter, kmer_size, maxDist, isContainment, threads);
	}
	else if(app.got_subcommand(dist)){
		cerr << "-----run the subcommand: dist" << endl;
		kssd_parameter_t kssd_parameter;
		if(!isSketchFile(refList) || !isSketchFile(queryList)){
			if(*dist_option_L){
				cerr << "---read the shuffle file: " << shuf_file << endl;
				shuffled_dim = read_shuffle_dim(shuf_file);
			}
			else{
				cerr << "---generate the shuffle file: " << endl;
				shuffled_dim = generate_shuffle_dim(half_subk);
			}
			kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);
		}
		command_dist(refList, queryList, outputFile, kssd_parameter, kmer_size, maxDist, isContainment, threads);
	}
	
	return 0;
}

//========================================================================================
//========================================================================================
//========================================================================================

#ifdef DIST_INDEX
	uint64_t hash_range = 1 << (4 * (half_k - drlevel));
	vector< vector<uint64_t> > hash_index_arr[threads];
	//uint64_t ** hash_index_arr = (uint64_t **)malloc(threads * sizeof(uint64_t *));
	for(int i = 0; i < threads; i++)
	{
		hash_index_arr[i].resize(hash_range);
	}
	
	//double t22 = get_sec();
	//cerr << "===================time of multiThread generate sketches and index  is: " << t22 - t2 << endl;

	vector< vector<uint64_t> > final_hash_index;
	final_hash_index.resize(hash_range);
	for(int i = 0; i < hash_range; i++)
	{
		for(int j = 0; j < threads; j++)
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
	for(int i = 0; i < threads; i++)
	{
		string tmpName = output + to_string(i);
		dist_file_list.push_back(tmpName);

		FILE * fp0 = fopen(tmpName.c_str(), "w");
		fpArr.push_back(fp0);
	}
	
	//FILE * fp = fopen(output.c_str(), "w");
	#pragma omp parallel for num_threads(threads) schedule(dynamic)
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
	for(int i = 0; i < threads; i++)
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





