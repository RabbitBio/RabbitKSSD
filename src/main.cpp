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
	CLI::App * setUnion = app.add_subcommand("union", "computing the set union from multiple sketches");
	CLI::App * sub = app.add_subcommand("sub", "subtracting the specific sketch from the query sketches");
	CLI::App * convert = app.add_subcommand("convert", "converting the sketches from Kssd format to RabbitKSSD format");
	CLI::App * merge = app.add_subcommand("merge", "merge multiple sketch files into one output sketch file");

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
	bool isDetail = false;
	int maxNeighbor = 1;
	int leastNumKmer = 1;

	auto sketch_option_i = sketch->add_option("-i, --input", refList, "list of input genome path, one genome per line");
	auto sketch_option_k = sketch->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto sketch_option_s = sketch->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto sketch_option_l = sketch->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto sketch_option_L = sketch->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto sketch_option_o = sketch->add_option("-o, --output", outputFile, "set the output file");
	auto sketch_option_t = sketch->add_option("-t, --threads", threads, "set the thread number, default all cpus of the platform");
	auto sketch_option_n = sketch->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	auto sketch_flag_q = sketch->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	sketch_option_k->excludes(sketch_option_L);
	sketch_option_s->excludes(sketch_option_L);
	sketch_option_l->excludes(sketch_option_L);
	sketch_option_i->required();
	sketch_option_o->required();


	bool to_Kssd_sketch = false;
	auto convert_option_i = convert->add_option("-i, --input", inputDir, "input Kssd sketches directory, including cofiles.stat, combco.index.0, combco.0");
	auto convert_option_o = convert->add_option("-o, --output", outputFile, "output sketches file in RabbitKSSD format");
	auto convert_option_t = convert->add_option("-t, --threads", threads, "set the thread number, default all cpus of the platform");
	auto convert_flag_q = convert->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	auto convert_flag_reverse = convert->add_flag("--reverse", to_Kssd_sketch, "convert sketch file from RabbitKSSD format to Kssd format");
	convert_option_i->required();
	convert_option_o->required();


	
	auto alldist_option_i = alldist->add_option("-i, --input", refList, "list of input genome path, one genome per line");
	auto alldist_option_D = alldist->add_option("-D, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	auto alldist_option_k = alldist->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto alldist_option_s = alldist->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto alldist_option_l = alldist->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto alldist_option_L = alldist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto alldist_option_o = alldist->add_option("-o, --output", outputFile, "set the output file");
	auto alldist_option_t = alldist->add_option("-t, --threads", threads, "set the thread number");
	auto alldist_option_M = alldist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	auto alldist_option_n = alldist->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	//auto alldist_option_N = alldist->add_option("-N, --neighborN_max", maxNeighbor, "maximum number of neighbor reference output");
	alldist_option_k->excludes(alldist_option_L);
	alldist_option_s->excludes(alldist_option_L);
	alldist_option_l->excludes(alldist_option_L);
	alldist_option_i->required();
	alldist_option_o->required();

	auto dist_option_r = dist->add_option("-r, --reference", refList, "list of reference genome path, one genome per line");
	auto dist_option_q = dist->add_option("-q, --query", queryList, "list of query genome path, one genome per line");
	auto dist_option_D = dist->add_option("-D, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	auto dist_option_k = dist->add_option("-k, --halfk", half_k, "the half length of kmer size");
	auto dist_option_s = dist->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto dist_option_l = dist->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto dist_option_L = dist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	auto dist_option_o = dist->add_option("-o, --output", outputFile, "set the output file");
	auto dist_option_t = dist->add_option("-t, --threads", threads, "set the thread number");
	auto dist_option_M = dist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	auto dist_option_N = dist->add_option("-N, --neighborN_max", maxNeighbor, "maximum number of neighbor reference output");
	auto dist_option_n = dist->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	dist_option_k->excludes(dist_option_L);
	dist_option_s->excludes(dist_option_L);
	dist_option_l->excludes(dist_option_L);
	dist_option_r->required();
	dist_option_q->required();
	dist_option_o->required();

	string sketchFile = "default";
	auto info_option_i = info->add_option("-i, --input", sketchFile, "input sketch file to get the infomation");
	auto info_option_o = info->add_option("-o, --output", outputFile, "the genome name and hash values in each sketch");
	auto info_flag_F = info->add_flag("-F, --Fined", isDetail, "output the detailed hash values of each sketch");
	info_option_i->required();
	info_option_o->required();

	auto union_option_i = setUnion->add_option("-i, --input", sketchFile, "sketch file including hash values from multi-sketches, which is saved from the distance computing");
	auto union_option_o = setUnion->add_option("-o, --output", outputFile, "result file for union hash values");
	auto union_option_t = setUnion->add_option("-t, --threads", threads, "set the thread number");
	union_option_i->required();
	union_option_o->required();

	auto merge_option_i = merge->add_option("-i, --input", refList, "list file of sketch files for merging, one file per line");
	auto merge_option_o = merge->add_option("-o, --output", outputFile, "result sketch file after merging");
	auto merge_option_t = merge->add_option("-t, --threads", threads, "set the thread number");
	merge_option_i->required();
	merge_option_o->required();

	string refSketchFile = "default";	
	string querySketchFile = "default";	
	auto sub_option_r = sub->add_option("--rs", refSketchFile, "the union sketches for reference sketches, which need to be subtracted from the query sketches");
	auto sub_option_q = sub->add_option("--qs", querySketchFile, "the query sketches file");
	auto sub_option_o = sub->add_option("-o, --output", outputFile, "result file for unioned hash values");
	auto sub_option_t = sub->add_option("-t, --threads", threads, "set the thread number");
	sub_option_r->required();
	sub_option_q->required();
	sub_option_o->required();
	

	CLI11_PARSE(app, argc, argv);

	if(app.got_subcommand(info)){
		cerr << "-----run the subcommand: info" << endl;
		command_info(sketchFile, isDetail, outputFile);
		return 0;
	}
	else if(app.got_subcommand(setUnion)){
		cerr << "-----run the subcommand: union" << endl;
		command_union(sketchFile, outputFile, threads);
		return 0;
	}
	else if(app.got_subcommand(sub)){
		cerr << "-----run the subcommand: sub" << endl;
		//command_sub(refSketchFile, querySketchFile, outputFile, threads);
		new_command_sub(refSketchFile, querySketchFile, outputFile, threads);
		return 0;
	}
	else if(app.got_subcommand(convert)){
		cerr << "-----run the subcommand: convert" << endl;
		command_convert(inputDir, to_Kssd_sketch, isQuery, outputFile, threads);
		return 0;
	}
	else if(app.got_subcommand(merge)){
		cerr << "-----run the subcommand: merge" << endl;
		command_merge(refList, outputFile, threads);
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
		if(isSketchFile(refList)){
			if(!isQuery){
				vector<sketch_t> sketches;
				sketchInfo_t info;
				readSketches(sketches, info, refList);
				string cmd0 = "cp " + refList + ' ' + outputFile;
				system(cmd0.c_str());
				//saveSketches(sketches, info, outputFile);
				double tstart = get_sec();
				string dictFile = outputFile + ".dict";
				string indexFile = outputFile + ".index";
				transSketches(sketches, info, dictFile, indexFile, threads);
				double tend = get_sec();
				cerr << "===============the time of transSketches is: " << tend - tstart << endl;
			}
			else{
				cerr << "input is a sketch file, do nothing" << endl;
			}
			return 0;
		}
		if(*sketch_option_L){
			cerr << "---read the shuffle file: " << shuf_file << endl;
			shuffled_dim = read_shuffle_dim(shuf_file);
		}
		else{
			cerr << "---generate the shuffle file: " << endl;
			shuffled_dim = generate_shuffle_dim(half_subk);
		}
		kssd_parameter_t kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);

		command_sketch(refList, isQuery, outputFile, kssd_parameter, leastNumKmer, threads);
		return 0;
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
		command_alldist(refList, outputFile, kssd_parameter, kmer_size, leastNumKmer, maxDist, isContainment, threads);
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
		bool isNeighbor = false;
		if(*dist_option_N){
			isNeighbor = true;
		}
		command_dist(refList, queryList, outputFile, kssd_parameter, kmer_size, leastNumKmer, maxDist, maxNeighbor, isNeighbor, isContainment, threads);
	}
	
	return 0;
}






