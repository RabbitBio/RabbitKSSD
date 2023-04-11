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

	CLI::App app{"rabbit_kssd: accelerating Kssd-based genome distance estimation on modern multi-core architectures"};
	app.require_subcommand(1);
	CLI::App * shuffle = app.add_subcommand("shuffle", "generate the shuffle file for sketching usage");
	CLI::App * sketch = app.add_subcommand("sketch", "compute sketches for the input genome list");//TODO: only support input list
	CLI::App * alldist = app.add_subcommand("alldist", "compute all-vs-all distances for one input dataset");
	CLI::App * dist = app.add_subcommand("dist", "compute the all-vs-all distances between reference genomes and query datasets");
	CLI::App * setUnion = app.add_subcommand("union", "compute the set union from multiple sketches");
	CLI::App * sub = app.add_subcommand("sub", "subtract the reference sketch from the query sketches");
	CLI::App * convert = app.add_subcommand("convert", "convert the sketches between Kssd format and RabbitKSSD format");
	CLI::App * merge = app.add_subcommand("merge", "merge multiple sketch files into one single sketch file");
	CLI::App * info = app.add_subcommand("info", "get the information of the sketch file");

	string refList = "default";
	string queryList = "default";
	double maxDist = 1.0;
	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	int threads = get_nprocs_conf();
	string shuf_file = "shuf_file/L3K10.shuf";
	string outputFile = "result.out";
	int isContainment = 0;
	bool isQuery = false;
	string inputDir = "default";
	bool isDetail = false;
	int maxNeighbor = 1;
	int leastNumKmer = 1;
	int leastQual = 0;

	auto shuffle_option_k = shuffle->add_option("-k, --halfk", half_k, "the half length of kmer size");
	shuffle->add_option("-s, --subk", half_subk, "the half length of substring space");
	auto shuffle_option_l = shuffle->add_option("-l, --reduction", drlevel, "the dimention reduction level");
	auto shuffle_option_o = shuffle->add_option("-o, --output", outputFile, "set the output shuffle file");
	shuffle_option_k->required();
	shuffle_option_l->required();
	shuffle_option_o->required();


	auto sketch_option_i = sketch->add_option("-i, --input", refList, "list of input genome path, one genome per line");//TODO: only support input list
	auto sketch_option_o = sketch->add_option("-o, --output", outputFile, "set the output file");
	sketch->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling");
	sketch->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	sketch->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	sketch->add_option("-Q, --leastQuality", leastQual, "Filter Kmer with lowest base quality < q in fastq file");
	sketch->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	sketch_option_i->required();
	sketch_option_o->required();


	bool to_Kssd_sketch = false;
	auto convert_option_i = convert->add_option("-i, --input", inputDir, "input Kssd sketches directory(including cofiles.stat, combco.index.0, combco.0) or RabbitKSSD sketches(with --reverse option)");
	auto convert_option_o = convert->add_option("-o, --output", outputFile, "output sketches file in RabbitKSSD format (in Kssd format with --reverse option)");
	convert->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	convert->add_flag("-q, --query", isQuery, "the input genomes are query genome, thus not generate the hash value index dictionary");
	convert->add_flag("--reverse", to_Kssd_sketch, "convert sketch file from RabbitKSSD format to Kssd format");
	convert_option_i->required();
	convert_option_o->required();


	
	auto alldist_option_i = alldist->add_option("-i, --input", refList, "input list of genome path or one sketch file(*.sketch format)");
	auto alldist_option_o = alldist->add_option("-o, --output", outputFile, "set the output file");
	alldist->add_option("-D, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	alldist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling, when input as list of genome path");
	alldist->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	alldist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	alldist->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	alldist->add_option("-Q, --leastQuality", leastQual, "Filter Kmer with lowest base quality < q in fastq file");
	//auto alldist_option_N = alldist->add_option("-N, --neighborN_max", maxNeighbor, "maximum number of neighbor reference output");
	alldist_option_i->required();
	alldist_option_o->required();

	auto dist_option_r = dist->add_option("-r, --reference", refList, "list of reference genome path or a reference sketch file(*.sketch format)");
	auto dist_option_q = dist->add_option("-q, --query", queryList, "list of query genome path or a query sketch file (*.sketch format)");
	auto dist_option_o = dist->add_option("-o, --output", outputFile, "set the output file");
	auto dist_option_N = dist->add_option("-N, --neighborN_max", maxNeighbor, "maximum number of neighbor reference output");
	dist->add_option("-D, --maxDist", maxDist, "maximum distance to save in the result, distances over the maximum distance are omitted");
	dist->add_option("-L", shuf_file, "load the existed shuffle file for Fisher_yates shuffling, when input as list of genome path");
	dist->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	dist->add_option("-M, --metric", isContainment, "output metric: 0, jaccard; 1, containment");
	dist->add_option("-n, --leastNumKmer", leastNumKmer, "specify the least kmer occurence in fastq file");
	dist->add_option("-Q, --leastQuality", leastQual, "Filter Kmer with lowest base quality < q in fastq file");
	dist_option_r->required();
	dist_option_q->required();
	dist_option_o->required();

	string sketchFile = "default";
	auto info_option_i = info->add_option("-i, --input", sketchFile, "input sketch file to get the infomation (*.sketch format)");
	auto info_option_o = info->add_option("-o, --output", outputFile, "set output file");
	info->add_flag("-F, --Fined", isDetail, "output the detailed hash values of each sketch");
	info_option_i->required();
	info_option_o->required();

	auto union_option_i = setUnion->add_option("-i, --input", sketchFile, "sketch file including hash values from multi-sketches (*.sketch format)");
	auto union_option_o = setUnion->add_option("-o, --output", outputFile, "result file for union hash values (*.sketch format)");
	setUnion->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	union_option_i->required();
	union_option_o->required();

	auto merge_option_i = merge->add_option("-i, --input", refList, "list of multiple sketch files list, one file per line");
	auto merge_option_o = merge->add_option("-o, --output", outputFile, "one single sketch file after merging (*.sketch format)");
	merge->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	merge_option_i->required();
	merge_option_o->required();

	string refSketchFile = "default";	
	string querySketchFile = "default";	
	auto sub_option_r = sub->add_option("--rs", refSketchFile, "reference sketch file (*.sketch format) to be subtracted from the query sketches");
	auto sub_option_q = sub->add_option("--qs", querySketchFile, "the query sketches file, (*.sketch format)");
	auto sub_option_o = sub->add_option("-o, --output", outputFile, "result sketch file (*.sketch format) after substraction");
	sub->add_option("-t, --threads", threads, "set the thread number, default all CPUs of the platform");
	sub_option_r->required();
	sub_option_q->required();
	sub_option_o->required();
	

	CLI11_PARSE(app, argc, argv);


	if(app.got_subcommand(shuffle)){
		cerr << "-----generate the shuffle file: " << outputFile << endl;
		dim_shuffle_stat_t shuffle_stat;
		shuffle_stat.k = half_k;
		shuffle_stat.subk = half_subk;
		shuffle_stat.drlevel = drlevel;
		write_shuffle_dim_file(&shuffle_stat, outputFile);
		return 0;
	}
	else if(app.got_subcommand(info)){
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
		command_sub(refSketchFile, querySketchFile, outputFile, threads);
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
	
	if(app.got_subcommand(sketch)){
		cerr << "-----run the subcommand: sketch" << endl;
		if(isSketchFile(refList)){
			if(!isQuery){
				vector<sketch_t> sketches;
				sketchInfo_t info;
				readSketches(sketches, info, refList);
				cerr << "input is a sketch file, rename the sketch file from: " << refList << " to: " << outputFile << endl;
				string cmd0 = "cp " + refList + ' ' + outputFile;
				int stop = system(cmd0.c_str());
				if(!stop)	cerr << cmd0 << endl;
				
				//saveSketches(sketches, info, outputFile);
				double tstart = get_sec();
				string dictFile = outputFile + ".dict";
				string indexFile = outputFile + ".index";
				transSketches(sketches, info, dictFile, indexFile, threads);
				double tend = get_sec();
				cerr << "===============the time of transSketches is: " << tend - tstart << endl;
			}
			else{
				//cerr << "input is a sketch file, do nothing" << endl;
				cerr << "input is a sketch file, rename the sketch file from: " << refList << " to: " << outputFile << endl;
				string cmd0 = "mv " + refList + ' ' + outputFile;
				int stop = system(cmd0.c_str());
				if(!stop)	cerr << cmd0 << endl;
			}
			return 0;
		}

		cerr << "---read the shuffle file: " << shuf_file << endl;
		dim_shuffle_t* shuffled_info = read_shuffle_dim(shuf_file);
		half_k = shuffled_info->dim_shuffle_stat.k;
		half_subk = shuffled_info->dim_shuffle_stat.subk;
		drlevel = shuffled_info->dim_shuffle_stat.drlevel;
		kssd_parameter_t kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_info->shuffled_dim);

		command_sketch(refList, isQuery, outputFile, kssd_parameter, leastQual, leastNumKmer, threads);
		return 0;
	}
	else if(app.got_subcommand(alldist)){
		cerr << "-----run the subcommand: alldist" << endl;
		kssd_parameter_t kssd_parameter;
		if(!isSketchFile(refList)){
			dim_shuffle_t* shuffled_info = read_shuffle_dim(shuf_file);
			half_k = shuffled_info->dim_shuffle_stat.k;
			half_subk = shuffled_info->dim_shuffle_stat.subk;
			drlevel = shuffled_info->dim_shuffle_stat.drlevel;
			kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_info->shuffled_dim);
		}
		command_alldist(refList, outputFile, kssd_parameter, leastQual, leastNumKmer, maxDist, isContainment, threads);
		return 0;
	}
	else if(app.got_subcommand(dist)){
		cerr << "-----run the subcommand: dist" << endl;
		kssd_parameter_t kssd_parameter;
		if(!isSketchFile(refList) || !isSketchFile(queryList)){
			dim_shuffle_t* shuffled_info = read_shuffle_dim(shuf_file);
			half_k = shuffled_info->dim_shuffle_stat.k;
			half_subk = shuffled_info->dim_shuffle_stat.subk;
			drlevel = shuffled_info->dim_shuffle_stat.drlevel;
			kssd_parameter = initParameter(half_k, half_subk, drlevel, shuffled_info->shuffled_dim);
		}
		bool isNeighbor = false;
		if(*dist_option_N){
			isNeighbor = true;
		}
		command_dist(refList, queryList, outputFile, kssd_parameter, leastQual, leastNumKmer, maxDist, maxNeighbor, isNeighbor, isContainment, threads);
		return 0;
	}
	
	return 0;
}






