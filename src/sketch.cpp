#include "sketch.h"
#include <err.h>
#include "kseq.h"
#include <zlib.h>
#include <omp.h>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include "robin_hood.h"
#include <sys/stat.h>

#include <unordered_set>
#include "Formater.h"
#include "RabbitFX.h"

KSEQ_INIT(gzFile, gzread);

using namespace std;

#define H1(K, HASH_SZ) ((K) % (HASH_SZ))
#define H2(K, HASH_SZ) (1 + (K) % ((HASH_SZ)-1))
#define HASH(K, I, HASH_SZ) ((H1(K, HASH_SZ) + I * H2(K, HASH_SZ)) % HASH_SZ)

#define PATHLEN 256

bool cmp(uint32_t a, uint32_t b){
	return a < b;
}

bool cmpSketch(sketch_t s1, sketch_t s2){
	return s1.id < s2.id;
}

bool cmpSketchName(sketch_t s1, sketch_t s2){
	return s1.fileName < s2.fileName;
}

bool cmpFile(fileInfo_t f1, fileInfo_t f2){
	return f1.fileSize > f2.fileSize;
}

bool existFile(string FileName){
	if(FILE * fp = fopen(FileName.c_str(), "r")){
		fclose(fp);
		return true;
	}
	else 
		return false;
}

bool isFasta(string inputFile){
	ifstream ifs(inputFile);
	string line;
	getline(ifs, line);
	if(line[0] == '>') return true;
	else return false;
}

bool isFastq(string inputFile){
	ifstream ifs(inputFile);
	string line;
	getline(ifs, line);
	if(line[0] == '@') return true;
	else return false;
}

bool isFastaList(string inputList){
	ifstream ifs(inputList);
	string line;
	bool res = true;
	while(getline(ifs, line)){
		if(!isFasta(line)){
			res = false;
			break;
		}
	}
	ifs.close();
	return res;
}

bool isFastqList(string inputList){
	ifstream ifs(inputList);
	string line;
	bool res = true;
	while(getline(ifs, line)){
		if(!isFastq(line)){
			res = false;
			break;
		}
	}
	ifs.close();
	return res;
}


bool isSketchFile(string inputFile){
	int startPos = inputFile.find_last_of('.');
	if(startPos == string::npos) return false;
	string suffixName = inputFile.substr(startPos+1);
	if(suffixName == "sketch")	return true;
	else return false;
}


//void consumer_fasta_task(FXReader<FA> &m_reader, kssd_parameter_t& parameter, robin_hood::unordered_map<uint32_t, int>& shuffled_map, vector<uint32_t>& hashArr){
void consumer_fasta_task(FXReader<FA> &m_reader, kssd_parameter_t& parameter, robin_hood::unordered_map<uint32_t, int>& shuffled_map, robin_hood::unordered_set<uint32_t>& hashValueSet){

	int half_k = parameter.half_k;
	int drlevel = parameter.drlevel;
	int rev_add_move = parameter.rev_add_move;
	int half_outctx_len = parameter.half_outctx_len;
	int * shuffled_dim = parameter.shuffled_dim;
	int dim_start = parameter.dim_start;
	int dim_end = parameter.dim_end;
	int kmer_size = parameter.kmer_size;
	uint32_t hashSize = parameter.hashSize;
	int hashLimit = parameter.hashLimit;
	uint64_t tupmask = parameter.tupmask;
	uint64_t undomask0 = parameter.undomask0;
	uint64_t undomask1 = parameter.undomask1;
	uint64_t domask = parameter.domask;

	//unordered_set<uint32_t> hashValueSet;
	while(1){
		auto data = m_reader.get_formated_reads();
		if(data.size() == 0) break;
		for(Reference &r: data){
			int length = r.length;
			string name = r.name;

			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int keyCount = 0;
			uint64_t base = 1;

			for(int i = 0; i < length; i++){
				char ch = r.seq[i];
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;
				}
				if(base > kmer_size)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
					if(shuffled_map.count(dim_id) == 0){
						continue;
					}
					pfilter = shuffled_map[dim_id];

					pfilter -= dim_start;
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 

					hashValueSet.insert(dr_tuple);
				}
			} //end for, of a sequence
		}//end of all sequence

	}//end while(1)

	//for(auto x : hashValueSet)
	//	hashArr.push_back(x);
}

void consumer_fastq_task(FXReader<FQ_SE> &m_reader, kssd_parameter_t& parameter, robin_hood::unordered_map<uint32_t, int>& shuffled_map, int& leastNumKmer, robin_hood::unordered_map<uint32_t, int>& hashValueMap){

	int half_k = parameter.half_k;
	int drlevel = parameter.drlevel;
	int rev_add_move = parameter.rev_add_move;
	int half_outctx_len = parameter.half_outctx_len;
	int * shuffled_dim = parameter.shuffled_dim;
	int dim_start = parameter.dim_start;
	int dim_end = parameter.dim_end;
	int kmer_size = parameter.kmer_size;
	uint32_t hashSize = parameter.hashSize;
	int hashLimit = parameter.hashLimit;
	uint64_t tupmask = parameter.tupmask;
	uint64_t undomask0 = parameter.undomask0;
	uint64_t undomask1 = parameter.undomask1;
	uint64_t domask = parameter.domask;
		
	//cerr << "the thread id is: " << id << endl;
	//string outFileName = "output" + to_string(id);
	//FILE * fpTmp = fopen(outFileName.c_str(), "w+");
	
	while(1){
		auto data = m_reader.get_formated_reads();
		if(data.size() == 0) break;
		for(Reference &r: data){
			int length = r.seq.length();
			string name = r.name;
			//cout << name << endl;
			//cout << r.seq << endl;
			//fprintf(fpTmp, "%s\n", r.seq.c_str());

			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int keyCount = 0;
			uint64_t base = 1;

			for(int i = 0; i < length; i++){
				char ch = r.seq[i];
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;
				}
				if(base > kmer_size)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
					if(shuffled_map.count(dim_id) == 0){
						continue;
					}
					pfilter = shuffled_map[dim_id];

					pfilter -= dim_start;
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 
					hashValueMap.insert({dr_tuple, 0});
					hashValueMap[dr_tuple]++;

				}
			} //end for, of a sequence
		}//end of all sequence

	}//end while(1)

}



bool sketchFastaFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, string outputFile){
	//cerr << "run the sketchFastaFile " << endl;
	int half_k = parameter.half_k;
	int half_subk = parameter.half_subk;
	int drlevel = parameter.drlevel;
	int rev_add_move = parameter.rev_add_move;
	int half_outctx_len = parameter.half_outctx_len;
	int * shuffled_dim = parameter.shuffled_dim;
	int dim_start = parameter.dim_start;
	int dim_end = parameter.dim_end;
	int kmer_size = parameter.kmer_size;
	uint32_t hashSize = parameter.hashSize;
	int hashLimit = parameter.hashLimit;
	uint64_t tupmask = parameter.tupmask;
	uint64_t undomask0 = parameter.undomask0;
	uint64_t undomask1 = parameter.undomask1;
	uint64_t domask = parameter.domask;

	int dim_size = 1 << 4 * (half_k - half_outctx_len);
	robin_hood::unordered_map<uint32_t, int> shuffled_map;
	for(int t = 0; t < dim_size; t++){
		if(shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start){
			shuffled_map.insert({t, shuffled_dim[t]});
		}
		//cout << shuffled_dim[t] << endl;
	}


	ifstream fs(inputFile);
	if(!fs){
		err(errno, "cannot open the inputFile: %s\n", inputFile.c_str());
	}
	vector<fileInfo_t> fileList;
	uint64_t totalSize = 0;
	string fileName;
	while(getline(fs, fileName)){
		struct stat cur_stat;
		stat(fileName.c_str(), &cur_stat);
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
		fileInfo_t tmpF;
		tmpF.fileName = fileName;
		tmpF.fileSize = curSize;
		fileList.push_back(tmpF);
	}
	std::sort(fileList.begin(), fileList.end(), cmpFile);
	uint64_t limitSize = totalSize / numThreads;
	vector<string> bigFastaArr;
	vector<string> smallFileArr;
	for(int i = 0; i < fileList.size(); i++){
		if(fileList[i].fileSize > limitSize)
			bigFastaArr.push_back(fileList[i].fileName);
		else
			smallFileArr.push_back(fileList[i].fileName);
	}
	int numBigFasta = bigFastaArr.size();
	cerr << "the total fileNumber is: " << fileList.size() << endl;
	cerr << "the big fasta file is: " << numBigFasta << endl;
	cerr << "the small fileNumber is: " << smallFileArr.size() << endl;
 
 	if(numBigFasta != 0){
		int numConsumer = 7;
		int numGroup = numThreads / (numConsumer + 1);
		if(numBigFasta < numGroup){
			numGroup = numBigFasta;
			numConsumer = numThreads / numGroup - 1;
		}
		cerr << "the numThreads is: " << numThreads << endl;
		cerr << "the numGroup for sketching fasta is: " << numGroup << endl;
		cerr << "the numConsumer for fasta is: " << numConsumer << endl;
		#pragma omp parallel for num_threads(numGroup)
		for(int i = 0; i < bigFastaArr.size(); i++){
			//vector<uint32_t> hashArr[numConsumer];
			robin_hood::unordered_set<uint32_t> hashValueSets[numConsumer];
			FXReader<FA> m_reader(bigFastaArr[i]);
			std::thread ** threads =new std::thread*[numConsumer];
			//std::thread ** threads =new std::thread*[numThreads];
			for(int t = 0; t < numConsumer; t++){
				threads[t] = new std::thread(std::bind(consumer_fasta_task, std::ref(m_reader), std::ref(parameter), std::ref(shuffled_map), std::ref(hashValueSets[t]))); 
			}
			m_reader.join_producer();
			for(int t = 0; t < numConsumer; t++){
				threads[t]->join();
			}
			sketch_t tmpSketch;
			unordered_set<uint32_t> finalSet;
			vector<uint32_t> finalHashArr;
			for(int i = 0; i < numConsumer; i++){
				for(auto x : hashValueSets[i]){
					finalSet.insert(x);	
				}
			}
			for(auto x : finalSet){
				finalHashArr.push_back(x);
			}

			//std::sort(finalHashArr.begin(), finalHashArr.end(), cmp);
			tmpSketch.fileName = bigFastaArr[i];
			tmpSketch.id = i;
			tmpSketch.hashSet = finalHashArr;
			#pragma omp critical
			{
			sketches.push_back(tmpSketch);
			}
		}
		cerr << "finished the sketching of the big fasta genomes" << endl;
	}

	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int t = 0; t <smallFileArr.size(); t++)
	{
		int tid = omp_get_thread_num();
		sketch_t tmpSketch;
		gzFile fp1;
		kseq_t * ks1;
		fp1 = gzopen(smallFileArr[t].c_str(), "r");
		if(fp1 == NULL){
			err(errno, "cannot open the genome file: %s\n", smallFileArr[t].c_str());
		}
		ks1 = kseq_init(fp1);
		uint64_t totalLength = 0;
		tmpSketch.fileName = smallFileArr[t];

		unordered_set<uint32_t> hashValueSet;
		//int time = 0;
		while(1)
		{
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;
			string name = ks1->name.s;
			string comment = ks1->comment.s;
			//char * sequence = ks1->seq.s;
			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int keyCount = 0;
			uint64_t base = 1;
			//cerr << "the length is: " << length << endl;
			//slide window to generate k-mers and get hashes
			for(int i = 0; i < length; i++)
			{
				//char ch = sequence[i];
				char ch = ks1->seq.s[i];
				int curCh = int(ch);
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;

				}
				if(base > kmer_size)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);

					//pfilter = shuffled_dim[dim_id];
					//if((pfilter >= dim_end) || (pfilter < dim_start)){
					//	continue;
					//}
					
					if(shuffled_map.count(dim_id) == 0){
						continue;
					}
					pfilter = shuffled_map[dim_id];

					pfilter -= dim_start;
					//dr_tuple = (((uni_tuple & undomask0) + ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) + pfilter; 
					////only when the dim_end is 4096(the pfilter is 12bit in binary and the lowerst 12bit is all 0 
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 
					
					hashValueSet.insert(dr_tuple);
				}//end if, i > kmer_size

			}//end for, of a sequence 
			//exit(0);

		}//end while, read the file
		//coList[t] = co;
		vector<uint32_t> hashArr;
		for(auto x : hashValueSet){
			hashArr.push_back(x);
		}

		tmpSketch.id = t + numBigFasta;
		//std::sort(hashArr.begin(), hashArr.end(), cmp);
		tmpSketch.hashSet = hashArr;

		gzclose(fp1);
		kseq_destroy(ks1);

		#pragma omp critical
		{
			sketches.push_back(tmpSketch);
			if(t % 10000 == 0){
				cerr << "finshed sketching: " << t << " genomes" << endl;
			}
		}

	}//end for, the fileList

	//std::sort(sketches.begin(), sketches.end(), cmpSketch);

	if(!isSketchFile(outputFile)){
		outputFile = outputFile + ".sketch";
	}
	sketchInfo_t info;
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	//info.genomeNumber = sketches.size();
	saveSketches(sketches, info, outputFile);
	cerr << "save the sketches into: " << outputFile << endl;

	if(!isQuery){
		double tstart = get_sec();
		string dictFile = outputFile + ".dict";
		string indexFile = outputFile + ".index";
		transSketches(sketches, info, dictFile, indexFile, numThreads);
		double tend = get_sec();
		cerr << "===============the time of transSketches is: " << tend - tstart << endl;
	}

	return true;

}


bool sketchFastqFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, int leastNumKmer, vector<sketch_t>& sketches, string outputFile){
	cerr << "run the sketchFastqFile " << endl;
	int half_k = parameter.half_k;
	int half_subk = parameter.half_subk;
	int drlevel = parameter.drlevel;
	int rev_add_move = parameter.rev_add_move;
	int half_outctx_len = parameter.half_outctx_len;
	int * shuffled_dim = parameter.shuffled_dim;
	int dim_start = parameter.dim_start;
	int dim_end = parameter.dim_end;
	int kmer_size = parameter.kmer_size;
	uint32_t hashSize = parameter.hashSize;
	int hashLimit = parameter.hashLimit;
	uint64_t tupmask = parameter.tupmask;
	uint64_t undomask0 = parameter.undomask0;
	uint64_t undomask1 = parameter.undomask1;
	uint64_t domask = parameter.domask;

	int dim_size = 1 << 4 * (half_k - half_outctx_len);
	robin_hood::unordered_map<uint32_t, int> shuffled_map;
	for(int t = 0; t < dim_size; t++){
		if(shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start){
			shuffled_map.insert({t, shuffled_dim[t]});
		}
		//cout << shuffled_dim[t] << endl;
	}


	ifstream fs(inputFile);
	if(!fs){
		err(errno, "cannot open the inputFile: %s\n", inputFile.c_str());
	}
	vector<fileInfo_t> fileList;
	uint64_t totalSize = 0;
	string fileName;
	while(getline(fs, fileName)){
		struct stat cur_stat;
		stat(fileName.c_str(), &cur_stat);
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
		fileInfo_t tmpF;
		tmpF.fileName = fileName;
		tmpF.fileSize = curSize;
		fileList.push_back(tmpF);
	}
	std::sort(fileList.begin(), fileList.end(), cmpFile);
	uint64_t limitSize = totalSize / numThreads;
	vector<string> bigFastqArr;
	vector<string> smallFileArr;
	for(int i = 0; i < fileList.size(); i++){
		if(fileList[i].fileSize > limitSize)
			bigFastqArr.push_back(fileList[i].fileName);
		else
			smallFileArr.push_back(fileList[i].fileName);
	}
	int numBigFastq = bigFastqArr.size();
	cerr << "the total fileNumber is: " << fileList.size() << endl;
	cerr << "the big fasta file is: " << numBigFastq << endl;
	cerr << "the small fileNumber is: " << smallFileArr.size() << endl;

 	if(numBigFastq != 0){
		int numConsumer = 7;
		int numGroup = numThreads / (numConsumer + 1);
		if(numBigFastq < numGroup){
			numGroup = numBigFastq;
			numConsumer = numThreads / numGroup - 1;
		}
		//numConsumer = 4;
		cerr << "the numThreads is: " << numThreads << endl;
		cerr << "the numGroup for sketching fastq is: " << numGroup << endl;
		cerr << "the numConsumer for fastq is: " << numConsumer << endl;
		#pragma omp parallel for num_threads(numGroup)
		for(int i = 0; i < bigFastqArr.size(); i++){
			//vector<uint32_t> hashArr[numConsumer];
			robin_hood::unordered_map<uint32_t, int> hashValueMaps[numConsumer];
			FXReader<FQ_SE> m_reader(bigFastqArr[i]);
			std::thread ** threads =new std::thread*[numConsumer];
			int thread_idArr[numConsumer];
			for(int t = 0; t < numConsumer; t++){
				thread_idArr[t] = t;
			}
			for(int t = 0; t < numConsumer; t++){
				threads[t] = new std::thread(std::bind(consumer_fastq_task, std::ref(m_reader), std::ref(parameter), std::ref(shuffled_map), std::ref(leastNumKmer), std::ref(hashValueMaps[t]))); 
			}
			m_reader.join_producer();
			for(int t = 0; t < numConsumer; t++){
				threads[t]->join();
			}
			
			sketch_t tmpSketch;
			robin_hood::unordered_map<uint32_t, int> finalHashValueMap;
			vector<uint32_t> finalHashArr;
			for(int i = 0; i < numConsumer; i++){
				for(auto x : hashValueMaps[i]){
					uint32_t key = x.first;
					finalHashValueMap.insert({key, 0});
					finalHashValueMap[key] += x.second;
					//finalMap.insert(x);	
				}
			}
			for(auto x : finalHashValueMap){
				if(x.second >= leastNumKmer){
					finalHashArr.push_back(x.first);
				}
			}
			//std::sort(finalHashArr.begin(), finalHashArr.end(), cmp);
			tmpSketch.fileName = bigFastqArr[i];
			tmpSketch.id = i;
			tmpSketch.hashSet = finalHashArr;
			#pragma omp critical
			{
			sketches.push_back(tmpSketch);
			}
		}
		cerr << "finished the sketching of the big fastq genomes" << endl;
	}
 

	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int t = 0; t <smallFileArr.size(); t++)
	{
		int tid = omp_get_thread_num();
		sketch_t tmpSketch;
		gzFile fp1;
		kseq_t * ks1;
		fp1 = gzopen(smallFileArr[t].c_str(), "r");
		if(fp1 == NULL){
			err(errno, "cannot open the genome file: %s\n", smallFileArr[t].c_str());
		}
		ks1 = kseq_init(fp1);
		uint64_t totalLength = 0;
		tmpSketch.fileName = smallFileArr[t];

		unordered_map<uint32_t, int> hashValueMap;
		//int time = 0;
		while(1)
		{
			int length = kseq_read(ks1);
			if(length < 0){
				break;
			}
			totalLength += length;
			string name = ks1->name.s;
			string comment = ks1->comment.s;
			//char * sequence = ks1->seq.s;
			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
			int keyCount = 0;
			uint64_t base = 1;
			//cerr << "the length is: " << length << endl;
			//slide window to generate k-mers and get hashes
			for(int i = 0; i < length; i++)
			{
				//char ch = sequence[i];
				char ch = ks1->seq.s[i];
				int curCh = int(ch);
				int basenum = BaseMap[(int)ch];
				if(basenum != -1)
				{
					tuple = ((tuple << 2) | basenum) & tupmask;
					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_add_move); 
					base++;
				}
				else{
					base = 1;

				}
				if(base > kmer_size)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);

					//pfilter = shuffled_dim[dim_id];
					//if((pfilter >= dim_end) || (pfilter < dim_start)){
					//	continue;
					//}
					
					if(shuffled_map.count(dim_id) == 0){
						continue;
					}
					pfilter = shuffled_map[dim_id];

					pfilter -= dim_start;
					//dr_tuple = (((uni_tuple & undomask0) + ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) + pfilter; 
					////only when the dim_end is 4096(the pfilter is 12bit in binary and the lowerst 12bit is all 0 
					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 
					
					//hashValueSet.insert(dr_tuple);
					hashValueMap.insert({dr_tuple, 0});
					hashValueMap[dr_tuple]++;
				}//end if, i > kmer_size

			}//end for, of a sequence 
			//exit(0);

		}//end while, read the file
		//coList[t] = co;
		vector<uint32_t> hashArr;
		for(auto x : hashValueMap){
			if(x.second >= leastNumKmer){
				hashArr.push_back(x.first);
			}
		}

		tmpSketch.id = t + numBigFastq;
		//std::sort(hashArr.begin(), hashArr.end(), cmp);
		tmpSketch.hashSet = hashArr;

		gzclose(fp1);
		kseq_destroy(ks1);

		#pragma omp critical
		{
			sketches.push_back(tmpSketch);
			if(t % 10000 == 0){
				cerr << "finshed sketching: " << t << " genomes" << endl;
			}
		}

	}//end for, the fileList

	//std::sort(sketches.begin(), sketches.end(), cmpSketch);

	if(!isSketchFile(outputFile)){
		outputFile = outputFile + ".sketch";
	}
	sketchInfo_t info;
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	//info.genomeNumber = sketches.size();
	saveSketches(sketches, info, outputFile);
	cerr << "save the sketches into: " << outputFile << endl;

	if(!isQuery){
		double tstart = get_sec();
		string dictFile = outputFile + ".dict";
		string indexFile = outputFile + ".index";
		transSketches(sketches, info, dictFile, indexFile, numThreads);
		double tend = get_sec();
		cerr << "===============the time of transSketches is: " << tend - tstart << endl;
	}

	return true;

}



void transSketches(vector<sketch_t>& sketches, sketchInfo_t info, string dictFile, string indexFile, int numThreads){
	double t0 = get_sec();
	int half_k = info.half_k;
	int drlevel = info.drlevel;
	size_t hashSize = 1 << (4 * (half_k - drlevel));
	vector<vector<int>> hashMapId;
	for(int i = 0; i < hashSize; i++){
		hashMapId.push_back(vector<int>());
	}
	int * offsetArr = (int*)calloc(hashSize, sizeof(int));

	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "start the hashMapId " << endl;
	cerr << "the thread number is: " << numThreads << endl;
	for(int i = 0; i < sketches.size(); i++){
		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			int hash = sketches[i].hashSet[j];
			hashMapId[hash].push_back(i);
		}
	}
	cerr << "finish the hashMapId " << endl;
	double tt0 = get_sec();
	cerr << "the time of generate the idx by multiple threads are: " << tt0 - t0 << endl;

	FILE * fp0 = fopen(dictFile.c_str(), "w+");
	uint64_t totalIndex = 0;
	for(int hash = 0; hash < hashSize; hash++){
		if(hashMapId[hash].size() != 0){
			fwrite(hashMapId[hash].data(), sizeof(int), hashMapId[hash].size(), fp0);
			totalIndex += hashMapId[hash].size();
			offsetArr[hash] = hashMapId[hash].size();
		}
	}
	fclose(fp0);
	
	double t1 = get_sec();
	cerr << "the time of merge multiple idx into final hashMap is: " << t1 - tt0 << endl;

	FILE * fp1 = fopen(indexFile.c_str(), "w+");
	fwrite(&hashSize, sizeof(size_t), 1, fp1);
	fwrite(&totalIndex, sizeof(uint64_t), 1, fp1);
	fwrite(offsetArr, sizeof(int), hashSize, fp1);
	double t2 = get_sec();
	fclose(fp1);
	cerr << "the time of write output file is: " << t2 - t1 << endl;
	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "the totalIndex is: " << totalIndex << endl;

}


void saveSketches(vector<sketch_t>& sketches, sketchInfo_t info, string outputFile){

	FILE * fp = fopen(outputFile.c_str(), "w+");
	int sketchNumber = sketches.size();
	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	//uint64_t totalNumber = 0;
	//uint64_t totalLength = 0;
	for(int i = 0; i < sketchNumber; i++){
		genomeNameSize[i] = sketches[i].fileName.length();
		hashSetSize[i] = sketches[i].hashSet.size();
		//totalNumber += hashSetSize[i];
		//totalLength += genomeNameSize[i];
	}
	//cerr << "the sketches size is: " << sketchNumber << endl;
	//cerr << "the total hash number is: " << totalNumber << endl;
	//cerr << "the total name length is: " << totalLength << endl;
	//fwrite(&parameter, sizeof(kssd_parameter_t), 1, fp);
	info.genomeNumber = sketchNumber;
	fwrite(&info, sizeof(sketchInfo_t), 1, fp);
	//fwrite(&sketchNumber, sizeof(int), 1, fp);
	fwrite(genomeNameSize, sizeof(int), sketchNumber, fp);
	fwrite(hashSetSize, sizeof(int), sketchNumber, fp);
	for(int i = 0; i < sketchNumber; i++){
		const char * namePoint = sketches[i].fileName.c_str();
		fwrite(namePoint, sizeof(char), genomeNameSize[i], fp);
		uint32_t * curPoint = sketches[i].hashSet.data();
		fwrite(curPoint, sizeof(uint32_t), hashSetSize[i], fp);
	}
	fclose(fp);
	delete genomeNameSize;
	delete hashSetSize;

}

void readSketches(vector<sketch_t>& sketches, sketchInfo_t& info, string inputFile){
	FILE * fp = fopen(inputFile.c_str(), "rb+");
	fread(&info, sizeof(sketchInfo_t), 1, fp);
	int sketchNumber = info.genomeNumber;
	//fread(&sketchNumber, sizeof(int), 1, fp);
	//cerr << "sketchNumber is: " << sketchNumber << endl;
	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	fread(genomeNameSize, sizeof(int), sketchNumber, fp);
	fread(hashSetSize, sizeof(int), sketchNumber, fp);

	//uint64_t totalNumber = 0;
	//uint64_t totalLength = 0;
	
	int maxNameLength = 1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	for(int i = 0; i < sketchNumber; i++){
		//read the genome name.
		int curLength = genomeNameSize[i];
		if(curLength > maxNameLength){
			maxNameLength = curLength;
			curName = new char[maxNameLength+1];
		}
		//totalLength += curLength;
		int nameLength = fread(curName, sizeof(char), curLength, fp);
		if(nameLength != curLength){
			cerr << "error: the read nameLength is not equal to the saved nameLength, exit!" << endl;
			exit(0);
		}
		string genomeName;
		genomeName.assign(curName, curName + curLength);

		//read the hash values in each sketch
		int curSize = hashSetSize[i];
		//totalNumber += curSize;
		if(curSize > maxHashSize){
			maxHashSize = curSize;
			curPoint = new uint32_t[maxHashSize];
		}
		int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp);
		if(hashSize != curSize){
			cerr << "error: the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
			exit(0);
		}
		vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
		sketch_t s;
		s.fileName = genomeName;
		s.id = i;
		s.hashSet=curHashSet;
		sketches.push_back(s);
	}
	delete [] curPoint;
	delete [] curName;
	fclose(fp);
	//cerr << "the total number of hash value of: " << inputFile << " is: " << totalNumber << endl;
	//cerr << "the total length of genome name of: " << inputFile << " is: " << totalLength << endl;

}

void printInfos(vector<sketch_t>& sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	fprintf(fp, "the number of sketches are: %d\n", sketches.size());
	for(int i = 0; i < sketches.size(); i++){
		fprintf(fp, "%s\t%d\n", sketches[i].fileName.c_str(), sketches[i].hashSet.size());
	}
	fclose(fp);
}

void printSketches(vector<sketch_t>& sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	fprintf(fp, "the number of sketches are: %d\n", sketches.size());
	for(int i = 0; i < sketches.size(); i++){
		fprintf(fp, "%s\t%d\n", sketches[i].fileName.c_str(), sketches[i].hashSet.size());
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			fprintf(fp, "%u\t", sketches[i].hashSet[j]);
			if(j % 10 == 9) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void convertSketch(vector<sketch_t>& sketches, sketchInfo_t& info, string inputDir, int numThreads){
	string stateFile = inputDir + '/' + "cofiles.stat";
	string indexFile = inputDir + '/' + "combco.index.0";
	string sketchFile = inputDir + '/' + "combco.0";

	FILE* fp_stat = fopen(stateFile.c_str(), "r");
	if(!fp_stat){
		cerr << "cannot open: " << stateFile << endl;
		exit(1);
	}
	co_dstat_t curStat;
	fread(&curStat, sizeof(co_dstat_t), 1, fp_stat);
	int shuf_id = curStat.shuf_id;
	int infile_num = curStat.infile_num;
	int kmerSize = curStat.kmerlen;
	int dim_rd_len = curStat.dim_rd_len;
	info.genomeNumber = infile_num;
	info.half_k = kmerSize / 2;
	info.half_subk = 6;
	info.drlevel = dim_rd_len / 2;
	uint64_t all_ctx_ct = curStat.all_ctx_ct;

	//cerr << "sizeof size_t is: " << sizeof(size_t) << endl;
	//cerr << "sizeof bool is: " << sizeof(bool) << endl;
	//cerr << "sizeof co_dstat_t is: " << sizeof(co_dstat_t) << endl;
	//cerr << "shuf_id is: " << shuf_id << endl;
	//cerr << "infile_num is: " << infile_num << endl;
	//cerr << "kmerSize is: " << kmerSize << endl;
	//cerr << "dim_rd_len is: " << dim_rd_len << endl;
	//cerr << "all_ctx_ct is: " << all_ctx_ct << endl;

	uint32_t* tmp_ctx_ct = (uint32_t*)malloc(sizeof(uint32_t) * infile_num);
	fread(tmp_ctx_ct, sizeof(uint32_t), infile_num, fp_stat);
	char ** tmpName = (char**)malloc(infile_num *sizeof(char*));
	for(int i = 0; i < infile_num; i++){
		tmpName[i] = (char*)malloc(PATHLEN * sizeof(char));
		fread(tmpName[i], sizeof(char), PATHLEN, fp_stat);
	}
	fclose(fp_stat);

	size_t* cbdcoindex = (size_t*)malloc(sizeof(size_t) * (infile_num+1));
	FILE* fp_index = fopen(indexFile.c_str(), "rb");
	if(!fp_index){
		cerr << "cannot open: " << indexFile << endl;
		exit(1);
	}
	fread(cbdcoindex, sizeof(size_t), infile_num + 1, fp_index);
	fclose(fp_index);

	FILE* fp_sketch = fopen(sketchFile.c_str(), "rb");
	if(!fp_sketch){
		cerr << "cannot open: " << sketchFile << endl;
		exit(1);
	}

	int fd;
	struct stat statbuf;
	stat(sketchFile.c_str(), &statbuf);
	size_t fileSize = statbuf.st_size;
	size_t hashNumber = fileSize / sizeof(uint32_t);
	if(hashNumber != all_ctx_ct){
		cerr << "the total hash number is not match to the state info, exit" << endl;
		exit(1);
	}
	//cerr << "the fileSize is: " << fileSize << endl;
	//cerr << "the hashNumber is: " << hashNumber << endl;

	int maxArrSize = 1 << 24;
	uint32_t* curSketchArr = (uint32_t*)malloc(maxArrSize * sizeof(uint32_t));
	for(int i = 0; i < infile_num; i++){
		string curFileName = tmpName[i];
		size_t curHashNumber = cbdcoindex[i+1] - cbdcoindex[i];
		if(curHashNumber > maxArrSize){
			maxArrSize = curHashNumber;
			//curSketchArr = (uint32_t*)realloc(curSketchArr, maxArrSize * sizeof(uint32_t));
			curSketchArr = (uint32_t*)malloc(maxArrSize * sizeof(uint32_t));
		}
		vector<uint32_t> curHashArr;
		fread(curSketchArr, sizeof(uint32_t), curHashNumber, fp_sketch);
		for(size_t k = 0; k < curHashNumber; k++){
			curHashArr.emplace_back(curSketchArr[k]);
		}
		sketch_t tmpSketch;
		tmpSketch.fileName = curFileName;
		tmpSketch.hashSet = curHashArr;
		tmpSketch.id = i;
		sketches.push_back(tmpSketch);

	}
	fclose(fp_sketch);
	free(curSketchArr);
	free(tmp_ctx_ct);
	free(tmpName);
	free(cbdcoindex);

}


void convert_from_RabbitKSSDSketch_to_KssdSketch(vector<sketch_t>& sketches, sketchInfo_t& info, string outputDir, int numThreads){
	string command0 = "mkdir -p " + outputDir;
	system(command0.c_str());
	string stateFile = outputDir + '/' + "cofiles.stat";
	string indexFile = outputDir + '/' + "combco.index.0";
	string sketchFile = outputDir + '/' + "combco.0";

	int genomeNumber = sketches.size();
	if(genomeNumber != info.genomeNumber){
		cerr << "mismatch of sketches sizes and info genome number " << endl;
		exit(1);
	}
	//for stateFile
	uint64_t all_ctx_ct = 0;
	uint32_t* tmp_ctx_ct = (uint32_t*)malloc(sizeof(uint32_t) * genomeNumber);
	char ** tmpName = (char**)malloc(genomeNumber * sizeof(char*));

	//for index file
	size_t * cbdcoindex = (size_t*)calloc(sizeof(size_t), genomeNumber+1);

	FILE* fp_sketch = fopen(sketchFile.c_str(), "w");
	if(!fp_sketch){
		cerr << "cannot open: " << sketchFile << endl;
		exit(1);
	}

	for(int i = 0; i < genomeNumber; i++){
		uint32_t* curPoint = sketches[i].hashSet.data();
		fwrite(curPoint, sizeof(uint32_t), sketches[i].hashSet.size(), fp_sketch);
		string curFileName = sketches[i].fileName;
		tmpName[i] = (char*)malloc(PATHLEN * sizeof(char));
		memcpy(tmpName[i], curFileName.data(), curFileName.length()+1);
		//printf("%s\n", tmpName);
		tmp_ctx_ct[i] = sketches[i].hashSet.size();
		all_ctx_ct += tmp_ctx_ct[i];
		cbdcoindex[i+1] = cbdcoindex[i] + sketches[i].hashSet.size();
	}
	fclose(fp_sketch);

	FILE* fp_index = fopen(indexFile.c_str(), "w+");
	if(!fp_index){
		cerr << "cannot open: " << indexFile << endl;
		exit(1);
	}
	fwrite(cbdcoindex, sizeof(size_t), genomeNumber+1, fp_index);
	fclose(fp_index);
	
	//for state file: 
	//1. co_dstat_t curStat;
	//2. uint32_t tmp_ctx_ct[infile_num]
	//3. char tmp_Name[infile_num][PATHLEN]
	
	FILE* fp_stat = fopen(stateFile.c_str(), "w+");
	if(!fp_stat){
		cerr << "cannot open: " << stateFile << endl;
		exit(1);
	}
	co_dstat_t curStat;
	curStat.shuf_id = 348842630;
	curStat.koc = false;
	curStat.kmerlen = info.half_k * 2;
	curStat.dim_rd_len = info.drlevel * 2;
	curStat.comp_num = 1;
	curStat.infile_num = info.genomeNumber;
	curStat.all_ctx_ct = all_ctx_ct; //need to update

	fwrite(&curStat, sizeof(co_dstat_t), 1, fp_stat);
	fwrite(tmp_ctx_ct, sizeof(uint32_t), genomeNumber, fp_stat);
	cerr << "start the tmpName fwrite" << endl;
	for(int i = 0; i < genomeNumber; i++){
		fwrite(tmpName[i], sizeof(char), PATHLEN, fp_stat);
	}
	cerr << "finish the tmpName fwrite" << endl;
	fclose(fp_stat);

}














	

