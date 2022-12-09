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

bool isSketchFile(string inputFile){
	int startPos = inputFile.find_last_of('.');
	if(startPos == string::npos) return false;
	string suffixName = inputFile.substr(startPos+1);
	if(suffixName == "sketch")	return true;
	else return false;
}


void consumer_fasta_task(FXReader<FA> &m_reader, kssd_parameter_t& parameter, robin_hood::unordered_map<uint32_t, int>& shuffled_map, vector<uint32_t>& hashArr){

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
	
	uint32_t* co = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
	memset(co, 0, hashSize * sizeof(uint32_t));
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

					uint32_t j, n;
					for(j = 0; j < hashSize; j++)
					{
						n = HASH(dr_tuple, j, hashSize);
						if(co[n] == 0)
						{
							co[n] = dr_tuple;
							keyCount++;
							if(keyCount > hashLimit){
								err(errno, "the context space is too crowd, try rerun the program using -k %d\n", half_k + 1);
							}
							break;
						}
						else if(co[n] == dr_tuple){
							break;
						}
					}//end for, hash computing 
				}
			} //end for, of a sequence
		}//end of all sequence

	}//end while(1)

	for(int k = 0; k < hashSize; k++){
		if(co[k] != 0){
			hashArr.push_back(co[k]);
		}
	}

	//std::sort(hashArr.begin(), hashArr.end(), cmp);
	free(co);
}

void consumer_fastq_task(FXReader<FQ_SE> &m_reader, kssd_parameter_t& parameter, robin_hood::unordered_map<uint32_t, int>& shuffled_map, vector<uint32_t>& hashArr){

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
	
	//cerr << "enter the consumer_fastq_task " << endl;
	uint32_t* co = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
	memset(co, 0, hashSize * sizeof(uint32_t));
	while(1){
		//cerr << "enter the consumer_fastq_task " << endl;
		auto data = m_reader.get_formated_reads();
		if(data.size() == 0) break;
		for(Reference &r: data){
			int length = r.seq.length();
			//cerr << "the length is: " << length << endl;
			//exit(0);
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

					uint32_t j, n;
					for(j = 0; j < hashSize; j++)
					{
						n = HASH(dr_tuple, j, hashSize);
						if(co[n] == 0)
						{
							co[n] = dr_tuple;
							keyCount++;
							if(keyCount > hashLimit){
								err(errno, "the context space is too crowd, try rerun the program using -k %d\n", half_k + 1);
							}
							break;
						}
						else if(co[n] == dr_tuple){
							break;
						}
					}//end for, hash computing 
				}
			} //end for, of a sequence
		}//end of all sequence

	}//end while(1)

	for(int k = 0; k < hashSize; k++){
		if(co[k] != 0){
			hashArr.push_back(co[k]);
		}
	}

	//std::sort(hashArr.begin(), hashArr.end(), cmp);
	free(co);
}


bool sketchFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, string outputFile){
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

	//cerr << "the size of shuffled_map is: " << shuffled_map.size() << endl;
	//exit(0);
	
	//cerr << "hash size is: " << hashSize << endl;
	//exit(0);
	//printf("tupmask: %lx\n", tupmask);
	//printf("undomask0: %lx\n", undomask0);
	//printf("undomask1: %lx\n", undomask1);
	//printf("domask: %lx\n", domask);
	//exit(0);

	ifstream fs(inputFile);
	if(!fs){
		err(errno, "cannot open the inputFile: %s\n", inputFile.c_str());
	}
	vector<fileInfo_t> fileList;
	//vector<uint64_t> fileSizeArr;
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
	vector<string> bigFileArr;
	vector<string> bigFastaArr;
	vector<string> bigFastqArr;
	vector<string> smallFileArr;
	for(int i = 0; i < fileList.size(); i++){
		if(fileList[i].fileSize > limitSize){
			bigFileArr.push_back(fileList[i].fileName);
			if(isFasta(fileList[i].fileName))
				bigFastaArr.push_back(fileList[i].fileName);
			else if(isFastq(fileList[i].fileName))
				bigFastqArr.push_back(fileList[i].fileName);
			else
				err(errno, "the file format of %s is not support, only support fasta and fastq format input\n", fileList[i].fileName.c_str());
		}
		else
			smallFileArr.push_back(fileList[i].fileName);
	}
	int numBigFile = bigFileArr.size();
	int numBigFasta = bigFastaArr.size();
	int numBigFastq = bigFastqArr.size();
	cerr << "the total fileNumber is: " << fileList.size() << endl;
	cerr << "the big fileNumber is: " << numBigFile << endl;
	cerr << "the big fasta file is: " << numBigFasta << endl;
	cerr << "the big fastq file is: " << numBigFastq << endl;
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
			vector<uint32_t> hashArr[numConsumer];
			FXReader<FA> m_reader(bigFastaArr[i]);
			std::thread ** threads =new std::thread*[numConsumer];
			//std::thread ** threads =new std::thread*[numThreads];
			for(int t = 0; t < numConsumer; t++){
				threads[t] = new std::thread(std::bind(consumer_fasta_task, std::ref(m_reader), std::ref(parameter), std::ref(shuffled_map), std::ref(hashArr[t]))); 
			}
			m_reader.join_producer();
			for(int t = 0; t < numConsumer; t++){
				threads[t]->join();
			}
			sketch_t tmpSketch;
			unordered_set<uint32_t> finalMap;
			vector<uint32_t> finalHashArr;
			for(int i = 0; i < numConsumer; i++){
				for(auto x : hashArr[i]){
					finalMap.insert(x);	
				}
			}
			for(auto x : finalMap){
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
 	if(numBigFastq != 0){
		int numConsumer = 7;
		int numGroup = numThreads / (numConsumer + 1);
		if(numBigFastq < numGroup){
			numGroup = numBigFastq;
			numConsumer = numThreads / numGroup - 1;
		}
		cerr << "the numThreads is: " << numThreads << endl;
		cerr << "the numGroup for sketching fastq is: " << numGroup << endl;
		cerr << "the numConsumer for fastq is: " << numConsumer << endl;
		#pragma omp parallel for num_threads(numGroup)
		for(int i = 0; i < bigFastqArr.size(); i++){
			vector<uint32_t> hashArr[numConsumer];
			FXReader<FQ_SE> m_reader(bigFastqArr[i]);
			std::thread ** threads =new std::thread*[numConsumer];
			//std::thread ** threads =new std::thread*[numThreads];
			for(int t = 0; t < numConsumer; t++){
				threads[t] = new std::thread(std::bind(consumer_fastq_task, std::ref(m_reader), std::ref(parameter), std::ref(shuffled_map), std::ref(hashArr[t]))); 
			}
			m_reader.join_producer();
			for(int t = 0; t < numConsumer; t++){
				threads[t]->join();
			}
			sketch_t tmpSketch;
			unordered_set<uint32_t> finalMap;
			vector<uint32_t> finalHashArr;
			for(int i = 0; i < numConsumer; i++){
				for(auto x : hashArr[i]){
					finalMap.insert(x);	
				}
			}
			for(auto x : finalMap){
				finalHashArr.push_back(x);
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

	uint32_t ** coArr = (uint32_t **)malloc(numThreads * sizeof(uint32_t*));
	for(int i = 0; i < numThreads; i++)
	{
		coArr[i] = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
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
		//uint64_t* co = (uint64_t*)malloc(hashSize * sizeof(uint64_t));
		uint32_t* co = coArr[tid];
		memset(co, 0LLU, hashSize * sizeof(uint32_t));
		//size_t foundIndex = smallList[t].find_last_of('/');
		//tmpSketch.fileName = smallList[t].substr(foundIndex+1);
		tmpSketch.fileName = smallFileArr[t];

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

					//printf("\n");
					//cout  << "the pfilter is: " << pfilter << endl;
					//cout << "the basenum is: " << basenum << endl;
					//printf("the base is: %ld\n", base);
					//printf("the unituple is: %lx\n", uni_tuple);
					//printf("the tuple is: %lx\n", tuple);
					//printf("the crvstuple is: %lx\n", rvs_tuple);
					//printf("the dim_tup is: %lx\n", dim_id);
					//printf("\n");
					//exit(0);

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
					//printf("\n");
					//printf("the base is: %ld\n", base);
					//printf("the pfilter is: %lx\n", pfilter);
					//printf("the drtuple is: %lx\n", dr_tuple);
					//printf("the test_drtuple is: %lx\n", test_dr_tuple);
					//time++;
					//if(time > 15)
					//	exit(0);
					//exit(0);
					
					uint32_t j, n;
					for(j = 0; j < hashSize; j++)
					{
						n = HASH(dr_tuple, j, hashSize);
						if(co[n] == 0)
						{
							co[n] = dr_tuple;
							keyCount++;
							if(keyCount > hashLimit){
								err(errno, "the context space is too crowd, try rerun the program using -k %d\n", half_k + 1);
							}
							break;
						}
						else if(co[n] == dr_tuple){
							break;
						}
					}//end for, hash computing 

				}//end if, i > kmer_size

			}//end for, of a sequence 
			//exit(0);

		}//end while, read the file
		//coList[t] = co;
		vector<uint32_t> hashArr;
		for(int k = 0; k < hashSize; k++){
			if(co[k] != 0){
				//cerr << co[k] << endl;
				hashArr.push_back(co[k]);
			}
		}
		tmpSketch.id = t + numBigFile;
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

void transSketches(vector<sketch_t> sketches, sketchInfo_t info, string dictFile, string indexFile, int numThreads){
	double t0 = get_sec();
	int half_k = info.half_k;
	int drlevel = info.drlevel;
	size_t hashSize = 1 << (4 * (half_k - drlevel));
	vector<vector<int>> hashMapId;
	hashMapId.resize(hashSize);
	int * offsetArr = (int*)malloc(hashSize * sizeof(int));

	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "start the hashMapId " << endl;
	for(int i = 0; i < sketches.size(); i++){
		#pragma omp parellel for num_threads(numThreads)
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			int hash = sketches[i].hashSet[j];
			//cerr << "the hash: " << hash << endl;
			hashMapId[hash].emplace_back(i);
		}
	}
	//cerr << "finish the hashMapId " << endl;
	
	double t1 = get_sec();
	cerr << "the time of transpose sketch is: " << t1 - t0 << endl;

	FILE * fp0 = fopen(dictFile.c_str(), "w+");
	uint64_t totalIndex = 0;
	for(size_t i = 0; i < hashSize; i++){
		offsetArr[i] = hashMapId[i].size();
		if(hashMapId[i].size() != 0){
			fwrite(hashMapId[i].data(), sizeof(int), hashMapId[i].size(), fp0);
			totalIndex += hashMapId[i].size();
		}
	}
	fclose(fp0);
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


void saveSketches(vector<sketch_t> sketches, sketchInfo_t info, string outputFile){

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
	for(int i = 0; i < sketchNumber; i++){
		//read the genome name.
		int curLength = genomeNameSize[i];
		//totalLength += curLength;
		char * curName = new char[curLength+1];
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
		uint32_t * curPoint = new uint32_t[curSize];
		int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp);
		if(hashSize != curSize){
			cerr << "error: the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
			exit(0);
		}
		vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
		delete [] curPoint;
		sketch_t s;
		s.fileName = genomeName;
		s.id = i;
		s.hashSet=curHashSet;
		sketches.push_back(s);
	}
	//cerr << "the total number of hash value of: " << inputFile << " is: " << totalNumber << endl;
	//cerr << "the total length of genome name of: " << inputFile << " is: " << totalLength << endl;

}

void printInfos(vector<sketch_t> sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	fprintf(fp, "the number of sketches are: %d\n", sketches.size());
	for(int i = 0; i < sketches.size(); i++){
		fprintf(fp, "%s\t%d\n", sketches[i].fileName.c_str(), sketches[i].hashSet.size());
	}
	fclose(fp);
}

void printSketches(vector<sketch_t> sketches, string outputFile){
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
	//cerr << "the file number is: " << infile_num << endl;
	//cerr << "finish read curStat " << endl;

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
	//cerr << "finish read cbdcoindex" << endl;

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
	//cerr << "the fileSize is: " << fileSize << endl;
	//cerr << "the hashNumber is: " << hashNumber << endl;

	uint32_t* totalSketchArr = (uint32_t*)malloc(hashNumber* sizeof(uint32_t));
	fread(totalSketchArr, sizeof(uint32_t), hashNumber, fp_sketch);
	fclose(fp_sketch);

	#pragma omp parallel for num_threads(numThreads)
	for(int i = 0; i < infile_num; i++){
		string curFileName = tmpName[i];
		vector<uint32_t> curHashArr;
		for(size_t k = cbdcoindex[i]; k < cbdcoindex[i+1]; k++){
			//cerr << k << endl;
			curHashArr.push_back(totalSketchArr[k]);
		}
		//std::sort(curHashArr.begin(), curHashArr.end());
		sketch_t tmpSketch;
		tmpSketch.fileName = curFileName;
		tmpSketch.hashSet = curHashArr;
		tmpSketch.id = i;
		#pragma omp critical
		{
		sketches.push_back(tmpSketch);
		}
	}

	free(totalSketchArr);
	free(tmp_ctx_ct);
	free(tmpName);
	free(cbdcoindex);

	//std::sort(sketches.begin(), sketches.end(), cmpSketch);
	//std::sort(sketches.begin(), sketches.end(), cmpSketchName);
}















	
