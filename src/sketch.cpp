#include "sketch.h"
#include <err.h>
#include "kseq.h"
#include <zlib.h>
#include <omp.h>
#include <fstream>

KSEQ_INIT(gzFile, gzread);

using namespace std;

#define H1(K, HASH_SZ) ((K) % (HASH_SZ))
#define H2(K, HASH_SZ) (1 + (K) % ((HASH_SZ)-1))
#define HASH(K, I, HASH_SZ) ((H1(K, HASH_SZ) + I * H2(K, HASH_SZ)) % HASH_SZ)

bool sketchFile(string inputFile, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches){
	int half_k = parameter.half_k;
	int drlevel = parameter.drlevel;
	int rev_add_move = parameter.rev_add_move;
	int half_outctx_len = parameter.half_outctx_len;
	int * shuffled_dim = parameter.shuffled_dim;
	int dim_start = parameter.dim_start;
	int dim_end = parameter.dim_end;
	int kmer_size = parameter.kmer_size;
	int hashSize = parameter.hashSize;
	int hashLimit = parameter.hashLimit;
	uint64_t tupmask = parameter.tupmask;
	uint64_t undomask0 = parameter.undomask0;
	uint64_t undomask1 = parameter.undomask1;
	uint64_t domask = parameter.domask;
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
	vector<string> fileList;
	string fileName;
	while(getline(fs, fileName)){
		fileList.push_back(fileName);
	}
	cerr << "the size of fileList is: " << fileList.size() << endl;
	cerr << "the numThreads is: " << numThreads << endl;

	uint64_t ** coArr = (uint64_t **)malloc(numThreads * sizeof(uint64_t*));
	for(int i = 0; i < numThreads; i++)
	{
		coArr[i] = (uint64_t*)malloc(hashSize * sizeof(uint64_t));
	}
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int t = 0; t <fileList.size(); t++)
	{
		int tid = omp_get_thread_num();
		sketch_t tmpSketch;
		gzFile fp1;
		kseq_t * ks1;
		fp1 = gzopen(fileList[t].c_str(), "r");
		if(fp1 == NULL){
			//fprintf(stderr, "cannot open the genome file: %s\n", fileList[t].c_str());
			//return 1;
			err(errno, "cannot open the genome file: %s\n", fileList[t].c_str());
		}
		ks1 = kseq_init(fp1);
		uint64_t totalLength = 0;
		//uint64_t* co = (uint64_t*)malloc(hashSize * sizeof(uint64_t));
		uint64_t* co = coArr[tid];
		memset(co, 0LLU, hashSize * sizeof(uint64_t));
		size_t foundIndex = fileList[t].find_last_of('/');
		//tmpSketch.fileName = fileList[t].substr(foundIndex+1);
		tmpSketch.fileName = fileList[t];

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
				if(i >= kmer_size-1)
				{
					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
					pfilter = shuffled_dim[dim_id];
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

					if((pfilter >= dim_end) || (pfilter < dim_start)){
						continue;
					}
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
		vector<uint64_t> hashArr;
		for(int k = 0; k < hashSize; k++){
			if(co[k] != 0){
				//cerr << co[k] << endl;
				hashArr.push_back(co[k]);
				#ifdef DIST_INDEX
				hash_index_arr[tid][co[k]].push_back(t);
				#endif
			}
		}
		tmpSketch.id = t;
		tmpSketch.hashSet = hashArr;

		gzclose(fp1);
		kseq_destroy(ks1);

		#pragma omp critical
		{
			sketches.push_back(tmpSketch);
		}

	}//end for, the fileList


	return true;

}

void saveSketches(vector<sketch_t> sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	for(int i = 0; i < sketches.size(); i++){
		for(int j = 0; j < sketches[i].hashSet.size(); j++){
			fprintf(fp, "%llu\t", sketches[i].hashSet[j]);
			if(j % 10 == 9) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		fprintf(fp, "\n");
	}
	fclose(fp);

}


	
