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

int u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3);
int u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, int &a, int &b);
setResult_t getJaccard(vector<uint32_t> list1, vector<uint32_t> list2);
setResult_t vgetJaccard(vector<uint32_t> list1, vector<uint32_t> list2);

void index_dist(vector<sketch_t> sketches, string refSketchOut, string outputFile, int kmer_size, double maxDist, int numThreads){
	string indexFile = refSketchOut + ".index";
	string dictFile = refSketchOut + ".dict";
	size_t hashSize;
	uint64_t totalIndex;
	FILE * fp0 = fopen(indexFile.c_str(), "rb");
	fread(&hashSize, sizeof(size_t), 1, fp0);
	fread(&totalIndex, sizeof(uint64_t), 1, fp0);
	int * sketchSizeArr = (int*)malloc(hashSize * sizeof(int));
	fread(sketchSizeArr, sizeof(int), hashSize, fp0);

	int * offset = (int*)malloc(hashSize * sizeof(int));
	uint64_t totalHashNumber = 0;
	for(int i = 0; i < hashSize; i++){
		totalHashNumber += sketchSizeArr[i];
		offset[i] = sketchSizeArr[i];
		if(i > 0) offset[i] += offset[i-1];
	}
	if(totalHashNumber != totalIndex){
		cerr << "error of the total hash number" << endl;
		exit(1);
	}
	cerr << "the hashSize is: " << hashSize << endl;
	cerr << "totalIndex is: " << totalIndex << endl;
	cerr << "totalHashNumber is: " << totalHashNumber << endl;
	cerr << "offset[n-1] is: " << offset[hashSize-1] << endl;;

	int * indexArr = (int*)malloc(totalHashNumber * sizeof(int));
	FILE * fp1 = fopen(dictFile.c_str(), "rb");
	fread(indexArr, sizeof(int), totalHashNumber, fp1);

	size_t page_size = sysconf(_SC_PAGESIZE);
	size_t num_pages = sysconf(_SC_PHYS_PAGES);
	size_t total_memory = page_size * num_pages;
	//size_t limit_memory = total_memory / 2;
	size_t betterSize = (size_t) 1 << 31;
	size_t limit_memory = std::min(betterSize, total_memory / 2);
	size_t limit_size = limit_memory / sizeof(int);

	double t0 = get_sec();
	size_t numRef = sketches.size();
	size_t numQuery = sketches.size();
	size_t distSize = numRef * numQuery;
	size_t num_parts = 1;
	size_t numQuery_eachPart = numQuery;
	size_t numQuery_remain = 0;
	if(distSize <= limit_size){
		num_parts = 1;
		numQuery_eachPart = numQuery;
		numQuery_remain = 0;
	}
	else{
		numQuery_eachPart = limit_size / numRef;
		num_parts	= numQuery / numQuery_eachPart;
		numQuery_remain = numQuery % numQuery_eachPart;
	}
	size_t partSize = numQuery_eachPart * numRef;
	cerr << "the numRef is: " << numRef << endl;
	cerr << "the distSize is: " << distSize << endl;
	cerr << "the limit_size is: " << limit_size << endl;
	cerr << "the partSize is: " << partSize << endl;
	cerr << "the num_query is: " << numQuery << endl;
	cerr << "the num_parts is: " << num_parts << endl;
	cerr << "the numQuery_eachPart is: " << numQuery_eachPart << endl;
	cerr << "the numQuery_remain is: " << numQuery_remain << endl;
	//exit(0);


	int * intersection_arr = (int*)malloc(partSize * sizeof(int));

	vector<FILE*> fpArr;
	vector<string> dist_file_list;
	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = outputFile + to_string(i);
		dist_file_list.push_back(tmpName);

		FILE * fp0 = fopen(tmpName.c_str(), "w");
		fpArr.push_back(fp0);
	}
	
	cerr << "before generate the interaection " << endl;
	double tt1 = get_sec();

	for(int it = 0; it < num_parts; it++){
		memset(intersection_arr, 0, partSize * sizeof(int));
		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(size_t i = 0; i < numQuery_eachPart; i++){
			size_t dist_offset = i * numRef;
			size_t qid = it * numQuery_eachPart + i;
			for(size_t j = 0; j < sketches[qid].hashSet.size(); j++){
				int hash = sketches[qid].hashSet[j];
				//if(hash > hashSize) cerr << "hash: " << hash << " hashSize: " << hashSize << endl;
				if(sketchSizeArr[hash] == 0) continue;
				int start = hash > 0 ? offset[hash-1] : 0;
				int end = offset[hash];
				//cout << hash << '\t' << start << '\t' << end << endl;;
				for(size_t k = start; k < end; k++){
					size_t curIndex = dist_offset + (size_t)indexArr[k];
					intersection_arr[curIndex]++;
				}
			}
		}


		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(size_t i = 0; i < numQuery_eachPart; i++){
			int tid = omp_get_thread_num();
			size_t dist_offset = i * numRef;
			size_t qid = it * numQuery_eachPart + i;
			for(size_t j = qid+1; j < numRef; j++){
				int common = intersection_arr[dist_offset + j];
				int size0 = sketches[qid].hashSet.size();
				int size1 = sketches[j].hashSet.size();
				int denom = size0 + size1 - common;
				double jaccard = (double)common / denom;
				double mashD;
				if(jaccard == 1.0)
					mashD = 0.0;
				else if(jaccard == 0.0)
					mashD = 1.0;
				else
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				if(mashD < maxDist)
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[qid].fileName.c_str(), sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
		}

	}//end for numParts

		memset(intersection_arr, 0, partSize * sizeof(int));
		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(size_t i = 0; i < numQuery_remain; i++){
			size_t dist_offset = i * numRef;
			size_t qid = num_parts * numQuery_eachPart + i;
			for(size_t j = 0; j < sketches[qid].hashSet.size(); j++){
				int hash = sketches[qid].hashSet[j];
				//if(hash > hashSize) cerr << "hash: " << hash << " hashSize: " << hashSize << endl;
				if(sketchSizeArr[hash] == 0) continue;
				int start = hash > 0 ? offset[hash-1] : 0;
				int end = offset[hash];
				//cout << hash << '\t' << start << '\t' << end << endl;;
				for(size_t k = start; k < end; k++){
					size_t curIndex = dist_offset + (size_t)indexArr[k];
					intersection_arr[curIndex]++;
				}
			}
		}

		#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
		for(size_t i = 0; i < numQuery_remain; i++){
			int tid = omp_get_thread_num();
			size_t dist_offset = i * numRef;
			size_t qid = num_parts * numQuery_eachPart + i;
			for(size_t j = qid+1; j < numRef; j++){
				int common = intersection_arr[dist_offset + j];
				int size0 = sketches[qid].hashSet.size();
				int size1 = sketches[j].hashSet.size();
				int denom = size0 + size1 - common;
				double jaccard = (double)common / denom;
				double mashD;
				if(jaccard == 1.0)
					mashD = 0.0;
				else if(jaccard == 0.0)
					mashD = 1.0;
				else
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				if(mashD < maxDist)
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[qid].fileName.c_str(), sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
		}

	cerr << "finished multithread computing" << endl;

	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
	}
	cerr << "finished fpArr fclose" << endl;
	free(intersection_arr);

	double t1 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

	FILE * cofp;
	FILE * com_cofp = fopen(outputFile.c_str(), "w");
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

void tri_dist(vector<sketch_t> sketches, string outputFile, int kmer_size, double maxDist, int numThreads){
	double t0 = get_sec();
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
	for(int i = 0; i < sketches.size(); i++)
	{
		int tid = omp_get_thread_num();
		for(int j = i+1; j < sketches.size(); j++)
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
				fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[i].fileName.c_str(), sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
		}
	}
	//close the subResultFile
	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
	}
	double t1 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

	struct stat cof_stat;
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
		//if((cofp = fopen(dist_file_list[i].c_str(), "r") == NULL)){
		//	err(errno,"cannot open %s\n", dist_file_list[i].c_str());
		//}

		//stat(dist_file_list[i].c_str(), &cof_stat);
		//uint64_t tmpSize = cof_stat.st_size;
		//cerr << "the size of " << dist_file_list[i] << " is: " << tmpSize << endl;
		//char * tmpco = (char*)malloc((tmpSize+1) * sizeof(char));
		while(1)
		{
			lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
			//cerr << "the lengthRead is: " << lengthRead << endl;
			fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
			if(lengthRead < bufSize) break;
		}
		fclose(cofp);

		remove(dist_file_list[i].c_str());

		//fread(tmpco, sizeof(char), tmpSize, cofp);
		//fwrite(tmpco, sizeof(char), tmpSize, com_cofp);
		//fclose(cofp);
		//free(tmpco);
	}

	free(bufRead);
	fclose(com_cofp);
	double t2 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t2 - t1 << endl;

}



void dist(vector<sketch_t> ref_sketches, vector<sketch_t> query_sketches, string outputFile, int kmer_size, double maxDist, int numThreads){
	double t0 = get_sec();
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
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", ref_sketches[i].fileName.c_str(), query_sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
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
					fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", ref_sketches[j].fileName.c_str(), query_sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
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

	struct stat cof_stat;
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
		//if((cofp = fopen(dist_file_list[i].c_str(), "r") == NULL)){
		//	err(errno,"cannot open %s\n", dist_file_list[i].c_str());
		//}

		//stat(dist_file_list[i].c_str(), &cof_stat);
		//uint64_t tmpSize = cof_stat.st_size;
		//cerr << "the size of " << dist_file_list[i] << " is: " << tmpSize << endl;
		//char * tmpco = (char*)malloc((tmpSize+1) * sizeof(char));
		while(1)
		{
			lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
			//cerr << "the lengthRead is: " << lengthRead << endl;
			fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
			if(lengthRead < bufSize) break;
		}
		fclose(cofp);

		remove(dist_file_list[i].c_str());

		//fread(tmpco, sizeof(char), tmpSize, cofp);
		//fwrite(tmpco, sizeof(char), tmpSize, com_cofp);
		//fclose(cofp);
		//free(tmpco);
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


int u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, int &a, int &b){
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
		int a = 0; 
		int b = 0;
		//int *i_a = &a;
		//int *i_b = &b;
		////int * i_a, * i_b;
    //*i_a = 0;
    //*i_b = 0;
    uint64_t st_a = (size1 / 8) * 8;
    uint64_t st_b = (size2 / 8) * 8;

    int i_a_s, i_b_s;

    if(size3 <= 16){
        count += u32_intersect_scalar_stop(list1, size1, list2, size2, size3, a, b);
        return count;
    }

    uint64_t stop = size3 - 16;
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












