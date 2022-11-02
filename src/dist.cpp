#include "dist.h"
#include <omp.h>
#include <cstdio>
#include <math.h>
#include <sys/stat.h>
#include <immintrin.h>

int u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3);
int u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, int &a, int &b);
setResult_t getJaccard(vector<uint32_t> list1, vector<uint32_t> list2);
setResult_t vgetJaccard(vector<uint32_t> list1, vector<uint32_t> list2);

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
	cerr << "the time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

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
	cerr << "the time of merge the subFiles into final files is: " << t2 - t1 << endl;

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
	cerr << "the time of multiple threads distance computing and save the subFile is: " << t1 - t0 << endl;

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
	cerr << "the time of merge the subFiles into final files is: " << t2 - t1 << endl;

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












