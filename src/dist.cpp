#include "dist.h"
#include <omp.h>
#include <cstdio>
#include <math.h>
#include <sys/stat.h>

void tri_dist(vector<sketch_t> sketches, string outputFile, int kmer_size, int numThreads){
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
			setResult_t tmpResult = getJaccard(sketches[i].hashSet, sketches[j].hashSet);
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
			if(mashD != 1.0)
				fprintf(fpArr[tid], " %s\t%s\t%d\t%d\t%d\t%lf\t%lf\n", sketches[i].fileName.c_str(), sketches[j].fileName.c_str(), common, size0, size1, jaccard, mashD);
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
	fprintf(com_cofp, " genome0\tgenome1\tcommon\tdenom\tjaccard\t%mashD\n");
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

	//for(int i = 0; i < numThreads; i++)
	//{
	//	remove(dist_file_list[i].c_str());
	//}
	free(bufRead);
	//fclose(cofp);
	fclose(com_cofp);
	double t2 = get_sec();
	cerr << "the time of merge the subFiles into final files is: " << t2 - t1 << endl;

}


setResult_t getJaccard(vector<uint64_t> list1, vector<uint64_t> list2)
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













