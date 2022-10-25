#include "shuffle.h"
#include <cstdio>
#include <err.h>


int * read_shuffle_dim(string shuffle_file){
	int * shuffled_dim;
	FILE * shuf_in = fopen(shuffle_file.c_str(), "rb");
	if(shuf_in == NULL){
		err(errno, "cannot read shuffle file: %s\n", shuffle_file.c_str());
	}
	dim_shuffle_t * dim_shuffle = (dim_shuffle_t *)malloc(sizeof(dim_shuffle_t));
	fread(&(dim_shuffle->dim_shuffle_stat), sizeof(dim_shuffle_stat_t), 1, shuf_in);
	int dim_size = 1 << 4* dim_shuffle->dim_shuffle_stat.subk;
	shuffled_dim = (int*) malloc(sizeof(int) * dim_size);
	fread(shuffled_dim, sizeof(int) * dim_size, 1, shuf_in);
	fprintf(stderr, "read from %s\n", shuffle_file.c_str());
	fprintf(stderr, "\tthe half_k is: %d\n", dim_shuffle->dim_shuffle_stat.k);
	fprintf(stderr, "\tthe half_subk is: %d\n", dim_shuffle->dim_shuffle_stat.subk);
	fprintf(stderr, "\tthe drlevel is: %d\n", dim_shuffle->dim_shuffle_stat.drlevel);

	//printf("print the shuffle_dim : \n");
	//for(int i = 0; i < dim_size; i++)
	//	printf("%lx\n", shuffled_dim[i]);
	//exit(0);

	return shuffled_dim;
}

int * generate_shuffle_dim(int half_subk){
	int dim_size = 1 << 4 * half_subk;
	int * shuffled_dim = shuffleN(dim_size, 0);
	shuffled_dim = shuffle(shuffled_dim, dim_size);

	//printf("print the shuffle_dim : \n");
	//for(int i = 0; i < dim_size; i++)
	//	printf("%lx\n", shuffled_dim[i]);
	//exit(0);

	return shuffled_dim;
}

int * shuffleN(int n, int base)
{
	int * arr;
	arr = (int* ) malloc(n * sizeof(int));
	for(int i = 0; i < n; i++){
		arr[i] = i + base;
	}
	
	return shuffle(arr, n);
}

int * shuffle(int arr[], int length)
{
	if(length > RAND_MAX){
		err(errno, "shuffling array length %d must be less than RAND_MAX: %d", length, RAND_MAX);
	}
	//srand(time(NULL));
	srand(23);
	int j, tmp;
	for(int i = length-1; i > 0; i--)
	{
		j = rand() % (i + 1);
		tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}
	
	return arr;
}







