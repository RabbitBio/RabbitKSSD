#ifndef _SHUFFLE_H_
#define _SHUFFLE_H_

#include <iostream>
#include <string>

using namespace std;

typedef struct dim_shuffle_stat
{
	int id;
	int k;
	int subk;
	int drlevel;
}	dim_shuffle_stat_t;

typedef struct dim_shuffle
{
	dim_shuffle_stat_t dim_shuffle_stat;
	int * shuffled_dim;
} dim_shuffle_t;

int * shuffleN(int n, int base);
int * shuffle(int arr[], int length);
int * read_shuffle_dim(string shuffle_file);
int * generate_shuffle_dim(int half_subk);


#endif
