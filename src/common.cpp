#include "common.h"
#include <err.h>

#define COMPONENT_SZ 7
#define MIN_SUBCTX_DIM_SMP_SZ 4096
#define _64MASK 0xffffffffffffffffLLU
#define CTX_SPC_USE_L 8
#define LD_FCTR 0.6

double get_sec()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000;
}

kssd_parameter_t initParameter(int half_k, int half_subk, int drlevel, int * shuffled_dim){
	// init kssd parameters
	kssd_parameter_t parameter;
	parameter.half_k = half_k;
	parameter.drlevel = drlevel;
	int half_outctx_len = half_k - half_subk;
	parameter.half_outctx_len = half_outctx_len;
	parameter.rev_add_move = 4 * half_k - 2;
	parameter.kmer_size = 2 * half_k;
	std::cerr << "ehhhhhhh" <<  BaseMap['A'] << std::endl;
	
	parameter.dim_start = 0;
	parameter.dim_end = MIN_SUBCTX_DIM_SMP_SZ;
	int hashSize = get_hashSize(half_k, drlevel);
	parameter.hashSize = hashSize;
	parameter.hashLimit = hashSize * LD_FCTR;
	parameter.shuffled_dim = shuffled_dim;

	//int component_num = half_k - drlevel > COMPONENT_SZ ? 1LU << 4 * (half_k - drlevel - COMPONENT_SZ) : 1;
	int comp_bittl = 64 - 4 * half_k;
	//cout << "the comp_bittl is: " << comp_bittl << endl;

	uint64_t tupmask = _64MASK >> comp_bittl;
	uint64_t domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
	uint64_t undomask = (tupmask ^ domask) & tupmask;
	uint64_t undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
	uint64_t undomask0 = undomask ^ undomask1;
	parameter.domask = domask;
	parameter.tupmask = tupmask;
	parameter.undomask1 = undomask1;
	parameter.undomask0 = undomask0;
	//printf("the tupmask is: %lx\n", tupmask);
	//printf("the domask is: %lx\n", domask);
	//printf("the undomask0 is: %lx\n", undomask0);
	//printf("the undomask1 is: %lx\n", undomask1);
	
	return parameter;
}


int get_hashSize(int half_k, int drlevel)
{
	int dim_reduce_rate = 1 << 4 * drlevel;
	uint64_t ctx_space_sz = 1LLU << 4 * (half_k - drlevel);
	int primer_id = 4 * (half_k - drlevel) - CTX_SPC_USE_L - 7;
	if(primer_id < 0 || primer_id > 24)
	{
		int k_add = primer_id < 0 ? (1 + (0 - primer_id) / 4) : -1 * (1 + (primer_id - 24) / 4);
		err(errno, "get_hashSize(): primer_id: %d out of range(0 ~ 24), by formula:\n"
								"int primer_id = 4 * (half_k - drlevel) - CTX_SPC_USE_L - 7;\n"
								"this might caused by too small or too large k\n"
								"half kmer length = %d\n"
								"dim reduction level %d\n"
								"ctx_space size = %llu\n"
								"try rerun the program with option -k = %d\n",
								primer_id, half_k, drlevel, ctx_space_sz, half_k + k_add);

	}
	int hashSize = primer[primer_id];
	fprintf(stderr, "dimension reduced %d\n"
									"ctx_space size  %llu\n"
									"half_k is: %d\n"
									"drlevel is: %d\n"
									"primer_id is: %d\n"
									"hashSize is: %u \n",
									dim_reduce_rate, ctx_space_sz, half_k, drlevel, primer_id, hashSize);
	
	return hashSize;
}

