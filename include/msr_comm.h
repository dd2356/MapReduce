#include "msr.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <mpi.h>
#include <stdlib.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration, int rank, int size, MPI_Offset file_size);

void displacement(int *src, int *dst, int size);
void define_pair_type(MPI_Datatype *type);
int sum(int *arr, int size);
int communicate(Pair *sendbuf, int *sendcounts, Pair **recvbuf, int n);
