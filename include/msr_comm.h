#include "msr.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <mpi.h>
#include <stdlib.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration, int buffer_idx, int rank, int size, 
    MPI_Offset file_size, MPI_Request *file_requests);

void displacement(int *src, int *dst, int size);
void define_pair_type(MPI_Datatype *type);
int sum(int *arr, int size);
int communicate(Pair **sendbuf, int **sendcounts, int **sdispl, int **recvcounts, int **recvdispls, 
    MPI_Request *requests, MPI_Request *all_to_all_requests, int idx, int all_to_all_idx, Pair **recvbuf, int n, bool send_counts, int rank);
