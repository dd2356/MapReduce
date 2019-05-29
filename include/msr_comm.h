#include "msr.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <mpi.h>
#include <stdlib.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration, int rank, int size, MPI_Offset file_size);
int communicate(Pair **sendbuf, int **sendcounts, int **sdispl, int **recvcounts, int **recvdispls, 
    MPI_Request *requests, MPI_Request *all_to_all_requests, int idx, int all_to_all_idx, Pair **recvbuf, int n, bool send_counts);
