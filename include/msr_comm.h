#include "msr.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <mpi.h>
#include <stdlib.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration);
int communicate(Pair *sendbuf, int *sendcounts, Pair *recvbuf, int n);
