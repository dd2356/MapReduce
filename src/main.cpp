#include "msr.h"
#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <map>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char **argv) {

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_File fh;

	MPI_File_open(MPI_COMM_WORLD, "dat/wiki_100k.txt", 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );

}