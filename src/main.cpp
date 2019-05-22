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
	MPI_Offset file_size, chunk_size, overlap;
	chunk_size = 64 << 20; // 64 kB
	overlap = 0 << 10; // 1kB
	char *buf = (char*) malloc((chunk_size + overlap + 1) * sizeof(char));

	MPI_File_open(MPI_COMM_WORLD, "dat/wiki_10GB.txt", 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );
	MPI_File_get_size(fh, &file_size);
	int loop_limit = file_size / chunk_size / size;

	for (int i = 0; i < loop_limit; i++) {
		read(&fh, buf, chunk_size, overlap, i);
		// printf("%s\n", buf);
	}

}