#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <map>
#include <stdlib.h>
#include <mpi.h>

void read(MPI_File *fh, char *buf, int iteration) {
	int size, rank;
	MPI_Offset offset, chunk_size, filesize, overlap;


	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	offset = (iteration * size + rank) * chunk_size;

	MPI_File_get_size(*fh, &filesize);
	filesize--;  /* get rid of text file eof */

	MPI_File_read_at_all(*fh, offset, buf,
                                chunk_size, MPI_CHAR,
                                MPI_STATUS_IGNORE);

	buf[chunk_size] = '\0';
	printf("chunk size: %d, %llu, %llu, %llu\n", 
		rank, chunk_size, offset, filesize);


}

void communicate(Pair *sendbuf, int *offsets, Pair *recvbuf) {

}
