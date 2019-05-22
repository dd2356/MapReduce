#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <mpi.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration) {
	
	int size, rank;
	MPI_Offset offset, filesize;

	// TODO: send these as arguments
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	offset = (iteration * size + rank) * chunk_size;

	MPI_File_get_size(*fh, &filesize);
	filesize--;  /* get rid of text file eof */

	if (offset + chunk_size + overlap > filesize) {
		chunk_size = filesize - offset;
		overlap = 0;
	}
	MPI_File_read_at_all(*fh, offset, buf, chunk_size + overlap, 
		MPI_CHAR, MPI_STATUS_IGNORE);

/*
	int start_offset = 0;
	int end_offset = chunk_size + overlap;
	if (offset > 0) {
		for (int i = 0; i < overlap; i++) {
			if (buf[i] == '\n') {
				start_offset = i;
			}
		}
	}
	for (int i = chunk_size; i < chunk_size + overlap; i++) {
		if (buf[i] == '\n') {
			end_offset = i;
		}
	}
*/
	buf[chunk_size + overlap] = '\0';
	// printf("chunk size: %d, %llu, %llu, %llu, %c, %c\n", 
		// rank, chunk_size, offset, filesize, buf[0], buf[chunk_size-1]);


}

void communicate(Pair *sendbuf, int *offsets, Pair *recvbuf) {

}
