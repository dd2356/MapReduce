#include "msr_comm.h"
#include <stdio.h>
#include <cstddef>
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

void displacement(int *src, int *dst, int size) {
    dst[0] = 0; 
    for(int i = 1; i < size; i++) { 
        dst[i]=dst[i-1]+src[i-1]; 
    }
}

int sum(int *arr, int size) { 
    int ans = 0; 
    for (int i = 0; i < size; i++) { 
        ans += arr[i]; 
    }
    return ans; 
}

void define_pair_type(MPI_Datatype *type)  { 
    MPI_Datatype pair_tmp_type; 
    MPI_Aint lb, extent; 

    MPI_Aint pair_displacements[2] = {offsetof(Pair, count), offsetof(Pair, word)};
    MPI_Datatype pair_types[2] = {MPI_LONG,MPI_CHAR}; 
    int blocklen[2] = {1, WORD_SIZE}; 
    
    MPI_Type_create_struct(2,blocklen, pair_displacements, pair_types, &pair_tmp_type); 
    MPI_Type_get_extent(pair_tmp_type, &lb, &extent); 
    MPI_Type_create_resized(pair_tmp_type, lb, extent, type); 
    MPI_Type_commit(type); 
}

int communicate(Pair *sendbuf, int *sendcounts, Pair **recvbuf, int n) {
    MPI_Comm comm = MPI_COMM_WORLD; 
    /* populate recvcounts*/ 
    int recvcounts[n]; 
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, comm); 
    /* calculate offsets */ 
    int sdispl[n], rdispl[n];
    displacement(sendcounts,     sdispl, n); 
    displacement(recvcounts,     rdispl, n); 
    /* allocate recvbuf */
    int size = sum(recvcounts, n); 
    (*recvbuf) = (Pair*)malloc(sizeof(Pair) * size);
    // printf("size: %d\n", size);
    /* create custom type */ 
    MPI_Datatype pair_type; 
    define_pair_type(&pair_type); 

    MPI_Alltoallv(sendbuf, sendcounts, sdispl, pair_type, 
        *recvbuf, recvcounts, rdispl, pair_type, comm);
    // printf("before loop\n");
    // for (int i = 0; i < size; i++) {
        // printf("word: %s -> %ld\n", (*recvbuf[i]).word, (*recvbuf[i]).count);
    // }
    // printf("after loop\n");
    return size; 
}
