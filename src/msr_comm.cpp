#include "msr_comm.h"
#include <stdio.h>
#include <cstddef>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <mpi.h>

void read(MPI_File *fh, char *buf, MPI_Offset chunk_size, 
	MPI_Offset overlap, int iteration, int rank, int size, 
    MPI_Offset file_size) {
	
	MPI_Offset offset = (iteration * size + rank) * chunk_size;

	file_size--;  // get rid of text file eof

    if (offset > file_size) {
        offset = file_size;
    }
	if (offset + chunk_size + overlap > file_size) {
		chunk_size = file_size - offset;
		overlap = 0;
	}
    printf("reading %lld chars on rank %d\n", chunk_size, rank);
    MPI_File_read_all(*fh, buf, chunk_size, MPI_CHAR, MPI_STATUS_IGNORE);
	buf[chunk_size + overlap] = '\0';
    char temp[41];
    memcpy(temp, buf, 40);
    temp[40] = '\0';
    printf("%s\n", temp);
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
    /* create custom type */ 
    MPI_Datatype pair_type; 
    define_pair_type(&pair_type); 

    MPI_Alltoallv(sendbuf, sendcounts, sdispl, pair_type, 
        *recvbuf, recvcounts, rdispl, pair_type, comm);
    return size; 
}
