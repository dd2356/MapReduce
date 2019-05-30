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
    MPI_File_read_all(*fh, buf, chunk_size, MPI_CHAR, MPI_STATUS_IGNORE);
	buf[chunk_size + overlap] = '\0';
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

int communicate(Pair **sendbuf, int **sendcounts, int **sdispl, 
    int **recvcounts, int **recvdispls, MPI_Request *requests, 
    MPI_Request *all_to_all_requests, int idx, int all_to_all_idx, 
    Pair **recvbuf, int n, bool send_counts, int rank) {

    MPI_Comm comm = MPI_COMM_WORLD; 
    if (send_counts) {
        printf("communicating sizes on %d from buffer %d ([%d, %d])\n", 
            rank, idx, sendcounts[idx][0], sendcounts[idx][1]);
        MPI_Ialltoall(sendcounts[idx], 1, MPI_INT, 
            recvcounts[idx], 1, MPI_INT, comm, &requests[idx]);
    }
    if (all_to_all_idx != -1) {
        int i = all_to_all_idx;
        MPI_Wait(&requests[i], MPI_STATUS_IGNORE);

        displacement(recvcounts[i], recvdispls[i], n);
        /* allocate recvbuf */
        int size = sum(recvcounts[i], n); 
        free(recvbuf[i]);
        recvbuf[i] = (Pair*)malloc(sizeof(Pair) * size);

        /* create custom type */ 
        MPI_Datatype pair_type; 
        define_pair_type(&pair_type); 
        int completed = 0;
        if (i > 0) {
            MPI_Test(
                &all_to_all_requests[i-1],
                &completed,
                MPI_STATUS_IGNORE
            );
        }
        printf("sending words on %d with buffer "
            "%d (%d, [%d, %d], [%d, %d]) ([%d, %d], [%d, %d]) %d\n", 
            rank, i, size, sendcounts[i][0], sendcounts[i][1], sdispl[i][0], sdispl[i][1],
            recvcounts[i][0], recvcounts[i][1], recvdispls[i][0], recvdispls[i][1], completed);
        MPI_Ialltoallv(sendbuf[i], sendcounts[i], sdispl[i], pair_type, 
            recvbuf[i], recvcounts[i], recvdispls[i], pair_type, 
            comm, &all_to_all_requests[i]);
        printf("communication done on %d\n", rank);
        return size;
    }
    return 0; 
}
