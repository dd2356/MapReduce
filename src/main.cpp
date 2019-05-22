#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv) {

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_File fh;
	MPI_Offset file_size, chunk_size, overlap, buffer_size;
	chunk_size = 64 << 10; // 64 kB
	overlap = 0 << 10; // 1kB
	buffer_size = chunk_size + overlap + 1;
	char *buf = (char*) malloc(buffer_size * sizeof(char));
	int *out_counts = (int*) malloc(size * sizeof(int));
	int *out_offsets = (int*) malloc(size * sizeof(int));

	MPI_File_open(MPI_COMM_WORLD, "dat/wiki_100k.txt", 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );
	MPI_File_get_size(fh, &file_size);
	int loop_limit = file_size / chunk_size / size;

	for (int i = 0; i < loop_limit; i++) {
		if (rank == 0) {
			printf("iteration: %d\n", i);
		}
		read(&fh, buf, chunk_size, overlap, i);
		std::unordered_map<Word,long> words;
		map(buf, buffer_size, words);
		Pair *out_data = (Pair*) malloc(words.size() * sizeof(Pair));
		shuffle(words, size, out_counts, out_offsets, out_data);
		// printf("%s\n", buf);
        int buff_size = communicate(out_data, out_counts, out_data, size); 
        printf("Recieved buffer of size %d",buff_size); 
	}

}
