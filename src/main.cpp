#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char **argv) {

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 2) {
		if (rank == 0) {
			printf("Usage: %s <input_filename>", argv[0]);
		}
		exit(1);
	}

	MPI_File fh;
	MPI_Offset file_size, chunk_size, overlap, buffer_size;
	chunk_size = 64 << 20; // 64 kB
	overlap = 0 << 10; // 1kB
	buffer_size = chunk_size + overlap + 1;
	char *buf = (char*) malloc(buffer_size * sizeof(char));
	int *out_counts = (int*) malloc(size * sizeof(int));
	int *out_offsets = (int*) malloc(size * sizeof(int));

	MPI_File_open(MPI_COMM_WORLD, argv[1], 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );
	MPI_File_get_size(fh, &file_size);
	int loop_limit = file_size / chunk_size / size;
	// loop_limit = 10;
	Pair *recvbuf;
	std::unordered_map<Word,long> process_map;

	clock_t start, end;
	double *times = (double*)calloc(5, sizeof(double));
	

	for (int i = 0; i < loop_limit; i++) {
		if (rank == 0) {
			printf("iteration: %d / %d\n", i, loop_limit);
		}
		start = clock();
		read(&fh, buf, chunk_size, overlap, i);
		end = clock(); times[0] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		std::unordered_map<Word,long> words;
		map(buf, buffer_size, words);
		end = clock(); times[1] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		Pair *out_data = (Pair*) malloc(words.size() * sizeof(Pair));
		shuffle(words, size, out_counts, out_offsets, out_data);
		// printf("%s\n", buf);
		end = clock(); times[2] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
	    int buff_size = communicate(out_data, out_counts, &recvbuf, size); 
        // printf("Recieved buffer of size %d\n",buff_size);
		end = clock(); times[3] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		reduce(recvbuf, buff_size, process_map);
		end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	Word max_word;
	long max_count = 0;
	for (auto& it: process_map) {
		if (it.second > max_count) {
			max_count = it.second;
			max_word = it.first;
		}
		// printf("%s -> %ld\n", it.first.word, it.second);
	}
	printf("process %d: %lu words (%s -> %ld)\n", 
		rank, process_map.size(), max_word.word, max_count);
	printf("times: %.2f, %.2f, %.2f, %.2f, %.2f\n", 
		times[0],
		times[1],
		times[2],
		times[3],
		times[4]
	);
}
