#include "msr_comm.h"
#include <stdio.h>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include <algorithm>

// descending sort
bool sort_words(Pair p1, Pair p2) {
	return p1.count > p2.count;
}

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
	chunk_size = 64 << 20; // 64 MB
	overlap = 0 << 20; // 2 MB
	buffer_size = chunk_size + overlap + 1;
	char *buf = (char*) malloc(buffer_size * sizeof(char));
	int *out_counts = (int*) malloc(size * sizeof(int));
	int *out_offsets = (int*) malloc(size * sizeof(int));

	MPI_File_open(MPI_COMM_WORLD, argv[1], 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );
	MPI_File_get_size(fh, &file_size);
	int loop_limit = file_size / chunk_size / size;
	loop_limit += (loop_limit == 0);
	Pair *recvbuf;
	std::unordered_map<Word,long> process_map;

	clock_t start, end;
	double *times = (double*)calloc(5, sizeof(double));
	// disable buffering for stdout
	setbuf(stdout, NULL);

	for (int i = 0; i < loop_limit; i++) {
		if (rank == 0) {
			printf("\riteration: %d / %d", i+1, loop_limit);
		}
		start = clock();
		read(&fh, buf, chunk_size, overlap, i);
		end = clock(); times[0] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		std::unordered_map<Word,long> words;
		map(buf, chunk_size, overlap, words);
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
	if (rank == 0) {
		printf("\n");
	}
	Word max_word;
	long max_count = 0;
	long total_chars = 0;
	std::vector<Pair> all_pairs;
	for (auto& it: process_map) {
		if (it.second > max_count) {
			max_count = it.second;
			max_word = it.first;
		}
		Pair p;
		memcpy(p.word, it.first.word, WORD_SIZE);
		p.count = it.second;
		all_pairs.push_back(p);
		int len = 0;
		while (it.first.word[++len] != '\0');
		total_chars += len * it.second;
		// printf("%s -> %ld\n", it.first.word, it.second);
	}
	// TODO: perform an allgather first, and have 
	// root process write to file after sorting words
	std::sort(all_pairs.begin(), all_pairs.end(), sort_words);
	usleep (10000 * rank);
	int top_10 = all_pairs.size() < 10 ? all_pairs.size() : 10;
	for (int i = 0; i < top_10; i++) {
		printf("Rank %d, top %d: %s -> %ld\n", 
			rank, i, all_pairs[i].word, all_pairs[i].count);
	}

	usleep(100000);
	usleep(10000 * rank);
	printf("times: read: %.2f, map: %.2f, shuffle: %.2f, "
		"communicate: %.2f, reduce: %.2f\n", 
		times[0],
		times[1],
		times[2],
		times[3],
		times[4]
	);
	printf("Total chars: %ld\n", total_chars);
	MPI_Finalize();
}
