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

void open_file_with_view(MPI_File &fh, int rank, int size, 
	MPI_Offset chunk_size, char *file_name, MPI_Offset &file_size) {

	MPI_Aint length = chunk_size * sizeof(char);
	MPI_Aint extent = size * length;
	MPI_Offset disp = rank * length;
	MPI_Datatype contig, filetype;
	MPI_Type_contiguous(chunk_size, MPI_CHAR, &contig);
	MPI_Type_create_resized(contig, 0, extent, &filetype);
	MPI_Type_commit(&filetype);
	MPI_File_open(MPI_COMM_WORLD, file_name, 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, disp, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
	MPI_File_get_size(fh, &file_size);
}

void setup_buffers_and_variables(MPI_Offset chunk_size, MPI_Offset overlap, 
	MPI_Offset &buffer_size, char ***buf, int ***out_counts, int ***out_offsets, 
	MPI_Request **requests, MPI_Request **all_to_all_requests, int size, 
	MPI_Offset file_size, int &loop_limit, int &buffers, double **times) {
	
	buffer_size = chunk_size + overlap + 1;
	// needs to be at least triple buffered
	buffers = 4;
	(*buf) = (char**)malloc(sizeof(char*) * buffers);
	(*out_counts) = (int**)malloc(sizeof(int*) * buffers);
	(*out_offsets) = (int**)malloc(sizeof(int*) * buffers);
	(*requests) = (MPI_Request*)calloc(sizeof(MPI_Request), buffers);
	(*all_to_all_requests) = (MPI_Request*)calloc(sizeof(MPI_Request), buffers);
	for (int i = 0; i < buffers; i++) {
		(*buf)[i] = (char*) malloc(buffer_size * sizeof(char));
		(*out_counts)[i] = (int*) malloc(size * sizeof(int));
		(*out_offsets)[i] = (int*) malloc(size * sizeof(int));
	}
	loop_limit = file_size / chunk_size / size;
	loop_limit += (loop_limit * chunk_size * size < file_size);
	(*times) = (double*)calloc(6, sizeof(double));
}

void cleanup_reduce(int loop_limit, int buffers, Pair **out_data,
	int **out_counts, int **out_offsets, int **recv_counts, int **recv_offsets, 
	MPI_Request *requests, MPI_Request *all_to_all_requests, 
	Pair **receive_buffer, int size, int *buff_sizes, double *times,
	std::unordered_map<Word,long> &process_map, int rank) {

	clock_t start, end;
	int idx = loop_limit % buffers;
	int last = (loop_limit-1) % buffers;
	int last_2 = (loop_limit-2) % buffers;

	start = clock();
	int buff_size = communicate(out_data, out_counts, out_offsets, 
		recv_counts, recv_offsets, requests, all_to_all_requests, 
		idx, last, receive_buffer, size, false, rank); 

	if (loop_limit > 0) {
		buff_sizes[last] = buff_size;
	}
	end = clock(); times[3] += ((double) (end - start)) / CLOCKS_PER_SEC;

	start = clock();
	MPI_Wait(&all_to_all_requests[last_2], MPI_STATUS_IGNORE);
	end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
	start = clock();
	reduce(receive_buffer[last_2], buff_sizes[last_2], process_map);
	end = clock(); times[5] += ((double) (end - start)) / CLOCKS_PER_SEC;

	start = clock();
	MPI_Wait(&all_to_all_requests[last], MPI_STATUS_IGNORE);
	end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
	start = clock();
	reduce(receive_buffer[last], buff_sizes[last], process_map);
	end = clock(); times[5] += ((double) (end - start)) / CLOCKS_PER_SEC;	
}

void mapreduce(int loop_limit, int rank, int size, MPI_File fh, char **buf, 
	MPI_Offset chunk_size, MPI_Offset overlap, MPI_Offset file_size,
	double *times, int **out_counts, int **out_offsets, MPI_Request *requests, 
	MPI_Request *all_to_all_requests, int buffers, 
	std::unordered_map<Word,long> &process_map) {
	
	// disable buffering for stdout
	setbuf(stdout, NULL);
	clock_t start, end;
	Pair **receive_buffer = (Pair**)malloc(buffers * sizeof(Pair*));
	Pair **out_data = (Pair**)malloc(buffers * sizeof(Pair*));
	int **recv_counts = (int**)malloc(buffers * sizeof(int*));
	int **recv_offsets = (int**)malloc(buffers * sizeof(int*));
	int *buff_sizes = (int*)malloc(buffers * sizeof(int));
	for (int i = 0; i < buffers; i++) {
		out_data[i] = NULL;
		receive_buffer[i] = NULL;
		recv_counts[i] = (int*)malloc(size * sizeof(int));
		recv_offsets[i] = (int*)malloc(size * sizeof(int));
	}

	int idx, last, last_2;
	for (int i = 0; i < loop_limit; i++) {
		if (rank == 0) {
			printf("\riteration: %d / %d\t\n", i+1, loop_limit);
		}

		idx = i % buffers;
		last = (i-1) % buffers;
		last_2 = (i-2) % buffers;
		// printf("read: %d\n", rank);
		start = clock();
		read(&fh, buf[idx], chunk_size, overlap, i, rank, size, file_size);
		// if (rank == 1) {
			// sleep(1);
		// }
		end = clock(); times[0] += ((double) (end - start)) / CLOCKS_PER_SEC;
		// printf("map\n");
		start = clock();
		std::unordered_map<Word,long> words;
		map(buf[idx], chunk_size, overlap, words);
		printf("mapped %lu words on process %d\n", words.size(), rank);
		end = clock(); times[1] += ((double) (end - start)) / CLOCKS_PER_SEC;
		// printf("shuffle\n");
		start = clock();
		free(out_data[idx]);
		out_data[idx] = (Pair*) malloc(words.size() * sizeof(Pair));
		printf("shuffling on rank %d using buffer %d\n", rank, idx);
		shuffle(words, size, out_counts[idx], out_offsets[idx], out_data[idx]);
		end = clock(); times[2] += ((double) (end - start)) / CLOCKS_PER_SEC;
		// printf("communicate\n");
		start = clock(); 
		int buff_size = communicate(out_data, out_counts, out_offsets, 
			recv_counts, recv_offsets, requests, all_to_all_requests, 
			idx, last, receive_buffer, size, true, rank); 

		if (i > 0) {
			buff_sizes[last] = buff_size;
		}
		end = clock(); times[3] += ((double) (end - start)) / CLOCKS_PER_SEC;

		if (i > 1) {
			start = clock();
			printf("waiting for MPI_Ialltoallv of buffer %d\n", last_2), 
			MPI_Wait(&all_to_all_requests[last_2], MPI_STATUS_IGNORE);
			end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
			printf("reducing %d with %d words\n", rank, buff_sizes[last_2]);
			start = clock();
			reduce(receive_buffer[last_2], buff_sizes[last_2], process_map);
			end = clock(); times[5] += ((double) (end - start)) / CLOCKS_PER_SEC;
		}
		printf("reached barrier on %d\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);
		printf("exited barrier on %d\n", rank);
		sleep(1);
		if (rank == 0) {
			printf("\n");
		}
	}

	cleanup_reduce(loop_limit, buffers, out_data, out_counts, out_offsets, 
		recv_counts, recv_offsets, requests, all_to_all_requests, 
		receive_buffer, size, buff_sizes, times, process_map, rank);

	if (rank == 0) {
		printf("\n");
	}

}

void recap(int rank, std::unordered_map<Word,long> process_map, double *times) {

	Word max_word;
	long max_count = 0;
	std::vector<Pair> all_pairs;

	for (std::pair<const Word, long>& wp: process_map) {
		if (wp.second > max_count) {
			max_count = wp.second;
			max_word = wp.first;
		}
		Pair p;
		memcpy(p.word, wp.first.word, WORD_SIZE);
		p.count = wp.second;
		all_pairs.push_back(p);
	}

	// TODO: perform an allgather first, and have 
	// root process write to file after sorting words
	std::sort(all_pairs.begin(), all_pairs.end(), sort_words);
	usleep (10000 * rank);
	int top_10 = all_pairs.size() < 10 ? all_pairs.size() : 10;
	for (int i = 0; i < top_10; i++) {
		printf("Rank %d, top %2d: %s -> %ld\n", 
			rank, i+1, all_pairs[i].word, all_pairs[i].count);
	}

	usleep(100000);
	usleep(10000 * rank);
	printf("times for process %d (%lu words)\nread: %.2f\tmap: %.2f\t"
		"shuffle: %.2f\tcommunicate: %.2f\treduce wait: %.2f\treduce: %.2f\n", 
		rank, all_pairs.size(),
		times[0], times[1], times[2], times[3], times[4], times[5]
	);
}

int main(int argc, char **argv) {

	int rank, size, loop_limit, buffers;
	MPI_File fh;
	MPI_Offset file_size, chunk_size, overlap, buffer_size;
	char **buf;
	int **out_counts, **out_offsets;
	MPI_Request *requests, *all_to_all_requests;
	std::unordered_map<Word,long> process_map;
	double *times;
	chunk_size = 64 << 20; // 64 MB
	overlap = 0 << 20; // 2 MB

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 2) {
		if (rank == 0) {
			printf("Usage: %s <input_filename>", argv[0]);
		}
		exit(1);
	}

	open_file_with_view(fh, rank, size, chunk_size, argv[1], file_size);

	setup_buffers_and_variables(chunk_size, overlap, buffer_size, 
		&buf, &out_counts, &out_offsets, &requests, &all_to_all_requests, 
		size, file_size, loop_limit, buffers, &times);

	mapreduce(loop_limit, rank, size, fh, buf, chunk_size, overlap, 
		file_size, times, out_counts, out_offsets, requests, 
		all_to_all_requests, buffers, process_map);	


	recap(rank, process_map, times);
	for (int i = 0; i < buffers; i++) {
		free(buf[i]);
		free(out_counts[i]);
		free(out_offsets[i]);		
	}
	free(buf);
	free(out_counts);
	free(out_offsets);
	MPI_Finalize();
}
