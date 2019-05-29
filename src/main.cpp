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

// #define DEBUG

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
#ifdef DEBUG
			printf("\riteration: %d / %d\t\n", i+1, loop_limit);
#else
			printf("\riteration: %d / %d\t", i+1, loop_limit);
#endif
		}

		idx = i % buffers;
		last = (i-1) % buffers;
		last_2 = (i-2) % buffers;
#ifdef DEBUG
		printf("read: %d\n", rank);
#endif
		start = clock();
		read(&fh, buf[idx], chunk_size, overlap, i, rank, size, file_size);
		end = clock(); times[0] += ((double) (end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG
		printf("map\n");
#endif
		start = clock();
		std::unordered_map<Word,long> words;
		map(buf[idx], chunk_size, overlap, words);
#ifdef DEBUG
		printf("mapped %lu words on process %d\n", words.size(), rank);
#endif
		end = clock(); times[1] += ((double) (end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG
		printf("shuffle\n");
#endif
		start = clock();
		free(out_data[idx]);
		out_data[idx] = (Pair*) malloc(words.size() * sizeof(Pair));
#ifdef DEBUG
		printf("shuffling on rank %d using buffer %d\n", rank, idx);
#endif
		shuffle(words, size, out_counts[idx], out_offsets[idx], out_data[idx]);
		end = clock(); times[2] += ((double) (end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG
		printf("communicate\n");
#endif
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
#ifdef DEBUG
			printf("waiting for MPI_Ialltoallv of buffer %d\n", last_2), 
#endif
			MPI_Wait(&all_to_all_requests[last_2], MPI_STATUS_IGNORE);
			end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG
			printf("reducing %d with %d words\n", rank, buff_sizes[last_2]);
#endif
			start = clock();
			reduce(receive_buffer[last_2], buff_sizes[last_2], process_map);
			end = clock(); times[5] += ((double) (end - start)) / CLOCKS_PER_SEC;
		}
#ifdef DEBUG
		if (rank == 0) {
			printf("\n");
		}
#endif
	}

	cleanup_reduce(loop_limit, buffers, out_data, out_counts, out_offsets, 
		recv_counts, recv_offsets, requests, all_to_all_requests, 
		receive_buffer, size, buff_sizes, times, process_map, rank);

	if (rank == 0) {
#ifdef DEBUG

		printf("\n");
#endif
	}

}

void recap(int rank, int world_size, std::unordered_map<Word,long> process_map, double *times) {
    printf("Enter: %d\n", rank);

	Word max_word;
	long max_count = 0;
	std::vector<Pair> local_pairs;

	for (auto& it: process_map) {
		if (it.second > max_count) {
			max_count = it.second;
			max_word = it.first;
		}
		Pair p;
		memcpy(p.word, it.first.word, WORD_SIZE);
		p.count = it.second;
		local_pairs.push_back(p);
	}

    // Gather the times
    double total_times[5]; 
    total_times[0] = total_times[1] = total_times[2] = total_times[3] = total_times[4] = 0; 
    double *t_recvbuf; 
    if (rank == 0) t_recvbuf = (double*) malloc(sizeof(double) * world_size * 5); 
    MPI_Gather(times, 5, MPI_DOUBLE, t_recvbuf, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    if (rank == 0) { 
        for (int i = 0; i < world_size * 5; i+=5) { 
            total_times[0] += t_recvbuf[i]; 
            total_times[1] += t_recvbuf[i+1]; 
            total_times[2] += t_recvbuf[i+2]; 
            total_times[3] += t_recvbuf[i+3]; 
            total_times[4] += t_recvbuf[i+4]; 
        }
    }
    if (rank == 0) free(t_recvbuf); 

    // Gather the total chars 
    long total_chars = 0; 
    long tc_recvbuf[world_size]; 
    MPI_Gather(&total_local_chars, 1, MPI_LONG, tc_recvbuf, 1, MPI_LONG, 0, MPI_COMM_WORLD); 
    
    if (rank == 0) { 
      for (int i = 0; i < world_size; i++) total_chars += tc_recvbuf[i];    
    }

    // Gather the size of the pair arrays
    int *recvcounts; 
    int local_size = local_pairs.size(); 
    if (rank == 0) {
    	recvcounts = (int*) malloc(sizeof(int) * world_size); 
    }
    MPI_Gather(&local_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    if (rank == 0) {
    	for (int i = 0; i < world_size; ++i) {
    		printf("sizes: %d\n", local_size);
    	}
    }

    // Gather all pairs
    int displs[world_size]; 
    displacement(recvcounts, displs, world_size);  
    int recvcount = sum(recvcounts, world_size); 
    MPI_Datatype datatype; 
    define_pair_type(&datatype);  
    Pair *recvbuf; 
    if (rank == 0) recvbuf = (Pair*) malloc(sizeof(Pair) * recvcount); 
    MPI_Gatherv(local_pairs.data(), local_size, datatype, recvbuf, recvcounts, displs, datatype, 0, MPI_COMM_WORLD); 
    
    if (rank == 0) {
    	free(recvcounts); 
    }
    if (rank == 0) { 
        // Init a vector from the buffer
        std::vector<Pair> all_pairs(recvbuf, recvbuf + recvcount); 
        std::sort(all_pairs.begin(), all_pairs.end(), sort_words);

        int top_10 = all_pairs.size() < 10 ? all_pairs.size() : 10;
        for (int i = 0; i < top_10; i++) {
            printf("Rank %d, top %d: %s -> %ld\n", 
                    rank, i, all_pairs[i].word, all_pairs[i].count);
        }

        printf("times: read: %.2f, map: %.2f, shuffle: %.2f, "
                "communicate: %.2f, reduce: %.2f\n", 
                total_times[0], total_times[1], total_times[2], total_times[3], total_times[4]
              );
        printf("Total chars: %ld\n", total_chars);
        free (recvbuf); 
    }
    printf("Done: %d\n", rank);
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
	chunk_size = 64 << 20	; // 64 MB
	overlap = 0 << 20; // 2 MB

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 2) {
		if (rank == 0) {
#ifdef DEBUG

			printf("Usage: %s <input_filename>", argv[0]);
#endif
		}
		exit(0);
	}

	open_file_with_view(fh, rank, size, chunk_size, argv[1], file_size);

	setup_buffers_and_variables(chunk_size, overlap, buffer_size, 
		&buf, &out_counts, &out_offsets, &requests, &all_to_all_requests, 
		size, file_size, loop_limit, buffers, &times);

	mapreduce(loop_limit, rank, size, fh, buf, chunk_size, overlap, 
		file_size, times, out_counts, out_offsets, requests, 
		all_to_all_requests, buffers, process_map);	


	recap(rank, size, process_map, times);
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
