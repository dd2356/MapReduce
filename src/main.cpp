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
	MPI_Offset &buffer_size, char **buf, int **out_counts, int **out_offsets, 
	int size, MPI_Offset file_size, int &loop_limit, double **times) {
	buffer_size = chunk_size + overlap + 1;
	(*buf) = (char*) malloc(buffer_size * sizeof(char));
	(*out_counts) = (int*) malloc(size * sizeof(int));
	(*out_offsets) = (int*) malloc(size * sizeof(int));
	loop_limit = file_size / chunk_size / size;
	loop_limit += (loop_limit * chunk_size * size < file_size);
	(*times) = (double*)calloc(5, sizeof(double));
}

void mapreduce(int loop_limit, int rank, int size, MPI_File fh, char *buf, 
	MPI_Offset chunk_size, MPI_Offset overlap, MPI_Offset file_size,
	double *times, int *out_counts, int *out_offsets,
	std::unordered_map<Word,long> &process_map) {
	// disable buffering for stdout
	setbuf(stdout, NULL);
	clock_t start, end;
	Pair *recvbuf;

	for (int i = 0; i < loop_limit; i++) {
		if (rank == 0) {
			printf("\riteration: %d / %d", i+1, loop_limit);
		}
		start = clock(); /*printf("map\n");*/
		read(&fh, buf, chunk_size, overlap, i, rank, size, file_size);
		end = clock(); times[0] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock(); /*printf("map\n");*/
		std::unordered_map<Word,long> words;
		map(buf, chunk_size, overlap, words);
		end = clock(); times[1] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock(); /*printf("map\n");*/
		Pair *out_data = (Pair*) malloc(words.size() * sizeof(Pair));
		shuffle(words, size, out_counts, out_offsets, out_data);
		end = clock(); times[2] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock(); /*printf("map\n");*/
	    int buff_size = communicate(out_data, out_counts, &recvbuf, size); 
		end = clock(); times[3] += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock(); /*printf("map\n");*/
		reduce(recvbuf, buff_size, process_map);
		end = clock(); times[4] += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	if (rank == 0) {
		printf("\n");
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

	int rank, size, loop_limit;
	MPI_File fh;
	MPI_Offset file_size, chunk_size, overlap, buffer_size;
	char *buf;
	int *out_counts, *out_offsets;
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
	// setup_buffers_and_variables();
	open_file_with_view(fh, rank, size, chunk_size, argv[1], file_size);
	setup_buffers_and_variables(chunk_size, overlap, buffer_size, 
		&buf, &out_counts, &out_offsets, 
	    size, file_size, loop_limit, &times);

	mapreduce(loop_limit, rank, size, fh, buf, chunk_size, overlap, 
		file_size, times, out_counts, out_offsets, process_map);	

	recap(rank, size, process_map, times);
	MPI_Finalize();
}
