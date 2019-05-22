#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

/* Test set_view with DISPLACEMENT_CURRENT */
int main( int argc, char *argv[] )
{
	int errs = 0, err;
	int size, rank, chunk_size;
	char *buf;
	MPI_Offset offset;
	MPI_File fh;
	MPI_Comm comm;
	MPI_Status status;

	chunk_size = 64 << 10;
	buf = malloc((chunk_size + 1) * sizeof(char));

	MPI_Init( &argc, &argv );
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* This test reads a header then sets the view to every "size" int,
		using set view and current displacement. The file is first written
		using a combination of collective and ordered writes */

	comm = MPI_COMM_WORLD;

	/* Reopen the file as sequential */
	err = MPI_File_open( comm, "../wiki_100k.txt", 
		MPI_MODE_RDONLY | MPI_MODE_SEQUENTIAL, MPI_INFO_NULL, &fh );
	if (err)
	{
		printf("failed to read file\n");
		MPI_Abort(MPI_COMM_WORLD, 911);
	}


	MPI_Barrier( comm );

	MPI_Offset filesize;
	/* All processes must provide the same file view for MODE_SEQUENTIAL */
	MPI_File_get_size(fh, &filesize);
	filesize--;  /* get rid of text file eof */
	int globalstart = rank * chunk_size;
	int globalend   = globalstart + chunk_size - 1;
	// if (rank == size-1) {
		// globalend = filesize-1;
	// }

	/* everyone reads in their part */
	err = MPI_File_read_at_all(fh, globalstart, buf, chunk_size, 
		MPI_CHAR, MPI_STATUS_IGNORE);
	// buf[chunk_size] = '\0';
	printf("chunk size: %d, %d, %d, %d, %llu, %d\n", 
		rank, chunk_size, globalstart, globalend, filesize, err);

	if (rank == 1) {
		printf("%s\n", buf);
	}
	free( buf );
	err = MPI_File_close( &fh );
	if (err) { errs++; }
	MPI_Finalize();
	return errs;
}