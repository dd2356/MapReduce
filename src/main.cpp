

int main(int argc, char **argv) {


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_File_open(MPI_COMM_WORLD, "dat/wiki_100k.txt", 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );

}