#include <stdio.h>
#include <math.h>
#include <map>
#include <mpi.h>
#include <stdlib.h>

struct Pair;
struct Word;

void read(MPI_File *fh, char *buf, int iteration);
void communicate(Pair *sendbuf, int *offsets, Pair *recvbuf);
