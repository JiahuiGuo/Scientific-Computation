#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv)
{

	int rank, size;
    int row, col;
    int N; // N = sqrt(P)
    MPI_Comm row_Comm, col_Comm;

	if(argc > 1)
		N = atoi(argv[1]);
	else
	{		
		printf("Usage: ./grid sqrt(P)");
		exit(1);
	}
	// Initial the MPI environment
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Determine Row and Col position
	row = rank / N;
	col = rank % N;

	// MPI_COMM_SPLIT( comm, color, key, newcomm )
	// Split comm into row and column comms
	MPI_Comm_split(MPI_COMM_WORLD, row, col, &row_Comm);
	MPI_Comm_split(MPI_COMM_WORLD, col, row, &col_Comm);

	printf("Processor %d's (row, col) are (%d, %d)\n", rank, row, col);

	MPI_Finalize();
	return 0;
}
