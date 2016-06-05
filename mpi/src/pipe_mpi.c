#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv)
{

	int rank, size;
    int row, col;
    int N; // N = sqrt(P)
    int token;
    MPI_Comm row_Comm, col_Comm;

	if(argc > 1)
		N = atoi(argv[1]);
	else
	{		
		printf("Usage: ./grid_mpi sqrt(P)");
		exit(1);
	}

	// Initial the MPI environment
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	// Determine Row and Col position
	row = rank / N;
	col = rank % N;

	token = size;

	// MPI_COMM_SPLIT(comm, color, key, newcomm)
	// Split comm into row and column comms
	MPI_Comm_split(MPI_COMM_WORLD, row, col, &row_Comm);
	MPI_Comm_split(MPI_COMM_WORLD, col, row, &col_Comm);

	// Obtain the rank in the row and col communicators
	int rank_row, rank_col;
	int size_row, size_col;
	MPI_Comm_rank(row_Comm, &rank_row);
	MPI_Comm_rank(col_Comm, &rank_col);
	MPI_Comm_size(row_Comm, &size_row);
	MPI_Comm_size(col_Comm, &size_col);

    printf("Processor %d's (row, col) are (%d, %d) and the rank in row/col comm is (%d, %d)\n", rank, row, col, rank_row, rank_col);

	// Pass the token
	if(rank_row != 0)
	{
		MPI_Recv(&token, 1, MPI_INT, rank_row - 1, 0, row_Comm, MPI_STATUS_IGNORE);
		printf("Processor %d received token %d from processor %d in Row Communicator!\n", rank_row, token, rank_row - 1);
	}

	MPI_Send(&token, 1, MPI_INT, (rank_row + 1) % size_row, 0, row_Comm);

	if(rank_col != 0)
	{
		MPI_Recv(&token, 1, MPI_INT, rank_col - 1, 0, col_Comm, MPI_STATUS_IGNORE);
		printf("Processor %d received token %d from processor %d in Col Communicator!\n", rank_col, token, rank_col - 1);
	}

	MPI_Send(&token, 1, MPI_INT, (rank_col + 1) % size_col, 0, col_Comm);

	MPI_Finalize();
	return 0;
}
