#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>


int main(int argc, char **argv)
{
	int rankCount, ID;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);




	MPI_Finalize();

	return 0;
}
