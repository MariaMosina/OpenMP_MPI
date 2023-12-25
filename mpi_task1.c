#include <mpi.h>
#include <stdio.h>
#include <math.h> 
#include <stdlib.h> 
int main(int argc, char* argv) {
	double x[100000], Min;//, ProcSum = 0.0; 
	int ProcRank, ProcNum, N=100000; 
	MPI_Status Status; 
	MPI_Init(NULL,NULL); 
	MPI_Comm_size(MPI_COMM_WORLD,&ProcNum); 
	MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank); 
	if ( ProcRank == 0 )
		for (int i = 0; i < N; i++) x[i] = N - i; 
	double t1 = MPI_Wtime();
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	int k = N / ProcNum; 
	int i1 = k * ProcRank; 
 	int i2 = k * ( ProcRank + 1 );
	float locMin = x[i1]; 
 	if ( ProcRank == ProcNum-1 ) i2 = N; 
 	for ( int i = i1; i < i2; i++ ) 
		if (x[i] < locMin) locMin = x[i];
	if ( ProcRank == 0 ) { 
 		Min = locMin; 
 		for ( int i=1; i < ProcNum; i++ ) { 
 			MPI_Recv(&locMin, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status); 
 			if (locMin < Min) Min = locMin; 
 		} 
 	} 
 	else 
 	MPI_Send(&locMin, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
    	if ( ProcRank == 0 ) 
 	printf("\nMin = %10.2f \n",Min);
	double t2 = MPI_Wtime(); 
	if (ProcRank == 0) printf("Time: %10.8f \n", t2-t1);
 	MPI_Finalize();
    return 0;
}

