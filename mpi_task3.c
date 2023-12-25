#include <mpi.h>
#include <stdio.h>
#include <math.h> 
#include <stdlib.h> 
int main(int argc, char* argv) {
	double a[100000], b[100000], a1[100000], b1[100000];//, ProcSum = 0.0; 
	int ProcRank, ProcNum, N=100000; 
	MPI_Status Status; 
	MPI_Init(NULL,NULL); 
	MPI_Comm_size(MPI_COMM_WORLD,&ProcNum); 
	MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank); 
	if ( ProcRank == 0 )
		for (int i = 0; i < N; i++) {
			a[i] = i;
			b[i] = i-1;
		}
	double t1 = MPI_Wtime();
	//MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int j = 0; j<100; j++){
	if ( ProcRank == 0 ) {  
 		MPI_Send(&a, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&a1, N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status); 
 	} 
 	else if (ProcRank == 1){
		MPI_Recv(&b1, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
 		MPI_Send(&b, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	} 
	}
	double t2 = MPI_Wtime(); 
	if (ProcRank == 0) printf("Time for %d: %10.8f \n", N, t2-t1);
 	MPI_Finalize();
    return 0;
}

