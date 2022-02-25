#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]){
    int id;
    int p;
    int sum[2]= {0, 0};
    int total_sum[2];

    MPI_Init(&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    double elapsedtime = - MPI_Wtime();

    //cyclic allocation
    for (int i = id; i < 100; i+=p){
        sum[0] += i;
    }

    //block allocation
    for (int i = (100/p)*id; i <(100/p)*(id+1); i++){
        sum[1] += i;
    }

    MPI_Reduce(&sum, &total_sum, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsedtime += MPI_Wtime();

    printf("\nLocal sum of %d = %d, %d %f\n", id, sum[0], sum[1], elapsedtime);
    fflush(stdout);

    MPI_Finalize();

    if (id == 0){
        printf("Total %d and %d\n", total_sum[0], total_sum[1]);
    }
    return 0;
}