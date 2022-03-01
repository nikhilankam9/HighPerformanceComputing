#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "decomposition.h"
#include <stdio.h>
#include <string.h>

ull minVal(ull a, ull b){
    if (a < b){
        return a;
    }
    return b;
}

int main(int argc, char *argv[]){
    int processors;
    int rank;
    double elapsed_time;
    ull n, low, high, size, prime, n1, lowVal;
    char *marked, *primes;
    int global_sum, local = 0;

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);

    n = atoll(argv[1]);
    n1 = (n%2 == 0) ? n/2 : (n/2 + 1);

    low = BLOCK_LOW(rank, processors, n1);
    high = BLOCK_HIGH(rank, processors, n1);
    size = BLOCK_SIZE(rank, processors, n1);

    // printf("processor: %d (%llu, %llu, %llu)\n", rank, low, high, size);
    // fflush(stdout);

    if ((n1)/processors < sqrt(n1)){
        if (rank == 0){
            printf("Too many processes\n");
        }
        MPI_Finalize();
        exit(1);
    }
    
    ull rootn = (ull) (sqrt(n));
    int primesSize = (rootn%2)? rootn/2 +1: rootn/2;
    primes = (char *)malloc(primesSize);
    memset(primes, 0, primesSize);
    primes[0] = 1; //marking 1 as not time (2i+1 i.e 2*0+1)
    for (ull p = 3; p <= rootn; p +=2){
        if (primes[(p-1)/2] == 0){
            for (ull itr = p*p; itr <= rootn; itr += p){
                if (itr %2 != 0){
                    primes[(itr-1)/2] = 1;
                }
            }
        }
    }

    marked = (char *) malloc(size);
    if (marked == NULL){
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    memset(marked, 0, size);
    if (!rank){
        marked[0] = 1;
        local = 1; //account 2 the only even prime number
    }

    int block = 1000000;
    for (int k = 0; k<size; k+=block){
        for (int itr = 1; itr < primesSize; itr++){
            if (primes[itr] == 0){
                prime = (ull) (2*itr +1);
                int start = k;
                lowVal = 2 * (low+k) + 1;
                ull rem = lowVal % prime;
                if (rem != 0){
                    ull startVal = lowVal + prime - rem;
                    if (!(startVal&1)){ //if startVal is even find the next time prime (e + o = o)
                        startVal += prime;
                    }
                    start =(startVal - 1)/2 - low ;
                }
                for (int i = start; i < minVal(k+block, size); i+= prime){
                    marked[i] = 1;
                }

                if (!rank){
                    int index = (prime - 1)/2;
                    marked[index] = 0;
                }
            }
        }
    }

    for (int i = 0; i <size; i++){
        if (marked[i] == 0){
            local += 1;
        }
    }
    MPI_Reduce(&local, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    MPI_Finalize();

    if (!rank){
        printf("\nTotal primes %d\n", global_sum);
        printf("Total time: %10.6f\n", elapsed_time);
    }

    return 0;
}
