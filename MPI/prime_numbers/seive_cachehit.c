#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "decomposition.h"
#include <stdio.h>
#include <string.h>

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
    
    ull rootn = (sqrt(n));
    int primesSize = (rootn%2)? rootn/2 +1: rootn/2;
    primes = (char *)malloc(primesSize);
    memset(primes, 0, primesSize);
    primes[0] = 1; //marking 1 as not time (2i+1 i.e 2*0+1)
    for (int p = 3; p <= (int) rootn; p +=2){
        if (primes[(p-1)/2] == 0){
            for (int itr = p*p; itr <= (int) rootn; itr += p){
                if (itr %2 != 0){
                    primes[(itr-1)/2] = 1;
                }
            }
        }
    }

    // printf("primes: ");
    // for (int p = 0; p < primesSize; p++){
    //     if (primes[p] == 0){
    //         printf("%d ", 2*p+1);
    //     }
    // }
    // printf("\n");
    // fflush(stdout);

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

    int b = 1000000;
    for (int k = 0; k<size; k+=b){
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
                // printf("prime = %lld %d  --- ",prime, start);
                for (int i = start; i < k+b; i+= prime){
                    // printf("%d,", 2*i+1);
                    marked[i] = 1;
                }
                // printf("\n");
            }

            if (!rank){
                int index = (prime - 1)/2;
                marked[index] = 0;
            }
        }
    }

    // ull itr = 0;
    // int block = 1000;
    // for (int i = 0; i< size; i+=block){
    //     for (int j = i; j< i+block; j++){
    //         ull val = 2*(j + low) + 1;
    //         // printf("%lld - ", val);
    //         for (int p = 1; p < primesSize; p++){
    //             if (primes[p] == 0){
    //                 ull pval = 2*p + 1;
    //                 if (j - pval >= i){
    //                     if (marked[j-pval] == 1){
    //                         marked[j] = 1;
    //                     }
    //                 }else if (val % pval == 0){
    //                     // printf(" %lld %lld ", pval, val);
    //                     if (val != pval){
    //                         marked[j] = 1;
    //                     }
    //                     break;
    //                 }
    //             }
    //         }
    //         // printf("\n");
    //     }
    // }


    // for (ull i = 0; i< size; i+=1){
    //     ull val = 2*(i + low) + 1;
    //     printf("%lld - ", val);
    //     for (int p = 1; p < primesSize; p++){
    //         if (primes[p] == 0){
    //             ull pval = 2*p + 1;
    //             printf("(prime = %lld)", pval);
    //             if (val != pval){
    //                 if (i > pval){
    //                     printf(" %lld -- %d ", i-pval , marked[i-pval]);
    //                     if (marked[i-pval] == 1  || (2*(i-pval + low) +1 == pval)){
    //                         marked[i] = 1;
    //                         break;
    //                     }
    //                 }
    //                 else if (val % pval == 0){
    //                     printf(" %lld, %lld ", pval, val);
    //                     marked[i] = 1;
    //                     break;
    //                 }
    //             }else{
    //                 break;
    //             }
                
    //         }
    //     }
    //     printf("  --------------- %d \n", marked[i]);
    // }



    // while (itr < size){
    //     ull val = 2*(itr + low) + 1;
    //     // printf("%lld - ", val);
    //     for (int p = 0; p < primesSize; p++){
    //         if (primes[p] == 0){
    //             ull pval = 2*p + 1;
    //             if (val % pval == 0){
    //                 // printf(" %lld ", pval);
    //                 if (val != pval){
    //                     marked[itr] = 1;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    //     // printf("\n");
    //     itr += 1;
    // }
    // fflush(stdout);


    // printf("res = ");
    for (int i = 0; i <size; i++){
        if (marked[i] == 0){
            // printf("%d, ", 2*i+1);
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