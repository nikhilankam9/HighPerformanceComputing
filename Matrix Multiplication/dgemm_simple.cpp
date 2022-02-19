#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>

class Timer{
public:
    void start(){
        err = clock_gettime(CLOCK_REALTIME, &start_time);
    }
    void end(){
        err = clock_gettime(CLOCK_REALTIME, &end_time);
    }
    double get(){
        return (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/(double)1000000000;
    }
private:
    int err = 0;
    struct timespec start_time, end_time;
};

int RandMatrixGen(double *M, int n){
    srand (time(NULL));
    for (int i = 0; i < n*n; i++)
      M[i] =  2.0*rand()/RAND_MAX - 1.0;
    return 0;
}

double verify(double *C0, double *C1, int n){
    double diff = fabs(C0[0] - C1[0]);
    for (int i = 0; i < n*n; i++)
      if (fabs(C0[i] - C1[i]) > diff) diff = fabs(C0[i] - C1[i]);
    return diff;
}

int dgemm(double *A, double *B, double *C, int n){
    // dgemm with loop order ijk 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_v(double *A, double *B, double *C, int n){
    // dgemm with variable reuse
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
          double r = C[i*n+j];
          for (int k = 0; k < n; k++)
            r += A[i*n+k] * B[k*n+j];
          C[i*n+j] = r;
      }
    return 0;
}

int main(int argc, char **argv){

    int n = atoi(argv[1]);
    double * A = (double *) malloc(n*n*sizeof(double));
    double * B = (double *) malloc(n*n*sizeof(double));
    double * C0 = (double *) malloc(n*n*sizeof(double));
    double * C1 = (double *) malloc(n*n*sizeof(double));

    Timer timer;
    // generate input matrix A and B
    RandMatrixGen(A, n);
    RandMatrixGen(B, n);

    // compute result using dgemm
    memset(C0, 0, sizeof(double)*n*n);
    timer.start();
    dgemm(A, B, C0, n);
    timer.end();
    std::cout << "dgemm execution time: " << timer.get() << "\n";

    // compute result using dgemm_v
    timer.start();
    memset(C1, 0, sizeof(double)*n*n);
    dgemm_v(A, B, C1, n);
    timer.end();
    std::cout << "dgemm_v execution time: " << timer.get() << "\n";

    double max_difference = verify(C0, C1, n);
    std::cout << "max difference = " << max_difference << std::endl;

    free(A);
    free(B);
    free(C0);
    free(C1);
    return 0;

}

