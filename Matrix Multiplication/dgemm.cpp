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

int dgemm_ijk(double *A, double *B, double *C, int n){
    // dgemm with loop order ijk 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_ikj(double *A, double *B, double *C, int n){
    // implement your dgemm with loop order ikj 
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_jik(double *A, double *B, double *C, int n){
    // implement your dgemm with loop order jik
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_jki(double *A, double *B, double *C, int n){
    // implement your dgemm with loop order jki
    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
            for (int i = 0; i < n; i++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_kij(double *A, double *B, double *C, int n){
    // implement your dgemm with loop order kij
    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_kji(double *A, double *B, double *C, int n){
    // implement your dgemm with loop order kji
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_register(double *A, double *B, double *C, int n){
    // BONUS example: dgemm with register reuse
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
          register double r = C[i*n+j];
          for (int k = 0; k < n; k++)
            r += A[i*n+k] * B[k*n+j];
          C[i*n+j] = r;
      }
    return 0;
}

int dgemm_register_blocked_2x2(double *A, double *B, double *C, int n){
    // BONUS: implement your dgemm with blocked cache/register reuse (block size = 2)
    for (int i = 0; i < n; i += 2)
        for (int j = 0; j < n; j += 2){
            register double c1 = C[i*n+j];
            register double c2 = C[i*n+j+1];
            register double c3 = C[(i+1)*n+j];
            register double c4 = C[(i+1)*n+j+1];
            for (int k = 0; k < n; k += 2){
                register double a1 = A[i*n+k];
                register double a2 = A[i*n+k+1];
                register double a3 = A[(i+1)*n+k];
                register double a4 = A[(i+1)*n+k+1];
                register double b1 = B[k*n+j];
                register double b2 = B[(k+1)*n+j];
                register double b3 = B[k*n+j+1];
                register double b4 = B[(k+1)*n+j+1];
                c1 += a1*b1 + a2*b2;
                c2 += a1*b3 + a2*b4;
                c3 += a3*b1 + a4*b2;
                c4 += a3*b3 + a4*b4;
            }
            C[i*n+j] = c1;
            C[i*n+j+1] = c2;
            C[(i+1)*n+j] = c3;
            C[(i+1)*n+j+1] = c4;
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
    // name of functions
    std::string function_names[8] = {"dgemm_ijk", "dgemm_ikj", "dgemm_jik", "dgemm_jki", "dgemm_kij", "dgemm_kji", "dgemm_register", "dgemm_register_blocked_2x2"};
    // function points
    std::vector<int (*)(double *, double *, double *, int)> functions;
    functions.push_back(dgemm_ijk);
    functions.push_back(dgemm_ikj);
    functions.push_back(dgemm_jik);
    functions.push_back(dgemm_jki);
    functions.push_back(dgemm_kij);
    functions.push_back(dgemm_kji);
    functions.push_back(dgemm_register);
    functions.push_back(dgemm_register_blocked_2x2);
    // generate input matrix A and B
    RandMatrixGen(A, n);
    RandMatrixGen(B, n);

    // compute result using dgemm_ijk
    memset(C0, 0, sizeof(double)*n*n);
    dgemm_ijk(A, B, C0, n);

    for(int i=0; i<functions.size(); i++){
        memset(C1, 0, sizeof(double)*n*n);
        timer.start();
        functions[i](A, B, C1, n);
        timer.end();
        double max_difference = verify(C0, C1, n);
        std::cout << function_names[i] << " execution time: " << timer.get() << ", ";
        std::cout << "max difference = " << max_difference << std::endl;
    }

    free(A);
    free(B);
    free(C0);
    free(C1);
    return 0;

}

