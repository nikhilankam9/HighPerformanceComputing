#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;

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
    double d = 0;
    for (int i = 1; i < n*n; i++) {
        d = fabs(C0[i] - C1[i]);        
        if (d > diff) diff = d;
    }
    return diff;
}

int dgemm(double *A, double *B, double *C, int n){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
}

int dgemm_block(double *A, double *B, double *C, int n, int block){
    for (int i = 0; i < n; i += block){
        for (int j = 0; j < n; j += block){
            for (int k = 0; k < n; k += block){
                for (int i1 = i; i1 < i + block; i1 += 1){
                    for (int j1 = j; j1 < j + block; j1 += 1){
                        for (int k1 = k; k1 < k + block; k1 += 1){
                            C[i1*n + j1] += A[i1*n + k1] * B[k1*n + j1];
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int main(int argc, char **argv){
    int n = 1024;
    double * A = (double *) malloc(n*n*sizeof(double));
    double * B = (double *) malloc(n*n*sizeof(double));
    double * C0 = (double *) malloc(n*n*sizeof(double));
    double * C1 = (double *) malloc(n*n*sizeof(double));

    Timer timer;
    RandMatrixGen(A, n);
    RandMatrixGen(B, n);

    memset(C0, 0, sizeof(double)*n*n);
    timer.start();
    dgemm(A, B, C0, n);
    timer.end();
    std::cout << "dgemm execution time: " << timer.get() << "\n";

    int blocks[6] = {2, 4, 8, 16, 32, 64};
    for (int itr = 0; itr < 6; itr += 1){
        memset(C1, 0, sizeof(double)*n*n);
        timer.start();
        dgemm_block(A, B, C1, n, blocks[itr]);
        timer.end();
        cout<<"dgemm block("<<blocks[itr]<<") time:"<<timer.get()<<endl;

        double max_difference = verify(C0, C1, n);
        std::cout << "max difference = " << max_difference << std::endl;
    }

    free(A);
    free(B);
    free(C0);
    free(C1);
    return 0;
}
