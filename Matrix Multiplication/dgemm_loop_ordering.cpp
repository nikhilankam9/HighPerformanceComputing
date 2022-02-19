#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
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
    // dgemm with loop order ijk 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
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

int dgemm_v(double *A, double *B, double *C, int n){
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
          double r = C[i*n+j];
          for (int k = 0; k < n; k++)
            r += A[i*n+k] * B[k*n+j];
          C[i*n+j] = r;
      }
    return 0;
}

int dgemm_register(double *A, double *B, double *C, int n){
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
          register double r = C[i*n+j];
          for (int k = 0; k < n; k++)
            r += A[i*n+k] * B[k*n+j];
          C[i*n+j] = r;
      }
    return 0;
}

int main(int argc, char **argv){
    double narray[5] = {64, 128, 256, 512, 1024};
    Timer timer;

    std::string function_names[8] = {"dgemm_ijk", "dgemm_ikj", "dgemm_jik", "dgemm_jki", "dgemm_kij", "dgemm_kji", "dgemm_register", "dgemm_register_blocked_2x2"};
    std::vector<int (*)(double *, double *, double *, int)> functions;
    functions.push_back(dgemm_ijk);
    functions.push_back(dgemm_ikj);
    functions.push_back(dgemm_jik);
    functions.push_back(dgemm_jki);
    functions.push_back(dgemm_kij);
    functions.push_back(dgemm_kji);
    // functions.push_back(dgemm_register);

    for (int itr = 0; itr < 5; itr += 1){
        double n = narray[itr];
        cout<<"n = "<<n<<endl;
        double * A = (double *) malloc(n*n*sizeof(double));
        double * B = (double *) malloc(n*n*sizeof(double));
        double * C0 = (double *) malloc(n*n*sizeof(double));
        double * C1 = (double *) malloc(n*n*sizeof(double));
        double dt;

        double no_ops = (n/10) * (n/10) * (n/10) * (3 + 1 + 1 + 1); //3 load, 1 store, 1 product, 1 sum

        RandMatrixGen(A, n);
        RandMatrixGen(B, n);

        memset(C0, 0, sizeof(double)*n*n);
        timer.start();
        dgemm(A, B, C0, n);
        timer.end();
        dt = timer.get();
        std::cout << "dgemm ijk"<<" GFLOPS: "<<(no_ops / dt)/1000000 <<", time: "<<dt<<endl;

        timer.start();
        memset(C1, 0, sizeof(double)*n*n);
        dgemm_v(A, B, C1, n);
        timer.end();
        dt = timer.get();
        double max_difference = verify(C0, C1, n);
        std::cout << "dgemm_v ijk"<<" time: "<<dt<<", error: "<<max_difference<<endl;

        for(int i=0; i<functions.size(); i++){
            memset(C1, 0, sizeof(double)*n*n);
            timer.start();
            functions[i](A, B, C1, n);
            timer.end();
            double max_difference = verify(C0, C1, n);
            dt = timer.get();
            cout<<function_names[i]<<" GFLOPS: "<<(no_ops / dt)/1000000 <<", time: "<<dt<<", error: "<<max_difference<<endl;
        }

        cout<<"-------------------------------"<<endl;

        free(A);
        free(B);
        free(C0);
        free(C1);
    }
    
    
    return 0;
}