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

int dgemm_v(double *A, double *B, double *C, int n){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++){
            double r = 0;
            for (int k = 0; k < n; k++)
                r += A[i*n+k] * B[k*n+j];
            C[i*n+j] = r;
        }
    return 0;
}

int dgemm_register_ijk(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i++){
        register int i_n = i * n;
        for (j = 0; j < n; j++){
            register double r = C[i_n + j];
            for (k = 0; k < n; k++)
                r += A[i_n + k] * B[k*n + j];
            C[i_n + j] = r;
        }
    }
    return 0;
}

int dgemm_register_ikj(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i++){
        register int i_n = i * n;
        for (k = 0; k < n; k++){
            register double r = A[i_n + k];
            register int k_n = k * n;
            for (j = 0; j < n; j++){
                C[i_n + j] += r * B[k_n + j];
            }
        }
    }
    return 0;
}

int dgemm_register_ijk_2x2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 2)
        for (j = 0; j < n; j += 2){
            register double c1 = C[i*n+j];
            register double c2 = C[i*n+j+1];
            register double c3 = C[(i+1)*n+j];
            register double c4 = C[(i+1)*n+j+1];
            for (k = 0; k < n; k += 2){
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

int dgemm_block_ikj_2(double *A, double *B, double *C, int n, int block){
    for (int i = 0; i < n; i += block){
        for (int k = 0; k < n; k += block){
            for (int j = 0; j < n; j += block){
                for (int i1 = i; i1 < i + block; i1 += 1){
                    for (int k1 = k; k1 < k + block; k1 += 1){
                        register double r = A[i1*n + k1];
                        for (int j1 = j; j1 < j + block; j1 += 1){
                            C[i1*n + j1] += r * B[k1*n + j1];
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int dgemm_block_ikj_2x2(double *A, double *B, double *C, int n, int block){
    for (int i = 0; i < n; i += block){
        for (int k = 0; k < n; k += block){
            for (int j = 0; j < n; j += block){
                for (int i1 = i; i1 < i + block; i1 += 2){
                    for (int k1 = k; k1 < k + block; k1 += 2){
                        register double a1 = A[i1*n + k1];
                        register double a2 = A[i1*n + k1 + 1];
                        register double a3 = A[(i1 + 1) * n + k1];
                        register double a4 = A[(i1 + 1) * n + k1 + 1];
                        for (int j1 = j; j1 < j + block; j1 += 2){
                            register double b1 = B[k1*n + j1];
                            register double b3 = B[k1*n + j1 + 1];
                            register double b2 = B[(k1 + 1) * n + j1];
                            register double b4 = B[(k1 + 1) * n + j1 + 1];
                            
                            C[i1*n + j1]             += a1*b1 + a2*b2;
                            C[i1*n + j1 + 1]         += a1*b3 + a2*b4;
                            C[(i1 + 1) * n + j1]     += a3*b1 + a4*b2;
                            C[(i1 + 1) * n + j1 + 1] += a3*b3 + a4*b4;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int dgemm_block_ikj_3x3(double *A, double *B, double *C, int n, int block){
    for (int i = 0; i < n; i += block){
        for (int k = 0; k < n; k += block){
            for (int j = 0; j < n; j += block){
                for (int i1 = i; i1 < i + block; i1 += 2){
                    for (int k1 = k; k1 < k + block; k1 += 2){
                        register double a1 = A[i1*n + k1];
                        register double a2 = A[i1*n + k1 + 1];
                        register double a3 = A[(i1 + 1) * n + k1];
                        register double a4 = A[(i1 + 1) * n + k1 + 1];
                        for (int j1 = j; j1 < j + block; j1 += 2){
                            register double b1 = B[k1*n + j1];
                            register double b3 = B[k1*n + j1 + 1];
                            register double b2 = B[(k1 + 1) * n + j1];
                            register double b4 = B[(k1 + 1) * n + j1 + 1];
                            
                            C[i1*n + j1]             += a1*b1 + a2*b2;
                            C[i1*n + j1 + 1]         += a1*b3 + a2*b4;
                            C[(i1 + 1) * n + j1]     += a3*b1 + a4*b2;
                            C[(i1 + 1) * n + j1 + 1] += a3*b3 + a4*b4;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int main(int argc, char **argv){
    int narray[3] = {1024, 2400, 4800};
    int blocks[6] = {2, 4, 8, 16, 32, 64};
    double max_difference = 0, dgemm_time, dgemm_v_time, tt;
    Timer timer;

    for (int itr = 0; itr < 1; itr += 1){
        int n = narray[itr];
        double * A = (double *) malloc(n*n*sizeof(double));
        double * B = (double *) malloc(n*n*sizeof(double));
        double * C0 = (double *) malloc(n*n*sizeof(double));
        double * C1 = (double *) malloc(n*n*sizeof(double));
        double * C2 = (double *) malloc(n*n*sizeof(double));
        
        RandMatrixGen(A, n);
        RandMatrixGen(B, n);

        timer.start();
        memset(C0, 0, sizeof(double)*n*n);
        dgemm_v(A, B, C0, n);
        timer.end();
        dgemm_v_time = timer.get();
        std::cout << "dgemm_v: " << dgemm_v_time << "\n";

        std::string function_names[14] = {"dgemm_block", "dgemm_block_ikj_2x2"};
        std::vector<int (*)(double *, double *, double *, int, int)> functions;
        functions.push_back(dgemm_block);
        functions.push_back(dgemm_block_ikj_2x2);

        for (int bitr = 0; bitr < 6; bitr+=1){
            cout<<"n = "<<n<<" b = "<<blocks[bitr]<<endl;
            for(int i=0; i<functions.size(); i++){
                memset(C1, 0, sizeof(double)*n*n);
                timer.start();
                functions[i](A, B, C1, n, blocks[bitr]);
                timer.end();
                max_difference = verify(C0, C1, n);
                tt = timer.get();
                cout <<function_names[i]<<"\ttime: "<< tt << "\tdgemm_v_time/time: "<<dgemm_v_time/tt<<"\terror:"<<max_difference<<endl;
            }
            cout<<"-------------------------------"<<endl;
        }

        free(A);
        free(B);
        free(C0);
        free(C1);
        free(C2);
    }
    return 0;
}