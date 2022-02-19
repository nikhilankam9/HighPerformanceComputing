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
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j];
    return 0;
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

int dgemm_register_ijk_3x2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 3)
        for (j = 0; j < n; j += 2){
            register double c1 = C[i*n+j];
            register double c2 = C[i*n+j+1];
            register double c3 = C[(i+1)*n+j];
            register double c4 = C[(i+1)*n+j+1];
            register double c5 = C[(i+2)*n+j];
            register double c6 = C[(i+2)*n+j+1];

            for (k = 0; k < n; k += 2){
                register double a1 = A[i*n+k];
                register double a2 = A[i*n+k+1];
                register double a3 = A[(i+1)*n+k];
                register double a4 = A[(i+1)*n+k+1];
                register double a5 = A[(i+2)*n+k];
                register double a6 = A[(i+2)*n+k+1];
                
                register double b1 = B[k*n+j];
                register double b2 = B[k*n+j+1];
                register double b3 = B[(k+1)*n+j];
                register double b4 = B[(k+1)*n+j+1];
                // register double b5 = B[(k+2)*n+j];  //not used
                // register double b6 = B[(k+2)*n+j+1]; //not used

                c1 += a1*b1 + a2*b3;
                c2 += a1*b2 + a2*b4;
                c3 += a3*b1 + a4*b3;
                c4 += a3*b2 + a4*b4;
                c5 += a5*b1 + a6*b3;
                c6 += a5*b2 + a6*b4;
            }
            C[i*n+j] = c1;
            C[i*n+j+1] = c2;
            C[(i+1)*n+j] = c3;
            C[(i+1)*n+j+1] = c4;
            C[(i+2)*n+j] = c5;
            C[(i+2)*n+j+1] = c6;
    }
    return 0;
}

int dgemm_register_ijk_2x2_2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for(i = 0; i < n; i += 2)
        for(j = 0; j < n; j += 2) {
            register int t = i*n+j; register int tt = t+n;
            register double c00 = C[t]; register double c01 = C[t+1]; register double c10 = C[tt]; register double c11 = C[tt+1];
            for(k = 0; k < n; k += 2) {
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                register double a00 = A[ta]; register double a10 = A[tta]; register double b00 = B[tb]; register double b01 = B[tb+1];
                
                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
                a00 = A[ta+1]; a10 = A[tta+1]; b00 = B[ttb]; b01 = B[ttb+1];
                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
            }

            C[t] = c00;
            C[t+1] = c01;
            C[tt] = c10;
            C[tt+1] = c11;
        }
    return 0;
}

int dgemm_register_ikj_2x2(double *A, double *B, double *C, int N){
    register int i, j, k;
    register int n = N;
    for (i = 0; i < n; i += 2){
        register int i_n = i * n;
        register int i_n_n = i_n + n;
        for (k = 0; k < n; k += 2){
            register double a1 = A[i_n + k];
            register double a2 = A[i_n + k + 1];
            register double a3 = A[i_n_n + k];
            register double a4 = A[i_n_n + k + 1];
            register int k_n = k * n;
            register int k_n_n = k * n + n;
            for (j = 0; j < n; j += 2){
                register double b1 = B[k_n + j];
                register double b2 = B[k_n_n + j];
                register double b3 = B[k_n + j + 1];
                register double b4 = B[k_n_n + j + 1];
                C[i_n + j]     += a1*b1 + a2*b2;
                C[i_n + j + 1] += a1*b3 + a2*b4;
                C[i_n_n + j]   += a3*b1 + a4*b2;
                C[i_n_n + j + 1] += a3*b3 + a4*b4;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_2x3(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 2){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            // register double a3 = A[i * n + k + 2]; //not used
            register double a4 = A[(i + 1) * n + k];
            register double a5 = A[(i + 1) * n + k + 1];
            // register double a6 = A[(i + 1) * n + k + 2]; //not used
            for (j = 0; j < n; j += 3){
                register double b1 = B[k * n + j];
                register double b4 = B[k * n + j + 1];
                register double b7 = B[k * n + j + 2];
                register double b2 = B[(k + 1)*n + j];
                register double b5 = B[(k + 1)*n + j + 1];
                register double b8 = B[(k + 1)*n + j + 2];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b4 + a2*b5;
                C[i * n + j + 2] += a1*b7 + a2*b8;

                C[(i + 1) * n + j]   += a4*b1 + a5*b2;
                C[(i + 1) * n + j + 1] += a4*b4 + a5*b5;
                C[(i + 1) * n + j + 2] += a4*b7 + a5*b8;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_2x4(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 2){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            // register double a3 = A[i * n + k + 2]; //not used
            // register double a4 = A[i * n + k + 3]; //not used
            register double a5 = A[(i + 1) * n + k];
            register double a6 = A[(i + 1) * n + k + 1];
            // register double a7 = A[(i + 1) * n + k + 2]; //not used
            // register double a8 = A[(i + 1) * n + k + 3]; //not used

            for (j = 0; j < n; j += 4){
                register double b1 = B[k * n + j];
                register double b3 = B[k * n + j + 1];
                register double b5 = B[k * n + j + 2];
                register double b7 = B[k * n + j + 3];
                register double b2 = B[(k + 1)*n + j];
                register double b4 = B[(k + 1)*n + j + 1];
                register double b6 = B[(k + 1)*n + j + 2];
                register double b8 = B[(k + 1)*n + j + 3];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b3 + a2*b4;
                C[i * n + j + 2] += a1*b5 + a2*b6;
                C[i * n + j + 3] += a1*b7 + a2*b8;

                C[(i + 1) * n + j]     += a5*b1 + a6*b2;
                C[(i + 1) * n + j + 1] += a5*b3 + a6*b4;
                C[(i + 1) * n + j + 2] += a5*b5 + a6*b6;
                C[(i + 1) * n + j + 3] += a5*b7 + a6*b8;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_2x5(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 2){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            // register double a3 = A[i * n + k + 2]; //not used
            // register double a4 = A[i * n + k + 3]; //not used
            // register double a5 = A[i * n + k + 4]; //not used
            register double a6 = A[(i + 1) * n + k];
            register double a7 = A[(i + 1) * n + k + 1];
            // register double a8 = A[(i + 1) * n + k + 2]; //not used
            // register double a9 = A[(i + 1) * n + k + 3]; //not used
            // register double a10 = A[(i + 1) * n + k + 4]; //not used

            for (j = 0; j < n; j += 5){
                register double b1 = B[k * n + j];
                register double b3 = B[k * n + j + 1];
                register double b5 = B[k * n + j + 2];
                register double b7 = B[k * n + j + 3];
                register double b9 = B[k * n + j + 4];
                register double b2 = B[(k + 1)*n + j];
                register double b4 = B[(k + 1)*n + j + 1];
                register double b6 = B[(k + 1)*n + j + 2];
                register double b8 = B[(k + 1)*n + j + 3];
                register double b10 = B[(k + 1)*n + j + 4];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b3 + a2*b4;
                C[i * n + j + 2] += a1*b5 + a2*b6;
                C[i * n + j + 3] += a1*b7 + a2*b8;
                C[i * n + j + 4] += a1*b9 + a2*b10;

                C[(i + 1) * n + j]     += a6*b1 + a7*b2;
                C[(i + 1) * n + j + 1] += a6*b3 + a7*b4;
                C[(i + 1) * n + j + 2] += a6*b5 + a7*b6;
                C[(i + 1) * n + j + 3] += a6*b7 + a7*b8;
                C[(i + 1) * n + j + 4] += a6*b9 + a7*b10;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_3x5(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 2){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            // register double a3 = A[i * n + k + 2]; //not used
            // register double a4 = A[i * n + k + 3]; //not used
            // register double a5 = A[i * n + k + 4]; //not used
            register double a6 = A[(i + 1) * n + k];
            register double a7 = A[(i + 1) * n + k + 1];
            // register double a8 = A[(i + 1) * n + k + 2]; //not used
            // register double a9 = A[(i + 1) * n + k + 3]; //not used
            // register double a10 = A[(i + 1) * n + k + 4]; //not used

            for (j = 0; j < n; j += 5){
                register double b1 = B[k * n + j];
                register double b3 = B[k * n + j + 1];
                register double b5 = B[k * n + j + 2];
                register double b7 = B[k * n + j + 3];
                register double b9 = B[k * n + j + 4];
                register double b2 = B[(k + 1)*n + j];
                register double b4 = B[(k + 1)*n + j + 1];
                register double b6 = B[(k + 1)*n + j + 2];
                register double b8 = B[(k + 1)*n + j + 3];
                register double b10 = B[(k + 1)*n + j + 4];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b3 + a2*b4;
                C[i * n + j + 2] += a1*b5 + a2*b6;
                C[i * n + j + 3] += a1*b7 + a2*b8;
                C[i * n + j + 4] += a1*b9 + a2*b10;

                C[(i + 1) * n + j]     += a6*b1 + a7*b2;
                C[(i + 1) * n + j + 1] += a6*b3 + a7*b4;
                C[(i + 1) * n + j + 2] += a6*b5 + a7*b6;
                C[(i + 1) * n + j + 3] += a6*b7 + a7*b8;
                C[(i + 1) * n + j + 4] += a6*b9 + a7*b10;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_3x2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 3){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            register double a4 = A[(i + 1) * n + k];
            register double a5 = A[(i + 1) * n + k + 1];
            register double a7 = A[(i + 2) * n + k];
            register double a8 = A[(i + 2) * n + k + 1];
            for (j = 0; j < n; j += 2){
                register double b1 = B[k * n + j];
                register double b4 = B[k * n + j + 1];
                register double b2 = B[(k + 1)*n + j];
                register double b5 = B[(k + 1)*n + j + 1];
                // register double b3 = B[(k + 2)*n + j]; //not used
                // register double b6 = B[(k + 2)*n + j + 1]; //not used

                register double c1 = C[i * n + j];
                register double c2 = C[i * n + j + 1];
                register double c3 = C[(i + 1) * n + j];
                register double c4 = C[(i + 1) * n + j + 1];
                register double c5 = C[(i + 2) * n + j];
                register double c6 = C[(i + 2) * n + j + 1];

                c1 += a1*b1 + a2*b2;
                c2 += a1*b4 + a2*b5;

                c3 += a4*b1 + a5*b2;
                c4 += a4*b4 + a5*b5;

                c5 += a7*b1 + a8*b2;
                c6 += a7*b4 + a8*b5;

                C[i * n + j] = c1;
                C[i * n + j + 1] = c2;
                C[(i + 1) * n + j] = c3;
                C[(i + 1) * n + j + 1] = c4;
                C[(i + 2) * n + j] = c5;
                C[(i + 2) * n + j + 1] = c6;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_4x2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 4){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            register double a3 = A[(i + 1) * n + k];
            register double a4 = A[(i + 1) * n + k + 1];
            register double a5 = A[(i + 2) * n + k];
            register double a6 = A[(i + 2) * n + k + 1];
            register double a7 = A[(i + 3) * n + k];
            register double a8 = A[(i + 3) * n + k + 1];

            for (j = 0; j < n; j += 2){
                register double b1 = B[k * n + j];
                register double b5 = B[k * n + j + 1];
                register double b2 = B[(k + 1)*n + j];
                register double b6 = B[(k + 1)*n + j + 1];
                // register double b3 = B[(k + 2)*n + j];
                // register double b7 = B[(k + 2)*n + j + 1];
                // register double b4 = B[(k + 3)*n + j];
                // register double b8 = B[(k + 3)*n + j + 1];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b5 + a2*b6;

                C[(i + 1) * n + j]     += a3*b1 + a4*b2;
                C[(i + 1) * n + j + 1] += a3*b5 + a4*b6;

                C[(i + 2) * n + j]     += a5*b1 + a6*b2;
                C[(i + 2) * n + j + 1] += a5*b5 + a6*b6;

                C[(i + 3) * n + j]     += a7*b1 + a8*b2;
                C[(i + 3) * n + j + 1] += a7*b5 + a8*b6;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_6x2(double *A, double *B, double *C, int n){
    register int i, j, k;
    for (i = 0; i < n; i += 6){
        for (k = 0; k < n; k += 2){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            register double a3 = A[(i + 1) * n + k];
            register double a4 = A[(i + 1) * n + k + 1];
            register double a5 = A[(i + 2) * n + k];
            register double a6 = A[(i + 2) * n + k + 1];
            register double a7 = A[(i + 3) * n + k];
            register double a8 = A[(i + 3) * n + k + 1];
            register double a9 = A[(i + 4) * n + k];
            register double a10 = A[(i + 4) * n + k + 1];
            register double a11 = A[(i + 5) * n + k];
            register double a12 = A[(i + 5) * n + k + 1];

            for (j = 0; j < n; j += 2){
                register double b1 = B[k * n + j];
                register double b5 = B[k * n + j + 1];
                register double b2 = B[(k + 1)*n + j];
                register double b6 = B[(k + 1)*n + j + 1];
                // register double b3 = B[(k + 2)*n + j];
                // register double b7 = B[(k + 2)*n + j + 1];
                // register double b4 = B[(k + 3)*n + j];
                // register double b8 = B[(k + 3)*n + j + 1];

                C[i * n + j]     += a1*b1 + a2*b2;
                C[i * n + j + 1] += a1*b5 + a2*b6;

                C[(i + 1) * n + j]     += a3*b1 + a4*b2;
                C[(i + 1) * n + j + 1] += a3*b5 + a4*b6;

                C[(i + 2) * n + j]     += a5*b1 + a6*b2;
                C[(i + 2) * n + j + 1] += a5*b5 + a6*b6;

                C[(i + 3) * n + j]     += a7*b1 + a8*b2;
                C[(i + 3) * n + j + 1] += a7*b5 + a8*b6;

                C[(i + 4) * n + j]     += a9*b1 + a10*b2;
                C[(i + 4) * n + j + 1] += a9*b5 + a10*b6;

                C[(i + 5) * n + j]     += a11*b1 + a12*b2;
                C[(i + 5) * n + j + 1] += a11*b5 + a12*b6;
            }
        }
    }
    return 0;
}

int dgemm_register_ikj_3x3(double *A, double *B, double *C, int n){
    int i, j, k;
    for (i = 0; i < n; i += 3){
        for (k = 0; k < n; k += 3){
            register double a1 = A[i * n + k];
            register double a2 = A[i * n + k + 1];
            register double a3 = A[i * n + k + 2];
            register double a4 = A[(i + 1) * n + k];
            register double a5 = A[(i + 1) * n + k + 1];
            register double a6 = A[(i + 1) * n + k + 2];
            register double a7 = A[(i + 2) * n + k];
            register double a8 = A[(i + 2) * n + k + 1];
            register double a9 = A[(i + 2) * n + k + 2];
            for (j = 0; j < n; j += 3){
                register double b1 = B[k * n + j];
                register double b4 = B[k * n + j + 1];
                register double b7 = B[k * n + j + 2];
                register double b2 = B[(k + 1)*n + j];
                register double b5 = B[(k + 1)*n + j + 1];
                register double b8 = B[(k + 1)*n + j + 2];
                register double b3 = B[(k + 2)*n + j];
                register double b6 = B[(k + 2)*n + j + 1];
                register double b9 = B[(k + 2)*n + j + 2];

                C[i * n + j]     += a1*b1 + a2*b2 + a3*b3;
                C[i * n + j + 1] += a1*b4 + a2*b5 + a3*b6;
                C[i * n + j + 2] += a1*b7 + a2*b8 + a3*b9;

                C[(i + 1) * n + j]   += a4*b1 + a5*b2 + a6*b3;
                C[(i + 1) * n + j + 1] += a4*b4 + a5*b5 + a6*b6;
                C[(i + 1) * n + j + 2] += a4*b7 + a5*b8 + a6*b9;

                C[(i + 2) * n + j]   += a7*b1 + a8*b2 + a9*b3;
                C[(i + 2) * n + j + 1] += a7*b4 + a8*b5 + a9*b6;
                C[(i + 2) * n + j + 2] += a7*b7 + a8*b8 + a9*b9;
            }
        }
    }
    return 0;
}

int main(int argc, char **argv){
    int narray[3] = {1200, 2400, 4800};
    double max_difference = 0, dgemm_time, dgemm_v_time, tt;
    Timer timer;

    for (int itr = 0; itr < 3; itr += 1){
        int n = narray[itr];
        cout<<"n = "<<n<<endl;
        double * A = (double *) malloc(n*n*sizeof(double));
        double * B = (double *) malloc(n*n*sizeof(double));
        double * C0 = (double *) malloc(n*n*sizeof(double));
        double * C1 = (double *) malloc(n*n*sizeof(double));
        double * C2 = (double *) malloc(n*n*sizeof(double));

        RandMatrixGen(A, n);
        RandMatrixGen(B, n);

        memset(C0, 0, sizeof(double)*n*n);
        timer.start();
        dgemm(A, B, C0, n);
        timer.end();
        dgemm_time = timer.get();
        std::cout << "dgemm: " << dgemm_time << "\n";

        timer.start();
        memset(C1, 0, sizeof(double)*n*n);
        dgemm_v(A, B, C1, n);
        timer.end();
        dgemm_v_time = timer.get();
        std::cout << "dgemm_v: " << dgemm_v_time << "\n";

        std::string function_names[14] = {"dgemm_register_ijk", "dgemm_register_ikj", "dgemm_register_ijk_2x2", 
                "dgemm_register_ijk_2x2_2", "dgemm_register_ikj_2x2", "dgemm_register_ikj_2x3", "dgemm_register_ikj_2x4",
                 "dgemm_register_ikj_2x5", "dgemm_register_ikj_3x2", "dgemm_register_ijk_3x2", "dgemm_register_ikj_4x2", 
                 "dgemm_register_ikj_6x2", "dgemm_register_ikj_3x3", "dgemm_register_ikj_3x5"};
        std::vector<int (*)(double *, double *, double *, int)> functions;
        functions.push_back(dgemm_register_ijk);
        functions.push_back(dgemm_register_ikj);
        functions.push_back(dgemm_register_ijk_2x2);
        functions.push_back(dgemm_register_ijk_2x2_2);
        functions.push_back(dgemm_register_ikj_2x2);
        functions.push_back(dgemm_register_ikj_2x3);
        functions.push_back(dgemm_register_ikj_2x4);
        functions.push_back(dgemm_register_ikj_2x5);
        functions.push_back(dgemm_register_ikj_3x2);
        functions.push_back(dgemm_register_ijk_3x2);
        functions.push_back(dgemm_register_ikj_4x2);
        functions.push_back(dgemm_register_ikj_6x2);
        functions.push_back(dgemm_register_ikj_3x3);
        functions.push_back(dgemm_register_ikj_3x5);

        for(int i=0; i<functions.size(); i++){
            memset(C1, 0, sizeof(double)*n*n);
            timer.start();
            functions[i](A, B, C1, n);
            timer.end();
            max_difference = verify(C0, C1, n);
            tt = timer.get();
            cout <<function_names[i]<<"\n\ttime: "<< tt << "\n\tdgemm_time/time: "<<dgemm_time / tt
                            <<"\n\tdgemm_v_time/time: "<<dgemm_v_time / tt<<"\n\terror:"<<max_difference<<endl<<endl;
        }
        cout<<"-------------------------------"<<endl;

        free(A);
        free(B);
        free(C0);
        free(C1);
        free(C2);
    }
    
    
    return 0;
}