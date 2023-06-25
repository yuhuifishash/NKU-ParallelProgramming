#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cuda_runtime.h>
#include <chrono>
#include"device_launch_parameters.h"
using namespace std;
const int MAXN = 4005, MAXR = 105;
class Complex
{
public:
    float Re;
    float Im;
    __device__ Complex operator-() {
        return Complex{ -Re,-Im };
    }
    __device__ Complex operator+(const Complex& b) {
        return Complex{ Re + b.Re,Im + b.Im };
    }
    __device__ Complex operator-(const Complex& b) {
        return Complex{ Re - b.Re,Im - b.Im };
    }
    __device__ Complex operator*(const Complex& b) {
        return Complex{ Re * b.Re - Im * b.Im,Im * b.Re + Re * b.Im };
    }
};
Complex* a, * Core;//b展开为一维之后的FFT结果
int n = 0, m = 0, t = 0, * c;
int l = 0, * r;
int limit = 1;
int blocks = 31;
//device
Complex* a_d, * Core_d;
int* r_d, * c_d;
__device__ void swap(Complex& a, Complex& b)
{
    Complex tmp = a;
    a = b;
    b = tmp;
}
__device__ void FFT(Complex* A, int type, int limit, int* r)
{
    const float PI = acos(-1.0);
    for (int i = 0; i < limit; ++i) {
        if (i < r[i]) {
            swap(A[i], A[r[i]]);
        }
    }
    for (int mid = 1; mid < limit; mid <<= 1) {
        Complex Wn{ cos(PI / mid),type * sin(PI / mid) };
        for (int R = mid << 1, j = 0; j < limit; j += R) {
            Complex w{ 1,0 };
            for (int k = 0; k < mid; ++k, w = w * Wn) {
                Complex x = A[j + k], y = w * A[j + mid + k];
                A[j + k] = x + y;
                A[j + mid + k] = x - y;
            }
        }
    }
    if (type == 1) { return; }
    for (int i = 0; i < limit; ++i) {
        A[i].Re = A[i].Re / limit;
        A[i].Im = A[i].Im / limit;
    }
}
void init()
{
    //cin>>n>>m>>t;
    n = 32, m = 32, t = 3;
    //n = 15,m = 15,t = 3;
    while (limit <= blocks + t - 2) {
        limit <<= 1; ++l;
    }
    for (int i = 0; i < limit; ++i) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            //cin>>a[i][j].Re;
            a[i * m + j].Re = rand() % 10 + 1;
            //cout << a[i*m+j].Re << " ";
        }//cout << "\n";
    }
    for (int i = 0; i < t; ++i) {
        for (int j = 0; j < t; ++j) {
            //cin>>b[i][j].Re;
            Core[i * limit + j].Re = rand() % 10 + 1;
        }
    }

}
__device__ void trans(Complex* t, int limit)
{
    for (int i = threadIdx.x; i < limit; i += limit) {
        for (int j = 0; j < i; ++j) {
            swap(t[i * limit + j], t[j * limit + i]);
        }
    }
}
__global__ void FFT_Core(Complex* Core, int* r, int limit)
{
    FFT(Core + threadIdx.x * limit, 1, limit, r);
    __syncthreads();
    trans(Core, limit);
    __syncthreads();
    FFT(Core + threadIdx.x * limit, 1, limit, r);
    __syncthreads();
    trans(Core, limit);
    __syncthreads();
}
__global__ void Conv_2D(Complex* a, Complex* Core, int* r, int* c, int limit, int blocks, int n, int m, int t)
{
    int i = blockIdx.x * blocks, j = blockIdx.y * blocks;
    __shared__ Complex t1[4096];
    for (int k = threadIdx.x; k < blocks + t; k += limit) {
        for (int s = 0; s < blocks + t; ++s) {
            t1[k * limit + s] = a[(k + i) * m + s + j];
            if (k + i > n || s + j > m) { t1[k * limit + s].Re = t1[k * limit + s].Im = 0; }
        }
    }
    __syncthreads();
    FFT(t1 + threadIdx.x * limit, 1, limit, r);
    __syncthreads();
    trans(t1, limit);
    __syncthreads();
    FFT(t1 + threadIdx.x * limit, 1, limit, r);
    __syncthreads();
    trans(t1, limit);
    __syncthreads();
    for (int k = threadIdx.x; k < limit; k += limit) {
        for (int s = 0; s < limit; ++s) {
            t1[k * limit + s] = t1[k * limit + s] * Core[k * limit + s];
        }
    }
    __syncthreads();
    FFT(t1 + threadIdx.x * limit, -1, limit, r);
    __syncthreads();
    trans(t1, limit);
    __syncthreads();
    FFT(t1 + threadIdx.x * limit, -1, limit, r);
    __syncthreads();
    trans(t1, limit);
    __syncthreads();

    /*if (blockIdx.x == 1 && blockIdx.y == 0 && threadIdx.x == 0) {
        for (int i = 0; i < limit; ++i) {
            for (int j = 0; j < limit; ++j) {
                printf("%.2f %.2fj ", t1[i * limit + j].Re, t1[i * limit + j].Im);
            }printf("\n");
        }
    }*/
    for (int k = threadIdx.x; k < blocks; k += limit) {
        for (int s = 0; s < blocks; ++s) {
            if (k + i >= n - t + 1 || s + j >= m - t + 1) { continue; }
            c[(k + i) * (m - t + 1) + s + j] = (int)t1[(k + t - 1) * limit + (s + t - 1)].Re;
        }
    }
}

int main()
{
    a = (Complex*)malloc(sizeof(Complex) * MAXN * MAXN);
    Core = (Complex*)malloc(sizeof(Complex) * MAXR * MAXR * 10);
    r = (int*)malloc(sizeof(int) * MAXN);
    c = (int*)malloc(sizeof(int) * MAXN * MAXN);
    memset(a, 0, sizeof(Complex) * MAXN * MAXN);
    memset(Core, 0, sizeof(Complex) * MAXR * MAXR * 10);
    memset(r, 0, sizeof(int) * MAXN);
    init();


    cudaMalloc((Complex**)&a_d, sizeof(Complex) * n * m);
    cudaMalloc((Complex**)&Core_d, sizeof(Complex) * limit * limit * 2);
    cudaMalloc((int**)&r_d, sizeof(int) * limit);
    cudaMalloc((int**)&c_d, sizeof(int) * n * m);

    cudaMemcpy(a_d, a, sizeof(Complex) * n * m, cudaMemcpyHostToDevice);
    cudaMemcpy(Core_d, Core, sizeof(Complex) * limit * limit * 2, cudaMemcpyHostToDevice);
    cudaMemcpy(r_d, r, sizeof(int) * limit, cudaMemcpyHostToDevice);

    dim3 FFT_grid(n / blocks + (bool)(n % blocks), m / blocks + (bool)(m % blocks));
    dim3 FFT_block(limit);

    auto start = std::chrono::high_resolution_clock::now();

    FFT_Core << <1, limit >> > (Core_d, r_d, limit);
    Conv_2D << <FFT_grid, FFT_block >> > (a_d, Core_d, r_d, c_d, limit, blocks, n, m, t);
    cudaMemcpy(c, c_d, sizeof(int) * n * m, cudaMemcpyDeviceToHost);

    /*for (int i = 0; i < n-t+1; ++i) {
        for (int j = 0; j < m-t+1; ++j) {
            cout << c[i * (m-t+1) + j] << " ";
        }
        cout << "\n";
    }*/
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n", elapsed.count());

    cudaFree(a_d);
    cudaFree(Core_d);
    cudaFree(r_d);
    cudaFree(c_d);

    free(a);
    free(Core);
    free(r);
    free(c);

    return 0;
}