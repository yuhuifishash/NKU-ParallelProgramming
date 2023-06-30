#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <immintrin.h>
#include <cuda_runtime.h>
#include"device_launch_parameters.h"
using namespace std;
const int MAXN = 2100, MAXR = 105;
float* A, *I, *R;
float* ans;
int n = 0, t = 0;
int cnt = -1;
float* I_d, *R_d,*ans_d;
void init()
{
    n = 512, t = 31;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = rand() % 10 + 1;
            //cout<<A[i*n+j]<<" ";
        }//cout<<"\n";
    }
    for (int i = 0; i < t; ++i) {
        for (int j = 0; j < t; ++j) {
            R[i * t + j] = rand() % 10 + 1;
        }
    }
}
void img2col()
{
    int col = t * t;
    int row = (n - t + 1) * (n - t + 1);

    for (int i = 0; i <= n - t; ++i) {
        for (int j = 0; j <= n - t; ++j) {
            ++cnt;
            for (int s1 = 0; s1 < t; ++s1) {
                for (int s2 = 0; s2 < t; ++s2) {
                    I[(s1 * t + s2) * row + cnt] = A[(i + s1) * n + j + s2];
                    //cout << I[cnt * col + (s1 * t + s2)] << " ";
                }
            }//cout << "\n";
        }
    }
    /*for (int i = 0; i < col; ++i) {
        for (int j = 0; j < row; ++j) {
            cout << I[i * row + j] << " ";
        }
        cout << "\n";
    }*/
}
//void GEMM()
//{
//    for (int i = 0; i <= cnt; ++i) {
//        for (int j = 0; j < t * t; ++j) {
//            ans[i] += I[i * (t * t) + j] * R[j];
//        }
//    }
//}
#define Block 512
#define tileX 16
__global__ void GEMM(float* A, float* B, float* C, int N, int M)
{
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    tx = tx * tileX;
    float B_t;
    for (int k = 0; k < M; ++k) {
        B_t = B[k];
        for (int i = 0; i < tileX; ++i) {
            C[tx + i] += B_t * A[ k * N + (tx + i)];
        }
    }

}
int main()
{

	A = (float*)malloc(sizeof(float) * MAXN * MAXN);
	R = (float*)malloc(sizeof(float) * MAXR * MAXR);
	I = (float*)malloc(sizeof(float) * MAXN * MAXN * 100);
	ans = (float*)malloc(sizeof(float) * MAXN * MAXN);

	cudaMalloc((float**)&I_d, sizeof(float) * MAXN * MAXN * 100);
	cudaMalloc((float**)&R_d, sizeof(float) * MAXR * MAXR);
	cudaMalloc((float**)&ans_d, sizeof(float) * MAXN * MAXN);
	init();

	auto start = std::chrono::high_resolution_clock::now();
	img2col();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
	printf("%f\n", elapsed.count());

	int col = (n - t + 1) * (n - t + 1);
	cudaMemcpy(I_d, I, sizeof(float) * col * t * t, cudaMemcpyHostToDevice);
	cudaMemcpy(R_d, R, sizeof(float) * t * t, cudaMemcpyHostToDevice);


	dim3 GEMM_grid(col / Block + (bool)(col % Block));
	dim3 GEMM_block(Block / tileX);

	start = std::chrono::high_resolution_clock::now();

	GEMM << <GEMM_grid, GEMM_block >> > (I_d, R_d, ans_d, col, t * t);

	cudaMemcpy(ans, ans_d, sizeof(float) * (n - t + 1) * (n - t + 1), cudaMemcpyDeviceToHost);

	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	printf("%f\n", elapsed.count());

	/*for(int i=0;i<n-t+1;++i){
		for(int j=0;j<n-t+1;++j){
			cout<<ans[i*(n-t+1)+j]<<" ";
		}
		cout<<"\n";
	}*/

    cudaFree(I_d);
    cudaFree(R_d);
    cudaFree(ans_d);

    free(A);
    free(R);
    free(I);
    free(ans);

	return 0;
}