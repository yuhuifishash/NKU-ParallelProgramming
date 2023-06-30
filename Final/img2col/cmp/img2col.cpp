#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <immintrin.h>
using namespace std;
const int MAXN = 2100, MAXR = 105;
float A[MAXN * MAXN], I[MAXN * MAXN * 100] = {},R[MAXR * MAXR] = {};
float ans[MAXN*MAXN]={};
int n = 0, t = 0;
int cnt = -1;
void init()
{
	n = 1024, t = 3;
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
    int col = t*t;
	for (int i = 0; i <= n - t; ++i) {
		for (int j = 0; j <= n - t; ++j) {
            ++cnt;
			for(int s1 = 0;s1 < t; ++s1){
                for(int s2 = 0; s2 < t; ++s2){
                    I[cnt*col + (s1*t+s2)] = A[(i+s1)*n+j+s2];
                }
            }
		}
	}
}
void GEMM()
{
    for(int i=0;i<=cnt;++i){
        for(int j=0;j<t*t;++j){
            ans[i]+=I[i*(t*t)+j]*R[j];
        }
    }
}
int main()
{
	init();

    auto start = std::chrono::high_resolution_clock::now();
    img2col();
    auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
    start = std::chrono::high_resolution_clock::now();
    GEMM();
    // for(int i=0;i<n-t+1;++i){
    //     for(int j=0;j<n-t+1;++j){
    //         cout<<ans[i*(n-t+1)+j]<<" ";
    //     }
    //     cout<<"\n";
    // }
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
    printf("%f\n",elapsed.count());
    return 0;
}