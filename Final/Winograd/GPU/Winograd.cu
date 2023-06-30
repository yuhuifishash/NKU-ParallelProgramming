#include<iostream>
#include<random>
#include<chrono>
#include <cuda_runtime.h>
#include"device_launch_parameters.h"
using namespace std;
const int MAXN = 2060, MAXR = 105;
float* a;
float* ans;
float* G, *D;

float* D_d;
float* ans_d;
float* G_d;
int n = 0, t = 0;
void init()
{
    //n = 16, t = 3;
    n = 1024, t = 3;
    a = (float*)malloc(sizeof(float) * n * n);
    G = (float*)malloc(sizeof(float) * t * t);
    ans = (float*)malloc(sizeof(float) * n * n);
    D = (float*)malloc(sizeof(float) * n * n * 8);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i * n + j] = rand() % 10 + 1;
            //cout << a[i * n + j] << " ";
        }//cout << "\n";
    }
    for (int i = 0; i < t; ++i) {
        for (int j = 0; j < t; ++j) {
            G[i * t + j] = rand() % 10 + 1;
        }
    }
}
__global__ void F_2x2_3x3(float* G, float* D, float* ans_d, int n,int t)
{
    int tx = threadIdx.x, ty = blockIdx.x;
    int cnt = ((n-4)/2+1)*tx+ty;
    
    //printf("%d %d %d\n", tx , ty , cnt);

    float g[16] = {};
    float tmp[16] = {};
    float d[16] = {};

    for (int i = 0; i < 9; ++i) {
        g[i] = G[i];
    }
    for (int i = 0; i < 16; ++i) {
        d[i] = D[cnt * 16 + i];
    }

    tmp[0] = d[0] - d[8];
    tmp[1] = d[1] - d[9];
    tmp[2] = d[2] - d[10];
    tmp[3] = d[3] - d[11];

    tmp[4] = d[4] + d[8];
    tmp[5] = d[5] + d[9];
    tmp[6] = d[6] + d[10];
    tmp[7] = d[7] + d[11];

    tmp[8] = -d[4] + d[8];
    tmp[9] = -d[5] + d[9];
    tmp[10] = -d[6] + d[10];
    tmp[11] = -d[7] + d[11];

    tmp[12] = d[4] - d[12];
    tmp[13] = d[5] - d[13];
    tmp[14] = d[6] - d[14];
    tmp[15] = d[7] - d[15];

    d[0] = tmp[0] - tmp[2];
    d[4] = tmp[4] - tmp[6];
    d[8] = tmp[8] - tmp[10];
    d[12] = tmp[12] - tmp[14];

    d[1] = tmp[1] + tmp[2];
    d[5] = tmp[5] + tmp[6];
    d[9] = tmp[9] + tmp[10];
    d[13] = tmp[13] + tmp[14];

    d[2] = -tmp[1] + tmp[2];
    d[6] = -tmp[5] + tmp[6];
    d[10] = -tmp[9] + tmp[10];
    d[14] = -tmp[13] + tmp[14];

    d[3] = tmp[1] - tmp[3];
    d[7] = tmp[5] - tmp[7];
    d[11] = tmp[9] - tmp[11];
    d[15] = tmp[13] - tmp[15];


    tmp[0] = g[0];
    tmp[1] = g[1];
    tmp[2] = g[2];

    tmp[3] = (g[0] + g[3] + g[6]) / 2;
    tmp[4] = (g[1] + g[4] + g[7]) / 2;
    tmp[5] = (g[2] + g[5] + g[8]) / 2;

    tmp[6] = (g[0] - g[3] + g[6]) / 2;
    tmp[7] = (g[1] - g[4] + g[7]) / 2;
    tmp[8] = (g[2] - g[5] + g[8]) / 2;

    tmp[9] = g[6];
    tmp[10] = g[7];
    tmp[11] = g[8];

    g[0] = tmp[0];
    g[1] = (tmp[0] + tmp[1] + tmp[2]) / 2;
    g[2] = (tmp[0] - tmp[1] + tmp[2]) / 2;
    g[3] = tmp[2];

    g[4] = tmp[3];
    g[5] = (tmp[3] + tmp[4] + tmp[5]) / 2;
    g[6] = (tmp[3] - tmp[4] + tmp[5]) / 2;
    g[7] = tmp[5];

    g[8] = tmp[6];
    g[9] = (tmp[6] + tmp[7] + tmp[8]) / 2;
    g[10] = (tmp[6] - tmp[7] + tmp[8]) / 2;
    g[11] = tmp[8];

    g[12] = tmp[9];
    g[13] = (tmp[9] + tmp[10] + tmp[11]) / 2;
    g[14] = (tmp[9] - tmp[10] + tmp[11]) / 2;
    g[15] = tmp[11];

    d[0] = d[0] * g[0];
    d[1] = d[1] * g[1];
    d[2] = d[2] * g[2];
    d[3] = d[3] * g[3];
    d[4] = d[4] * g[4];
    d[5] = d[5] * g[5];
    d[6] = d[6] * g[6];
    d[7] = d[7] * g[7];
    d[8] = d[8] * g[8];
    d[9] = d[9] * g[9];
    d[10] = d[10] * g[10];
    d[11] = d[11] * g[11];
    d[12] = d[12] * g[12];
    d[13] = d[13] * g[13];
    d[14] = d[14] * g[14];
    d[15] = d[15] * g[15];

    tmp[0] = d[0] + d[4] + d[8];
    tmp[1] = d[1] + d[5] + d[9];
    tmp[2] = d[2] + d[6] + d[10];
    tmp[3] = d[3] + d[7] + d[11];

    tmp[4] = d[4] - d[8] - d[12];
    tmp[5] = d[5] - d[9] - d[13];
    tmp[6] = d[6] - d[10] - d[14];
    tmp[7] = d[7] - d[11] - d[15];

    tx *= 2; 
    ty *= 2;
    ans_d[tx*(n-t+1)+ty] = (tmp[0] + tmp[1] + tmp[2]);
    ans_d[(tx+1)*(n-t+1)+ty] = (tmp[4] + tmp[5] + tmp[6]);
    ans_d[tx*(n-t+1)+ty+1] = (tmp[1] - tmp[2] - tmp[3]);
    ans_d[(tx+1)*(n-t+1)+ty+1] = (tmp[5] - tmp[6] - tmp[7]);
}
void img2col()
{
    int cnt = -1;
    for (int i = 0; i <= n - 4; i += 2) {
        for (int j = 0; j <= n - 4; j += 2) {
            ++cnt;
            //printf("%d %d %d\n", i/2, j/2, cnt);
            D[cnt*16] = a[i * n + j]; D[cnt*16+1] = a[i * n + j + 1];
            D[cnt*16+2] = a[i * n + j + 2]; D[cnt*16+3] = a[i * n + j + 3];

            D[cnt*16+4] = a[(i + 1) * n + j]; D[cnt*16+5] = a[(i + 1) * n + j + 1];
            D[cnt*16+6] = a[(i + 1) * n + j + 2]; D[cnt*16+7] = a[(i + 1) * n + j + 3];

            D[cnt*16+8] = a[(i + 2) * n + j]; D[cnt*16+9] = a[(i + 2) * n + j + 1];
            D[cnt*16+10] = a[(i + 2) * n + j + 2]; D[cnt*16+11] = a[(i + 2) * n + j + 3];

            D[cnt*16+12] = a[(i + 3) * n + j]; D[cnt*16+13] = a[(i + 3) * n + j + 1];
            D[cnt*16+14] = a[(i + 3) * n + j + 2]; D[cnt*16+15] = a[(i + 3) * n + j + 3];
        }
    }
}
void Wiongrad()
{
    init();
    img2col();

    cudaMalloc((float**)&D_d, sizeof(float) * n * n * 8);
    cudaMalloc((float**)&G_d, sizeof(float) * t * t);
    cudaMalloc((float**)&ans_d, sizeof(float) * n * n);

    cudaMemcpy(D_d, D, sizeof(float) * n * n * 8, cudaMemcpyHostToDevice);
    cudaMemcpy(G_d, G, sizeof(float) * t * t, cudaMemcpyHostToDevice);

    auto start = std::chrono::high_resolution_clock::now();

    dim3 Winograd_grid((n - 4) / 2 + 1);
    dim3 Winograd_block((n - 4) / 2 + 1);

    F_2x2_3x3 << <Winograd_grid, Winograd_block >> > (G_d, D_d, ans_d, n, t);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    printf("%f\n", elapsed.count());

    cudaMemcpy(ans, ans_d, sizeof(int) * n * n, cudaMemcpyDeviceToHost);

    cudaFree(D_d);
    cudaFree(G_d);
    cudaFree(ans_d);

    free(a);
    free(G);
    free(D);
}

int main()
{
    ios::sync_with_stdio(false);
	cin.tie(0);

	Wiongrad();

	/*for (int i = 0; i < n - t + 1; ++i) {
		for (int j = 0; j < n - t + 1; ++j) {
			cout << ans[i * (n - t + 1) + j] << " ";
		}cout << "\n";
	}*/

	free(ans);
	return 0;
}