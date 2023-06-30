#include<iostream>
#include<random>
#include<chrono>
using namespace std;
const int MAXN = 2060,MAXR = 105;
float a[MAXN*MAXN]={},R[MAXR*MAXR]={};
float ans[MAXN*MAXN]={};
float g[16]={},d[MAXN*MAXN/4][16]={},P[16]={};
float re_g[16]={},re_d[16]={},result[4]={};
int n = 0,t = 0;
void init()
{
    n = 1024, t = 3;
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            a[i*n+j] = rand()%10+1;
        }
    }
    for(int i=0;i<t;++i){
        for(int j=0;j<t;++j){
            g[i*t+j] = rand()%10+1;
        }
    }
}
void F_2x2_3x3(float* g,float* d,float* re_d,float* re_g,float* P,float* result)
{
    float tmp[16]={};//BT x d
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

    re_d[0] = tmp[0] - tmp[2];
    re_d[4] = tmp[4] - tmp[6];
    re_d[8] = tmp[8] - tmp[10];
    re_d[12] = tmp[12] - tmp[14];

    re_d[1] = tmp[1] + tmp[2];
    re_d[5] = tmp[5] + tmp[6];
    re_d[9] = tmp[9] + tmp[10];
    re_d[13] = tmp[13] + tmp[14];

    re_d[2] = -tmp[1] + tmp[2];
    re_d[6] = -tmp[5] + tmp[6];
    re_d[10] = -tmp[9] + tmp[10];
    re_d[14] = -tmp[13] + tmp[14];

    re_d[3] = tmp[1] - tmp[3];
    re_d[7] = tmp[5] - tmp[7];
    re_d[11] = tmp[9] - tmp[11];
    re_d[15] = tmp[13] - tmp[15];


    tmp[0] = g[0];
    tmp[1] = g[1];
    tmp[2] = g[2];

    tmp[3] = (g[0] + g[3] + g[6])/2;
    tmp[4] = (g[1] + g[4] + g[7])/2;
    tmp[5] = (g[2] + g[5] + g[8])/2;

    tmp[6] = (g[0] - g[3] + g[6])/2;
    tmp[7] = (g[1] - g[4] + g[7])/2;
    tmp[8] = (g[2] - g[5] + g[8])/2;

    tmp[9] = g[6];
    tmp[10] = g[7];
    tmp[11] = g[8];

    re_g[0] = tmp[0];
    re_g[1] = (tmp[0] + tmp[1] + tmp[2])/2;
    re_g[2] = (tmp[0] - tmp[1] + tmp[2])/2;
    re_g[3] = tmp[2];

    re_g[4] = tmp[3];
    re_g[5] = (tmp[3] + tmp[4] + tmp[5])/2;
    re_g[6] = (tmp[3] - tmp[4] + tmp[5])/2;
    re_g[7] = tmp[5];

    re_g[8] = tmp[6];
    re_g[9] = (tmp[6] + tmp[7] + tmp[8])/2;
    re_g[10] = (tmp[6] - tmp[7] + tmp[8])/2;
    re_g[11] = tmp[8];

    re_g[12] = tmp[9];
    re_g[13] = (tmp[9] + tmp[10] + tmp[11])/2;
    re_g[14] = (tmp[9] - tmp[10] + tmp[11])/2;
    re_g[15] = tmp[11];

    P[0] = re_d[0]*re_g[0];
    P[1] = re_d[1]*re_g[1];
    P[2] = re_d[2]*re_g[2];
    P[3] = re_d[3]*re_g[3];
    P[4] = re_d[4]*re_g[4];
    P[5] = re_d[5]*re_g[5];
    P[6] = re_d[6]*re_g[6];
    P[7] = re_d[7]*re_g[7];
    P[8] = re_d[8]*re_g[8];
    P[9] = re_d[9]*re_g[9];
    P[10] = re_d[10]*re_g[10];
    P[11] = re_d[11]*re_g[11];
    P[12] = re_d[12]*re_g[12];
    P[13] = re_d[13]*re_g[13];
    P[14] = re_d[14]*re_g[14];
    P[15] = re_d[15]*re_g[15];

    tmp[0] = P[0] + P[4] + P[8];
    tmp[1] = P[1] + P[5] + P[9];
    tmp[2] = P[2] + P[6] + P[10];
    tmp[3] = P[3] + P[7] + P[11];
    
    tmp[4] = P[4] - P[8] - P[12];
    tmp[5] = P[5] - P[9] - P[13];
    tmp[6] = P[6] - P[10] - P[14];
    tmp[7] = P[7] - P[11] - P[15];

    result[0] = (tmp[0] + tmp[1] + tmp[2]);
    result[2] = (tmp[4] + tmp[5] + tmp[6]);
    result[1] = (tmp[1] - tmp[2] - tmp[3]);
    result[3] = (tmp[5] - tmp[6] - tmp[7]);
}
void img2col()
{
    int cnt = -1;
    for(int i=0;i<=n-4;i+=2){
        for(int j=0;j<=n-4;j+=2){
            ++cnt;
            d[cnt][0] = a[i*n+j];d[cnt][1] = a[i*n+j+1];
            d[cnt][2] = a[i*n+j+2];d[cnt][3] = a[i*n+j+3];

            d[cnt][4] = a[(i+1)*n+j];d[cnt][5] = a[(i+1)*n+j+1];
            d[cnt][6] = a[(i+1)*n+j+2];d[cnt][7] = a[(i+1)*n+j+3];

            d[cnt][8] = a[(i+2)*n+j];d[cnt][9] = a[(i+2)*n+j+1];
            d[cnt][10] = a[(i+2)*n+j+2];d[cnt][11] = a[(i+2)*n+j+3];

            d[cnt][12] = a[(i+3)*n+j];d[cnt][13] = a[(i+3)*n+j+1];
            d[cnt][14] = a[(i+3)*n+j+2];d[cnt][15] = a[(i+3)*n+j+3];
        }
    }
}
void Wiongrad()
{
    init();
    img2col();
    auto start = std::chrono::high_resolution_clock::now();
    int cnt = -1;
    for(int i=0;i<=n-4;i+=2){
        for(int j=0;j<=n-4;j+=2){
            ++cnt;
            F_2x2_3x3(g,d[cnt],re_d,re_g,P,result);
            
            ans[i*(n-t+1)+j] = result[0];
            ans[i*(n-t+1)+j+1] = result[1];
            ans[(i+1)*(n-t+1)+j] = result[2];
            ans[(i+1)*(n-t+1)+j+1] = result[3];
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);

    Wiongrad();

    // for(int i=0;i<n-t+1;++i){
    //     for(int j=0;j<n-t+1;++j){
    //         cout<<ans[i*(n-t+1)+j]<<" ";
    //     }cout<<"\n";
    // }
    return 0;
}