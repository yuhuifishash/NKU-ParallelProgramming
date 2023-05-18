//g++ -march=armv8-a fft_blocks/cmp1/fft.cpp -o test -g -std=c++11 -fopenmp
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<iostream>
#include<omp.h>
#include<chrono>
#include<cstring>
using namespace std;
const int MAXN=2e7+5;
const float PI=acos(-1.0);
const int thread_count = 8;
class Complex
{
public:
    float Re;
    float Im;
    Complex operator-(){
        return Complex{-Re,-Im};
    }
    Complex operator+(const Complex& b){
        return Complex{Re+b.Re,Im+b.Im};
    }
    Complex operator-(const Complex& b){
        return Complex{Re-b.Re,Im-b.Im};
    }
    Complex operator*(const Complex& b){
        return Complex{Re*b.Re-Im*b.Im,Im*b.Re+Re*b.Im};
    }
};
void V_complex_mul(Complex* x, Complex* y,Complex* c,float* neg)//利用neon指令集同时进行四个复数乘法,结果存储在C中
{
    asm volatile(
        "ldr q0,[%0]\n"
        "ldr q2,[%1]\n"
        "ldr q3,[%3]\n"
        "trn2 v1.4s,v0.4s,v0.4s\n"
        "rev64 v4.4s,v2.4s\n"
        "trn1 v0.4s,v0.4s,v0.4s\n"
        "fmul v1.4s,v1.4s,v4.4s\n"
        "fmul v0.4s,v0.4s,v2.4s\n"
        "fmla v0.4s,v1.4s,v3.4s\n"
        "str q0,[%2]\n"
        :"+r"((float*)x),"+r"((float*)y),"+r"((float*)c),"+r"(neg)
        :
        : "v0","v1","v2","v3","q0","q2","q3","memory","cc"
    );
}
Complex a[MAXN]={},b[205]={},ans[MAXN]={};
Complex W[2][20][100005]={};//预处理单位根
int n=0,m=0,l=0,r[MAXN]={},limit=1;
float neg[4]={-1,1,-1,1};
void init()
{
    //n=50,m=7;
    //n=8388608,m=7;
    limit = 1,l = 0;
    for(int i=0;i<n;++i){
        a[i].Re=rand()%10+1;
    }
    for(int i=0;i<m;++i){
        b[i].Re=rand()%10+1;
    }
}
void FFT(Complex* A,int type)
{
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    for(int s=0;(1<<s)<limit;++s){
        int mid = 1<<s;//log2(mid)
        int R=mid<<1;
        for(int j=0;j<limit;j+=R){
            for(int k=0;k<mid;++k){
                Complex x=A[j+k],y=W[type == 1 ? 0:1][s][k]*A[j+mid+k];
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
        }
    }
    if(type == 1){return;}
    for(int i=0;i<=limit;++i){
        A[i].Re/=limit;
        A[i].Im/=limit;
    }
}
int blocks = 24;
void Conv_blocks(Complex* A,Complex* B,Complex* C)
{
    FFT(A,1);
    for(int i=0;i<=limit;i+=2){
        V_complex_mul(A+i,B+i,C+i,neg);
    }
    FFT(C,-1);

}
void init_root(int limit)
{
    for(int s=0;(1<<s)<limit;++s){
        int mid = (1<<s);
        Complex Wn{cos(PI/mid),sin(PI/mid)};
        Complex w{1,0};
        for(int k=0;k<mid;++k,w=w*Wn){
            W[0][s][k] = w;
        }
    }

    for(int s=0;(1<<s)<limit;++s){
        int mid = (1<<s);
        Complex Wn{cos(PI/mid),-sin(PI/mid)};
        Complex w{1,0};
        for(int k=0;k<mid;++k,w=w*Wn){
            W[1][s][k] = w;
        }
    }
}
void Conv(Complex* A,Complex* B,Complex* C)
{
    init();

    auto start = std::chrono::high_resolution_clock::now();

    while(limit<=blocks+m-2){limit<<=1;++l;}
    for(int i=0;i<limit;++i){r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));}
    init_root(limit);

    FFT(B,1);

    for(int i=0;i<n;i+=blocks){
        Complex* ta = new Complex[limit+5];
        Complex* tc = new Complex[limit+5];

        memcpy(ta,A+i,sizeof(Complex)*(blocks+m));
        Conv_blocks(ta,B,tc);
        memcpy(C+i,tc+m-1,sizeof(Complex)*blocks);

        delete []ta;
        delete []tc;
    }

    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}
int main()
{
    int t[8]={3,5,7,9,11,31,101};
    for(int i=0;i<7;++i){
        n = 8388608,m = t[i];
        blocks = 253 - m;
        Conv(a,b,ans);
    }
    return 0;
}