//g++ -march=armv8-a fft_final_optimize_openmp/fft.cpp -o test -g -std=c++11 -fopenmp
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<iostream>
#include<omp.h>
#include<chrono>
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
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={},W[MAXN]={};
int n=0,m=0,l=0,r[MAXN]={},limit=1;
float neg[4]={-1,1,-1,1};
void init()
{
    //n=10,m=10;
    n=5000000,m=5000000;
    for(int i=0;i<n;++i){
        a[i].Re=rand()%10+1;
        //scanf("%lf",&a[i].Re);
    }
    for(int i=0;i<m;++i){
        b[i].Re=rand()%10+1;
        //scanf("%lf",&b[i].Re);
    }
}
void FFT(Complex* A,int type)
{
    #pragma omp single
    {   
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    }
    #pragma omp barrier
    for(int mid=1;mid<limit;mid<<=1){
        Complex Wn{cos(PI/mid),type*sin(PI/mid)};
        int R=mid<<1;
        if(mid >= 4){
            #pragma omp single
            {
            Complex Wn_2 = Wn*Wn;
            Complex Mul_epi[2]={Wn_2,Wn_2};
            W[0] = Complex{1,0};
            W[1] = Complex{1,0}*Wn;
            for(int k=2;k<mid;k+=2){
                V_complex_mul(W+k-2,Mul_epi,W+k,neg);
            }
            }
            if(limit/(mid<<1) < thread_count ){
                #pragma omp for collapse(2)
                for(int j=0;j<limit;j+=R){
                    for(int k=0;k<mid;k+=2){
                        V_complex_mul(A+j+k+mid,W+k,A+j+k+mid,neg);
                        float* A1 = (float*)(A+j+k);
                        float* A2 = (float*)(A+j+k+mid);
                        asm volatile(
                            "ldr q0,[%0]\n"
                            "ldr q1,[%1]\n"
                            "fadd v2.4s,v1.4s,v0.4s\n"
                            "fsub v0.4s,v0.4s,v1.4s\n"
                            "str q2,[%0]\n"
                            "str q0,[%1]\n"
                            :"+r"(A1),"+r"(A2)
                            :
                            :"v0","v1","v2","q0","q1","q2","memory","cc"
                        );
                    }   
                }
            }
            else{
                #pragma omp for
                for(int j=0;j<limit;j+=R){
                    for(int k=0;k<mid;k+=2){
                        V_complex_mul(A+j+k+mid,W+k,A+j+k+mid,neg);
                        float* A1 = (float*)(A+j+k);
                        float* A2 = (float*)(A+j+k+mid);
                        asm volatile(
                            "ldr q0,[%0]\n"
                            "ldr q1,[%1]\n"
                            "fadd v2.4s,v1.4s,v0.4s\n"
                            "fsub v0.4s,v0.4s,v1.4s\n"
                            "str q2,[%0]\n"
                            "str q0,[%1]\n"
                            :"+r"(A1),"+r"(A2)
                            :
                            :"v0","v1","v2","q0","q1","q2","memory","cc"
                        );
                    }   
                }
            }
        }
        else if(mid < limit>>1){
            #pragma omp for
            for(int j=0;j<limit;j+=2*R){
                Complex w{1,0};
                Complex A_t[2],B_t[2];
                for(int k=0;k<mid;++k,w=w*Wn){
                    A_t[0] = A[j+k];A_t[1] = A[j+k+R];
                    B_t[0] = w*A[j+k+mid];B_t[1] = w*A[j+k+R+mid];
                    float* A_tt = (float*)A_t;
                    float* B_tt = (float*)B_t;
                    asm volatile(
                        "ldr q0,[%0]\n"
                        "ldr q1,[%1]\n"
                        "fadd v2.4s,v1.4s,v0.4s\n"
                        "fsub v0.4s,v0.4s,v1.4s\n"
                        "str q2,[%0]\n"
                        "str q0,[%1]\n"
                        :"+r"(A_tt),"+r"(B_tt)
                        :
                        :"v0","v1","v2","q0","q1","q2","memory","cc"
                    );
                    A[j+k] = A_t[0];A[j+k+R] = A_t[1];
                    A[j+k+mid] = B_t[0];A[j+k+R+mid] = B_t[1];
                }
            }
        }
        else{
            #pragma omp for
            for(int j=0;j<limit;j+=R){
                Complex w{1,0};
                for(int k=0;k<mid;++k,w=w*Wn){
                    Complex x=A[j+k],y=w*A[j+mid+k];
                    A[j+k]=x+y;
                    A[j+mid+k]=x-y;
                }
            }
        }
    }
    if(type == 1){return;}
    float LIM[4]={(float)limit,(float)limit,(float)limit,(float)limit};
    float* LIMt = (float*)LIM;
    #pragma omp for
    for(int i=0;i<=limit;i+=2){
        float* A_t = (float*)(A+i);
        asm volatile(
            "ldr q0,[%0]\n"
            "ldr q1,[%1]\n"
            "fdiv v0.4s,v0.4s,v1.4s\n"
            "str q0,[%0]\n"
            :"+r"(A_t),"+r"(LIMt)
            :
            :"v0","v1","q0","q1","memory","cc"
        );
    }
    
}
void Conv(Complex* A,Complex* B,Complex* C)
{
    init();

    auto start = std::chrono::high_resolution_clock::now();
    while(limit<=n+m-2){limit<<=1;++l;}
    for(int i=0;i<limit;++i){r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));}
    #pragma omp parallel num_threads(thread_count)
    {
        FFT(A,1);
        FFT(B,1);

        #pragma omp for
        for(int i=0;i<=limit;++i){
            C[i]=A[i]*B[i];
        }

        FFT(C,-1);
    }
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}
int main()
{
    Conv(a,b,c);
    for(int i=0;i<n+m-1;++i){
        //printf("%d ",(int)(c[i].Re+0.5));
    }
    return 0;
}