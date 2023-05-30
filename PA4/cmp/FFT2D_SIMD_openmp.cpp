//g++ cmp/FFT2D_SIMD_openmp.cpp -o test -g -std=c++11 -fopenmp
#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>
#include<omp.h>
using namespace std;
const int MAXN=6005;
const float PI=acos(-1.0);
int thread_count = 8;
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
Complex a[MAXN][MAXN]={},b[MAXN][MAXN]={},c[MAXN][MAXN]={};
Complex rev_col[MAXN]={};
int n=0,m=0,s=0,t=0;
int l=0,r[MAXN]={};
int limit=1;
Complex W[MAXN]={};float neg[4]={-1,1,-1,1};
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
void init()
{
    n = 2000,m = 2000,s = 2000,t = 2000;
    //n = 10,m = 10,s = 10,t = 10;
    for(int i=0;i<n;++i){
        for(int j=0;j<m;++j){
            //cin>>a[i][j].Re;
            a[i][j].Re = rand()%10+1;
        }
    }
    for(int i=0;i<s;++i){
        for(int j=0;j<t;++j){
            //cin>>b[i][j].Re;
            b[i][j].Re = rand()%10+1;
        }
    }

}
void FFT_2D(Complex A[][MAXN],int n,int m,int type)
{
    #pragma omp parallel num_threads(thread_count)
    {
        for(int i=0;i<n;++i){
            FFT(A[i],type);
        }
    //for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<A[i][j].Re<<" "<<A[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";
        for(int j=0;j<m;++j){
            #pragma omp for
            for(int i=0;i<n;++i){
                rev_col[i] = A[i][j];
            }
            FFT(rev_col,type);
            #pragma omp for
            for(int i=0;i<n;++i){
                A[i][j] = rev_col[i];
            }
        }
    }
}
void CONV_2D()
{
    while(limit<=max(n+s-2,m+t-2)){
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
    auto start = std::chrono::high_resolution_clock::now();
    FFT_2D(a,limit,limit,1);
    //for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<a[i][j].Re<<" "<<a[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";
    FFT_2D(b,limit,limit,1);
    for(int i = 0;i<=limit;++i){
        for(int j = 0;j<=limit;++j){
            c[i][j] = a[i][j]*b[i][j];
        }
    }
    FFT_2D(c,limit,limit,-1);
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    init();
    CONV_2D();   
    // for(int i = 0;i<=n+s-2;++i){
    //     for(int j = 0;j<=m+t-2;++j){
    //         cout<<c[i][j].Re<<" ";
    //     }
    //     cout<<"\n";
    // }
    return 0;
}