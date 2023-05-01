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
struct thread_id
{
    int id;
    int type;
    Complex* arr;
};
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={};
int n=0,m=0,l=0,r[MAXN]={},limit=1;
void init()
{
    scanf("%d%d",&n,&m);
    for(int i=0;i<n;++i){
        a[i].Re=rand()%10+1;
        //scanf("%lf",&a[i].Re);
    }
    for(int i=0;i<m;++i){
        b[i].Re=rand()%10+1;
        //scanf("%lf",&b[i].Re);
    }
    while(limit<=n+m-2){
        limit<<=1;
        ++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
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
        int R=(mid<<1);

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

    if(type == -1){
        #pragma omp for
        for(int i=0;i<=limit;++i){
            A[i].Re=(int)(A[i].Re/limit+0.5);
        }
    }
}
int main()
{
    init();

    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel num_threads(thread_count)
    {
        FFT(a,1);
        FFT(b,1);

        #pragma omp for
        for(int i=0;i<=limit;++i){
            c[i]=a[i]*b[i];
        }

        FFT(c,-1);
    }
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
    for(int i=0;i<=n+m-2;++i){
        //printf("%d ",(int)c[i].Re);
    }
    return 0;
}