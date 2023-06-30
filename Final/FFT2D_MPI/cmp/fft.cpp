//mpic++ FFT_blocks/cmp/fft.cpp -o test -g -std=c++11 -fopenmp
#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>
#include<mpi.h>
using namespace std;
const int MAXN=4005,MAXR = 105;
const double PI=acos(-1.0);
class Complex
{
public:
    double Re;
    double Im;
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
Complex a[MAXN*MAXN]={},b[MAXR][MAXR]={},c[MAXN][MAXN]={};
Complex Core[MAXR*MAXR*10]={};//b展开为一维之后的FFT结果
Complex W[2][15][MAXN]={};
int n=0,m=0,t=0;
int l=0,r[MAXN]={};
int limit=1;
int blocks = 24;
void FFT(Complex* A,int type)
{
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    for(int s=0;(1<<s)<limit;++s){
        int mid = 1<<s;
        for(int R=mid<<1,j=0;j<limit;j+=R){
            for(int k=0;k<mid;++k){
                Complex x=A[j+k],y=W[type == 1 ? 0:1][s][k]*A[j+mid+k];
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
        }
    }
    if(type == 1){return;}
    for(int i=0;i<limit;++i){
        A[i].Re=A[i].Re/limit;
        A[i].Im=A[i].Im/limit;
    }
}
void init()
{
    //cin>>n>>m>>t;
    n = 2000,m = 2000,t = 3;
    //n = 15,m = 15,t = 7;
    for(int i=0;i<n;++i){
        for(int j=0;j<m;++j){
            //cin>>a[i][j].Re;
            a[i*m+j].Re = rand()%10+1;
        }
    }
    for(int i=0;i<t;++i){
        for(int j=0;j<t;++j){
            //cin>>b[i][j].Re;
            b[i][j].Re = rand()%10+1;
        }
    }
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
void FFT_2D_row(Complex* A,int type)
{
    for(int i=0;i<limit;++i){
        FFT(A+i*limit,type);
    }
}
void FFT_2D_col(Complex* A,int type)
{
    Complex* Trans = new Complex[limit*limit+limit];
    for(int i=0;i<limit;++i){
        for(int j=0;j<limit;++j){
            Trans[i*limit+j] = A[j*limit+i];
        }
    }
    FFT_2D_row(Trans,type);
    for(int i=0;i<limit;++i){
        for(int j=0;j<limit;++j){
            A[i*limit+j] = Trans[j*limit+i];
        }
    }
}
void FFT_2D(Complex* A,int type)
{
    FFT_2D_row(A,type);
    FFT_2D_col(A,type);
}
void CONV_2D(Complex* A,Complex* B,Complex* C)
{
    FFT_2D(A,1);
    for(int i = 0;i<limit;++i){
        for(int j = 0;j<limit;++j){
            C[i*limit+j] = A[i*limit+j]*B[i*limit+j];
        }
    }
    FFT_2D(C,-1);
}
void CONV_blocks()
{
    init();
    auto start = std::chrono::high_resolution_clock::now();

    while(limit<=blocks+t-2){
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
    init_root(limit);
    for(int i=0;i<limit;++i){
        for(int j=0;j<limit;++j){
            Core[i*limit+j] = b[i][j];
        }
    }
    FFT_2D(Core,1);

    for(int i=0;i<n;i+=blocks){
        for(int j=0;j<m;j+=blocks){
            Complex* t1 = new Complex[limit*limit+limit];
            Complex* t2 = new Complex[limit*limit+limit];
            for(int k=0;k<blocks+t;++k){
                for(int s=0;s<blocks+t;++s){
                    t1[k*limit+s] = a[(k+i)*m+s+j];
                }
            }
            CONV_2D(t1,Core,t2);

            for(int k=0;k<blocks;++k){
                for(int s=0;s<blocks;++s){
                    c[k+i][s+j] = t2[(k+t-1)*limit+(s+t-1)];
                }
            }
            delete []t1;
            delete []t2;
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
    CONV_blocks();   
    for(int i = 0;i<n-t+1;++i){
        for(int j = 0;j<m-t+1;++j){
            //cout<<c[i][j].Re<<" ";
        }
        //cout<<"\n";
    }
    return 0;
}