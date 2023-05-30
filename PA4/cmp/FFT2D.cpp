//g++ cmp/FFT2D.cpp -o test -g -std=c++11 
#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>
#include<mpi.h>
using namespace std;
const int MAXN=6005;
const float PI=acos(-1.0);
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
void FFT(Complex* A,int type)
{
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    for(int mid=1;mid<limit;mid<<=1){
        Complex Wn{cos(PI/mid),type*sin(PI/mid)};
        for(int R=mid<<1,j=0;j<limit;j+=R){
            Complex w{1,0};
            for(int k=0;k<mid;++k,w=w*Wn){
                Complex x=A[j+k],y=w*A[j+mid+k];
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
        }
    }
    if(type == 1){return;}
    for(int i=0;i<=limit;++i){
        A[i].Re=A[i].Re/limit;
        A[i].Im=A[i].Im/limit;
    }
}
void init()
{
    n = 2000,m = 2000,s = 2000,t = 2000;
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
    for(int i=0;i<n;++i){
        FFT(A[i],type);
    }
    //for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<A[i][j].Re<<" "<<A[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";
    for(int j=0;j<m;++j){
        for(int i=0;i<n;++i){
            rev_col[i] = A[i][j];
        }
        FFT(rev_col,type);
        for(int i=0;i<n;++i){
            A[i][j] = rev_col[i];
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
    for(int i = 0;i<=n+s-2;++i){
        for(int j = 0;j<=m+t-2;++j){
            //cout<<c[i][j].Re<<" ";
        }
        //cout<<"\n";
    }
    return 0;
}