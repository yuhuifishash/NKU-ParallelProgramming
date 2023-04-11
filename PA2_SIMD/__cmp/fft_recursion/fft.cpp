#include<iostream>
#include<cmath>
#include<chrono>
using namespace std;
const int MAXN=20000005;
const double PI=acos(-1.0);
int N=0,M=0,limit=1;
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
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={};
void FFT(Complex* A,int n,int type)
{
    if(n == 1){return;}
    int m = n>>1;
    Complex* A_t1 = new Complex[m];
    Complex* A_t2 = new Complex[m];
    for(int i = 0;i < m;++i){
        A_t1[i] = A[2*i];
        A_t2[i] = A[2*i+1];
    }
    FFT(A_t1,m,type);
    FFT(A_t2,m,type);
    Complex Wn{cos(2*PI/n),type*sin(2*PI/n)};
    Complex w{1,0};
    for(int i = 0;i < m;++i,w=w*Wn){
        Complex x = A_t1[i],y = w*A_t2[i];
        A[i] = x+y;
        A[i+m] = x-y;
    }
    delete []A_t1;
    delete []A_t2;
}
void Conv(Complex* A,Complex* B,Complex* ans)
{
    cin>>N>>M;
    for(int i=0;i<N;++i){
        a[i].Re=rand()%10+1;
        //cin>>A[i].Re;
    }
    //输入b为卷积核的翻转。
    for(int i=0;i<M;++i){
        b[i].Re=rand()%10+1;
        //cin>>B[i].Re;
    }

    auto start = std::chrono::high_resolution_clock::now();
    while(limit<=N+M-2){
        limit<<=1;
    }
    FFT(A,limit,1);
    FFT(B,limit,1);

    for(int i=0;i<=limit;++i){
        ans[i]=A[i]*B[i];
    }

    FFT(ans,limit,-1);
    for(int i=0;i<=limit;++i){
        ans[i].Re/=limit;
        ans[i].Im/=limit;
    }

    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    Conv(a,b,c);
    for(int i=0;i<N+M-1;++i){
        //cout<<(int)c[i].Re<<" ";
    }//cout<<"\n";
    return 0;
}