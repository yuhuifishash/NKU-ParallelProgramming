#include<iostream>
#include<cmath>
#include<chrono>
using namespace std;
const int MAXN=20000005;
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
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={};
int N=0,M=0;
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
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
    
    FFT(A,1);
    FFT(B,1);

    for(int i=0;i<=limit;++i){
        ans[i]=A[i]*B[i];
    }
    FFT(ans,-1);
    for(int i=0;i<=N+M-2;++i){
        ans[i].Re=(int)(ans[i].Re/limit+0.5);
    }
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
}
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    Conv(a,b,c);
    for(int i=0;i<N+M-1;++i){
        //cout<<(int)c[i].Re<<" ";
    }
    return 0;
}