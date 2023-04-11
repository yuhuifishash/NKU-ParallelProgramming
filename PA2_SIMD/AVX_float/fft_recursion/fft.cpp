#include<iostream>
#include<cmath>
#include<chrono>
#include<immintrin.h>
using namespace std;
const int MAXN=20000005;
const float PI=acos(-1.0);
int N=0,M=0,limit=1;
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
void FFT(Complex* A,int n,int type)
{
    if(n == 1){return;}
    int m = n>>1;
    Complex* A_t1 = new Complex[m+2];
    Complex* A_t2 = new Complex[m+2];
    for(int i = 0;i < m;++i){
        A_t1[i] = A[i<<1];
        A_t2[i] = A[(i<<1)+1];
    }
    FFT(A_t1,m,type);
    FFT(A_t2,m,type);
    Complex Wn{cos(2*PI/n),type*sin(2*PI/n)};
    Complex w{1,0};
    __m256 t1,t2,t3;
    if(m >= 4){
        for(int i = 0;i < m;i+=4,w=w*Wn){
            A_t2[i] = w*A_t2[i];
            w = w*Wn; 
            A_t2[i+1] = w*A_t2[i+1];
            w = w*Wn;
            A_t2[i+2] = w*A_t2[i+2];
            w = w*Wn;
            A_t2[i+3] = w*A_t2[i+3];
            t1 = _mm256_loadu_ps((float*)(A_t1+i));
            t2 = _mm256_loadu_ps((float*)(A_t2+i));
            t3 = _mm256_add_ps(t1,t2);
            t2 = _mm256_sub_ps(t1,t2);
            _mm256_storeu_ps((float*)(A+i),t3);
            _mm256_storeu_ps((float*)(A+i+m),t2);
        }
    }
    else{
        for(int i = 0;i < m;++i,w=w*Wn){
            Complex x = A_t1[i],y = w*A_t2[i];
            A[i] = x+y;
            A[i+m] = x-y;
        }
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

    __m256 t2 = _mm256_set1_ps(limit);
    for(int i=0;i<=limit;i+=4){
        __m256 t1 = _mm256_loadu_ps((float*)(ans+i));
        t1 = _mm256_div_ps(t1,t2);
        _mm256_storeu_ps((float*)(ans+i),t1);
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
    }//cout<<"\n";
    return 0;
}