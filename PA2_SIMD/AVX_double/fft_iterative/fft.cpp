#include<iostream>
#include<cmath>
#include<immintrin.h>
#include<chrono>
using namespace std;
const int MAXN=20000005;
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
        if(mid < limit>>1){
            for(int R=mid<<1,j=0;j<limit;j+=2*R){
                Complex w{1,0};
                Complex A_t[2],B_t[2];
                for(int k=0;k<mid;++k,w=w*Wn){
                    A_t[0] = A[j+k];B_t[0] = w*A[j+k+mid];
                    A_t[1] = A[j+R+k];B_t[1] = w*A[j+k+R+mid];
                    __m256d t1 = _mm256_loadu_pd((double*)A_t);
                    __m256d t2 = _mm256_loadu_pd((double*)B_t);
                    __m256d t3 = _mm256_add_pd(t1,t2);
                    t2 = _mm256_sub_pd(t1,t2);
                    _mm256_storeu_pd((double*)A_t,t3);
                    _mm256_storeu_pd((double*)B_t,t2);
                    A[j+k] = A_t[0];A[j+k+mid] = B_t[0];
                    A[j+R+k] = A_t[1];A[j+k+R+mid] = B_t[1];
                }
            }
        }
        else{
            Complex w{1,0};
            for(int k=0;k<mid;++k,w=w*Wn){
                Complex x=A[k],y=w*A[mid+k];
                A[k]=x+y;
                A[mid+k]=x-y;
            }
        }
    }
    if(type == 1){return;}
    __m256d t2 = _mm256_set1_pd(limit);
    for(int i=0;i<=limit;i+=2){
        __m256d t1 = _mm256_loadu_pd((double*)(A+i));
        t1 = _mm256_div_pd(t1,t2);
        _mm256_storeu_pd((double*)(A+i),t1);
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