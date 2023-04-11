#include<iostream>
#include<cmath>
#include<immintrin.h>
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
        if(mid < limit>>2){
            for(int R=mid<<1,j=0;j<limit;j+=4*R){
                Complex w{1,0};
                Complex A_t[4],B_t[4];
                __m256 t1,t2,t3;
                for(int k=0;k<mid;++k,w=w*Wn){
                    A_t[0] = A[j+k];A_t[1] = A[j+k+R];
                    A_t[2] = A[j+k+2*R];A_t[3] = A[j+k+3*R];
                    B_t[0] = w*A[j+k+mid];B_t[1] = w*A[j+k+R+mid];
                    B_t[2] = w*A[j+k+2*R+mid];B_t[3] = w*A[j+k+3*R+mid];
                    t1 = _mm256_loadu_ps((float*)A_t);
                    t2 = _mm256_loadu_ps((float*)B_t);
                    t3 = _mm256_add_ps(t1,t2);
                    t2 = _mm256_sub_ps(t1,t2);
                    _mm256_storeu_ps((float*)A_t,t3);
                    _mm256_storeu_ps((float*)B_t,t2);
                    A[j+k] = A_t[0];A[j+k+R] = A_t[1];
                    A[j+k+2*R] = A_t[2];A[j+k+3*R] = A_t[3];
                    A[j+k+mid] = B_t[0];A[j+k+R+mid] = B_t[1];
                    A[j+k+2*R+mid] = B_t[2];A[j+k+3*R+mid] = B_t[3];
                }
            }
        }
        else{
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
    if(type == 1){return;}
    __m256 t2 = _mm256_set1_ps(limit);
    for(int i=0;i<=limit;i+=4){
        __m256 t1 = _mm256_loadu_ps((float*)(A+i));
        t1 = _mm256_div_ps(t1,t2);
        _mm256_storeu_ps((float*)(A+i),t1);
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