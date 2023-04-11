#include<iostream>
#include<cmath>
#include<iomanip>
#include<cstring>
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
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={},W[MAXN]={};
int N=0,M=0;
int l=0,r[MAXN]={};
int limit=1;
// 15=1111b   b[0] b[0] b[1] b[1]
// 5=0101b    d[0] c[0] d[1] c[1]   
// r  b[0]*d[0]  b[0]*c[0]  b[1]*d[1]  b[1]*c[1]
//    a[0] a[0] a[1] a[1]
//(a+bi)(c+di)=(ac-bd)+(ad+bc)i
// fmaddsub   a[0]*c[0]-b[0]*d[0] a[0]*d[0]+b[0]*c[0] a[1]*c[1]-d[1]*d[1] a[1]*d[1]-b[1]*c[1]
// fmaddsub(a,b,c)  res = a*b +/- c
__m128 _mm_complex_mul(Complex* x, Complex* y)//利用SSE指令集同时进行两个复数乘法
{
    __m128 a_vec,b_vec,r;
    a_vec = _mm_loadu_ps((float*)(x));//a[0] b[0] a[1] b[1]
    b_vec = _mm_loadu_ps((float*)(y));//c[0] d[0] c[1] d[1]    
    r = _mm_permute_ps(a_vec,245)*_mm_permute_ps(b_vec,177);
    r = _mm_fmaddsub_ps(b_vec,_mm_permute_ps(a_vec,160),r);
    return r;
}
/*Wn = cos(PI/mid) +/- sin(PI/mid)i
  [1,Wn,WN^2,Wn^3,...,Wn^(mid-1)]
  [1,Wn] *Wn^2-> [Wn^2,Wn^3] *Wn^2-> [Wn^4,Wn^5] *Wn^2-> ...

  Vec_A = [A[j],A[j+1],...,A[j+mid-1]]
  Vec_B = [A[j+mid],Wn*A[j+mid+1],(Wn^2)*A[j+mid+2],...,(Wn^(mid-1))*A[j+2*mid-1]]
  Vec_A = Vec_A + Vec_B
  Vec_B = Vec_A - Vec_B
*/
void FFT(Complex* A,int type)
{
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    for(int mid=1;mid<limit;mid<<=1){
        Complex Wn{cos(PI/mid),type*sin(PI/mid)};
        if(mid >= 4){//预处理单位根
            Complex Wn_2 = Wn*Wn;
            Complex Mul_epi[2]={Wn_2,Wn_2};
            W[0] = Complex{1,0};
            W[1] = Complex{1,0}*Wn;
            for(int k=2;k<mid;k+=2){
                __m128 t = _mm_complex_mul(W+k-2,Mul_epi);
                _mm_storeu_ps((float*)(W+k),t);
            }
            for(int R=mid<<1,j=0;j<limit;j+=R){
                Complex w{1,0};
                for(int k=0;k<mid;k+=2){
                    __m128 t1 = _mm_loadu_ps((float*)(A+j+k));
                    __m128 t2 = _mm_complex_mul(A+j+k+mid,W+k);
                    __m128 t3 = _mm_add_ps(t1,t2);
                    t2 = _mm_sub_ps(t1,t2);
                    _mm_storeu_ps((float*)(A+j+k),t3);
                    _mm_storeu_ps((float*)(A+j+k+mid),t2);
                }
            }
        }
        else if(mid < limit>>1){
            for(int R=mid<<1,j=0;j<limit;j+=2*R){
                Complex w{1,0};
                Complex A_t[2],B_t[2];
                __m128 t1,t2,t3;
                for(int k=0;k<mid;++k,w=w*Wn){
                    A_t[0] = A[j+k];B_t[0] = w*A[j+k+mid];
                    A_t[1] = A[j+R+k];B_t[1] = w*A[j+k+R+mid];
                    t1 = _mm_loadu_ps((float*)A_t);
                    t2 = _mm_loadu_ps((float*)B_t);
                    t3 = _mm_add_ps(t1,t2);
                    t2 = _mm_sub_ps(t1,t2);
                    _mm_storeu_ps((float*)A_t,t3);
                    _mm_storeu_ps((float*)B_t,t2);
                    A[j+k] = A_t[0];A[j+k+mid] = B_t[0];
                    A[j+R+k] = A_t[1];A[j+k+R+mid] = B_t[1];
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
    __m128 t2 = _mm_set1_ps(limit);
    for(int i=0;i<=limit;i+=2){
        __m128 t1 = _mm_loadu_ps((float*)(A+i));
        t1 = _mm_div_ps(t1,t2);
        _mm_storeu_ps((float*)(A+i),t1);
    }
    
}
void Conv(Complex* A,Complex* B,Complex* ans)
{
    cin>>N>>M;
    for(int i=0;i<N;++i){
        A[i].Re=rand()%10+1;
        //cin>>A[i].Re;
    }
    //输入b为卷积核的翻转。
    for(int i=0;i<M;++i){
        B[i].Re=rand()%10+1;
        //cin>>B[i].Re;
    }
    auto start = std::chrono::high_resolution_clock::now();
    while(limit<=N+M-2){limit<<=1;++l;}
    for(int i=0;i<limit;++i){r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));}

    FFT(A,1);
    FFT(B,1);
    
    for(int i=0;i<=limit;i+=2){
        __m128 t = _mm_complex_mul(A+i,B+i);
        _mm_storeu_ps((float*)(ans+i),t);
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
        //cout<<(int)(c[i].Re+0.5)<<" ";
    }//cout<<"\n";
    return 0;
}