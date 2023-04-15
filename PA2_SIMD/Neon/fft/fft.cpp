//g++ -march=native test.cpp -o test -g -std=c++11
#include<iostream>
#include<cmath>
#include<iomanip>
#include<cstring>
#include<arm_neon.h>
#include<chrono>
using namespace std;
const int MAXN=30000005;
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
float neg[4]={-1,1,-1,1};
void V_complex_mul(Complex* x, Complex* y,Complex* c,float* neg)//利用neon指令集同时进行两个复数乘法,结果存储在C中
{
    // float32x4_t a = vld1q_f32((float*)x);
    // float32x4_t b = vld1q_f32((float*)y);
    // float32x4x2_t TMP = vtrnq_f32(a,a);
    // float32x4_t a1 = TMP.val[1]; float32x4_t a2 = TMP.val[0];
    // float32x4_t b1 = vrev64q_f32(b);
    // float32x4_t r = vmulq_f32(a1,b1);
    // r = vmulq_f32(r,*(float32x4_t*)neg);
    // r = vmlaq_f32(r,a2,b);
    // vst1q_f32((float*)c,r);
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
        if(mid >= 4){
            Complex Wn_2 = Wn*Wn;
            Complex Mul_epi[2]={Wn_2,Wn_2};
            W[0] = Complex{1,0};
            W[1] = Complex{1,0}*Wn;
            for(int k=2;k<mid;k+=2){
                V_complex_mul(W+k-2,Mul_epi,W+k,neg);
            }
            for(int R=mid<<1,j=0;j<limit;j+=R){
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
                    // float32x4_t t1 = vld1q_f32((float*)(A+j+k));
                    // float32x4_t t2 = vld1q_f32((float*)(A+j+k+mid));
                    // float32x4_t t3 = vaddq_f32(t1,t2);
                    // t2 = vsubq_f32(t1,t2);
                    // vst1q_f32((float*)(A+j+k),t3);
                    // vst1q_f32((float*)(A+j+k+mid),t2);
                }   
            }
        }
        else if(mid < limit>>1){
            for(int R=mid<<1,j=0;j<limit;j+=2*R){
                Complex w{1,0};
                Complex A_t[2],B_t[2];
                float32x4_t t1,t2,t3;
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
                    // t1 = vld1q_f32((float*)A_t);
                    // t2 = vld1q_f32((float*)B_t);
                    // t3 = vaddq_f32(t1,t2);
                    // t2 = vsubq_f32(t1,t2);
                    // vst1q_f32((float*)A_t,t3);
                    // vst1q_f32((float*)B_t,t2);
                    A[j+k] = A_t[0];A[j+k+R] = A_t[1];
                    A[j+k+mid] = B_t[0];A[j+k+R+mid] = B_t[1];
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
    float LIM[4]={(float)limit,(float)limit,(float)limit,(float)limit};
    float* LIMt = (float*)LIM;
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
        // float32x4_t t1 = vld1q_f32(A_t);
        // float32x4_t t2 = vld1q_f32((float*)LIM);
        // t1 = vdivq_f32(t1,t2);
        // vst1q_f32(A_t,t1);
    }
    
}
void Conv(Complex* A,Complex* B,Complex* ans,int N,int M)
{
    //cin>>N>>M;
    //N = 5000000;M = 5000000;
    //N=10,M=10;
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
        V_complex_mul(A+i,B+i,ans+i,neg);
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
    // Conv(a,b,c);
    // for(int i=0;i<N+M-1;++i){
    //     //cout<<(int)(c[i].Re+0.5)<<" ";
    // }//cout<<"\n";
    for(int i=0;i<=23;++i){
        Conv(a,b,c,1<<i,1<<i);
    }
    return 0;
}