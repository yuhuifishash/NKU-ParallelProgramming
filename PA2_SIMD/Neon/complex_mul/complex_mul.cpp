#include<iostream>
#include<cmath>
#include<arm_neon.h>
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
//g++ -march=armv8-a test.cpp -o test -g -std=c++11
Complex A[MAXN],B[MAXN],C[MAXN];
float neg[4]={-1,1,-1,1};
void V_complex_mul(Complex* x, Complex* y,Complex* c,float* neg)//利用neon指令集同时进行两个复数乘法,结果存储在C中
{
    float32x4_t a = vld1q_f32((float*)x);
    float32x4_t b = vld1q_f32((float*)y);
    float32x4x2_t TMP = vtrnq_f32(a,a);
    float32x4_t a1 = TMP.val[1]; float32x4_t a2 = TMP.val[0];
    float32x4_t b1 = vrev64q_f32(b);
    float32x4_t r = vmulq_f32(a1,b1);
    r = vmulq_f32(r,*(float32x4_t*)neg);
    r = vmlaq_f32(r,a2,b);
    vst1q_f32((float*)c,r);
    // asm volatile(
    //     "ldr q0,[%0]\n"
    //     "ldr q2,[%1]\n"
    //     "ldr q3,[%3]\n"
    //     "trn2 v1.4s,v0.4s,v0.4s\n"
    //     "rev64 v4.4s,v2.4s\n"
    //     "trn1 v0.4s,v0.4s,v0.4s\n"
    //     "fmul v1.4s,v1.4s,v4.4s\n"
    //     "fmul v0.4s,v0.4s,v2.4s\n"
    //     "fmla v0.4s,v1.4s,v3.4s\n"
    //     "str q0,[%2]\n"
    //     :"+r"((float*)x),"+r"((float*)y),"+r"((float*)c),"+r"(neg)
    //     :
    //     : "v0","v1","v2","v3","q0","q2","q3","memory","cc"
    // );
}
int main()
{
    for(int i=0;i<5000000;++i){
        A[i].Re = 3.0+0.1*i; A[i].Im = 4.0+i;
        B[i].Re = -3.0-i; B[i].Im = -12-3*i;
    }
    auto start = std::chrono::high_resolution_clock::now();
    
    for(int i=1;i<=5000000;++i){
        C[i]=A[i]*B[i];
    }
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish-start;
    std::cout<<1000*elapsed.count()<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<5000000;i+=2){
        V_complex_mul(A+i,B+i,C+i,neg);
    }
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()<<"\n";
    return 0;
    // for(int i = 0;i < 4;++i){
    //     A[i].Re = 3 + i; A[i].Im = -4 + i;
    //     B[i].Re = 4 + i; B[i].Im = 1 + i;
    // }
    // V_complex_mul(A,B,C,neg);
    // for(int i = 0;i < 2;++i){
    //     cout<<C[i].Re<<" "<<C[i].Im<<" ";
    // }
    // return 0;
}