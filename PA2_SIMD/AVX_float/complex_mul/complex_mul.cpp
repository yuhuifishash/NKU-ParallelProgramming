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
Complex A[MAXN],B[MAXN],C[MAXN];
// float tmp[MAXN];
// void print()
// {
//     cout<<tmp[0]<<" "<<tmp[1]<<" "<<tmp[2]<<" "<<tmp[3]<<" "<<tmp[4]<<" "<<tmp[5]<<" "<<tmp[6]<<" "<<tmp[7]<<"\n";;
// }
__m256 _mm256_complex_mul(Complex* x, Complex* y)//利用AVX指令集同时进行两个复数乘法
{
    __m256 a_vec,b_vec,r;
    a_vec = _mm256_loadu_ps((float*)(x));
    //_mm256_storeu_ps(tmp,a_vec);print();
    b_vec = _mm256_loadu_ps((float*)(y));
    //_mm256_storeu_ps(tmp,b_vec);print();
    r = _mm256_permute_ps(a_vec,245)*_mm256_permute_ps(b_vec,177);
    //_mm256_storeu_ps(tmp,_mm256_permute_ps(a_vec,160));print();
    //_mm256_storeu_ps(tmp,_mm256_permute_ps(b_vec,5));print();
    r = _mm256_fmaddsub_ps(b_vec,_mm256_permute_ps(a_vec,160),r);
    return r;
}
int main()
{
    // for(int i=0;i<5000000;++i){
    //     A[i].Re = 3.0+0.1*i; A[i].Im = 4.0+i;
    //     B[i].Re = -3.0-i; B[i].Im = -12-3*i;
    // }
    // auto start = std::chrono::high_resolution_clock::now();
    // for(int i=0;i<5000000;++i){
    //     C[i] = A[i]*B[i];
    // }
    // auto finish = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = finish -start;
	// std::cout<<1000*elapsed.count()<<"\n";

    // start = std::chrono::high_resolution_clock::now();
    // for(int i=0;i<5000000;i+=4){
    //     __m256 ans = _mm256_complex_mul(A+i,B+i);
    //     _mm256_storeu_ps((float*)(C+i),ans);
    // }
    // finish = std::chrono::high_resolution_clock::now();
	// elapsed = finish -start;
	// std::cout<<1000*elapsed.count()<<"\n";
    // return 0;
    for(int i = 0;i < 4;++i){
        A[i].Re = 3 + i; A[i].Im = -4 + i;
        B[i].Re = 4 + i; B[i].Im = 1 + i;
    }
    __m256 ans = _mm256_complex_mul(A,B);
    _mm256_storeu_ps((float*)(C),ans);
    for(int i = 0; i < 4;++i){
        cout<<C[i].Re<<" "<<C[i].Im<<"\n";
    }
    return 0;
}