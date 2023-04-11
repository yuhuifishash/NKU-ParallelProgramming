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
Complex A[MAXN],B[MAXN],C[MAXN];
double tmp[MAXN];
// void print()
// {
//     cout<<tmp[0]<<" "<<tmp[1]<<" "<<tmp[2]<<" "<<tmp[3]<<"\n";
// }
__m256d _mm256_complex_mul(Complex* x, Complex* y)//利用AVX指令集同时进行两个复数乘法
{
    __m256d a_vec,b_vec,r;
    a_vec = _mm256_loadu_pd((double*)(x));//a[0] b[0] a[1] b[1] 1 2 2 2
    //_mm256_storeu_pd(tmp,a_vec);print();
    b_vec = _mm256_loadu_pd((double*)(y));//c[0] d[0] c[1] d[1] 3 4 1 2  
    //_mm256_storeu_pd(tmp,b_vec);print();
    //_mm256_storeu_pd(tmp,_mm256_permute_pd(b_vec,5));print();  
    r = _mm256_permute_pd(a_vec,15)*_mm256_permute_pd(b_vec,5);// 8 6 4 2
    //_mm256_storeu_pd(tmp,r);print();
    r = _mm256_fmaddsub_pd(b_vec,_mm256_movedup_pd(a_vec),r);
    //_mm256_storeu_pd(tmp,r);print();
    return r;
}
int main()
{
    for(int i=0;i<5000000;++i){
        A[i].Re = 3.0+0.1*i; A[i].Im = 4.0+i;
        B[i].Re = -3.0-i; B[i].Im = -12-3*i;
    }
    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<5000000;++i){
        C[i] = A[i]*B[i];
    }
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish -start;
	std::cout<<1000*elapsed.count()<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<5000000;i+=2){
        __m256d ans = _mm256_complex_mul(A+i,B+i);
        _mm256_storeu_pd((double*)(C+i),ans);
    }
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()<<"\n";
    return 0;
}