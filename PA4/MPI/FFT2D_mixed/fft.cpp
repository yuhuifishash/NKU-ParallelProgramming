//mpic++ MPI/FFT2D_mixed/fft.cpp -o test -g -std=c++11 -fopenmp
#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>
#include<mpi.h>
using namespace std;
const int MAXN=5005;
const float PI=acos(-1.0);
int thread_count = 8;
int pro_size;
int pro_id;
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
Complex a[MAXN][MAXN]={},b[MAXN][MAXN]={},c[MAXN*MAXN]={};
Complex Cal[MAXN*MAXN]={},Cal_rev[MAXN*MAXN]={},Ta[MAXN*MAXN]={},Tb[MAXN*MAXN]={};
Complex W[2][20][MAXN]={};float neg[4]={-1,1,-1,1};
int n=0,m=0,s=0,t=0;
int l=0,r[MAXN]={};
int limit=1;
void V_complex_mul(Complex* x, Complex* y,Complex* c,float* neg)//利用neon指令集同时进行四个复数乘法,结果存储在C中
{
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
void init_root(int limit)
{
    for(int s=0;(1<<s)<limit;++s){
        int mid = (1<<s);
        Complex Wn{cos(PI/mid),sin(PI/mid)};
        Complex w{1,0};
        for(int k=0;k<mid;++k,w=w*Wn){
            W[0][s][k] = w;
        }
    }

    for(int s=0;(1<<s)<limit;++s){
        int mid = (1<<s);
        Complex Wn{cos(PI/mid),-sin(PI/mid)};
        Complex w{1,0};
        for(int k=0;k<mid;++k,w=w*Wn){
            W[1][s][k] = w;
        }
    }
}
void FFT(Complex* A,int type)
{ 
    for(int i=0;i<limit;++i){
        if(i<r[i]){
            swap(A[i],A[r[i]]);
        }
    }
    for(int s=0;(1<<s)<limit;++s){
        int mid = 1<<s;
        if(mid >= 4){
            for(int R=mid<<1,j=0;j<limit;j+=R){
                for(int k=0;k<mid;k+=2){
                    V_complex_mul(A+j+k+mid,W[type == 1 ? 0:1][s]+k,A+j+k+mid,neg);
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
                }   
            }
        }
        else if(mid < limit>>1){
            for(int R=mid<<1,j=0;j<limit;j+=2*R){
                Complex A_t[2],B_t[2];
                for(int k=0;k<mid;++k){
                    A_t[0] = A[j+k];A_t[1] = A[j+k+R];
                    B_t[0] = W[type == 1 ? 0:1][s][k]*A[j+k+mid];
                    B_t[1] = W[type == 1 ? 0:1][s][k]*A[j+k+R+mid];
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
                    A[j+k] = A_t[0];A[j+k+R] = A_t[1];
                    A[j+k+mid] = B_t[0];A[j+k+R+mid] = B_t[1];
                }
            }
        }
        else{
            for(int R=mid<<1,j=0;j<limit;j+=R){
                for(int k=0;k<mid;++k){
                    Complex x=A[j+k],y=W[type == 1 ? 0:1][s][k]*A[j+mid+k];
                    A[j+k]=x+y;
                    A[j+mid+k]=x-y;
                }
            }
        }
    }
    if(type == 1){return;}
    float LIM[4]={(float)limit,(float)limit,(float)limit,(float)limit};
    float* LIMt = (float*)LIM;
    for(int i=0;i<limit;i+=2){
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
    }
}
void init()
{
    //n = 10,m = 10,s = 10,t = 10;
    //n = 2000,m = 2000,s = 2000,t = 2000;
    if(pro_id == 0){
        //cin>>n>>m>>s>>t;
        for(int i=0;i<n;++i){
            for(int j=0;j<m;++j){
                //cin>>a[i][j].Re;
                a[i][j].Re = rand()%10+1;
            }
        }
        for(int i=0;i<s;++i){
            for(int j=0;j<t;++j){
                //cin>>b[i][j].Re;
                b[i][j].Re = rand()%10+1;
            }
        }
    }
}
void FFT_2D_row(Complex *T,int type,int typecl)
{
    int blocks = limit/pro_size;
    if(typecl == 0 && pro_id == 0){
        #pragma omp for
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                Cal[i*limit+j] = T[i*limit+j];
            }
        }
    }
    #pragma omp single
    {
    MPI_Scatter((double*)Cal,limit*blocks,MPI_DOUBLE,(double*)Cal_rev,limit*blocks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    
    #pragma omp for
    for(int i=0;i<blocks;++i){
            FFT(Cal_rev+i*limit,type);
    }
    #pragma omp single
    {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(Cal_rev,limit*blocks,MPI_DOUBLE,Cal,limit*blocks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    if(typecl == 0 && pro_id == 0){
        #pragma omp for
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                T[i*limit+j] = Cal[i*limit+j];  
            }
        }
    }

}
void FFT_2D_col(Complex *T,int type)
{
    if(pro_id == 0){
        #pragma omp for
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                Cal[i*limit+j] = T[j*limit+i];
            }
        }
    }
    FFT_2D_row(T,type,1);
    if(pro_id == 0){
        #pragma omp for
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                T[i*limit+j] = Cal[j*limit+i];
            }
        }
    }

}
void FFT_2D(Complex *T,int type)
{
    FFT_2D_row(T,type,0);
    #pragma omp single
    {
    MPI_Barrier(MPI_COMM_WORLD);
    }

    FFT_2D_col(T,type);

    #pragma omp single
    {
    MPI_Barrier(MPI_COMM_WORLD);
    }
}
void CONV_2D()
{
    while(limit<=max(n+s-2,m+t-2)){
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
    init_root(limit);
    #pragma omp parallel num_threads(thread_count)
    {
        #pragma omp for
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                Ta[i*limit+j] = a[i][j];
                Tb[i*limit+j] = b[i][j];
            }
        }

        FFT_2D(Ta,1);
        FFT_2D(Tb,1);
        if(pro_id == 0){
            #pragma omp for
            for(int i = 0;i<limit;++i){
                for(int j = 0;j<limit;++j){
                    c[i*limit+j] = Ta[i*limit+j]*Tb[i*limit+j];
                }
            }
        }
        FFT_2D(c,-1);
    }
}
int main()
{
    
    int provide;
    MPI_Init_thread(NULL,NULL,MPI_THREAD_SERIALIZED,&provide);
    MPI_Comm_rank(MPI_COMM_WORLD,&pro_id);
    MPI_Comm_size(MPI_COMM_WORLD,&pro_size);

    // double T[20]={0.56,1.1,3.5,13,49,0.27,1.38,5.72};
    // int thr[20]={14,30,62,126,254,510,1022,2046};
    int thr[20]={1,2,4,8,16,32};
    for(int i=0;i<6;++i){
        n = m = s = t = 2000;
        thread_count = thr[i];
        init();
        double start = MPI_Wtime();
        CONV_2D();
        double end = MPI_Wtime();
        if(pro_id==0){cout<<fixed<<setprecision(2)<<end-start<<"s/"<<5.72/(end-start)<<" & ";}
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // if(pro_id == 0){
    //     for(int i = 0;i<=n+s-2;++i){
    //         for(int j = 0;j<=m+t-2;++j){
    //             cout<<(int)(c[i*limit+j].Re+0.5)<<" ";
    //         }
    //         cout<<"\n";
    //     }
    // }
    
    MPI_Finalize();
    return 0;
}