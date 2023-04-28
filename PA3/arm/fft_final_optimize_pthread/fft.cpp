//g++ -march=armv8-a fft_final_optimize_pthread/fft.cpp -o test -g -std=c++11 -lpthread
#include<stdio.h>
#include<cmath>
#include<iomanip>
#include<cstring>
#include<arm_neon.h>
#include<semaphore.h>
#include<pthread.h>
#include<chrono>
using namespace std;
const int MAXN=20000005;
const float PI=acos(-1.0);
const int thread_count = 4;
pthread_barrier_t barr_merge;
pthread_barrier_t barr_root;
sem_t sem_expand;
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
struct fft_thread
{
    int id;
    int type;
    Complex* arr;
};
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={},W[MAXN]={};
int N=0,M=0;
int l=0,r[MAXN]={};
int limit=1;
float neg[4]={-1,1,-1,1};
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
void* FFT(void* rank)
{
    fft_thread* t_id = (fft_thread*)rank;
    int id = t_id -> id;
    int type = t_id -> type;
    Complex* A = t_id -> arr;

    if(id == 0){
        for(int i=0;i<limit;++i){
            if(i<r[i]){
                swap(A[i],A[r[i]]);
            }
        }
    }
    else{
        sem_wait(&sem_expand);
    }
    if(id == 0){
        for(int i=0;i<thread_count-1;++i){
            sem_post(&sem_expand);
        }
    }

    for(int mid=1;mid<limit;mid<<=1){
        Complex Wn{cos(PI/mid),type*sin(PI/mid)};
        if(mid >= 4){
            int t1 = limit/(mid<<1),count1,count2;
            if(t1 < thread_count){
                count1 = t1;
                count2 = thread_count/t1;
            }
            else{t1 = thread_count;}

            if(id == 0){
                Complex Wn_2 = Wn*Wn;
                Complex Mul_epi[2]={Wn_2,Wn_2};
                W[0] = Complex{1,0};
                W[1] = Complex{1,0}*Wn;
                for(int k=2;k<mid;k+=2){
                    V_complex_mul(W+k-2,Mul_epi,W+k,neg);
                }
            }
            pthread_barrier_wait(&barr_root);
            
            if(t1 == thread_count){
                int blocks = max(1,limit/(mid<<1)/thread_count);
                int begin = min(limit,id*blocks*(mid<<1)),end = min(limit,(id+1)*blocks*(mid<<1));

                for(int R=mid<<1,j=begin;j<end;j+=R){
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
                    }   
                }
            }
            else{
                int b1 = max(1,limit/(mid<<1)/count1),b2 = max(1,mid/2/count2);

                int begin1 = min(limit,id/count2*b1*(mid<<1)),end1 = min(limit,(id/count2+1)*b1*(mid<<1));

                int begin2 = min(mid,id%count2*b2*2),end2 = min(mid,(id%count2+1)*b2*2);

                for(int R=mid<<1,j=begin1;j<end1;j+=R){
                    for(int k=begin2;k<end2;k+=2){
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
                    }   
                }
            }
        }
        else if(mid < limit>>1){
            int blocks = max(1,limit/2/(mid<<1)/thread_count);
            int begin = min(limit,id*blocks*(mid<<1)*2),end = min(limit,(id+1)*blocks*(mid<<1)*2);
            for(int R=mid<<1,j=begin;j<end;j+=2*R){
                Complex w{1,0};
                Complex A_t[2],B_t[2];
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
                    A[j+k] = A_t[0];A[j+k+R] = A_t[1];
                    A[j+k+mid] = B_t[0];A[j+k+R+mid] = B_t[1];
                }
            }
        }
        else{
            int blocks = max(1,limit/(mid<<1)/thread_count);
            int begin = min(limit,id*blocks*(mid<<1)),end = min(limit,(id+1)*blocks*(mid<<1));
            for(int R=mid<<1,j=begin;j<end;j+=R){
                Complex w{1,0};
                for(int k=0;k<mid;++k,w=w*Wn){
                    Complex x=A[j+k],y=w*A[j+mid+k];
                    A[j+k]=x+y;
                    A[j+mid+k]=x-y;
                }
            }
        }
        pthread_barrier_wait(&barr_merge);
    }
    if(type == 1){return NULL;}

    float LIM[4]={(float)limit,(float)limit,(float)limit,(float)limit};
    float* LIMt = (float*)LIM;
    int blocks = max(1,limit/2/thread_count);
    int begin = min(limit,id*2*blocks) , end = min(limit,(id+1)*2*blocks);
    for(int i=begin;i<end;i+=2){
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
    return NULL;
}
void Conv(Complex* A,Complex* B,Complex* ans)
{
    //N=20;M=20;
    N=5000000,M=5000000;
    for(int i=0;i<N;++i){
        A[i].Re=rand()%10+1;
        //cin>>A[i].Re;
    }
    //输入b为卷积核的翻转。
    for(int i=0;i<M;++i){
        B[i].Re=rand()%10+1;
        //cin>>B[i].Re;
    }
    pthread_t* thread_handles = new pthread_t[thread_count];
    fft_thread* id = new fft_thread[thread_count];

    auto start = std::chrono::high_resolution_clock::now();

    while(limit<=N+M-2){limit<<=1;++l;}
    for(int i=0;i<limit;++i){r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));}
    for(int i = 0 ; i < thread_count ; ++i){//FFT(A,1)
        id[i].arr=A ; id[i].id=i ; id[i].type=1;
        pthread_create(&thread_handles[i],NULL,FFT,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }

    for(int i = 0 ; i < thread_count ; ++i){//FFT(B,1)
        id[i].arr=B ; id[i].id=i ; id[i].type=1;
        pthread_create(&thread_handles[i],NULL,FFT,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }
    
    for(int i=0;i<=limit;i+=2){
        V_complex_mul(A+i,B+i,ans+i,neg);
    }

    for(int i = 0 ; i < thread_count ; ++i){//FFT(ans,-1)
        id[i].arr=ans ; id[i].id=i ; id[i].type=-1;
        pthread_create(&thread_handles[i],NULL,FFT,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }

    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
    printf("%f\n",elapsed.count());

    delete []thread_handles;
    delete []id;
}
int main()
{

    sem_init(&sem_expand,0,0);
    pthread_barrier_init(&barr_merge,NULL,thread_count);
    pthread_barrier_init(&barr_root,NULL,thread_count);

    Conv(a,b,c);
    for(int i=0;i<N+M-1;++i){
        //printf("%d ",(int)(c[i].Re+0.5));
    }//cout<<"\n";

    sem_destroy(&sem_expand);
    pthread_barrier_destroy(&barr_merge);
    pthread_barrier_destroy(&barr_root);
    return 0;
}