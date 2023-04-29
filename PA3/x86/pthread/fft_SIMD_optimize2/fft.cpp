#include<stdio.h>
#include<cmath>
#include<iomanip>
#include<cstring>
#include<immintrin.h>
#include<semaphore.h>
#include<pthread.h>
#include<chrono>
using namespace std;
const int MAXN=20000005;
const float PI=acos(-1.0);
const int thread_count = 8;
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
__m256 _mm256_complex_mul(Complex* x, Complex* y)//利用AVX指令集同时进行四个复数乘法
{
    __m256 a_vec,b_vec,r;
    a_vec = _mm256_loadu_ps((float*)(x));
    b_vec = _mm256_loadu_ps((float*)(y));   
    r = _mm256_permute_ps(a_vec,245)*_mm256_permute_ps(b_vec,177);
    r = _mm256_fmaddsub_ps(b_vec,_mm256_permute_ps(a_vec,160),r);
    return r;
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
        if(mid >= 8){
            if(id == 0){
                Complex Wn_4 = Wn*Wn*Wn*Wn;
                Complex Mul_epi[4]={Wn_4,Wn_4,Wn_4,Wn_4};
                W[0] = Complex{1,0};W[1] = W[0]*Wn;
                W[2] = W[1]*Wn;W[3] = W[2]*Wn;
                for(int k=4;k<mid;k+=4){
                    __m256 t = _mm256_complex_mul(W+k-4,Mul_epi);
                    _mm256_storeu_ps((float*)(W+k),t);
                }
            }
            pthread_barrier_wait(&barr_root);
            for(int R=mid<<1,j=id*(mid<<1);j<limit;j+=thread_count*R){
                for(int k=0;k<mid;k+=4){
                    __m256 t1 = _mm256_loadu_ps((float*)(A+j+k));
                    __m256 t2 = _mm256_complex_mul(A+j+k+mid,W+k);
                    __m256 t3 = _mm256_add_ps(t1,t2);
                    t2 = _mm256_sub_ps(t1,t2);
                    _mm256_storeu_ps((float*)(A+j+k),t3);
                    _mm256_storeu_ps((float*)(A+j+k+mid),t2);
                }   
            }
        }
        else if(mid < limit>>2){
            for(int R=mid<<1,j=0;j<limit;j+=4*thread_count*R){
                Complex w{1,0};
                Complex A_t[4],B_t[4];
                __m256 t1,t2,t3;
                for(int k=0;k<mid;++k,w=w*Wn){
                    A_t[0] = A[j+k+4*id*R];A_t[1] = A[j+k+R+4*id*R];
                    A_t[2] = A[j+k+2*R+4*id*R];A_t[3] = A[j+k+3*R+4*id*R];
                    B_t[0] = A[j+k+mid+4*id*R];B_t[1] = A[j+k+R+mid+4*id*R];
                    B_t[2] = A[j+k+2*R+mid+4*id*R];B_t[3] = A[j+k+3*R+mid+4*id*R];
                    Complex W_t[4]={w,w,w,w};
                    t1 = _mm256_loadu_ps((float*)A_t);
                    t2 = _mm256_complex_mul(W_t,B_t);
                    t3 = _mm256_add_ps(t1,t2);
                    t2 = _mm256_sub_ps(t1,t2);
                    _mm256_storeu_ps((float*)A_t,t3);
                    _mm256_storeu_ps((float*)B_t,t2);
                    A[j+k+4*id*R] = A_t[0];A[j+k+R+4*id*R] = A_t[1];
                    A[j+k+2*R+4*id*R] = A_t[2];A[j+k+3*R+4*id*R] = A_t[3];
                    A[j+k+mid+4*id*R] = B_t[0];A[j+k+R+mid+4*id*R] = B_t[1];
                    A[j+k+2*R+mid+4*id*R] = B_t[2];A[j+k+3*R+mid+4*id*R] = B_t[3];
                }
            }
        }
        else{
            for(int R=mid<<1,j=id*(mid<<1);j<limit;j+=thread_count*R){
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

    __m256 t2 = _mm256_set1_ps(limit);
    for(int i=0;i<limit;i+=4*thread_count){
        __m256 t1 = _mm256_loadu_ps((float*)(A+i+4*id));
        t1 = _mm256_div_ps(t1,t2);
        _mm256_storeu_ps((float*)(A+i+4*id),t1);
    }
    return NULL;
}
void Conv(Complex* A,Complex* B,Complex* ans)
{
    scanf("%d%d",&N,&M);
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
    
    for(int i=0;i<=limit;i+=4){
        __m256 t = _mm256_complex_mul(A+i,B+i);
        _mm256_storeu_ps((float*)(ans+i),t);
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