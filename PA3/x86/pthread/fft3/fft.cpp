#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<semaphore.h>
#include<iostream>
#include<chrono>
#include<cmath>
using namespace std;
const int MAXN=2e7+5;
const float PI=acos(-1.0);
const int thread_count = 8;
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
struct thread_id
{
    int id;
    int type;
    Complex* arr;
};
Complex a[MAXN]={},b[MAXN]={},c[MAXN]={},W[MAXN]={};
int n=0,m=0,l=0,r[MAXN]={},limit=1;
void init()
{
    scanf("%d%d",&n,&m);
    for(int i=0;i<n;++i){
        a[i].Re=rand()%10+1;
        //scanf("%lf",&a[i].Re);
    }
    for(int i=0;i<m;++i){
        b[i].Re=rand()%10+1;
        //scanf("%lf",&b[i].Re);
    }
    while(limit<=n+m-2){
        limit<<=1;
        ++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
}
pthread_barrier_t barr_merge;
pthread_barrier_t barr_fft;
sem_t sem_expand;
void* fft_dft(void* rank)
{
    thread_id* t_id = (thread_id*)rank;
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
        for(int i = 0 ; i < thread_count-1 ;++i){
            sem_post(&sem_expand);
        }
    }
    for(int mid=1;mid<limit;mid<<=1){

        if( limit/(mid<<1) > thread_count){
            Complex Wn{cos(PI/mid),type*sin(PI/mid)};
            for(int j=id*(mid<<1);j<limit;j+=thread_count*(mid<<1)){
                Complex w{1,0}; 
                for(int k=0;k<mid;++k,w=w*Wn){
                    Complex x=A[j+k],y=w*A[j+mid+k];
                    A[j+k]=x+y;
                    A[j+mid+k]=x-y;
                }
            } 
        }
        else{
            if(id == 0){
                Complex Wn{cos(PI/mid),type*sin(PI/mid)};
                W[0] = Complex{1,0};
                for(int k=1;k<mid;++k){
                    W[k] = W[k-1]*Wn;
                }
            }
            pthread_barrier_wait(&barr_fft);

            for(int j=0;j<limit;j+=(mid<<1)){
                for(int k=id;k<mid;k+=thread_count){
                    Complex x=A[j+k],y=W[k]*A[j+mid+k];
                    A[j+k]=x+y;
                    A[j+mid+k]=x-y;
                }
            }
        }

        pthread_barrier_wait(&barr_merge);
    }
    if(type == -1){
        for(int i=id;i<=limit;i+=thread_count){
            A[i].Re=(int)(A[i].Re/limit+0.5);
        }
    }
    return NULL;
}
int main()
{
    init();
    
    auto start = std::chrono::high_resolution_clock::now();

    sem_init(&sem_expand,0,0);
    pthread_barrier_init(&barr_merge,NULL,thread_count);
    pthread_barrier_init(&barr_fft,NULL,thread_count);

    pthread_t* thread_handles = new pthread_t[thread_count];
    thread_id* id = new thread_id[thread_count];
    for(int i = 0 ; i < thread_count ; ++i){
        id[i].arr=a ; id[i].id=i ; id[i].type=1;
        pthread_create(&thread_handles[i],NULL,fft_dft,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }

    for(int i = 0 ; i < thread_count ; ++i){
        id[i].arr=b ; id[i].id=i ; id[i].type=1;
        pthread_create(&thread_handles[i],NULL,fft_dft,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }

    for(int i=0;i<=limit;++i){
        c[i]=a[i]*b[i];
    }

    for(int i = 0 ; i < thread_count ; ++i){
        id[i].arr=c ; id[i].id=i ; id[i].type=-1;
        pthread_create(&thread_handles[i],NULL,fft_dft,(void*)(id+i));
    }
    for(int i = 0 ; i < thread_count ; ++i){
        pthread_join(thread_handles[i],NULL);
    }

    sem_destroy(&sem_expand);
    pthread_barrier_destroy(&barr_merge);
    pthread_barrier_destroy(&barr_fft);
    delete []thread_handles;
    delete []id;

    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = finish - start;
    printf("%f\n",elapsed.count());
    for(int i=0;i<n+m-1;++i){
        //printf("%d ",(int)c[i].Re);
    }
    return 0;
}