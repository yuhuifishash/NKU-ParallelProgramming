//mpic++ MPI/FFT2D_cache_optimize/fft.cpp -o test -g -std=c++11 -fopenmp
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
int n=0,m=0,s=0,t=0;
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
        for(int R=mid<<1,j=0;j<limit;j+=R){
            Complex w{1,0};
            for(int k=0;k<mid;++k,w=w*Wn){
                Complex x=A[j+k],y=w*A[j+mid+k];
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
        }
    }
    if(type == 1){return;}
    for(int i=0;i<limit;++i){
        A[i].Re=A[i].Re/limit;
        A[i].Im=A[i].Im/limit;
    }
}
void init()
{
    //n = 10,m = 10,s = 10,t = 10;
    n = 2000,m = 2000,s = 2000,t = 2000;
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
void FFT_2D_row(Complex* T,int type,int typecl)
{
    int blocks = limit/pro_size;

    if(typecl == 0 && pro_id == 0){
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                Cal[i*limit+j] = T[i*limit+j];
            }
        }
    }
    MPI_Scatter((double*)Cal,limit*blocks,MPI_DOUBLE,(double*)Cal_rev,limit*blocks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //if(pro_id == 0 && type == -1){for(int i=0;i<limit;++i){for(int j=0;j<limit;++j){cout<<Cal_rev[i*limit+j].Re<<" "<<Cal_rev[i*limit+j].Im<<"j ";}cout<<"\n";}cout<<"\n\n";}
    #pragma omp parallel num_threads(thread_count)
    {
        #pragma omp for
        for(int i=0;i<blocks;++i){
            FFT(Cal_rev+i*limit,type);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(Cal_rev,limit*blocks,MPI_DOUBLE,Cal,limit*blocks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(typecl == 0 && pro_id == 0){
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                T[i*limit+j] = Cal[i*limit+j];  
            }
        }
    }

}
void transform(Complex* x,Complex* y)
{
    for(int i=0;i<limit;i+=8){
        for(int j=0;j<limit;j+=8){
            for(int k=i;k<i+8;++k){
                Complex a_0 = x[k*limit+j];
                Complex a_1 = x[k*limit+j+1];
                Complex a_2 = x[k*limit+j+2];
                Complex a_3 = x[k*limit+j+3];
                Complex a_4 = x[k*limit+j+4];
                Complex a_5 = x[k*limit+j+5];
                Complex a_6 = x[k*limit+j+6];
                Complex a_7 = x[k*limit+j+7];
                y[j*limit+k] = a_0;
                y[(j+1)*limit+k] = a_1;
                y[(j+2)*limit+k] = a_2;
                y[(j+3)*limit+k] = a_3;
                y[(j+4)*limit+k] = a_4;
                y[(j+5)*limit+k] = a_5;
                y[(j+6)*limit+k] = a_6;
                y[(j+7)*limit+k] = a_7;
            }
        }
    }
}
void FFT_2D_col(Complex* T,int type)
{
    if(pro_id == 0){
        transform(T,Cal);
    }
    FFT_2D_row(T,type,1);
    if(pro_id == 0){
        transform(Cal,T);
    }

}
void FFT_2D(Complex* T,int type)
{
    FFT_2D_row(T,type,0);
    MPI_Barrier(MPI_COMM_WORLD);

    FFT_2D_col(T,type);
    MPI_Barrier(MPI_COMM_WORLD);

    //if(pro_id == 0){for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<A[i][j].Re<<" "<<A[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n\n";}
}
void CONV_2D()
{
    while(limit<=max(n+s-2,m+t-2)){
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }
    for(int i=0;i<limit;++i){
        for(int j=0;j<limit;++j){
            Ta[i*limit+j] = a[i][j];
            Tb[i*limit+j] = b[i][j];
        }
    }
    FFT_2D(Ta,1);
    FFT_2D(Tb,1);
    //if(pro_id == 0){for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<b[i][j].Re<<" "<<b[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";}
    if(pro_id == 0){
        for(int i = 0;i<limit;++i){
            for(int j = 0;j<limit;++j){
                c[i*limit+j] = Ta[i*limit+j]*Tb[i*limit+j];
            }
        }
    }
    //if(pro_id == 0){for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<c[i][j].Re<<" "<<c[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";}
    FFT_2D(c,-1);
}
int main()
{
    
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&pro_id);
    MPI_Comm_size(MPI_COMM_WORLD,&pro_size);

    init();
    double start = MPI_Wtime();
    CONV_2D();
    double end = MPI_Wtime();

    if(pro_id==0){cout<<fixed<<setprecision(6)<<end-start<<"\n";}
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