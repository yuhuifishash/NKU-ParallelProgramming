//mpic++ MPI/FFT2D_1/fft.cpp -o test -g -std=c++11
#include<iostream>
#include<cmath>
#include<chrono>
#include<iomanip>
#include<mpi.h>
using namespace std;
const int MAXN=5005;
const float PI=acos(-1.0);
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
Complex a[MAXN][MAXN]={},b[MAXN][MAXN]={},c[MAXN][MAXN]={};
Complex Tran[MAXN][MAXN]={};//矩阵转置
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
    for(int i=0;i<=limit;++i){
        A[i].Re=A[i].Re/limit;
        A[i].Im=A[i].Im/limit;
    }
}
void init()
{
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
void FFT_2D_row(Complex A[][MAXN],int n,int m,int type)
{
    int blocks = limit/pro_size;
    if(pro_id == 0){//传递数据
        for(int dest=1;dest<pro_size;++dest){
            int begin = dest*blocks,end = (dest+1)*blocks;
            for(int i=begin;i<end;++i){
                MPI_Send(A[i],2*limit,MPI_DOUBLE,dest,i,MPI_COMM_WORLD);
            }
        }
    }
    else{
        int begin = pro_id*blocks,end = (pro_id+1)*blocks;
        for(int i=begin;i<end;++i){
            MPI_Recv(A[i],2*limit,MPI_DOUBLE,0,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    //if(pro_id == 1){for(int i=0;i<limit;++i){for(int j=0;j<limit;++j){cout<<A[i][j].Re<<" ";}cout<<"\n";}cout<<"\n\n";}
    int begin = pro_id*blocks,end = (pro_id+1)*blocks;
    for(int i=begin;i<end;++i){
        FFT(A[i],type);
    }

    if(pro_id != 0){
        int begin = pro_id*blocks,end = (pro_id+1)*blocks;
        for(int i=begin;i<end;++i){
            MPI_Send(A[i],2*limit,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
        }
    }
    else{
        for(int source=1;source<pro_size;++source){
            int begin = source*blocks,end = (source+1)*blocks;
            for(int i=begin;i<end;++i){
                MPI_Recv(A[i],2*limit,MPI_DOUBLE,source,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
    }

}
void FFT_2D_col(Complex A[][MAXN],int n,int m,int type)
{
    if(pro_id == 0){
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                Tran[i][j] = A[j][i];
            }
        }
    }
    FFT_2D_row(Tran,limit,limit,type);
    if(pro_id == 0){
        for(int i=0;i<limit;++i){
            for(int j=0;j<limit;++j){
                A[i][j] = Tran[j][i];
            }
        }
    }

}
void FFT_2D(Complex A[][MAXN],int n,int m,int type)
{
    FFT_2D_row(A,limit,limit,type);
    MPI_Barrier(MPI_COMM_WORLD);
    //if(pro_id == 0){for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<a[i][j].Re<<" "<<a[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n\n";}

    FFT_2D_col(A,limit,limit,type);
    MPI_Barrier(MPI_COMM_WORLD);
}
void CONV_2D()
{
    while(limit<=max(n+s-2,m+t-2)){
        limit<<=1;++l;
    }
    for(int i=0;i<limit;++i){
        r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    }

    FFT_2D(a,limit,limit,1);
    //if(pro_id == 0){for(int i=0;i<n;++i){for(int j=0;j<m;++j){cout<<fixed<<setprecision(2)<<a[i][j].Re<<" "<<a[i][j].Im<<"j ";}cout<<"\n";}cout<<"\n";}
    FFT_2D(b,limit,limit,1);
    if(pro_id == 0){
        for(int i = 0;i<=limit;++i){
            for(int j = 0;j<=limit;++j){
                c[i][j] = a[i][j]*b[i][j];
            }
        }
    }
    FFT_2D(c,limit,limit,-1);
}
int main()
{
    

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&pro_id);
    MPI_Comm_size(MPI_COMM_WORLD,&pro_size);


    double start = MPI_Wtime();
    init();
    CONV_2D();
    double end = MPI_Wtime();
    if(pro_id==0){cout<<fixed<<setprecision(6)<<end-start<<"\n";}
    // if(pro_id == 0){
    //     for(int i = 0;i<=n+s-2;++i){
    //         for(int j = 0;j<=m+t-2;++j){
    //             cout<<(int)(c[i][j].Re+0.5)<<" ";
    //         }
    //         cout<<"\n";
    //     }
    // }
    

    MPI_Finalize();
    return 0;
}