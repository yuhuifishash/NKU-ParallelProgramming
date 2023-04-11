#include<iostream>
#include<chrono>
#include<sys/time.h>
using namespace std;
const int MAXN = 3e7+5;
int a[MAXN]={},n=0;
void init()
{
    for(int i = 1;i <= n; ++i){
        a[i] = 1;
    }
}
void sum_list()
{
    int sum = 0;
    for(int i=1 ; i<=n;++i){
        sum+=a[i];
    }
    int ans = sum;
}
void sum_2()
{
    int sum1 = 0, sum2 = 0;
    for(int i=1 ;i<=n ;i+=2){
        sum1+=a[i];
        sum2+=a[i+1];
    }
    int ans = sum1 + sum2;
}
void sum_4()
{
    int sum[4]={};
    for(int i=1 ;i<=n;i+=4){
        sum[0]+=a[i];
        sum[1]+=a[i+1];
        sum[2]+=a[i+2];
        sum[3]+=a[i+3];
    }
    int ans = sum[0]+sum[1]+sum[2]+sum[3];
}

void sum_8()
{
    int sum[8]={};
    for(int i=1 ;i<=n;i+=8){
        sum[0]+=a[i];
        sum[1]+=a[i+1];
        sum[2]+=a[i+2];
        sum[3]+=a[i+3];
        sum[4]+=a[i+4];
        sum[5]+=a[i+5];
        sum[6]+=a[i+6];
        sum[7]+=a[i+7];
    }
    int ans = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7];
}

void sum_16()
{
    int sum[16]={};
    for(int i=1 ;i<=n;i+=16){
        sum[0]+=a[i];
        sum[1]+=a[i+1];
        sum[2]+=a[i+2];
        sum[3]+=a[i+3];
        sum[4]+=a[i+4];
        sum[5]+=a[i+5];
        sum[6]+=a[i+6];
        sum[7]+=a[i+7];
        sum[8]+=a[i+8];
        sum[9]+=a[i+9];
        sum[10]+=a[i+10];
        sum[11]+=a[i+11];
        sum[12]+=a[i+12];
        sum[13]+=a[i+13];
        sum[14]+=a[i+14];
        sum[15]+=a[i+15];
    }
    int ans = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]
    +sum[8]+sum[9]+sum[10]+sum[11]+sum[12]+sum[13]+sum[14]+sum[15];
}
void sum_divide()
{
    for(int i = n; i > 1 ; i>>=1){
        for(int j = 1; j <= i/2 ; ++j){
            a[j] = a[2*j-1]+a[2*j];
        }
    }
    int ans = a[1];
}
int main()
{
    int cnt = 100;
    double st = 0, ed = 0;
    cin>>n>>cnt;
    init();

    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_list();}      
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_2();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_4();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_8();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_16();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){sum_divide();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    return 0;
}
