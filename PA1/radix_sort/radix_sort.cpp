#include<iostream>
#include<random>
#include<cmath>
#include<sys/time.h>
using namespace std;
typedef unsigned int uint;
int n = 0,d = 0,k = 0,s = 0;
int cnt[12][65540]={};
uint a[200000005],b[200000005];
double get_time()
{
    struct timeval tv;
    double timenow;
    gettimeofday(&tv,NULL);
    timenow = (double)tv.tv_sec*1000000;
    timenow +=tv.tv_usec;
    return timenow;
}
void init()
{
    for(int i = 0;i < n;++i){
        a[i]=rand()*rand();
    }
}
void sort()
{
    for(int t = 0;t < k;++t){
        for(int i = 0;i < d;++i){
            cnt[t][i] = 0;
        }
    }
    for(int i = 0;i < n;++i){
        uint val = a[i];
        for(int t = 0;t < k;++t){
            ++cnt[t][(val>>(t*s))&(d-1)];
        }
    }
    for(int t = 0;t < k;++t){ 
        for(int i = 0;i < d;++i){
            cnt[t][i] += cnt[t][i-1];
        }
    }
    for(int t = 0;t < k;++t){
        if(t&1){
            for(int i = n-1 ;i >= 0; --i){
                a[--cnt[t][(b[i]>>(t*s))&(d-1)]]=b[i];
            }
        }
        if(!(t&1)){
            for(int i = n-1 ;i >= 0; --i){
                b[--cnt[t][(a[i]>>(t*s))&(d-1)]]=a[i];
            }
        }
    }
}
int main()
{
    srand(time(0));
    int cnt = 0;
    double st = 0, ed = 0;
    cin>>n>>d>>cnt;
    switch(d){
        case 16 : k = 8,s = 4;break;
        case 32 : k = 7,s = 5;break;
        case 64 : k = 6,s = 6;break;
        case 128 : k = 5,s = 7;break;
        case 256 : k = 4,s = 8;break;
        case 1024 : k = 4,s = 10;break;
        case 2048 : k = 3,s = 11;break;
        case 4096 : k = 3,s = 12;break;
        case 8192 : k = 3,s = 13;break;
        case 32768 : k = 3,s = 15;break;
        case 65536 : k = 2,s = 16;break;
    } 
    init();

    st = get_time();
    for(int i = 1;i <= cnt;++i){sort();}      
    ed = get_time();
    cout<<(ed-st)/1000.0/(double)cnt<<"\n";
    return 0;
}
