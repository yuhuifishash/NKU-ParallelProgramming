#include<iostream>
#include<sys/time.h>
using namespace std;
typedef unsigned int uint;
int n = 0;
uint a[200000005],b[200000005];
int cnt0[256],cnt1[256],cnt2[256],cnt3[256];
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
        a[i] = 5e8 - i;
        if(i&1){a[i]+=7*i+1e8;}
    }
}
void sort()
{
    for(int i = 0;i <= 255;++i){
        cnt0[i] = 0; cnt1[i] = 0;cnt2[i] = 0; cnt3[i] = 0;
    }
    for(int i = 0;i < n;++i){
        uint t = a[i];
        ++cnt0[t&255];++cnt1[(t>>8)&255];
		    ++cnt2[(t>>16)&255];++cnt3[t>>24];
    }
    for(int i = 1;i <= 255;++i){
        cnt0[i]+=cnt0[i-1];cnt1[i]+=cnt1[i-1];
        cnt2[i]+=cnt2[i-1];cnt3[i]+=cnt3[i-1];
    }
    for(int i = n-1;i >= 0;--i){b[--cnt0[a[i]&255]]=a[i];}
    for(int i = n-1;i >= 0;--i){a[--cnt1[(b[i]>>8)&255]]=b[i];}
    for(int i = n-1;i >= 0;--i){b[--cnt2[(a[i]>>16)&255]]=a[i];}
    for(int i = n-1;i >= 0;--i){a[--cnt3[b[i]>>24]]=b[i];}
}
int main()
{
    int cnt = 0;
    double st = 0, ed = 0;
    cin>>n>>cnt;
    init();

    st = get_time();
    for(int i = 1;i <= cnt;++i){sort();}      
    ed = get_time();
    cout<<(ed-st)/1000.0/(double)cnt<<"\n";
    return 0;
}
