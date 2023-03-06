#include<iostream>
#include<sys/time.h>
using namespace std;
int n = 0;
int M[10005][10005]={},a[10005]={},ans[10005]={};
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
    for(int i=1;i<=n;++i){
        a[i] = 2;
        for(int j=1;j<=n;++j){
            M[i][j] = 1;
        }
    }
}

void dot_col()//列访问
{
    for(int i=1;i<=n;++i){
        ans[i] = 0;
        for(int j=1;j<=n;++j){
            ans[i] += a[j]*M[j][i];
        }
    }
}
void dot_row()//行访问
{
    for(int i=1;i<=n;++i){ans[i] = 0;}
    for(int j=1;j<=n;++j){
        for(int i=1;i<=n;++i){
            ans[i] += a[j]*M[j][i];
        }
    }
}
void dot_row2()//行访问+局部变量优化算法
{
    for(int i=1;i<=n;++i){ans[i] = 0;}
    for(int j=1;j<=n;++j){
        int t = a[j];
        for(int i=1;i<=n;++i){
            ans[i] += t*M[j][i];
        }
    }
}
int main()
{
    int cnt = 0;
    double st = 0, ed = 0;
    cin>>n>>cnt;
    init();

    st = get_time();
    for(int i = 1;i <= cnt;++i){dot_col();}      
    ed = get_time();
    cout<<(ed-st)/1000.0/(double)cnt<<"\n";

    st = get_time();
    for(int i = 1;i <= cnt;++i){dot_row();}      
    ed = get_time();
    cout<<(ed-st)/1000.0/(double)cnt<<"\n";

    st = get_time();
    for(int i = 1;i <= cnt;++i){dot_row2();}      
    ed = get_time();
    cout<<(ed-st)/1000.0/(double)cnt<<"\n";
}
