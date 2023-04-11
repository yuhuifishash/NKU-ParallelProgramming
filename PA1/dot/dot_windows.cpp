#include<iostream>
#include<chrono>
#include<sys/time.h>
using namespace std;
int n = 0;
int M[10005][10005]={},a[10005]={},ans[10005]={};
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
    int cnt = 1;
    double st = 0, ed = 0;
    cin>>n;
    init();

    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){dot_col();}      
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){dot_row();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";

    start = std::chrono::high_resolution_clock::now();
    for(int i = 1;i <= cnt;++i){dot_row2();}      
    finish = std::chrono::high_resolution_clock::now();
	elapsed = finish -start;
	std::cout<<1000*elapsed.count()/cnt<<"\n";
}
