#include<iostream>


using namespace std;
int main()
{
    const long MAX=2147483647;
    long  reverse(long);
    bool judge(long);
    long n,s;
    int count=0;
    cout<<"请输入一个大于10的正整数:";
    cin>>n;
    if(n>=MAX||n<10)
    {
        cout<<"输入错误，请重新输入";
    }
    cout<<"获得的回文数的过程如下："<<endl;
    for(;;)
    {
        count++;
        s=n+reverse(n);
        cout<<"["<<count<<"]:"<<n<<"+"<<reverse(n)<<"="<<s<<endl;
        if(judge(s))
        {
            cout<<"获得的回文数是："<<s<<endl;
            break;
        }
        else
        {
            n=s;
        }
        
    }
    return 0;
}
long reverse(long n)
{
    long r;
    for(r=0;n>0;n/=10)
    {
        r=r*10+n%10;
    }
    return r;
}
bool judge(long s)
{
    bool success=false;
    if(reverse(s)==s) success=true;
    return success;
}