#include<iostream>
#include<string>
#include "stdlib.h"
#define Max 100000
char num1[Max];   //num1存放原数
char num2[Max];   //num2存放逆序数
char sum[Max];   //sum存放和的数组
char sumout[Max]; //对应于sum的输出格式转换数组
void inttochar(int number); //输入的整形转换成字符型函数
bool Cmp();    //回文比较函数
void add();  //两数相加函数
void adjust(); //数组调整hanshu 
int length();  //计算长度函数
int times=0;  

using namespace std;
int main()
{
    bool equal=false;
    memset(num1,'\0',sizeof(num1));   //使数组1d的元素初始化为0
    memset(num2,'\0',sizeof(num2));
    memset(sum,'\0',sizeof(sum));
    int number;
    cout<<"number:";
    cin>>number;
    inttochar(number);
    int Length=length();     //求数组num1中有多少个‘\0’
    int tmp=0;
    for(int i=Max-1;i>=Length;i--)    //将num2与num1的数字位置相反及为逆序数
    {
        num2[Length+tmp]=num1[i];
        tmp++;
    }
    while(false==equal)       //判断回文数的循环
    {
        ++times;
        add();
        adjust();
        cout<<"["<<times<<"]:";
        cout<<sumout<<endl;
        equal=Cmp();
    }
    cout<<"total times is "<<times<<endl;
}
    void inttochar(int number)
    {
        char tmp[100],Length;
        memset(tmp,'\0',sizeof(tmp));
        itoa(number,tmp,10);   //整型数字转换为字符型数组，10代表10进制，包含于stdlib.h中
        for(int i=0;i<100;i++)
        {
            if(tmp[i]=='\0')
            {
                Length=i;
                break;
            }
        }
        for(int i=0;i<Length;i++)    //把num1初始化为“\0\0\0\0.....\0196”的形式
        {
            num1[Max-Length+i]=tmp[i];
        }
    }
    bool Cmp()
    {
        int Length,tmp,tmp1;
        Length=length();
        tmp1=(Max-Length)/2;
        tmp=0;
        for(int i=Length;i<Max;i++)
        {
            if(sum[i]!=sum[Max-1-tmp])
            {
                return false;
            }
            tmp++;
            if(tmp>tmp1)
            {
                return true;
            }
        }
        return true;
    }
    void add()
    {
        char tmpsum;
        int Length,i,cf;    //cf是进位
        cf=0;
        Length=length();
        for(i=Max-1;i>+Length;i--)
        {
            tmpsum=num1[i]+num2[i]+cf-48;
            if(tmpsum>57)
            {
                tmpsum-=10;
                cf=1;
            }
            else
            {
                {
                    cf=0;
                }
            }
            sum[i]=tmpsum;
            
        }
        if(cf==1)
        {
            sum[i]='1';
        }

    }
    void adjust()
    {
        int Length,tmp;
        tmp=0;
        for(int i=0;i<Max;i++)
        {
            num1[i]=sum[i];

        }
        Length=length();
        for(int i=Max-1;i>=Length;i--)
        {
            num2[Length+tmp]=num1[i];
            tmp++;
        }
        for(int i=Max-1;i>=Length;i--)
        {
            sumout[i-Length]=num1[i];
        }
    }
    int length()
    {
        int tmp;
        for(int i=0;i<Max;i++)
        {
            if(num1[i]!='\0')
            {
                tmp=i;
                break;
            }
        }
        if(0==tmp)
        {
            cout<<"The number is so long,The longest length is  "<<Max<<endl;
            cout<<"Total times is  "<<--times<<endl;
            exit(0);
        }
        return tmp;
    }
