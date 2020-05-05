#include<iostream>
#include<string>

using namespace std;
void DtoB(double);   //将输入的十进制的整数部分和小数部分分开
void DtoBI(int);      //将整数部分转换为十进制
void DtoBF(double);   //将小数部分转换为十进制
char bin[50];    //十进制整数数组
char binf[50];  //十进制小数数组
int main()

{
    double num,fnum;
    int inum;
    cout<<endl;
    cout<<"**************十进制向二进制转换器*************"<<endl;
    cout<<"请输入一个十进制数:";
    cin>>num;
    cout<<endl;
    while(num>10000000)
    {
        cout<<"数字太大，请重新输入"<<endl;
        cin>>num;
    }
    while(num<0)
    {
        cout<<"请输入正数:";
        cin>>num;
    }
    DtoB(num);
}
void DtoB(double num)
{
    int inum;
    double fnum;
    inum=int(num);
    fnum=num-inum;
    DtoBI(inum);
    if(fnum==0)
    {
        cout<<"十进制数"<<num<<"转换成二进制为:"<<bin<<endl;
    }
    else
    {
        DtoBF(fnum);
        cout<<"十进制数"<<num<<"转换成二进制为："<<bin<<"."<<binf<<endl;
    }
    
}
void DtoBI(int inum)
{
    int str[50];
    memset(bin,'\0',sizeof(bin));
    memset(str,'\0',sizeof(str));
    int i=0;
    while(1){
        if(inum==0)break;
        str[i]=inum%2;
        inum/=2;
        i++;
    }
    for(int j=0;j<i+1;j++)
    {
        bin[i-j-1]=str[j]+48;

    }
}
void DtoBF(double fnum)
{
    int str[50];
    memset(binf,'\0',sizeof(binf));
    memset(str,'\0',sizeof(bin));
    int i=0;
    while(fnum)
    {
        fnum*=2;
        str[i]=int(fnum)+48;
        i++;
        if(fnum>=1)
        {
            fnum--;

        }
        if(i==5)
        {
            break;
        }
    }
    strcpy(binf,str);
}
