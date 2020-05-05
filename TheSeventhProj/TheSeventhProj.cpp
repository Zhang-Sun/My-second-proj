#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>

using namespace std;
typedef char const *charpointer;
charpointer a[10][10]={
        {"太原","北京","天津","重庆","广州","香港","澳门","南宁","台北","上海"},
        {"商丘","郑州","平顶山","洛阳","驻马店","济源","漯河","信阳","开封","南阳"},
        {"石家庄","保定","雄安","张家口","秦皇岛","唐山","廊坊","衡水","巨鹿","遵化"},
        {"运城","霍州","侯马","临汾","武乡","长治","平遥","榆次","大同","朔州"},
        {"景德镇","崇义","瑞景","赣州","井冈山","广昌","德兴","武宁","庐山","南昌"},
        {"铜陵","太湖","马鞍山","南陵","合肥","芜湖","桐城","黄山","宿州","安庆"},
        {"新丰","乐昌","隆化","汕头","珠海","佛山","广宁","雷州","东莞","深圳"},
        {"武汉","武昌","襄樊","襄阳","鄂州","广水","黄冈","咸宁","佳木斯","宜昌"},
        {"神农架","荆门","巴东","牡丹江","玉门","敦煌","泰安","平凉","舟曲","碌曲"},
        {"玛曲","武都","兰州","哈尔滨","齐齐哈尔","公安","漠河","大庆","嫩江","鸡西"}
    };
void Jiami()
{
    int i,x,y;
    int code;
    char name[20];
    ofstream outfile("密文.txt",ios::out);
    while(1)
    {
        cout<<"请输入您的用户名:";
        cin>>name;
        if(strlen(name)>20)
        {
            cout<<"输入的用户名过长,请重新输入:"<<endl;
            continue;
        }
        else if(!(('A'<=name[i]&&name[i]<='Z')||('a'<=name[i]&&name[i]<='z')))
        {
            cout<<"用户名只接受字母,请重新输入:"<<endl;
            continue;
        }
        cout<<"请输入你的二位数密码:";
        cin>>code;
        if(code>100||code<0)
        {
            cout<<"输入的密码不合格，请重新输入"<<endl;
            continue;
        }
        break;
    }
    //开始加密
    for(unsigned int i=0;i<strlen(name);i++)
    {
        if(name[i]=='Z'||name[i]=='z')
        {
            name[i]-=25;
        }
        if(name[i]=='A'||name[i]=='a')
        {
            name[i]+=25;
        }
        else
        {
            name[i]+=1;
        }
        
    }
    x=code/10;
    y=code%10;
    outfile<<name<<" ";
    outfile<<a[x][y];
    outfile.close();
    cout<<"加密成功"<<endl;
}
void cryptograph()   //密文显示
{
    char temp[20];
    ifstream infile("密文.txt",ios::in);
    infile>>temp;
    cout<<endl<<"用户名的密文为:"<<temp;
    infile>>temp;
    cout<<endl<<"密码的密文为："<<temp;
}
void decode()
{
    char name[20];
    int code;
    char pt[2];
    ifstream infile("密文.txt",ios::in);
    infile>>name;
    for(unsigned int i=0;i<strlen(name);i++)
    {
        if(name[i]=='a'||name[i]=='A')
        {
            name[i]+=25;
        }
        else if(name[i]=='Z'||name[i]=='z')
        {
            name[i]-=25;
        }
        else
        {
            name[i]-=1;
        }
        cout<<name[i];
    }
    char temp[20];
    infile>>temp;
    int x,y;
    for(int m=0;m<10;m++)
    {
        for(int n=0;n<10;n++)
        {
            if(strcmp(a[m][n],temp)==0)
            {
                x=m;
                y=n;
            }
        }
    }
    pt[0]=x+48;
    pt[1]=y+48;
    cout<<"\n经解密,您的密码为:"<<pt[0]<<pt[1]<<endl;

}
int main()
{
    void jiami();
    void cryptograph();
    while(1)
    {
        int k=0;
        cout<<"请选择要执行的任务:1.加密   2.密文  3.解密  4.推出系统"<<endl;
        cin>>k;
        while(k<1||k>4)
        {
            cout<<"请重新输入:";
            cin>>k;
        }
        switch(k)
        {
            case 1:Jiami();break;
            case 2:decode();break;
            case 3:cryptograph();break;
            case 4:break;
        }
        return 0;
        
    }   

}