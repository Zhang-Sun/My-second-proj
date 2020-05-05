#include<iostream>
#include<fstream>
#include<string>
using namespace std;

int main()
{
    string filename;
    cout<<"请输入文件的名称：";
    cin>>filename;
    ifstream file(filename);
    string str;
    string buff[14][5];
    int i=0,j=0;
    while(file>>str)
    {
        if(j<4)
        {
            buff[i][j]=str;
            j++;
        }
        else
        {
            i++;
            j=0;
        }
    }
}

