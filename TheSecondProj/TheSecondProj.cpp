#include<iostream>
#include<time.h>
#include<stdlib.h>

using namespace std;
//用户所猜数字：stunum

class Guess
{
    private:
    int value;
    int CompareTimes;
    public:
    Guess();
    int Compare(int InputValue)
    {
        CompareTimes++;
        return InputValue-value;
    }
    int GetCompareTimes()
    {
        return CompareTimes;
    }

};
Guess::Guess()
{
    CompareTimes=0;
    srand((unsigned)time(NULL));
    value=rand()%100;
}
int main()
{
    int InputValue;
    cout<<"\n**************欢迎使用本程序******************\n";
    for(;;)
    {
        char Select;
        Guess guessobj;
        for(;;)
        {
            int CompareResult;
    
            cout<<"我已经想好数字啦(0-99)，请猜猜吧！"<<endl;
          
            cout<<"我想的是:";
            cin>>InputValue;
            CompareResult=guessobj.Compare(InputValue);
            if(CompareResult==0)
            {
                int GuessTimes=guessobj.GetCompareTimes();
                cout<<"\n恭喜你猜对了！"<<endl<<"您一共猜了"<<GuessTimes<<"次！"<<endl;
                break;
            }
            else if(CompareResult>0)
            {
                cout<<"\n哎呀，猜的数大了\n"<<endl;
            }
            else
            {
                cout<<"\n哎呀，猜的数小了\n"<<endl;                
            }
            

        }
        cout<<"\n您还想再玩吗？('n'=No,Other=Yes)\n"<<endl;
        cin>>Select;
        if(Select=='n')
        {
            break;
        }


    }
    cout<<"*********感谢您的使用********"<<endl;
    return 0;
}