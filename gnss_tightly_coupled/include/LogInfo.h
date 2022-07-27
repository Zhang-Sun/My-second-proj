//
// Created by wlzhang on 2022/7/5.
//

#ifndef FAST_LIO_IGCOUPLED_LOGINFO_H
#define FAST_LIO_IGCOUPLED_LOGINFO_H

#include "easylogging++.h"

using namespace std;

namespace IGCoupled{

#ifdef WIN32
    #define LOGINI_PATH "..\\conf\\log.ini"
#else
    #define LOGINI_PATH "../conf/log.ini"
#endif

    string SetLogConfPath(string path);
    int SetLogLevel(int level);
    void InitLog(int argc, char** argv, string path, int level);


}

#endif //FAST_LIO_IGCOUPLED_LOGINFO_H
