### FAST_LIO_GNSS_TIGHTLY_COUPLED/SRC
```
此文件夹包含了fast_lio_gnss_tightly_coupled的主要src文件，对应于头文件../include，可根据需求进行增添删除

```

* 在日常维护中需要注意以下问题：
  1. ReadFile.cc 中可根据需要对读取固定格式文件部分代码进行修改。因为产品或文件的命名规则在不断变化，或者使其适用于自己的文件。
  2. InsFunc.cc 中目前只支持零偏的估计，暂不支持杆臂、安装误差等的估计。bool类型的ErrModel为留的接口。
  3. InsFunc.cc 中的cInsAlign部分暂时也没有使用。
  4. InsFunc.cc 中的状态转移矩阵构建部分，可根据需求进行修改。
  5. DecodeRaw.cc 中可根据需求添加需要Decode的原始文件。
  6. CmnFunc.cc 中可根据需求增添常用的函数。
  7. 日志打印使用C++第三方头文件easylogging++.h 与 easylogging++.cc，当使用ceres时，会与ceres中goole log 产生冲突，此时需修改easylogging++.cc中的函数。


* 待解决问题：
  1. 
