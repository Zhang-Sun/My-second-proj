## FAST_LIO_GNSS_TIGHTLY_COUPLED/INCLUDE
```
此文件夹包含了Fast_lio_gnss_tightly_coupled所涉及的全部头文件，对应gnss_tightly_coupled中不同的板块构建，可根据需要进行增添删除

包括：文件读取、GNSS误差构建、估计器、结果输出、处理器等。
```

* 在日常维护中，需要注意以下几个部分：
  1. Constant.h 中需要维护卫星个数、频率信息、观测值类型等。一般卫星个数需要根据观测值中的最大观测去定义，可以超过现有卫星个数，但不能少于。
  2. Constant.h 中kUtcLeapSec需要根据UTC跳秒进行更新，一般为每2-3年更新一次。
  3. Solver.h 中定义了GNSS处理的Solver。可以在此扩充Solver处理类型。
  4. ReadFile.h 中可以扩充需要读取文件的类型.

* 待解决的问题:
  1. InsFuc.h中#include<ceres/ceres.h>需要确定电脑ceres的安装路径，或者可将ceres放入3rdparty中。
  2. ceres中的glog 与 easylogging++中的LOG冲突，需要将代码库中涉及easylogging++中的LOG使用原名。
  
