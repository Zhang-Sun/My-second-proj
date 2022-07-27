## GNSS-LIO tightly coupled

### ENGLISH


This folder contains all code of fast_lio_gnss_tightly_coupled module. It can be an independent system, or it can be connected with fast_lio other modules.

**This module contains the following functions(at present):**
1. GNSS GPS + BDS + GLONASS + GALIEO  1/2/3 frequency PPP/SPP/PPK
2. GNSS/IMU loosely/ tightly coupled(PPP/INS,SPP/INS, PPK/INS)
3. Estimator contains lsq(spp)/EKF/FGO
4. IMU mechanical arrangement and pre-integration

**To be added**
1. GNSS/INS FGO tightly coupled
2. Lidar or combine with other modules of fast_lio

**How To Use**
```
* bulid:
    mkdir build
    cd build
    cmake ..
    make

* before run:
- FILE PREPARATION:
    Rover observation file(o-RINEX file)  - for SPP PPP PPK 
	Base observaion file (o-RINEX file)   - for PPK
	Broadcast ephemeris("ftp://igs.gnsswhu.cn/pub/gps/data/daily/<YEAR>/<DOY>/<YY>/")    - for SPP PPP PPK
	Differetial code bias("ftp://igs.ign.fr/pub/igs/products/mgex/dcb/<YEAR>/CAS0MGXRAP_<YYYY><DOY>0000_01D_01D_DCB.BSX.gz")   - for SPP PPP PPK 
	GNSS orbit sp3 file("ftp://igs.gnsswhu.cn/pub/whu/phasepias/<YEAR>/orbit/")(3-days)     - for PPP
	GNSS clock file ("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/clock/")                - for PPP PPK 
	Earth Rotation file("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/orbit/")             - for PPP PPK 
	Phase Bias file ("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/bias/")                 - for PPP
	Atenna information file("http://files.igs.org/pub/station/general/igs14.atx")           - for PPP PPK

- CONFIG FILE PREPARATION:
    prc_date (YEAR/MONTH/DAY)
	data_dir(absolute path)
	site_name(up to your observation file, name as '<site_name><doy0>.xxo')
	imu(name as '*.imu')
	gsof/gnss_sol (GNSS solution file , only use for INS/SOL loosely coupled)
	rslt (result file name , only for INS/SOL loosely coupled)
	estimator (up to your choice)
	filter_type (F/B/C)
	fgo_window (sliding window size , only for fgo estimator)
	integration_len (integration length , only for fgo estimator)

	frq_opt (up to your choice  1/ 2 / 3 frequency)
	ac_opt (up to your accuracy product organization , ftps above are wum)
	eph_opt (precise for PPP ， broadcast for SPP or PPK )
	ion_opt (ion-free for PPP , klobuchar for SPP , off for PPK)
	trp_opt (est for PPP , saas for SPP,  off for PPK)
	tid_opt (on for PPP,  off for SPP/PPK)
	sat_pcv / rec_ant (on for PPP, off for PPK / SPP)
	ar_mode (ppp-ar for PPP, fix and hold for PPK, off for SPP)
	ar_prod (ftp above is osb_whu)

	IMU_type (up to your IMU file)
	coord_type (RFU and FRD )
	data_format (increment or rate)
	gyro_val_format (degree or rad)
	sample_hz (IMU sample hz)
	ins_align (gnss sol for  pos/gsof  observation for other)

	sol_fmt (Igcoupled format and rtklib format)
	out_ins_mech_frq (base on imu)

- RUNNING
	
	 **run following at bin/**

	./IGCoupledMain -C ../example/conf/PPP/conf_file  -M PPP-KINE -S G -L 1

	-C  configuration file path
	-M  processing mode(SPP-KINE,PPP-KINE,PPK-KINE,IGLC-GSOF, IGLC-SOL, IGLC-SPP ,IGLC-PPP,IGLC-PPK, IGTC-SPP, IGTC-PPP,IGTC-PPK,INS)
	-S  selected GNSS(G,B,E,R,GBJ,GR,..)
	-L  debug level(1:Debug 32:Warning 128:Info)

- RESULT
	
	Result directory is in your dataset directory; 
```


**3rd Party and version:**
```
Eigen : 3.3.x
Ceres : <= 3.1.0 
```

**Tips**

* This module can be used in Windows/Linux
* The executable files are generated in /bin folder after building
* The result is in you data file
* The module needs configuration file

**PROJECT FRAME**
.
├── CMakeLists.txt
├── README.md
├── example
│   ├── conf
│   │   ├── LC
│   │   │   ├── IGCoupled_LC_GSOF_M39.ini
│   │   │   ├── IGCoupled_LC_PPK_CPT.ini
│   │   │   ├── IGCoupled_LC_PPP_CPT.ini
│   │   │   └── IGCoupled_LC_SOL_CPT.ini
│   │   ├── PPK
│   │   │   ├── IGCP_PPK_CPT.ini
│   │   │   └── myeasylog.log
│   │   ├── PPP
│   │   │   ├── IGCP_PPP_CPT.ini
│   │   │   ├── IGCP_PPP_FCB.ini
│   │   │   ├── IGCP_PPP_IRC.ini
│   │   │   └── myeasylog.log
│   │   ├── SPP
│   │   │   ├── IGCP_SPP_M39.ini
│   │   │   └── myeasylog.log
│   │   ├── TC
│   │   │   ├── IGCP_TC_PPK_CPT.ini
│   │   │   ├── IGCP_TC_PPK_POS.ini
│   │   │   └── IGCP_TC_PPP_CPT.ini
│   │   └── log.ini
│   └── dataset
│       └── README.md
├── exe
│   ├── CMakeLists.txt
│   └── IGCoupledMain.cc
├── include
│   ├── AdjFunc.h
│   ├── CmnFunc.h
│   ├── Constant.h
│   ├── DecodeRaw.h
│   ├── GnssAR.h
│   ├── GnssErrorModel.h
│   ├── GnssFunc.h
│   ├── InsFunc.h
│   ├── LogInfo.h
│   ├── OutSol.h
│   ├── README.md
│   ├── ReadFiles.h
│   ├── Solver.h
│   ├── dirent.h
│   ├── easylogging++.h
│   └── unistd.h
├── source
│   ├── 1.gif
│   ├── FGO.gif
│   ├── Framge1.png
│   └── Framge2.png
└── src
    ├── AdjFunc
    │   ├── AdjFunc.cc
    │   └── CMakeLists.txt
    ├── CMakeLists.txt
    ├── CmnFunc
    │   ├── CMakeLists.txt
    │   └── CmnFunc.cc
    ├── DecodeRaw
    │   ├── CMakeLists.txt
    │   └── DecodeRaw.cc
    ├── GnssFunc
    │   ├── CMakeLists.txt
    │   ├── GnssAR.cc
    │   ├── GnssErrorModel.cc
    │   └── GnssFunc.cc
    ├── InsFunc
    │   ├── CMakeLists.txt
    │   └── InsFunc.cc
    ├── LogInfo
    │   ├── CMakeLists.txt
    │   ├── LogInfo.cc
    │   └── easylogging++.cc
    ├── OutSol
    │   ├── CMakeLists.txt
    │   └── OutSol.cc
    ├── README.md
    ├── ReadFiles
    │   ├── CMakeLists.txt
    │   └── ReadFile.cc
    └── Solver
        ├── CMakeLists.txt
        └── Solver.cc

**Code Frame**
![image](gnss_tightly_coupled/source/Framge1.png)
![image](gnss_tightly_coupled/source/Framge2.png)


**FGO Loosely Coupled**
![image](gnss_tightly_coupled/source/FGO.png)

**FGO Dynamic Demonstration**
![image](gnss_tightly_coupled/source/1.gif)

**Update Logs:**
* Add time outage simulation.(2022/07/08 by wlzhang)
* Add F/B EKF of Fusion Solver.(2022/07/08 by wlzhang)
* Add INS Pre-Integration Factor(2022/07/11 by wlzhang)
* Add Mariginalization / Residuals / Pose local parameterization Factor(2022/07/12 by wlzhang)
* Add FGO Solver (2022/07/12 by wlzhang)
* Add WUM OSB File Reading and correction(2022/07/19 by wlzhang)




----------------------------------------border line-------------------------------------------------------------


### 中文


此文件夹包含fast_lio_gnss_tightly_coupled板块的所有代码。可以作为单独的系统进行解算，也可以将其与fast_lio其他板块进行融合(待完成)

**此模块包含以下功能：**
1. GNSS GPS + BDS + GLONASS + GALIEO（多系统）- 单频、双频、三频（多频） - PPP/SPP/PPK
2. GNSS/IMU 松、紧组合（PPP/INS、SPP/INS、PPK/INS）
3. 包含的估计器：lsq（最小二乘）、EKF、FGO（因子图优化）
4. IMU机械编排与预积分

**待添加：**
1. GNSS/INS 因子图优化紧组合
2. 加入Lidar观测值或与fast_lio其他模块融合
3. 加入视觉


**如何使用**
```
* 编译:
    mkdir build
    cd build
    cmake ..
    make

* 运行前准备:
- 文件准备:
    移动站观测文件(o-RINEX file)  - for SPP PPP PPK 
	基准站观测文件 (o-RINEX file)   - for PPK
	广播星历("ftp://igs.gnsswhu.cn/pub/gps/data/daily/<YEAR>/<DOY>/<YY>/")    - for SPP PPP PPK
	DCB文件("ftp://igs.ign.fr/pub/igs/products/mgex/dcb/<YEAR>/CAS0MGXRAP_<YYYY><DOY>0000_01D_01D_DCB.BSX.gz")   - for SPP PPP PPK 
	卫星轨道文件("ftp://igs.gnsswhu.cn/pub/whu/phasepias/<YEAR>/orbit/")(3-days)     - for PPP
	卫星钟差文件 ("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/clock/")                - for PPP PPK 
	地球旋转文件("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/orbit/")             - for PPP PPK 
	相位偏差文件 ("ftp://igs.gnsswhu.cn/pub/whu/phasebias/<YEAR>/bias/")                 - for PPP
	天线信息文件("http://files.igs.org/pub/station/general/igs14.atx")           - for PPP PPK

- 配置文件:
    prc_date (YEAR/MONTH/DAY)
	data_dir(absolute path)
	site_name(up to your observation file, name as '<site_name><doy0>.xxo')
	imu(name as '*.imu')
	gsof/gnss_sol (GNSS solution file , only use for INS/SOL loosely coupled)
	rslt (result file name , only for INS/SOL loosely coupled)
	estimator (up to your choice)
	filter_type (F/B/C)
	fgo_window (sliding window size , only for fgo estimator)
	integration_len (integration length , only for fgo estimator)

	frq_opt (up to your choice  1/ 2 / 3 frequency)
	ac_opt (up to your accuracy product organization , ftps above are wum)
	eph_opt (precise for PPP ， broadcast for SPP or PPK )
	ion_opt (ion-free for PPP , klobuchar for SPP , off for PPK)
	trp_opt (est for PPP , saas for SPP,  off for PPK)
	tid_opt (on for PPP,  off for SPP/PPK)
	sat_pcv / rec_ant (on for PPP, off for PPK / SPP)
	ar_mode (ppp-ar for PPP, fix and hold for PPK, off for SPP)
	ar_prod (ftp above is osb_whu)

	IMU_type (up to your IMU file)
	coord_type (RFU and FRD )
	data_format (increment or rate)
	gyro_val_format (degree or rad)
	sample_hz (IMU sample hz)
	ins_align (gnss sol for  pos/gsof  observation for other)

	sol_fmt (Igcoupled format and rtklib format)
	out_ins_mech_frq (base on imu)

- 运行
	
	 **在bin/ 下运行下列命令**

	./IGCoupledMain -C ../example/conf/PPP/conf_file  -M PPP-KINE -S G -L 1

	-C  configuration file path
	-M  processing mode(SPP-KINE,PPP-KINE,PPK-KINE,IGLC-GSOF, IGLC-SOL, IGLC-SPP ,IGLC-PPP,IGLC-PPK, IGTC-SPP, IGTC-PPP,IGTC-PPK,INS)
	-S  selected GNSS(G,B,E,R,GBJ,GR,..)
	-L  debug level(1:Debug 32:Warning 128:Info)

- 结果
	
	结果文件夹在数据集文件内; 
```


**第三方库及其版本**
```
Eigen : 3.3.x
Ceres : <= 3.1.0
```

**注意：**
* 此模块可在Windows/Linux系统中编译
* 可执行文件在本文件夹的bin/目录下
* 结果文件在数据文件夹中
* 运行需要通过配置文件，配置文件为程序解算提供必要的参数支持，后续可根据需要删除配置文件部分

**工程结构**
.
├── CMakeLists.txt
├── README.md
├── example
│   ├── conf
│   │   ├── LC
│   │   │   ├── IGCoupled_LC_GSOF_M39.ini
│   │   │   ├── IGCoupled_LC_PPK_CPT.ini
│   │   │   ├── IGCoupled_LC_PPP_CPT.ini
│   │   │   └── IGCoupled_LC_SOL_CPT.ini
│   │   ├── PPK
│   │   │   ├── IGCP_PPK_CPT.ini
│   │   │   └── myeasylog.log
│   │   ├── PPP
│   │   │   ├── IGCP_PPP_CPT.ini
│   │   │   ├── IGCP_PPP_FCB.ini
│   │   │   ├── IGCP_PPP_IRC.ini
│   │   │   └── myeasylog.log
│   │   ├── SPP
│   │   │   ├── IGCP_SPP_M39.ini
│   │   │   └── myeasylog.log
│   │   ├── TC
│   │   │   ├── IGCP_TC_PPK_CPT.ini
│   │   │   ├── IGCP_TC_PPK_POS.ini
│   │   │   └── IGCP_TC_PPP_CPT.ini
│   │   └── log.ini
│   └── dataset
│       └── README.md
├── exe
│   ├── CMakeLists.txt
│   └── IGCoupledMain.cc
├── include
│   ├── AdjFunc.h
│   ├── CmnFunc.h
│   ├── Constant.h
│   ├── DecodeRaw.h
│   ├── GnssAR.h
│   ├── GnssErrorModel.h
│   ├── GnssFunc.h
│   ├── InsFunc.h
│   ├── LogInfo.h
│   ├── OutSol.h
│   ├── README.md
│   ├── ReadFiles.h
│   ├── Solver.h
│   ├── dirent.h
│   ├── easylogging++.h
│   └── unistd.h
├── source
│   ├── 1.gif
│   ├── FGO.gif
│   ├── Framge1.png
│   └── Framge2.png
└── src
    ├── AdjFunc
    │   ├── AdjFunc.cc
    │   └── CMakeLists.txt
    ├── CMakeLists.txt
    ├── CmnFunc
    │   ├── CMakeLists.txt
    │   └── CmnFunc.cc
    ├── DecodeRaw
    │   ├── CMakeLists.txt
    │   └── DecodeRaw.cc
    ├── GnssFunc
    │   ├── CMakeLists.txt
    │   ├── GnssAR.cc
    │   ├── GnssErrorModel.cc
    │   └── GnssFunc.cc
    ├── InsFunc
    │   ├── CMakeLists.txt
    │   └── InsFunc.cc
    ├── LogInfo
    │   ├── CMakeLists.txt
    │   ├── LogInfo.cc
    │   └── easylogging++.cc
    ├── OutSol
    │   ├── CMakeLists.txt
    │   └── OutSol.cc
    ├── README.md
    ├── ReadFiles
    │   ├── CMakeLists.txt
    │   └── ReadFile.cc
    └── Solver
        ├── CMakeLists.txt
        └── Solver.cc


**代码框架**
![image](gnss_tightly_coupled/source/Framge1.png)
![image](gnss_tightly_coupled/source/Framge2.png)


**FGO 松组合流程**
![image](gnss_tightly_coupled/source/FGO.png)

**FGO动态演示**
![image](gnss_tightly_coupled/source/1.gif)

**更新日志：**
* 加入gnss缺失时间模拟。(2022/07/08/ by wlzhang)
* 加入双向滤波算法.(2022/07/08 by wlzhang)
* 加入IMU预积分及其残差（2022/07/11 by wlzhang）
* 加入边际化、残差块、位姿局部参数化等模块（2022/07/12 by wlzhang）
* 加入FGO处理器 （2022/07/12 by wlzhang）
* GNSS模块加入WHU的OSB改正观测值从而进行模糊度固定的方法。（2022/07/19 by wlzhang）
