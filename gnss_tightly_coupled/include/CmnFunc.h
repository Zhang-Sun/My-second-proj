// some common function / struct / class  in GNSS and INS calculate
// Created by wlzhang on 6/28/22.
//

#ifndef FAST_LIO_IGCOUPLED_CMNFUNC_H
#define FAST_LIO_IGCOUPLED_CMNFUNC_H

#include <Eigen/Dense>
#include <iomanip>
#include <regex>
#include "LogInfo.h"
#include "Constant.h"

#ifdef WIN32
#define FILEPATHSEP '\\'
#include "dirent.h"
#include "unistd.h"
#else
#define FILEPATHSEP '/'
#include <unistd.h>
#include <dirent.h>
#endif

using namespace std;
using namespace Eigen;

namespace IGCoupled{
    void CoutVector(VectorXd& vec,int p,int q,string s,bool trans); // output vector
    void CoutMatrix(MatrixXd& mat,int p,int q,string s,bool trans);

    int Round(double d);                                            
    double VectorMean(vector<double>& seri);
    double NormDistribution(const double u);                         // Normalization 
    double ReNorm(double p);                                         // 


    // vector dot multiplication
    template <typename Iter1,typename Iter2>
    double Dot(const Iter1 VecA,const Iter2 VecB,int SizeVec){
        double dInn=0.0;

        while (--SizeVec>=0){
            dInn+=VecA[SizeVec]*VecB[SizeVec];
        }
        return dInn;
    }

    // Vector module length
    template <typename Iter>
    double Norm(const Iter VecA,int SizeVec){
        return sqrt(Dot(VecA,VecA,SizeVec));
    }

    // calculate unit vector
    template <typename Iter1,typename Iter2>
    int NormV3(const Iter1 vec1,Iter2 vec2){
        double r;
        if((r=Norm(vec1,3))<=0.0) return 0;
        vec2[0]=vec1[0]/r;
        vec2[1]=vec1[1]/r;
        vec2[2]=vec1[2]/r;
        return 1;
    }

    // vector cross multiply
    template <typename Iter1,typename Iter2, typename Iter3>
    void CrossVec3(const Iter1 vec1,const Iter2 vec2,Iter3 vec3){
        vec3[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
        vec3[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
        vec3[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
    }


    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec);           // 
    extern void MatMul(const char *tr, int n, int k, int m, double alpha,
                       const double *A, const double *B, double beta, double *C);  // 
    extern int MatInv(double *A, int n);            // 

    string Doul2Str(int str_len, int dec_len, const string str_filler, const double src_num, string &dst_str);  // double to str
    string Int2Str(int str_len, const string str_filler, const int src_num, string &dst_str); // int to str
    int Str2Double(string src_str,double &dst_num); // str to double 
    int Str2Int(const string src_str,int &dst_num); // str to int 
    string StrTrim(string s);                       // 
    void SplitString(const string& s, vector<string>& v,string c); // 
    vector<string> MultiSplitStr(const string &s, const string &seperator); // 
    vector<string> TextSplit(const string &in, const string &delim); //
    extern void CreateDir(const char *path);   //
    extern int ExPath(char *path,char *paths[],int nmax);
    typedef struct {
        time_t long_time;
        double sec;
    }tTime;                        

    typedef struct {
        long sn;
        double tos;
    }tSod;

    typedef struct {
        long day;
        tSod sod;
    }tMjd;

    class cTime {
    public:
        cTime();
        cTime(const double *ep);
        cTime(string str_time);
        cTime operator=(const tTime t);
        cTime operator+(double sec);
        void operator+=(double sec);
        ~cTime();

    public:
        double* GetEpoch();
        string GetTimeStr(int n);
        tMjd* GetMjd();
        int GetDoy();
        double TimeDiff(tTime t1);
        int Str2Time(string s);

        cTime *Epoch2Time(const double *ep);
        void Time2Epoch();
        double Time2Gpst(int* week,int* day, int sys);
        cTime* Gpst2Time(int week, double wos, int sys);
        cTime Utc2Gpst();
        cTime Gpst2Utc();
        cTime Gpst2Bdst();
        cTime Bdst2Gpst();
        double Utc2Gmst(double ut1_utc);
        cTime* AdjWeek(cTime t);
        cTime* AdjDay(cTime t);
        double Time2Doy();

    private:
        string Time2Str(int n);
        void Time2Mjd();
        double Time2Sec(cTime& day);

    private:
        double epoch_[6];
        string time_str_;
        int doy_;
        tMjd  mjd_;

    public:
        tTime t_;
    };

#if 0
    class cCoord{
    public:
        cCoord();
        cCoord(const Vector3d& coord, const COORDINATE_TYPE coord_type);
        cCoord(const double *coord, const COORDINATE_TYPE coord_type);
        ~cCoord();

    public:
        Vector3d GetCoordEnu(cCoord ref_coord);
        Vector3d GetCoordXyz();
        Vector3d GetCoordBlh();
        Vector3d GetCoordNed(cCoord ref_coord);
        Matrix3d GetCne(const COORDINATE_TYPE type);

    private:
        void Xyz2Blh();
        void Blh2Xyz();
        void Xyz2Enu(cCoord ref_coord);
        void Enu2Xyz(cCoord ref_coord);
        void Enu2Ned();
        void Ned2Enu();
        void CalcCne(const COORDINATE_TYPE type);
        double CalcLat();

    private:
        Vector3d coord_XYZ_;
        Vector3d coord_BLH_;
        Vector3d coord_ENU_;
        Vector3d coord_NED_;
        Matrix3d Cne_;
    };
#endif
    Vector3d Xyz2Blh(Vector3d& coord_xyz);
    Vector3d Blh2Xyz(Vector3d& coord_blh);
    Vector3d Xyz2Enu(Vector3d& coord_blh,Vector3d& coord_xyz);
    Vector3d Enu2Xyz(Vector3d& coord_blh,Vector3d& coord_enu);
    Vector3d Enu2Ned(Vector3d& coord_enu);
    Vector3d Ned2Enu(Vector3d& coord_ned);
    Matrix3d CalcCen(Vector3d& coord_blh,COORDINATE_TYPE lf_type);

    typedef struct{
        int restart_gap;
        int nav_sys;
        double sample_rate;
        double ele_min;
        GNSS_FRQ_OPT frq_opt;
        int gnss_frq[NSYS+1][MAX_GNSS_USED_FRQ_NUM];
        bool use_bd3;
        bool est_bd3_isb;
        bool adj_obs;
        bool csc;
        bool use_doppler;
        double code_phase_ratio;
        double meas_err_factor[5]; // [1-3] error factor a/b/c/of phase
        GNSS_AC_OPT  ac_opt;
        GNSS_EPH_OPT eph_opt;
        GNSS_ION_OPT ion_opt;
        GNSS_TRP_OPT trp_opt;
        GNSS_TID_OPT tid_opt;
        GLO_IFCB_OPT glo_ifcb_opt;
        bool sat_pcv;
        bool rec_ant;
        double cs_thres[2];   //mw and gf
        double max_pdop;
        double max_prior;
        int max_out;
        double max_inno;
        double ait_psd[3];
        bool check_dual_phase;
        GNSS_AR_MODE ar_mode;
        GNSS_AR_PROD ar_prod;
        GLO_AR_MODE glo_ar_mode;
        bool bds_ar_mode;
        bool gal_ar_mode;
        double ar_thres[8];
        double ar_el_mask;
        int min_sat_num2fix;
        int min_sat_num2drop;
        int min_lock2fix;
        double hold_er_mask;
        int min_sat_num2hold;
        int min_fix2hold;
        bool ar_filter;
        bool par_ar;
        bool res_qc;
        Vector3d rb;
		bool use_outage;
		double outage_time;
		double outage_len;
		double outage_period;
    }tGnssConf;

    typedef struct{
        IMU_TYPE imu_type;
        IMU_COORD_TYPE coord_type;
        IMU_DATA_FORMAT data_format;
        GYRO_DATA_FORMAT gyro_val_format;
        INS_ALIGN ins_align;
        double sample_rate;
        Vector3d lever;
        double correction_time_ba;
        double correction_time_bg;
        double init_pos_unc;
        double init_vel_unc;
        double init_att_unc;
        double init_ba_unc;
        double init_bg_unc;
        double init_sa_unc;
        double init_sg_unc;
        double init_ra_unc;
        double init_rg_unc;
        double init_lever_unc;
		double init_odo_unc;
        double psd_ba;
        double psd_bg;
        double psd_acce;
        double psd_gyro;
        bool err_model;
        bool est_sa;
        bool est_sg;
        bool est_ra;
        bool est_rg;
        bool est_level;
		bool use_odo;
		Vector3d abv; // b 系与v系的安装角 rad
		Vector3d lodo; // b系下的历程计
		Vector3d odo_std; // 历程计std m/s
		double odo_srw; // 历程计比例因子随机游走
    }tInsConf;

    typedef struct{
        string rover;
        string base;
        string brd;
        string cbias;
        string cod_dcb;
        string fcb;
		string osb;
        string clk;
        string pr_clk;
        string sp3[3];
        string erp;
        string atx;
        string gim;
        string blq;
        string imu;
        string gsof;
        string kine_ref;
        string snx;
        string sol;
        string sol_stat;
    }tFileConf;

    typedef struct{
        int sol_fmt;     //0: igcoupled fmt  1: rtklib fmt
        bool out_sol;
        bool out_head;
        bool out_vel;
        bool out_trp;
        bool out_att;
        bool out_ba;
        bool out_bg;
        bool out_err_fmt;
        bool out_stat;
        int out_ins_mech_frq;
        COORDINATE_TYPE sol_coord;
    }tSolConf;

	 typedef struct tFGOState {
        cTime t_tag;

        Vector3d p{0, 0, 0};
        Quaterniond q{0, 0, 0, 0};
        Vector3d v{0, 0, 0};

        Vector3d bg{0, 0, 0};
        Vector3d ba{0, 0, 0};

        Vector3d s{0, 0, 0};
        double sodo{0};
        Vector2d abv{0, 0};

        Vector3d sg{0, 0, 0};
        Vector3d sa{0, 0, 0};

        double pt[NSYS+1];//伪距钟
        double lt[NSYS];//相位钟
        double ifb[NSYS];
        double gloIfcb[NUM_GLO_SAT];
        double gloIfpb[2];


        double trp[6]; //对流层 （ppp ： 1 + 2  ppk： 1 + 2  1 + 2）


        double ion[MAX_SAT_NUM];
        double amb[MAX_SAT_NUM * 3]; // 模糊度 共AX_SAT_NUM*GetGnssUsedFrqs()个



    } tFGOState;

    typedef struct{
        cTime prc_date;
        string data_dir;
        string site_name;
        int site_idx;
        bool use_custom_dir;
        IGCOUPLED_MODE mode;
        IGCOUPLED_MODE_OPT mode_opt;
        int dynamic;
        SOLVE_ESTIMATOR estimator;
        FILTER_TYPE filter_type;
		int fgo_window;
		double integration_len;
		bool use_image;
        tGnssConf gnssC;
        tInsConf  insC;
        tFileConf fileC;
        tSolConf  solC;
    }tIGCOUPLEDConf;

    typedef struct {
        cTime t_tag;
        int epoch_idx;
        SOL_STAT stat;
        SOL_INS_STAT ins_stat;
        int observed_sat_num;
        int valid_sat_num;
        Vector3d pos;
        Vector3d vel;
        double q_pos[6];
        double q_vel[6];
        double clk_error[NSYS+1]; //clock+isb
        double rec_dcb[NSYS+1];
        double rec_ifcb[NSYS+1];   // inter-frequency code bias
        double zenith_trp_delay[4]; // ZHD ZWD ge gn
        double dops[4];
        float ratio;
        double age;
        int num_ar_sat;
        double sigma;

        Vector3d att{0,0,0};  //roll pitch yaw
        double q_att[6];
        Vector3d gyro_bias{0,0,0};
        Vector3d accl_bias{0,0,0};
    }tSolInfoUnit;

    class cParSetting{
    public:
        cParSetting();
        cParSetting(tIGCOUPLEDConf conf);
        ~cParSetting();

    public:
        int GetGnssUsedFrqs();
        int GetNumObsType();
        int GetInsTransParNum(tIGCOUPLEDConf C);
        int GetIGCOUPLEDPar(tIGCOUPLEDConf C);
        int GetRealFixParNum(tIGCOUPLEDConf C);

        int NumPos();
        int NumVel();

        // para of INS
        int NumAtt();
        int NumBa();
        int NumBg();
        int NumSa();
        int NumSg();
        int NumRa();
        int NumRg();
        int NumLever();

        // para of GNSS
        int NumClPar();
        int NumClk();
        int NumPhaseClk();
        int NumClkDrift();
        int NumDcb();
        int NumIfb();
        int NumGloIfcb();
        int NumGloIfpb();
        int NumTrp();
        int NumIon();
        int NumAmb();

        int IndexPos();
        int IndexVel();
        int IndexAtt();
        int IndexBa();
        int IndexBg();
        int IndexSa();
        int IndexSg();
        int IndexRa();
        int IndexRg();
        int IndexLever();

        int IndexClk(int sys_index);
        int IndexPhaseClk(int sys_index);
        int IndexClkDritf(int sys_index);
        int IndexDcb(int sys_index);
        int IndexIfb(int sys_index);
        int IndexGloIfcb(int i);
        int IndexGloIfpb();
        int IndexTrp();
        int IndexIon(int sat_no);
        int IndexAmb(int f,int sat_no);

    public:
        tIGCOUPLEDConf IGCOUPLEDC_;
    };

    class Config {
    public:
        using Ptr_=std::shared_ptr<Config>;

    private:
        Config() = default;
        Config(Config &&) = delete;
        Config(const Config &)= delete;
        Config &operator=(Config &&)= delete;
        Config &operator=(const Config &)= delete;
        static Ptr_ config_info_;
        std::map<std::string, std::string> storage_;

    public:
        ~Config() = default;

    public:
        static Ptr_ GetInstance();
        bool Open(std::string config_file);
        template <typename T>
        T Get(std::string key){
            transform(key.begin(),key.end(),key.begin(),::tolower);
            if(storage_.count(key)>0){
                try{
                    double value=stod(storage_[key]);
                    return static_cast<T>(value);
                }
                catch (const std::exception &e){
                    std::cerr<<e.what()<<'\n';
                }
            }
            else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
                return T(0x0);
            }
        }
        template <typename  T>
        std::vector<T> GetArray(std::string key){
            std::vector<T> data;
            transform(key.begin(),key.end(),key.begin(),::tolower);
            if(storage_.count(key)>0){
                try{
                    auto text=TextSplit(storage_[key],",");
                    for(auto index:text){
                        double value=stod(index);
                        data.emplace_back(static_cast<T>(value));
                    }
                }
                catch(const std::exception &e){
                    std::cerr<<e.what()<<'\n';
                }
            }
            else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
            }
            return data;
        }
    };

    template <>
    inline std::string Config::Get<std::string>(std::string key){
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                return std::string(storage_[key]);
            }
            catch (const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
//            getchar();
            return "";
        }
    }

    template <>
    inline std::vector<std::string> Config::GetArray<std::string>(std::string key){
        std::vector<std::string> data;
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                data=TextSplit(storage_[key],",");
            }
            catch(const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
            getchar();
        }
        return data;
    }




}


#endif //FAST_LIO_IGCOUPLED_CMNFUNC_H
