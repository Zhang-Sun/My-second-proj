// some GNSS common function / struct / class  statement 
// Created by wlzhang on 6/30/22.
// 

#ifndef FAST_LIO_IGCOUPLED_GNSSFUNC_H
#define FAST_LIO_IGCOUPLED_GNSSFUNC_H

#include "CmnFunc.h"

namespace IGCoupled{

    typedef struct{
        int no;   
        int sys;
        int prn;
        string id;
        int sys_idx;
    }tSat; // struct of Satellite number

    class cSat{
    public:
        cSat();
        cSat(int sat_no);
        cSat(int sat_sys,int sat_prn);
        cSat(string sat_id);
        ~cSat();

    public:
        void SatPrn2No();
        void SatNo2Prn();
        void SatId2No();
        void SatNo2Id();

    public:
        tSat sat_={0};
    };// class of Satellite 

    typedef struct{
        cSat sat;
        int iode,iodc;
        int sva,svh;
        int week,code;
        cTime toe,toc,ttr;
        double A,e,i0,Omg0,omg,M0,deln,Omgd,idot;
        double crc,crs,cuc,cus,cic,cis;
        double toes;
        double f0,f1,f2;
        Vector4d tgd;
        double Adot,ndot;
    }tBrdEphUnit;// struct of Broadcast ephemeris

    typedef struct{
        cSat sat;
        int iode,frq;
        int svh,sva,age;
        cTime toe,tof;
        double taun,gamn,dtaun;
        Vector3d pos;
        Vector3d vel;
        Vector3d acc;
    }tBrdGloEphUnit; // struct of Glonass Broadcast ephemeris

    typedef struct{
        cTime t_tag;
        Vector4d pos[MAX_SAT_NUM];
        Vector4d vel[MAX_SAT_NUM];
        Vector4d std_pos[MAX_SAT_NUM];
        Vector4d std_vel[MAX_SAT_NUM];
        double clk[MAX_SAT_NUM];
        double std_clk[MAX_SAT_NUM];
    }tPreOrbUnit; // //struct of Precise orbit(orb) product Unit

    typedef struct{
        cTime t_tag;
        double clk[MAX_SAT_NUM];
        float  std[MAX_SAT_NUM];
    }tPreClkUnit; // struct of Precise clock(clk) product Unit

    typedef struct{
        double mjd;
        double xp,yp;
        double xpr,ypr;
        double ut1_utc;
        double lod;
    }tErpUnit;  // struct of Earth Rotation Product(ERP) Unit

    typedef struct{
        cSat  sat;
        cTime ts,te;
        string ant_type;
        string ser_code;
        Vector3d pco[MAX_GNSS_FRQ_NUM*NSYS];
        double   pcv[MAX_GNSS_FRQ_NUM*NSYS][80*30];
        double dazi;
        double zen1,zen2,dzen;
        Vector3d rec_ant_del[2];
    }tAntUnit;//  struct of antenna(Ant) Unit

    typedef struct{
        cTime t_tag;
        int ndata[3]{};
        double re;
        double lats[3]{};
        double lons[3]{};
        double hgts[3]{};
        vector<double> data;
        vector<float>  rms;
    }tTecUnit;  // struct of Ionospheric electron content Unit

    typedef struct {
        string name;
        string marker;
        string ant_desc;
        string ant_seri;
        string rec_type;
        string firm_ver;
        Vector3d del;
        Vector3d apr_pos;
        double ant_hgt;
    }tStaInfoUnit; // struct of Station infomation Unit

    typedef struct {
        cTime ts,te;
        double nl_fcb[MAX_SAT_NUM];
    }tFcbUnit;  // struct of fractional code bias (fcb) Unit 

    typedef struct {
        cTime ts,te;
        double code_osb[MAX_GNSS_CODE_TYPE];
        double phase_osb[MAX_GNSS_CODE_TYPE];
    }tOsbUnit;  // struct of Observation bias (OSB) Unit

    typedef struct{
        vector<tBrdEphUnit>    brd_eph;
        vector<tBrdGloEphUnit> brd_glo_eph;
        vector<tPreOrbUnit> pre_eph;
        vector<tPreClkUnit> pre_clk;
        vector<tPreClkUnit> pre_pr_clk;
        vector<tErpUnit>erp_paras;
        vector<tTecUnit>tec_paras;

        tStaInfoUnit sta_paras[2];
        tAntUnit sat_ant[MAX_SAT_NUM];
        tAntUnit rec_ant[2];

        int glo_frq_num[GLO_MAX_PRN+1];
        double glo_cp_bias[4];
        int leaps;
        double ion_para[NSYS][8];
        double utc_para[NSYS][4];
        double code_bias[MAX_SAT_NUM][MAX_GNSS_CODE_BIAS_PAIRS];
        double ocean_paras[2][6*11];

        double wide_line_bias[MAX_SAT_NUM];
        vector<tFcbUnit> nl_fcbs;
        vector<tOsbUnit> obs_osbs;

		//osb for wum
        double code_osb[MAX_SAT_NUM][MAX_GNSS_CODE_TYPE];
        double phase_osb[MAX_SAT_NUM][MAX_GNSS_CODE_TYPE];
    }tNav; // struct of Navigation information Unit

    typedef struct{
        int n;
        int frq[MAX_GNSS_OBS_TYPE];
        int pos[MAX_GNSS_OBS_TYPE];
        unsigned char pri[MAX_GNSS_OBS_TYPE];
        unsigned char type[MAX_GNSS_OBS_TYPE];
        unsigned char code[MAX_GNSS_OBS_TYPE];
        double shift[MAX_GNSS_OBS_TYPE];
    }tGnssSignal; // struct of GNSS signal information

    typedef struct{
        cSat sat;
        double P[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char code[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double L[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        float D[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char SNR[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char LLI[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double frq[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double CSC_P[MAX_GNSS_FRQ_NUM];
    }tSatObsUnit; // struct of Satellite observation Unit 

    typedef struct{
        cTime obs_time;
        int sat_num;
        vector<tSatObsUnit> epoch_data;
    }tEpochSatUnit; // struct of satellite observation of one epoch Unit

    class cGnssObs{
    public:
        cGnssObs();
        ~cGnssObs();

    public:
        void SetTimeSpan(cTime* ts, cTime* te);
        cTime* GetStartTime();
        cTime* GetEndTime();
        void SetRcvIdx(RECEIVER_INDEX rcv);
        vector<tEpochSatUnit>& GetGnssObs();
        tStaInfoUnit* GetStation();
        void ResetGnssObs();
        void CarrierSmoothCode(int smooth_length);
		void ReverseObs();

    public:
        int epoch_num;
        double sample_;
        RECEIVER_INDEX rcv_idx_;
        tGnssSignal *signal_[NSYS];

    private:
        cTime ts_,te_;
        vector<tEpochSatUnit> obs_;
        tStaInfoUnit station_;
    }; // calss of GNSS Observation 

    typedef struct {
        cTime t_tag[4];          // last epoch time
        int n[4];                // number of epochs
        double lc_amb[4];        // linear combination average */
        double var_amb[4];       // linear combination variance */
        int fix_cnt;             // fix count */
        char flags[MAX_SAT_NUM]; // fix flags */
    }tLcAmb;  // struct of linear combination ambiguity 

    typedef struct {
        int ref_sat,sat;
        int f;
        int inherit_flag;
    }tDdAmb;                     // struct of double difference ambiguity

    typedef struct {
        int ref_sat_idx,sat_idx;
        int f;
    }tSdAmb;                     // struct of  single difference ambiguity

    typedef struct{
        cSat sat;
        cTime t_tag;
        cTime t_trans;
        GNSS_SAT_STAT stat;
        int brd_eph_index;

        unsigned char P_code[MAX_GNSS_USED_FRQ_NUM];
        double raw_P[MAX_GNSS_USED_FRQ_NUM];        // L1 L2 L5
        double raw_L[MAX_GNSS_USED_FRQ_NUM];
        double raw_D[MAX_GNSS_USED_FRQ_NUM];
        unsigned char raw_S[MAX_GNSS_USED_FRQ_NUM];
        unsigned char LLI[MAX_GNSS_USED_FRQ_NUM];
        double csc_P[MAX_GNSS_USED_FRQ_NUM];
        double cor_P[MAX_GNSS_USED_FRQ_NUM];        // corrected code bias BDS satellite-specific multipath
        double cor_L[MAX_GNSS_USED_FRQ_NUM];        // corrected phase bias for PPP-AR
        double cor_D[MAX_GNSS_USED_FRQ_NUM];        // corrected doppler 
        double lam[MAX_GNSS_USED_FRQ_NUM];
        double frq[MAX_GNSS_USED_FRQ_NUM];

        Vector3d brd_pos;
        Vector3d brd_vel;
        Vector2d brd_clk;  // clk, clk-rate
        Vector3d pre_pos;
        Vector3d pre_vel;
        Vector2d pre_clk;
        Vector2d pre_pr_clk;

        Vector2d trp_dry_delay; // slant_dry,map_dry
        Vector4d trp_wet_delay; // slant_wet,map_wet,grand_e,grand_n
        Vector2d ion_delay; // L1_slant_ion, map_ion;
        double clk_rel;
        double sagnac;
        double shapiro;

        double code_bias[MAX_GNSS_USED_FRQ_NUM];
        double cp_bias[MAX_GNSS_USED_FRQ_NUM];
        double bd2_mp[3];
        double phase_wp;
        double float_amb[MAX_GNSS_USED_FRQ_NUM]; //L1,L2,L5
        double fix_amb[MAX_GNSS_USED_FRQ_NUM];

        double omc_P[MAX_GNSS_USED_FRQ_NUM];
        double omc_L[MAX_GNSS_USED_FRQ_NUM];

        double prior_res_P[MAX_GNSS_USED_FRQ_NUM];
        double prior_res_L[MAX_GNSS_USED_FRQ_NUM];
        double post_res_P[MAX_GNSS_USED_FRQ_NUM];
        double post_res_L[MAX_GNSS_USED_FRQ_NUM];

        Vector3d sig_vec;
        Vector2d el_az;
        Vector2d raw_mw;  //L1_L2, L1_L5
        Vector2d sm_mw;   //L1_L2, L1_L5
        Vector2d var_mw;
        int mw_idx[2];
        Vector2d gf;          //L1_L2, L1_L5
        Vector2d multipath_comb;
        double glo_bias[MAX_GNSS_USED_FRQ_NUM];

        double tdcp;
        double cmc_P;             // code minus carrier
        double cor_if_P[2];     // SF-IF or DF-IF or TF-IF1/TF-IF2
        double cor_if_L[2];     // 0 :L1-L2-L5 1: L1-L2 2:L1-L5
        int tf_if_idx[3];
        tLcAmb lc_amb;

        int svh;
        double brd_eph_var;
        double pre_eph_var;
        double ion_var;
        double trp_var;
        double c_var_factor[MAX_GNSS_USED_FRQ_NUM];
        double p_var_factor[MAX_GNSS_USED_FRQ_NUM];
        int outc[MAX_GNSS_USED_FRQ_NUM];
        int lock[MAX_GNSS_USED_FRQ_NUM];
        unsigned char rejc[MAX_GNSS_USED_FRQ_NUM]; // reject flag for exclude satellite
        unsigned char rejc_phase[MAX_GNSS_USED_FRQ_NUM]; // reject flag for exclude satellite due to phase residual
        unsigned char slip[MAX_GNSS_USED_FRQ_NUM];
        unsigned char fix[MAX_GNSS_USED_FRQ_NUM];
        unsigned char vsat[MAX_GNSS_USED_FRQ_NUM];
        unsigned char res_idx[MAX_GNSS_USED_FRQ_NUM];

        double res_wl;  // for ppp-ar, cycle;
        double res_nl;
    }tSatInfoUnit; // struct of satellite information Unit

    class cGnssObsOperator {
    public:
        cGnssObsOperator();
        ~cGnssObsOperator();


    public:
        void ReAlignObs(tIGCOUPLEDConf C,tSatInfoUnit& sat_info, tSatObsUnit sat_obs,int f,int raw_idx,int frq_idx,int *glo_frq);

        double GnssObsIfComb(double obs1,double obs2,double f1,double f2);
        double GnssObsMwComb(double obs_P1,double obs_P2,double obs_L1,double obs_L2,double lam1,double lam2);
        double GnssObsGfComb(double obs1,double obs2,double lam1,double lam2);
        double GnssObsCmcComb(double obs_P,double obs_L1,double obs_L2,double f1,double f2);
        double GnssObsTdComb(double cur_obs,double pre_obs);
        void MakeGnssObsComb(tIGCOUPLEDConf C,GNSS_OBS_COMB type,tSatInfoUnit* sat_info,const tSatInfoUnit previous_sat_info);
        double GnssSdObs(const tSatInfoUnit& sat1,const tSatInfoUnit& sat2,int f,GNSS_OBS type);
        double GnssObsLinearComb(tIGCOUPLEDConf C,int i,int j,int k,tSatInfoUnit& sat_info, GNSS_OBS type,double* var);
        double LinearCombLam(int i,int j,int k,tSatInfoUnit& sat_info);

        bool CodeMinuPhase(tIGCOUPLEDConf C, tSatInfoUnit& sat_info,tSatInfoUnit& pre_sat_info);
        bool LliCycleSlip(tIGCOUPLEDConf C, tSatInfoUnit& sat_info,int nf,double tt,RECEIVER_INDEX rcv, unsigned char* s);
        bool LliCycleSlip(tIGCOUPLEDConf C, tSatInfoUnit& sat_info, tSatInfoUnit& pre_sat_info,int nf, double tt, RECEIVER_INDEX rcv);
        void MwCycleSlip(tIGCOUPLEDConf C,double sample_dt,double dt,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat,tTime last_time);
        void GfCycleSlip(tIGCOUPLEDConf C,double sample_dt,double dt,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat);
        void SmoothMw(tIGCOUPLEDConf C,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat);

    }; // class of GNSS Observation operator function

    typedef struct {
        cTime t_tag;
        int num_use_sat;
        int rcv_status;
        int no;

        int sol_level;
        int vel_flag;
        Vector3d pos;
        Vector3d llh;
        Vector3d vel;
        Vector3d rb_delta;
        double dop[4];
        float cov[8];
        float sig[6]; //enu
    }tGsofUnit; // struct of GSOF Unit 

    typedef struct {
        int n,nmax;
        tGsofUnit *data;
    }tGsofs;

    double GeoDist(Vector3d sat_pos,Vector3d rec_pos,Vector3d& sig_vec); // Distance between station and satellite
    double SatElAz(Vector3d rec_blh,Vector3d sig_vec,Vector2d& el_az); // Satellite azimuth and elevation 
    double GnssMeasVar(tIGCOUPLEDConf C, GNSS_OBS obs_type,tSatInfoUnit sat_info); // GNSS measurement variance
    double SunMoonPos(cTime ut1t,const double *erp_val,Vector3d& sun_pos, Vector3d& moon_pos); // Sun and Moon position
}



#endif //FAST_LIO_IGCOUPLED_GNSSFUNC_H
