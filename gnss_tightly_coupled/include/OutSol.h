//The class of outputing soulution file
// Created by wlzhang on 7/1/22.
//

#ifndef FAST_LIO_IGCOUPLED_OUTSOL_H
#define FAST_LIO_IGCOUPLED_OUTSOL_H

#include "CmnFunc.h"
#include "GnssFunc.h"
#include "InsFunc.h"

#define COMMENTH   "%"

namespace IGCoupled {
    class cOutSol{
    public:
        cOutSol();
        cOutSol(tIGCOUPLEDConf C);
        cOutSol(tIGCOUPLEDConf C,vector<tSolInfoUnit>& ref_sols);
        ~cOutSol();

    private:
        int MatchRefSol(cTime sol_time);
        tSolInfoUnit CompareSol(tSolInfoUnit& sol,tSolInfoUnit& ref_sol);

        int OutEcef(unsigned char* buff,const char *s,tSolInfoUnit& sol);
        int OutEnu(unsigned char* buff,const char *s,tSolInfoUnit& sol);
        int OutBlh(unsigned char* buff,const char *s,tSolInfoUnit& sol);
        int OutSolStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos,char *buff);
        int IGCoupledOut(tSolInfoUnit *sol,vector<tSatInfoUnit>& epoch_sat_info);
        void IGCoupledOutSat(vector<tSatInfoUnit>& epoch_sat_info);

    public:
        bool InitOutSol(tIGCOUPLEDConf C,string file);
        void WriteHead();

        void WriteSol(tSolInfoUnit sol,vector<tSatInfoUnit>&epoch_sat_info);
        void WriteSatStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos);

        void WriteImuHead();
        void WriteImuObs();

    public:
        tSolInfoUnit ref_sol_;

    private:
        int ref_index_=0;
        tIGCOUPLEDConf C_;
    public:
        string buff_;
        FILE *fout_;
        FILE *fout_stat_;
        ofstream  igcoupled_out_;

    public:
        vector<tSolInfoUnit> ref_sols_;
        cImuData *imus= nullptr;
    };
}


#endif //PPPLIB_OUTSOL_H
