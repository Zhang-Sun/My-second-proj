//
// Created by wlzhang on 7/5/22.
// 

#include "OutSol.h"

using namespace IGCoupled;

namespace IGCoupled{

    /*
    * Function : The constructor and destructor of class cOutSol
    */
    cOutSol::cOutSol() {}

    cOutSol::cOutSol(tIGCOUPLEDConf C) {
        C_=C;
    }

    cOutSol::cOutSol(tIGCOUPLEDConf C,vector<tSolInfoUnit>& ref_sols) {C_=C;ref_sols_=ref_sols;}

    cOutSol::~cOutSol() {}

    /*
    * Function : Match sol time with reference sol
    */
    int cOutSol::MatchRefSol(cTime sol_time) {
        double sow1,sow2;
        int i=0;
        bool stat=false;

        for(i=ref_index_-100<0?0:ref_index_-10;i<ref_sols_.size();i++){
            sow1=ref_sols_.at(i).t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            sow2=sol_time.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<DTTOL){
                ref_index_=i;stat=true;break;
            }
            else if((sow1-sow2)>2.0*DTTOL){
                stat=false;break;
            }
        }
        return stat;
    }

    /*
    * Function : Return the difference of position/velocity between sol and reference sol
    */
    // 返回位置和速度的差值
    tSolInfoUnit cOutSol::CompareSol(tSolInfoUnit &sol, tSolInfoUnit &ref_sol) {
        tSolInfoUnit dsol;
        Vector3d ref_blh=Xyz2Blh(ref_sol.pos);
        Vector3d dxyz=sol.pos-ref_sol.pos;
        if(C_.solC.sol_coord==COORD_ENU){
            dsol.pos=Xyz2Enu(ref_blh,dxyz);
        }
        dsol.vel=sol.vel-ref_sol.vel;

        return dsol;
    }

    static double SqrtCovVar(double cov_var){
        return cov_var<0.0?-sqrt(-cov_var):sqrt(cov_var);
    }

    /*
    * Function : The state of out solution
    */
    int cOutSol::OutSolStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos, char *buff) {
        if(sol->stat<=SOL_NONE) return 0;

//        if(C_.mode==MODE_PPP||C_.mode_opt==MODE_OPT_PPP){
//            return 0;
//        }

        char *p=buff;

        int week;
        double tow;
        tow=sol->t_tag.Time2Gpst(&week, nullptr,SYS_GPS);

        // position
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   sol->stat,sol->pos[0],sol->pos[1],sol->pos[2],
                   0.0,0.0,0.0);

        // clock
        p+=sprintf(p,"$CLK,%d,%.3f,%d,%.4f %.4f %.4f %.4f %.4f\n",week,tow,sol->stat,sol->clk_error[0],sol->clk_error[1],sol->clk_error[2],sol->clk_error[3],sol->clk_error[4]);

        // ifb
        if(C_.gnssC.frq_opt==FRQ_TRIPLE){
            p+=sprintf(p,"$IFB,%d,%.3f,%d,%.4f %.4f %.4f %.4f %.4f\n",week,tow,sol->stat,sol->rec_ifcb[0],sol->rec_ifcb[1],sol->rec_ifcb[2],sol->rec_ifcb[3],sol->rec_ifcb[4]);
        }

        // dcb

        // trp
        if(C_.gnssC.trp_opt>=TRP_EST_WET){
            p+=sprintf(p,"$TRP,%d,%.3f,%d,%.4f,%.4f\n",week,tow,sol->stat,sol->zenith_trp_delay[0],sol->zenith_trp_delay[1]);
        }

        return (int)(p-buff);
    }

    /*
    * Function : Write the state of satellite
    */
    void cOutSol::WriteSatStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos) {
        char buff[8191+1];

        int n=OutSolStat(sol,sat_infos,buff);

        fputs(buff,fout_stat_);
        buff[n]='\0';
        if(sol->stat<=SOL_NONE) return;

        int ion_flag=0;
        double ion;
        if(C_.gnssC.ion_opt>=ION_EST){
            ion_flag=1;
        }

        tSatInfoUnit *sat_info= nullptr;
        int nf=C_.gnssC.frq_opt+1;// the number of frq 
        if(C_.gnssC.ion_opt==ION_IF&&C_.gnssC.frq_opt<=FRQ_DUAL) nf=1;
        else if(C_.gnssC.ion_opt==ION_IF&&C_.gnssC.frq_opt==FRQ_TRIPLE) nf=2;
        if(C_.gnssC.ion_opt==ION_IF_DUAL&&C_.gnssC.frq_opt==FRQ_TRIPLE) nf=2;

        int i,j,week;
        double tow,amb;
        tow=sol->t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
        for(i=0;i<MAX_SAT_NUM;i++){// every satellite
            sat_info=&sat_infos[i];
            if(!sat_info->vsat[0]) continue;
            //for(j=0;j<nf;j++){// every frq 
                for (j=0;j<nf;j++) { // every frq
                    fprintf(fout_stat_,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f",
                               week,tow,sat_info->sat.sat_.id.c_str(),j+1,sat_info->el_az[1]*R2D,sat_info->el_az[0]*R2D,
                               sat_info->post_res_P[j],sat_info->post_res_L[j],sat_info->vsat[j],sat_info->raw_S[j]*0.25,
                               sat_info->fix[j],sat_info->slip[j]&3,sat_info->lock[j],sat_info->outc[j],
                               0,sat_info->rejc[j],sat_info->float_amb[j],sat_info->fix_amb[j],sat_info->raw_mw[0],sat_info->sm_mw[0],sat_info->res_wl,sat_info->res_nl);
                    if(ion_flag){
                        fprintf(fout_stat_,",%.4f",sat_info->ion_delay[0]);
                    }
                    fprintf(fout_stat_,"\n");
                }
           // }
        }
    }

    /*
    * 
    */
    int cOutSol::OutEcef(unsigned char *buff, const char *s, tSolInfoUnit &sol) {
        const char *sep=" ";
        char *p=(char *)buff;
        Vector3d blh=Xyz2Blh(sol.pos);

        p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f%s%6.2f",
                   s,sep,sol.pos[0],sep,sol.pos[1],sep,sol.pos[2],sep,sol.stat,sep,
                   sol.valid_sat_num,sep,SQRT(sol.q_pos[0]),sep,SQRT(sol.q_pos[1]),sep,SQRT(sol.q_pos[2]),
                   sep,SqrtCovVar(sol.q_pos[3]),sep,SqrtCovVar(sol.q_pos[4]),sep,SQRT(sol.q_pos[5]),
                   sep,sol.age,sep,sol.ratio,sep,sol.sigma);
        if(C_.mode>=MODE_INS){
            p+=sprintf(p,"%s%4d",sep,sol.ins_stat);
        }

        if(C_.solC.out_vel){
            p+=sprintf(p,"%s%10.5f%s%10.5f%s%10.5f%s%9.5f%s%8.5f%s%8.5f%s%8.5f%s%8.5f%s%8.5f",
                       sep,sol.vel[0],sep,sol.vel[1],sep,sol.vel[2],sep,
                       SQRT(sol.q_vel[0]),sep,SQRT(sol.q_vel[1]),sep,SQRT(sol.q_vel[2]),
                       sep,SqrtCovVar(sol.q_vel[3]),sep,SqrtCovVar(sol.q_vel[4]),sep,
                       SqrtCovVar(sol.q_vel[5]));
        }
        if(C_.solC.out_att&&C_.mode>=MODE_INS){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf",
                       sep,sol.att[0]*R2D,sep,sol.att[1]*R2D,sep,sol.att[2]*R2D<0.0?sol.att[2]*R2D+360:sol.att[2]*R2D,
                       sep,SQRT(sol.q_att[0])*R2D,sep,SQRT(sol.q_att[1])*R2D,sep,SQRT(sol.q_att[2])*R2D);
        }
        if(C_.solC.out_ba&&C_.mode>=MODE_IGLC){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf",sep,sol.accl_bias[0],sep,sol.accl_bias[1],sep,sol.accl_bias[2]);

        }
        if(C_.solC.out_bg&&C_.mode>=MODE_IGLC){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf",sep,sol.gyro_bias[0],sep,sol.gyro_bias[1],sep,sol.gyro_bias[2]);
        }

        p+=sprintf(p,"\n");
        return p-(char*)buff;
    }

    /*
    * Function : Output IGCoupled solution
    * -Args :
    *       tSolInfoUnit sol               I        solution information
    *       tSatInfoUnit epoch_sat_info    I        the satellite information vector
    * -Returns :
    *       
    */
    int cOutSol::IGCoupledOut(tSolInfoUnit *sol, vector<tSatInfoUnit>& epoch_sat_info) {
        Vector3d pre_pos(0,0,0),out_pos(0,0,0),dr(0,0,0);
        bool dgnss=C_.mode==MODE_PPK||C_.mode_opt==MODE_OPT_PPK;

        
        Vector3d blh(0,0,0);
        if(C_.solC.sol_coord==COORD_ENU){ //输出enu下的坐标
            if(dgnss){
                blh=Xyz2Blh(C_.gnssC.rb);
                dr=sol->pos-C_.gnssC.rb;
                out_pos=Xyz2Enu(blh,dr);
                if(C_.solC.out_err_fmt){
                    Vector3d ref_dr=ref_sol_.pos-C_.gnssC.rb;
                    Vector3d ref_pos_enu=Xyz2Enu(blh,ref_dr);
                    out_pos-=ref_pos_enu;// PPK 
                }
            }
            else{
                blh=Xyz2Blh(ref_sol_.pos);
                dr=sol->pos-ref_sol_.pos;
                out_pos=Xyz2Enu(blh,dr);  //PPP 
            }
        }
        else if(C_.solC.sol_coord==COORD_BLH){


        }
        else if(C_.solC.sol_coord==COORD_XYZ){ //输出xyz坐标
            blh=Xyz2Blh(ref_sol_.pos);
            out_pos=sol->pos;
            if(C_.solC.out_err_fmt){
                dr=out_pos-ref_sol_.pos;
                out_pos=Xyz2Enu(blh,dr);
            }
        }

        double wos=0;
        int week=0;

        wos=sol->t_tag.Time2Gpst(&week,nullptr,SYS_GPS);
        char epoch[1024]={'\0'};
        sprintf(epoch,"$EPOCH    %12s %4d %10.2f %5d %d %3d %3d %5.1f",
                sol->t_tag.GetTimeStr(1).c_str(),week,wos,sol->epoch_idx,sol->stat,sol->observed_sat_num,sol->valid_sat_num,sol->dops[1]);
        igcoupled_out_<<epoch<<endl;

        if(sol->stat>SOL_NONE){
            char pos[1024]={'\0'};
            int d=12;
            if(C_.solC.out_err_fmt) d=5;
            sprintf(pos,"$POS   %15.3f %15.3f %15.3f\n",out_pos[0],out_pos[1],out_pos[2]);
            char clk[1024]={'\0'};
            sprintf(clk,"$CLK   %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f\n",sol->clk_error[0],sol->clk_error[1],sol->clk_error[2],sol->clk_error[3],sol->clk_error[4],sol->clk_error[5]);
            char trp[1024]={'\0'};
            sprintf(trp,"$TRP   %15.3f %15.3f\n",sol->zenith_trp_delay[0],sol->zenith_trp_delay[1]);
            igcoupled_out_<<pos;
            igcoupled_out_<<clk;
            igcoupled_out_<<trp;
            if(C_.solC.out_stat){
             IGCoupledOutSat(epoch_sat_info);
            }
            igcoupled_out_<<endl;
        }
        else{
            if(C_.solC.out_stat){
                IGCoupledOutSat(epoch_sat_info);
            }
            igcoupled_out_<<endl;
        }
    }

    /*
    * Function : The state of IGCoupled solution 
    */
    void cOutSol::IGCoupledOutSat(vector<tSatInfoUnit>& epoch_sat_info) {
        int i,sat_no=0;
        tSatInfoUnit sat_info;

        for(i=0;i<epoch_sat_info.size();i++){
            sat_info=epoch_sat_info[i];
            if(sat_info.stat==SAT_NO_USE||sat_info.stat>SAT_USED){ //如果卫星没用，残差设置为0
                for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++){
                    sat_info.post_res_L[j]=sat_info.post_res_P[j]=0.0;
                }
            }

            char info[1024]={'\0'};
            char res[1024]={'\0'};
            char amb[1024]={'\0'};
            char trp[1024]={'\0'};
            char ion[1024]={'\0'};
            char loc[1024]={'\0'};
            // info: id, stat, az, el;
            sprintf(info,"$SAT     %4s %d %6.1f %6.1f",sat_info.sat.sat_.id.c_str(),sat_info.stat,sat_info.el_az[1]*R2D,sat_info.el_az[0]*R2D);

            // res: P1(IF_P),P2,P3(m), L1(IF_L),L2,L3(mm)
            sprintf(res,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",sat_info.post_res_P[0],sat_info.post_res_P[1],
                    sat_info.post_res_P[2],sat_info.post_res_L[0],sat_info.post_res_L[1],sat_info.post_res_L[2]);

            // amb: L1(IF_L), L2, L3, MW, SMW
            sprintf(amb,"%8.3f %8.3f %8.3f %8.3f %8.3f",sat_info.float_amb[0],sat_info.float_amb[1],sat_info.float_amb[2],sat_info.raw_mw[0],sat_info.sm_mw[0]);

            // trp: strp_h,strp_w,map_h,map_w
            sprintf(trp,"%6.3f %6.3f %6.3f %6.3f",sat_info.trp_dry_delay[0],sat_info.trp_wet_delay[0],sat_info.trp_dry_delay[1],sat_info.trp_wet_delay[1]);

            // ion: sion
            double fact=40.30E16/SQR(sat_info.frq[0]);
            sprintf(ion,"%6.3f",sat_info.ion_delay[0]*fact);

            // lock and outc
            sprintf(loc,"%5d %5d",sat_info.lock[0],sat_info.outc[0]);

            igcoupled_out_<<info<<" "<<res<<" "<<amb<<" "<<trp<<" "<<ion<<" "<<loc<<endl;
        }
    }

    /*
    * Function : Initialize the output solution ,open the output file
    */
    bool cOutSol::InitOutSol(tIGCOUPLEDConf C, string file) {
        C_=C;
        if(!C_.solC.sol_fmt){
            igcoupled_out_.open(C_.fileC.sol,ios::out);
        }
        else{
            fout_=fopen(file.c_str(),"w");
            if(C_.solC.out_stat&&!C_.fileC.sol_stat.empty()){
                fout_stat_=fopen(C_.fileC.sol_stat.c_str(),"w");
            }
        }

    }

    /*
    * Function : Write the head of output file
    */
    void cOutSol::WriteHead() {
        if(!C_.solC.out_head) return;

        if(C_.solC.sol_fmt){//rtklib formate sol file
            const char *sep=" ";
            fprintf(fout_,"%s  %-*s%s%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s%s%6s",
                    COMMENTH,20,"GPST",sep,"x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
                    "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,
                    "sdyz(m)",sep,"sdzx(m)",sep,"age(s)",sep,"ratio",sep,"sigma");

            if(C_.mode>=MODE_INS){
                fprintf(fout_,"%s%4s",sep,"QINS");
            }
            if(C_.solC.out_vel){
                fprintf(fout_,"%s%10s%s%10s%s%10s%s%9s%s%8s%s%8s%s%8s%s%8s%s%8s",
                        sep,"vn(m/s)",sep,"ve(m/s)",sep,"vu(m/s)",sep,"sdvn",sep,
                        "sdve",sep,"sdvu",sep,"sdvne",sep,"sdveu",sep,"sdvun");
            }
            if(C_.solC.out_att&&C_.mode>=MODE_INS){
                fprintf(fout_,"%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s",
                        sep,"roll(deg)",sep,"pitch(deg)",sep,"yaw(deg)",
                        sep,"sdroll",sep,"sdpitch",sep,"sdyaw");
            }
            if(C_.solC.out_ba&&C_.mode>=MODE_IGLC){
                fprintf(fout_,"%s%10s%s%10s%s%10s",sep,"bax(m/s^2)",sep,"bay(m/s^2)",sep,"baz(m/s^2)");
            }
            if(C_.solC.out_bg&&C_.mode>=MODE_IGLC){
                fprintf(fout_,"%s%10s%s%10s%s%10s",sep,"bgx(rad/s)",sep,"bgy(rad/s)",sep,"bgz(rad/s)");
            }
            fprintf(fout_,"\n");
        }
        else{// IGCoupled sol  format 
            string comm="%";
            string use_sys,clk;

            if(C_.gnssC.nav_sys&SYS_GPS){
                use_sys+="GPS ";
            }
            if(C_.gnssC.nav_sys&SYS_BDS){
                if(C_.gnssC.est_bd3_isb){
                    use_sys+="BD2 BD3 ";
                }
                else{
                    use_sys+="BDS ";
                }
            }
            if(C_.gnssC.nav_sys&SYS_GAL){
                use_sys+="GAL ";
            }
            if(C_.gnssC.nav_sys&SYS_GLO){
                use_sys+="GLO ";
            }
            if(C_.gnssC.nav_sys&SYS_QZS){
                use_sys+="QZS ";
            }

            igcoupled_out_<<"+ IGCoupled HEADER"<<endl;
            igcoupled_out_<<"  Station : "<<endl;
            igcoupled_out_<<"  Rec Type: "<<endl;
            igcoupled_out_<<"  Ant Type: "<<endl;
            if(C_.mode!=MODE_PPK){
                igcoupled_out_<<"  Sta  Pos: "<<endl;
            }
            igcoupled_out_<<"  Nav  Sys: "<<use_sys<<endl;
            igcoupled_out_<<"  Trp Mode: "<<endl;
            igcoupled_out_<<"  Ion Mode: "<<endl;
            igcoupled_out_<<"  Sat  Eph: "<<endl;
            igcoupled_out_<<"  $EPOCH  : "<<"Y/M/D H:M:S | WEEK | WOS | Epoch | SolStat | ObsNum | ValidSat | PDOP | SIGMA0 "<<endl;
            igcoupled_out_<<"  $POS    : "<<endl;
            igcoupled_out_<<"  $CLK    : "<<"Clock | ISB"<<endl;
            igcoupled_out_<<"  $TRP    : "<<"Hydrostatic | Wet "<<endl;
            igcoupled_out_<<"  $SAT    : "<<"ID | Stat | Obs | Az | El | Res_P1 | Res_P2 | Res_P3 | Res_L1 | Res_L2 | Res_L3 | Amb1 | Amb2 | Amb3 | MwAmb | sMwAmb | sTrp_h | sTrp_w | map_h | map_w | sIon | Lock | Outc"<<endl;
//            if(C_.gnssC.ion_opt!=ION_IF) igcoupled_out_<<"| Ion "<<endl;
//            else igcoupled_out_<<endl;
            igcoupled_out_<<"- IGCoupled HEADER"<<endl;
            igcoupled_out_<<endl;
        }
    }

    /*
    * Function : Write the solution in output file 
    */
    void cOutSol::WriteSol(tSolInfoUnit sol,vector<tSatInfoUnit>& epoch_sat_info) {
        tSolInfoUnit dsol;
        dsol=sol;

        unsigned char p[8191+1];
        char s[256];
        int n=0;

        if(C_.solC.sol_fmt){ // rtklib format
            switch(C_.solC.sol_coord){
                case COORD_BLH:
                case COORD_XYZ:
                case COORD_ENU:
                    n=OutEcef(p, sol.t_tag.GetTimeStr(3).c_str(),dsol);
            }
            fwrite(p,n,1,fout_);
        }
        else{  // igcoupled format
            IGCoupledOut(&sol,epoch_sat_info);
        }
    }

    /*
    * Function : Write the head of IMU obs file 
    */
    void cOutSol::WriteImuHead() {
        if(!C_.solC.out_head) return;
        const char *sep=" ";
        fprintf(fout_,"%s %4s%s%14s%s%14s%s%14s%s%14s%s%14s%s%14s%s%14s",
                COMMENTH,"GPSW",sep,"SOW",sep,"x-acce(m/s^2)",sep,"y-acce(m/s^2)",sep,"z-acce(m/s^2)",sep,"x-gyro(rad/s)",sep,"y-gyro(rad/s)",sep,"z-gyro(rad/s)");
        fprintf(fout_,"\n");
    }

    /*
    * Function : Write the body of IMU obs file 
    */
    void cOutSol::WriteImuObs() {
        int n=imus->data_.size();
        cTime t_tag;
        tImuDataUnit data;
        const char *sep=" ";
        int week;
        double sow;

        for(int i=0;i<n;i++){
            data=imus->data_[i];
            sow=data.t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
            fprintf(fout_,"%6d%s%14.4f%s%14.5f%s%14.5f%s%14.5f%s%14.5f%s%14.5f%s%14.5f\n",
                    week,sep,sow,sep,data.acce[0],sep,data.acce[1],sep,data.acce[2],sep,data.gyro[0],sep,data.gyro[1],sep,data.gyro[2]);
        }
        fclose(fout_);
    }
}
