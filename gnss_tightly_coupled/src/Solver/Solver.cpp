//
// Created by wlzhang on 7/5/22.
// 

#include "Solver.h"
#include "ReadFiles.h"
#include "DecodeRaw.h"

namespace IGCoupled{
    cSolver::cSolver() {
        epoch_idx_=0;
    }

    cSolver::~cSolver() {}

    int cSolver::GnssObsRes(int post, tIGCOUPLEDConf C, double *x) {
        // ????
    }

    /*
    * Function : Combine forward and backward solution
    */
    void cSolver::CombFbSol(tIGCOUPLEDConf C) {
        int i,j,k,num_solb=solb_.size(),num_solf=solf_.size();
        double tt;
        tSolInfoUnit com_sol={0};

        Matrix3d Qf,Qb;
        Qf=Matrix3d::Zero(),Qb=Matrix3d::Zero();
        for(i=0,j=num_solb-1;i<num_solf&&j>=0;i++,j--){// for every forward and backward solution 
            if((tt=solf_[i].t_tag.TimeDiff(solb_[j].t_tag.t_))<-DTTOL){
                com_sol=solf_[i];
                j++;
            }
            else if(tt>DTTOL){
                com_sol=solb_[j];
                i--;
            }
            else if(solf_[i].stat<solb_[j].stat){
                com_sol=solb_[j];
            }
            else{
                com_sol=solf_[i];
                com_sol.t_tag+=-tt/2.0;

                for(k=0;k<3;k++){
                    Qf(k,k)=solf_[i].q_pos[k];
                    Qb(k,k)=solb_[j].q_pos[k];
                }
                Qf(1,0)=Qf(0,1)=solf_[i].q_pos[3];
                Qf(2,1)=Qf(1,2)=solf_[i].q_pos[4];
                Qf(2,0)=Qf(0,2)=solf_[i].q_pos[5];
                Qb(1,0)=Qb(0,1)=solb_[i].q_pos[3];
                Qb(2,1)=Qb(1,2)=solb_[i].q_pos[4];
                Qb(2,0)=Qb(0,2)=solb_[i].q_pos[5];

                VectorXd pos=com_sol.pos;
                MatrixXd Qs;
                if(Smoother(solf_[i].pos,Qf,solb_[j].pos,Qb,pos,Qs,3)) continue; // 前后向滤波平滑
                for(int n=0;n<3;n++) com_sol.pos[n]=pos[n];
                com_sol.q_pos[0]=Qs(0,0);
                com_sol.q_pos[1]=Qs(1,1);
                com_sol.q_pos[2]=Qs(2,2);
                com_sol.q_pos[3]=Qs(1,0);
                com_sol.q_pos[4]=Qs(2,1);
                com_sol.q_pos[5]=Qs(2,0);

                if(C.solC.out_sol) out_->WriteSol(com_sol,epoch_sat_info_collect_);
                if(C.solC.out_stat&&C.solC.sol_fmt) out_->WriteSatStat(&com_sol,previous_sat_info_);
            }
        }
    }

    /*
    * Function : Forward and backward filtering and smoothing 
    * -Args :
    *       VectorXd xf                       I             forward filtering position solution
    *       MatrixXd Qf                       I             forward filtering cofactor matrix
    *       VectorXd xb                       I             backward filtering position solution
    *       MatrixXd Qb                       I             backward filtering cofactor matrix
    *       VectorXd xs                       O             the position after smoothing 
    *       MatrixXd Qs                       O             the cofactor after smoothing 
    *       int      n                        I             the number of position solution
    * -Returns :
    *       int info (0:success      -1 : failed)
    */
    int cSolver::Smoother(const VectorXd xf,const MatrixXd Qf,const VectorXd xb,const MatrixXd Qb,VectorXd& xs,MatrixXd& Qs,int n) {
        MatrixXd Qf_inv,Qb_inv;
        Qf_inv=MatrixXd::Zero(n,n),Qb_inv=MatrixXd::Zero(n,n);
        Qf_inv=Qf;Qb_inv=Qb;
        int info=-1;

        if(!MatInv(Qf_inv.data(),n)&&!MatInv(Qb_inv.data(),n)){
            Qs=Qf_inv+Qb_inv;// Qs = Qf^ + Qb ^
            if(!(info=MatInv(Qs.data(),n))){
                xs=Qs*(Qf_inv*xf+Qb_inv*xb); 
            }
        }

        return info;
    }

    /*
    * Function : Get GNSS singal index
    * -Args : 
    *        tSatObsUnit sat_obs            I        The satellite observation
    *        int         *f                 IO       The index of frq
    * -Returns:
    *        1 : ok  0: failed
    */
    int cSolver::GetSingalInd(tSatObsUnit& sat_obs,int *f){//这里是一个指针变量，当匹配成功后，*f的频率值变为移动站该频率的地址
        const tGnssSignal *ps= nullptr;
        int n,p;

        if(sat_obs.L[*f]==0.0||sat_obs.P[*f]==0.0){
            ps=rover_obs_.signal_[sat_obs.sat.sat_.sys_idx];
            for(n=0,p=-1;n<ps->n;n++){//移动站观测值卫星系统的所有GNSS信号
                if(ps->frq[n]==*f+1&&ps->pri[n]&&(p<0||ps->pri[n]>ps->pri[p])&&sat_obs.L[ps->pos[n]]&&sat_obs.P[ps->pos[n]]){
                    p=n;
                }
            }
            if(p<0) return 0;
            *f=ps->pos[p];

            return 1;
        }
        else return 0;
    }

    /*
    *  Function : Adjust GNSS obervation frequency
    */
    int cSolver::AdjustObsInd(tSatObsUnit& sat_obs,int *i, int *j, int *k) {// i、j、k为三个频率的地址
        int info=0;
        info|=i&&GetSingalInd(sat_obs,i); // 做或位运算，并把值付给info。即有一个频率匹配成功info即不为0
        info|=j&&GetSingalInd(sat_obs,j);
        info|=k&&GetSingalInd(sat_obs,k);
        return info;
    }

    /*
    * Function : Initialize satellites informations for all epoch
    */
    void cSolver::InitEpochSatInfo(vector<tSatInfoUnit> &sat_infos) {
        for(int j=0;j<MAX_SAT_NUM;j++){//遍历每颗卫星
            for(int k=0;k<MAX_GNSS_USED_FRQ_NUM;k++){//遍历每个使用的频率
                if(epoch_idx_==1) {//如果是第一个历元，初始化前一个历元的信息
                    previous_sat_info_[j].stat=SAT_NO_USE;
                    previous_sat_info_[j].rejc[k]=0;
                    previous_sat_info_[j].rejc_phase[k]=0;
                }
                previous_sat_info_[j].vsat[k]=0;
                previous_sat_info_[j].outc[k]++;
                previous_sat_info_[j].p_var_factor[k]=previous_sat_info_[j].c_var_factor[k]=1.0;
                if(k<2&&epoch_idx_==1){
                    previous_sat_info_[j].sm_mw[k]=0.0;
                    previous_sat_info_[j].mw_idx[k]=0.0;
                    previous_sat_info_[j].var_mw[k]=0.0;
                    previous_sat_info_[j].gf[k]=0.0;
                    previous_sat_info_[j].phase_wp=0.0;
                    previous_sat_info_[j].lc_amb={0};
                }
            }
        }

        for(int j=0;j<sat_infos.size();j++){//遍历每个历元
            for(int k=0;k<MAX_GNSS_USED_FRQ_NUM;k++){//遍历每个频率
//                sat_infos.at(j).slip[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].slip[k];
                sat_infos.at(j).outc[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].outc[k];
                sat_infos.at(j).lock[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].lock[k];
                sat_infos.at(j).rejc[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].rejc[k];
                sat_infos.at(j).rejc_phase[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].rejc_phase[k];
                if(previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].rejc[k]>2) sat_infos.at(j).rejc[k]=0;
                if(previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].rejc_phase[k]>2) sat_infos.at(j).rejc_phase[k]=0;

            }
            for(int k=0;k<2;k++){
                sat_infos.at(j).sm_mw[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].sm_mw[k];
                sat_infos.at(j).mw_idx[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].mw_idx[k];
                sat_infos.at(j).var_mw[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].var_mw[k];
                sat_infos.at(j).gf[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].gf[k];
                sat_infos.at(j).phase_wp=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].phase_wp;
            }
            sat_infos.at(j).lc_amb=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].lc_amb;
        }
    }

    /*
    * Function : Update satellite info for epochs 
    */
    void cSolver::UpdateSatInfo(vector<tSatInfoUnit> &sat_infos) {
        for(int i=0;i<sat_infos.size();i++){//遍历每个历元
            for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++) sat_infos.at(i).slip[j]=previous_sat_info_[sat_infos.at(i).sat.sat_.no-1].slip[j];
            previous_sat_info_[sat_infos.at(i).sat.sat_.no-1]=sat_infos.at(i);//sat 信息更新
        }
    }

    /*
    *  Function : Update gnss obervation informations
    * -Args :
    *       tIGCOUPLEDConf      C                  I           configure informations
    *       tEpochSatUnit  epoch_sat_obs      I           satellite observation for a epoch
    *       RECEIVER_INDEX rec                I           receiver index
    * -Return :
    *       None
    */
    void cSolver::UpdateGnssObs(tIGCOUPLEDConf C,tEpochSatUnit& epoch_sat_obs,RECEIVER_INDEX rec) {
        int sys,i,j,k;
        int *frqs;
        int f1,f2,f3;
        bool ppp=C.mode==MODE_PPP||C.mode_opt==MODE_OPT_PPP;

        for(i=0;i<epoch_sat_obs.sat_num;i++){//该历元的所有卫星遍历

            tSatInfoUnit sat_info={0};
            sys=epoch_sat_obs.epoch_data.at(i).sat.sat_.sys;//GNSS系统
            sat_info.t_tag=epoch_sat_obs.obs_time;//历元时刻
            sat_info.sat=epoch_sat_obs.epoch_data.at(i).sat;
            if(sat_info.sat.sat_.no==4) continue;
            if(sat_info.sat.sat_.no==50) continue;
            if(sat_info.sat.sat_.no>=33&&sat_info.sat.sat_.no<=37) continue;
            if(!C.gnssC.use_bd3&&sat_info.sat.sat_.sys==SYS_BDS&&sat_info.sat.sat_.prn>18) continue;
            sat_info.stat=SAT_USED;
            switch(sys){
                case SYS_BDS:
                    if(sat_info.sat.sat_.prn>18){
                        frqs=C.gnssC.gnss_frq[NSYS];break;//BDS-3
                    }
                    frqs=C.gnssC.gnss_frq[SYS_INDEX_BDS];break;//BDS-2
                case SYS_GAL: frqs=C.gnssC.gnss_frq[SYS_INDEX_GAL];break;
                case SYS_GLO: frqs=C.gnssC.gnss_frq[SYS_INDEX_GLO];break;
                case SYS_QZS: frqs=C.gnssC.gnss_frq[SYS_INDEX_QZS];break;
                default:
                    frqs=C.gnssC.gnss_frq[SYS_INDEX_GPS];break;//默认GPS
            }

            //first frequency
            f1=frqs[0];f2=frqs[1];f3=frqs[2];
            if(C.gnssC.adj_obs) AdjustObsInd(epoch_sat_obs.epoch_data.at(i),&f1,&f2,&f3);//调整观测值频率index

            if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),0,frqs[0],f1,nav_.glo_frq_num);
            else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),0,frqs[0],f1,nullptr);

            sat_info.c_var_factor[0]=1.0;
            sat_info.p_var_factor[0]=1.0;

            //second frequency
            if(C.gnssC.frq_opt==FRQ_DUAL){
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nullptr);
                sat_info.c_var_factor[1]=1.0;
                sat_info.p_var_factor[1]=1.0;
                if(ppp){//如果进行PPP
                    if(sat_info.raw_L[0]*sat_info.raw_L[1]==0.0) sat_info.stat=SAT_NO_USE;
                    if(fabs(sat_info.raw_P[0]-sat_info.raw_P[1])>=200.0) sat_info.stat=SAT_NO_USE;
                }
            }
            // third frequency
            else if(C.gnssC.frq_opt==FRQ_TRIPLE){
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nullptr);
                sat_info.c_var_factor[1]=1.0;
                sat_info.p_var_factor[1]=1.0;
                if(sat_info.raw_L[0]*sat_info.raw_L[1]==0.0) sat_info.stat=SAT_NO_USE;
                if(fabs(sat_info.raw_P[0]-sat_info.raw_P[1])>=200.0) sat_info.stat=SAT_NO_USE;

                f3=frqs[2];
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),2,frqs[2],f3,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),2,frqs[2],f3,nullptr);
                sat_info.c_var_factor[2]=1.0;
                sat_info.p_var_factor[2]=1.0;
            }
            //判断接收机是否是移动站
            rec==REC_ROVER?epoch_sat_info_collect_.push_back(sat_info):base_sat_info_collect_.push_back(sat_info);
        }
    }

    /*
    * Function : Correct GNSS observations
    * -Args :
    *       tIGCOUPLEDConf C                    I              configure informations
    *       Vector3d  rr                   I              The position of receiver
    * -Return:
    *       None
    */
    void cSolver::CorrGnssObs(tIGCOUPLEDConf C, Vector3d& rr) {

        tSatInfoUnit* sat_info = nullptr;
        double pcv_dants[MAX_GNSS_USED_FRQ_NUM] = { 0 }, dantr[MAX_GNSS_USED_FRQ_NUM] = { 0 };
        for (int i = 0; i < epoch_sat_info_collect_.size(); i++) {//遍历该历元所有卫星
            sat_info = &epoch_sat_info_collect_.at(i);

            gnss_err_corr_.cbia_model_.InitSatInfo(sat_info, nullptr);//初始化
            gnss_err_corr_.cbia_model_.GetCodeBias();//获取码偏差
            gnss_err_corr_.cbia_model_.UpdateSatInfo();//更新信息

            if (C.mode == MODE_PPP) {//如果进行PPP
                if (sat_info->el_az[0] * R2D < C.gnssC.ele_min) continue;
                if (C.gnssC.sat_pcv) gnss_err_corr_.ant_model_.SatPcvCorr(sat_info, rr, pcv_dants);//卫星端PCV 改正
                if (C.gnssC.rec_ant) gnss_err_corr_.ant_model_.RecAntCorr(sat_info, dantr, REC_ROVER); //接收机端天线改正
                gnss_err_corr_.PhaseWindUp(*sat_info, rr); // 卫星相位缠绕改正
            }


            for (int j = 0; j < MAX_GNSS_USED_FRQ_NUM; j++) {//遍历每个频率
                if (sat_info->raw_P[j] == 0.0) continue;
                sat_info->cor_P[j] = sat_info->raw_P[j] - sat_info->code_bias[j] - dantr[j] - pcv_dants[j];//伪距观测值改正
                if (sat_info->raw_L[j] == 0.0) continue;
                sat_info->cor_L[j] = sat_info->raw_L[j] * sat_info->lam[j] - dantr[j] - pcv_dants[j] - sat_info->phase_wp * sat_info->lam[j];//相位观测值改正
            }

            //BDS-2 PPP
            if (sat_info->sat.sat_.sys == SYS_BDS && sat_info->sat.sat_.prn <= 18 && (C.mode == MODE_PPP || C.mode_opt == MODE_OPT_PPP)) {
                gnss_err_corr_.BD2MultipathModel(sat_info);//BDS-2 多路径
                for (int j = 0; j < 3; j++) {
                    if (sat_info->raw_P[j] == 0.0) continue;
                    sat_info->cor_P[j] -= sat_info->bd2_mp[j];//改正多路径
                }
            }
            //双频 无电离层组合
            if (C.gnssC.ion_opt == ION_IF || C.gnssC.ion_opt == ION_IF_DUAL)
                gnss_obs_operator_.MakeGnssObsComb(C, COMB_IF, sat_info, previous_sat_info_[sat_info->sat.sat_.no - 1]);//组成无电离层观测值
        }
    }
    
    /*
    * Function : Read files(obs, brd, dcb)
    */
    bool cSolver::InitReader(tIGCOUPLEDConf C) {
		if((C.mode == MODE_PPP|| C.mode_opt == MODE_OPT_PPP) && C.gnssC.ar_mode==AR_PPP_AR && C.gnssC.ar_prod == AR_PROD_OSB_WUH){//读取 WUM osb文件 zwl 2022.4.7
            if(C.fileC.osb.empty()){
                CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<"NO WHU OSB FILE"<<endl;
                return false;
            }
            cReadOsb *osb_reader;
            osb_reader = new cReadOsb(C.fileC.osb,nav_);
            osb_reader->Reading();
        }

        if(!C.fileC.rover.empty()){// 移动站观测文件存在
            rover_obs_.ResetGnssObs();//重置rover_obs_
            cReadGnssObs *rover_reader;
            rover_reader=new cReadGnssObs(C.fileC.rover,nav_,rover_obs_,REC_ROVER);//读取移动站的观测文件
            rover_reader->SetGnssSysMask(C.gnssC.nav_sys);
            rover_reader->Reading();
        }

        if(!C.fileC.brd.empty()&&C.site_idx==1){//广播星历文件存在
            cReadGnssBrdEph brd_reader(C.fileC.brd, nav_);//读取广播星历
            brd_reader.SetGnssSysMask(C.gnssC.nav_sys);
            brd_reader.Reading();
        }

        // code dcb
		if((C.mode == MODE_PPP|| C.mode_opt == MODE_OPT_PPP) && C.gnssC.ar_mode==AR_PPP_AR && C.gnssC.ar_prod == AR_PROD_OSB_WUH){

        }else {
        	bool no_code_dcb_flag=true;
       		if(!C.fileC.cod_dcb.empty()&&C.site_idx==1){//dcb文件存在
            	int i,j,n,mon=C.prc_date.GetEpoch()[1],m;
            	char *ex_files[36]={nullptr};
            	string file_name,sep;
#if WIN32
            	sep="//.";
#else
            	sep="/.";
#endif
            	for(i=0;i<36;i++){
                	if(!(ex_files[i]=(char *)malloc(1024))){
                    	for(i--;i>=0;i--) free(ex_files[i]);
                    	return false;
                	}
            	}
            	n=ExPath((char *)C.fileC.cod_dcb.c_str(),ex_files,36);

            	if(n>0){
                	cReadGnssCodeBias dcb_reader(C.fileC.cod_dcb, nav_,0);//读取dcb文件
                	for(i=0;i<n;i++){
#if WIN32
                    	vector<string> splits=MultiSplitStr(ex_files[i],sep);
                	string a=((splits[0]).substr(splits[0].length()-2,2));
                	m=atoi(a.c_str());
                	if(mon!=m) continue;
#else
                    	vector<string> splits=MultiSplitStr(ex_files[i],sep);
                    	m=atoi(((splits.end()-2)->substr(6,2)).c_str());
                    	if(mon!=m) continue;
#endif
                    	no_code_dcb_flag= false;
                    	dcb_reader.file_=ex_files[i];
                    	dcb_reader.Reading();
                	}

                	for(i=0;i<36;i++) free(ex_files[i]);
            	}
        	}

        	// cas dcb
        	if((!C.fileC.cbias.empty()&&C.site_idx==1)){
            	cReadGnssCodeBias cbias_reader(C.fileC.cbias, nav_,1);
            	cbias_reader.no_code_dcb_=no_code_dcb_flag;
            	cbias_reader.Reading();
        	}
		}
		CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INIT FILE READER SUCCESS" <<endl;
        return true;
   	}


    void cSolver::InitSolver(tIGCOUPLEDConf C) {
        para_=cParSetting(C);
    }

    bool cSolver::SolverProcess(tIGCOUPLEDConf C,int idx) {}

    bool cSolver::SolverStart(int i,int idx) {}

    bool cSolver::SolverEpoch() {}

    bool cSolver::Estimator(tIGCOUPLEDConf C) {}

    bool cSolver::SolutionUpdate() {}

    void cSolver::ReinitSolver(tIGCOUPLEDConf C) {
        para_=cParSetting(C);
        num_full_x_=para_.GetIGCOUPLEDPar(C);//获取待估参数个数
        full_x_=VectorXd::Zero(num_full_x_);//待估参数向量
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);//待估参数权矩阵

        num_real_x_fix_=para_.GetRealFixParNum(C);
        real_x_fix_=VectorXd::Zero(num_real_x_fix_);
        real_Px_fix_=MatrixXd::Zero(num_real_x_fix_,num_real_x_fix_);
    }

    /*
    *  Function : Initialize P-matrix using configure file
    * -Args :
    *        double   unc                I                the value from configure file
    *        double   unc0               I                value (default )
    *        int      row_s              I                the start row index  of p-matrix
    *        int      col_s              I                the start column index of p-matrix
    *        MatrixXd P                  O                P-matrix after initializing
    *        int      is                 I                the start index of parameter for estimated
    *        int      ns                 I                the number of this parameter
    *        int      nx                 I                the number of total parameters
    * Return :
    *        None
    */
    static void InitP(double unc,double unc0,int row_s,int col_s,MatrixXd& P,int is,int ns,int nx){
        double q=unc==0.0?SQR(unc0):SQR(unc);//unc是配置文件读出值，unc0是默认值。如果读出值为0则用默认值，否则用读出值
//        Vector3d vec(q,q,q);
//        P.block<3,3>(row_s,col_s)=vec.asDiagonal();
        int i,j;
        for(i=is;i<is+ns;i++) for(j=0;j<nx;j++){
            if(j==i) P.data()[j+i*nx]=q;//对角线
            else P.data()[j+i*nx]=P.data()[i+j*nx]=0.0;//非对角线为0
        }
    }

    /*
    *  Function : Initialize Px-matrix using function : InitP()
    */
    void cSolver::InitInsPx(tIGCOUPLEDConf C,int nx,MatrixXd& Px) {
        Px=MatrixXd::Zero(nx,nx);
        Vector3d vec(0,0,0);
        //the start index of parameters
        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();
        int isa=para_.IndexSa();
        int isg=para_.IndexSg();
        int ira=para_.IndexRa();
        int irg=para_.IndexRg();
        int ilev=para_.IndexLever();

        if(para_.NumPos()>0){
            InitP(C.insC.init_pos_unc,UNC_POS,ip,ip,Px,ip,3,nx);//用配置文件初始化权矩阵
        }
        if(para_.NumVel()>0){
            InitP(C.insC.init_vel_unc,UNC_VEL,iv,iv,Px,iv,3,nx);
        }
        if(para_.NumAtt()>0){
            InitP(C.insC.init_att_unc,UNC_ATT,ia,ia,Px,ia,3,nx);
        }
        if(para_.NumBa()>0){
            InitP(C.insC.init_ba_unc,UNC_BA,iba,iba,Px,iba,3,nx);
        }
        if(para_.NumBg()>0){
            InitP(C.insC.init_bg_unc,UNC_BG,ibg,ibg,Px,ibg,3,nx);
        }
        if(para_.NumSa()>0){
            InitP(C.insC.init_sa_unc,UNC_SA,isa,isa,Px,isa,3,nx);
        }
        if(para_.NumSg()>0){
            InitP(C.insC.init_sg_unc,UNC_SG,isg,isg,Px,isg,3,nx);
        }
        if(para_.NumRa()>0){
            InitP(C.insC.init_ra_unc,UNC_RA,ira,ira,Px,ira,6,nx);
        }
        if(para_.NumRg()>0){
            InitP(C.insC.init_rg_unc,UNC_RA,irg,irg,Px,irg,6,nx);
        }
        if(para_.NumLever()>0){
            InitP(C.insC.init_lever_unc,UNC_LEVER,ilev,ilev,Px,ilev,3,nx);
        }
    }

    /*
    * Function : Initialize the paramters vector
    * -Args :
    *       double xi                  I          the initial value of parameters
    *       double var                 I          the initial value of variance
    *       int    idx                 I          the start index of parameter
    *       double *x                  O          the parameters vector after initializing
    *       double *p                  O          the co-variance matrix after initializing
    */
    void cSolver::InitX(double xi, double var, int idx,double *x,double *p) {
        x[idx]=xi;
        for(int j=0;j<num_full_x_;j++){
//            full_Px_(idx,j)=full_Px_(j,idx)=idx==j?var:0.0;
            p[idx+j*num_full_x_]=p[j+idx*num_full_x_]=idx==j?var:0.0;
        }
    }

    /*
    * Function : Set the Q-matrix using psd(power spectrum density)
    * -Args : 
    *        double   psd                 I             power spectrum density
    *        double   dt                  I             time interval
    *        int      row_s               I             the start row index of this parameter
    *        int      col_s               I             the start column index of this paramter
    *        MatrixXd &Q                  Q             the Q-matrix (noise matrix) after initializing 
    * -Return :
    *        None
    */
    static void SetPsd(double psd,double dt,int row_s,int col_s,MatrixXd& Q){
        Vector3d vec((psd*dt),(psd*dt),(psd*dt));
        Q.block<3,3>(row_s,col_s)=vec.asDiagonal();
    }

    /*
    * Function : Initialize Q-matrix
    */
    Eigen::MatrixXd cSolver::InitQ(tIGCOUPLEDConf C,double dt,int nx) {
        Eigen::MatrixXd Q;
        Q=MatrixXd::Zero(nx,nx);

        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();

        SetPsd(C.insC.psd_acce,dt,iv,iv,Q);
        SetPsd(C.insC.psd_gyro*D2R,dt,ia,ia,Q);
        SetPsd(C.insC.psd_ba,dt,iba,iba,Q);
        SetPsd(C.insC.psd_bg*D2R,dt,ibg,ibg,Q);

        return Q;
    }

    /*
    * 
    */
    Eigen::MatrixXd cSolver::InitPrecQ(tIGCOUPLEDConf C, double dt, int nx,Matrix3d Cbe) {
        MatrixXd Q, G;
        int nprn=15;
        Q=MatrixXd::Zero(nprn,nprn);
        G=MatrixXd::Zero(nx,nprn);

        int iv=3;
        int ia=6;
        int iba=9;
        int ibg=12;
        int i;

        for(i=iv;i<iv+3;i++) Q.data()[i+i*nprn]=C.insC.psd_acce*fabs(dt);
        for(i=ia;i<ia+3;i++) Q.data()[i+i*nprn]=C.insC.psd_gyro*fabs(dt);
        for(i=iba;i<iba+3;i++) Q.data()[i+i*nprn]=C.insC.psd_ba*fabs(dt);
        for(i=ibg;i<ibg+3;i++) Q.data()[i+i*nprn]=C.insC.psd_bg*fabs(dt);

        G.block<3,3>(iv,iv)=-Cbe;
        G.block<3,3>(ia,ia)=Cbe;
        G.block<3,3>(iba,iba)=Matrix3d::Identity();
        G.block<3,3>(ibg,ibg)=Matrix3d::Identity();

//        cout<<G.transpose()<<endl<<endl;

//        cout<<Q<<endl<<endl;

        Q=G*Q*G.transpose();

//        cout<<Q<<endl<<endl;
        return Q;
    }

    /*
    * Function : Remove lever between ins and gnss
    * -Args :
    *       tImuInfoUnit imu_info            I              imu informations
    *       Vector3d     lever               I              the distance between gnss antenna and imu center
    *       Vector3d     gnss_re             O              corrected imu position
    *       Vector3d     gnss_ve             O              corrected imu velocity 
    * -Return :
    *       None
    */
    void cSolver::RemoveLever(const tImuInfoUnit &imu_info, Vector3d &lever, Vector3d &gnss_re,
                                    Vector3d &gnss_ve) {
        Matrix3d Cbe=imu_info.Cbe;
        Vector3d wiee(0.0,0.0,OMGE_GPS);
        Vector3d T=Cbe*lever;

        // position correction
        gnss_re=imu_info.re+T;

        //velocity correction
        Vector3d omge=imu_info.cor_gyro,W;
        T=Cbe*VectorSkew(omge)*lever;
        W=VectorSkew(wiee)*Cbe*lever;
        gnss_ve=imu_info.ve+T-W;
    }

    /*
    * Function : Close loop state correction 
    * -Args : 
    *       VectorXd     x              I            Estimated parameters 
    *       tImuInfoUnit imu_info_corr  IO           The imu informations after close loop correction
    * -Return :
    *        None
    */
    void cSolver::CloseLoopState(VectorXd& x,tImuInfoUnit* imu_info_corr) {//闭环修正
        char buff[1024]={'\0'};
        int ip=para_.IndexPos();
        imu_info_corr->re[0]-=x[ip+0];
        imu_info_corr->re[1]-=x[ip+1];
        imu_info_corr->re[2]-=x[ip+2];
        imu_info_corr->rn=Xyz2Blh(imu_info_corr->re);
        sprintf(buff,"%10.3f - %5.3f, %10.3f - %5.3f  %10.3f - %5.3f",
                imu_info_corr->re[0],x[ip],imu_info_corr->re[1],x[ip+1],imu_info_corr->re[2],x[ip+2]);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<imu_info_corr->t_tag.GetTimeStr(1)<<" CLOSE LOOP STATE: ";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"POSITION: "<<buff;
        buff[0]='\0';

        int iv=para_.IndexVel();
        imu_info_corr->ve[0]-=x[iv+0];
        imu_info_corr->ve[1]-=x[iv+1];
        imu_info_corr->ve[2]-=x[iv+2];
        Matrix3d Cne=CalcCen(imu_info_corr->rn,COORD_NED).transpose();
        imu_info_corr->vn=Cne.transpose()*imu_info_corr->ve;
        sprintf(buff,"%10.3f - %5.3f, %10.3f - %5.3f  %10.3f - %5.3f",
                imu_info_corr->ve[0],x[iv],imu_info_corr->ve[1],x[iv+1],imu_info_corr->ve[2],x[iv+2]);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"VELOCITY: "<<buff;
        buff[0]='\0';

        int ia=para_.IndexAtt();
        if(x[ia]!=DIS_FLAG){
            Vector3d att(x[ia],x[ia+1],x[ia+2]);
            Matrix3d T=Matrix3d::Identity()-VectorSkew(att);//平台失准角误差
            imu_info_corr->Cbe=T*imu_info_corr->Cbe;
        }

        int iba=para_.IndexBa();
        if(x[iba]!=DIS_FLAG){
            imu_info_corr->ba[0]-=x[iba+0];
            imu_info_corr->ba[1]-=x[iba+1];
            imu_info_corr->ba[2]-=x[iba+2];
            sprintf(buff,"%10.3f - %5.3f, %10.3f - %5.3f  %10.3f - %5.3f",
                    imu_info_corr->ba[0],x[iba],imu_info_corr->ba[1],x[iba+1],imu_info_corr->ba[2],x[iba+2]);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"Ba: "<<buff;
            buff[0]='\0';
        }

        int ibg=para_.IndexBg();
        if(x[ibg]!=DIS_FLAG){
            imu_info_corr->bg[0]+=x[ibg+0];
            imu_info_corr->bg[1]+=x[ibg+1];
            imu_info_corr->bg[2]+=x[ibg+2];
            sprintf(buff,"%10.3f - %5.3f, %10.3f - %5.3f  %10.3f - %5.3f",
                    imu_info_corr->bg[0],x[ibg],imu_info_corr->bg[1],x[ibg+1],imu_info_corr->bg[2],x[ibg+2]);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"Bg: "<<buff;
            buff[0]='\0';
        }
    }

    /*
    * Function : The constructors and detructor function of class cSppSolver 
    */
    cSppSolver::cSppSolver() {
        spp_conf_.mode=MODE_SPP;
        spp_conf_.gnssC.ele_min=10.0;
        spp_conf_.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf_.gnssC.trp_opt=TRP_SAAS;
        spp_conf_.gnssC.ion_opt=ION_KLB;
        spp_conf_.estimator=SOLVE_LSQ;

        para_=cParSetting();
        num_full_x_=para_.GetIGCOUPLEDPar(spp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    cSppSolver::cSppSolver(tIGCOUPLEDConf conf) {
        spp_conf_=conf;
        spp_conf_.mode=MODE_SPP;
        spp_conf_.mode_opt=MODE_OPT_KINEMATIC;
        para_=cParSetting(spp_conf_);
        num_full_x_=para_.GetIGCOUPLEDPar(spp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    cSppSolver::~cSppSolver() {}

    /*
    *  Function : Correction of doppler observations
    * -Args :
    *       tSatInfoUnit sat_info              IO            satellite informations
    *       Vector3d     rover_xyz              I            the rover receiver position
    *       int          f                      I            the frq of doppler 
    * -Return :
    *       None
    */
    void cSppSolver::CorrDoppler(tSatInfoUnit &sat_info, Vector3d &rover_xyz,int f) {
        tSatInfoUnit sat_info0;
        Vector3d rover_blh=Xyz2Blh(rover_xyz);// 移动站位置
        Vector3d sat_pos=sat_info.brd_pos,sat_vel=sat_info.brd_vel;//广播星历卫星位置和速度
        sat_info0.brd_pos=sat_pos-sat_vel*1.0;
        Vector3d sig_vec0;
        GeoDist(sat_info0.brd_pos,rover_xyz,sig_vec0);// 卫地几何距离
        SatElAz(rover_blh,sig_vec0,sat_info0.el_az);//计算卫星高度角
        //计算大气误差
        gnss_err_corr_.trp_model_.InitSatInfo(&sat_info0,&rover_blh);
        gnss_err_corr_.trp_model_.GetSaasTrp(0.7, nullptr, nullptr);
        gnss_err_corr_.trp_model_.UpdateSatInfo();
        gnss_err_corr_.ion_model_.InitSatInfo(&sat_info0,&rover_blh);
        gnss_err_corr_.ion_model_.UpdateSatInfo();
        double trp_delta=((sat_info.trp_dry_delay[0]+sat_info.trp_wet_delay[0])-(sat_info0.trp_dry_delay[0]+sat_info0.trp_wet_delay[0])/1.0);//对流层改正值
        double ion_delta=((sat_info.ion_delay[0]-sat_info0.ion_delay[0])/1.0);//电离层改正值
        sat_info.cor_D[f]=sat_info.raw_D[f]-trp_delta+ion_delta;//多普勒观测值改正
    }

    /*
    * Function : Correct the Sagnac effect of satellites velocity
    * -Args :
    *       Vector3d sat_vel                I       The velocity of satellites
    *       double   tau                    I       The time between satellite and receiver
    * -Return :
    *       Vector3d                                Corrected velocity
    */
    Vector3d cSppSolver::SatVelSagnacCorr(const Vector3d &sat_vel, const double tau) {
        Vector3d vel_corr(0,0,0);
        double a=OMGE_GPS*tau;

        vel_corr[0]=sat_vel[0]+a*sat_vel[1];
        vel_corr[1]=sat_vel[1]-a*sat_vel[0];
        vel_corr[2]=sat_vel[2];
        return vel_corr;
    }

    /*
    * Function : The residual of doppler observations
    * -Args :
    *       tIGCOUPLEDConf C                        I            Configure file
    *       MatrixXd  H_mat                    I            Design matrix
    *       MatrixXd  R_mat                    I            noise matrix
    *       VectorXd  x                        I            parameters to be estimated
    *       Vector3d  rover_xyz                I            the position of rover 
    * -Return :
    *       int num_doppler                                 the number of doppler 
    */
    int cSppSolver::DopplerRes(tIGCOUPLEDConf C,MatrixXd& H_mat, MatrixXd& R_mat,VectorXd& L,VectorXd& x,Vector3d rover_xyz) {
        vector<double>H,omcs,meas_var_vec;
        int sys,sys_mask[NSYS]={0},num_doppler=0;
        double omc,r,meas_var;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;
        Matrix3d Cen=CalcCen(rover_blh,COORD_ENU);

        int num_used_frq=para_.GetGnssUsedFrqs();//用了多少频率

        for(int i=0;i<epoch_sat_info_collect_.size();i++){//逐历元遍历
            sat_info=&epoch_sat_info_collect_.at(i);//第i历元的卫星信息
            if(sat_info->sat.sat_.no==4) continue;
            if(sat_info->sat.sat_.no==55) continue;
            sys=sat_info->sat.sat_.sys;
//            if(sys==SYS_BDS&&(1<=sat_info->sat.sat_.prn&&sat_info->sat.sat_.prn<=5||sat_info->sat.sat_.prn>18)) continue;
//            if(sys!=SYS_GPS) continue; // only use GPS dopppler for velocity estimate
            if(sat_info->stat!=SAT_USED) continue;

            CLOG_IF(i==0,DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"DOPPLER RESIDUAL"<<"("<<epoch_idx_<<")"<<": "<<sat_info->t_tag.GetTimeStr(1)
                              <<setprecision(12)<<" "<<"VX="<<x[0]<<" "<<"VY="<<x[1]<<" "<<"VZ="<<x[2];

            for(int f=0;f<num_used_frq;f++){//逐频率遍历
                double doppler=sat_info->raw_D[f];//原始多普勒观测值
                double lam=sat_info->lam[f];//波长
                if(doppler==0.0||lam==0.0) continue;
                CorrDoppler(*sat_info,rover_xyz,f);//改正大气误差的多普勒观测值
                Vector3d sat_pos=sat_info->brd_pos;//广播星历的卫星位置和速度
                Vector3d sat_vel=sat_info->brd_vel;
                Vector3d s_r=(sat_pos-rover_xyz);//站星距离
                double tau=s_r.norm()/CLIGHT;//站星距离除以光速
                Vector3d sat_vel_corr=SatVelSagnacCorr(sat_vel,tau);//卫星速度萨格纳克改正
                Vector3d rover_vel(x[0],x[1],x[2]);
                Vector3d e=s_r/s_r.norm();//站星方向单位向量
                Vector3d relative_vel=sat_vel_corr-rover_vel;//卫星和移动站的相对速度
                double rate=e.dot(relative_vel);

#if 0
                double cos_el=cos(sat_info->el_az[0]);
                Vector3d a(0.0,0.0,0.0);
                a[0]=sin(sat_info->el_az[1])*cos_el;
                a[1]=cos(sat_info->el_az[1])*cos_el;
                a[2]=sin(sat_info->el_az[0]);
                Vector3d e=Cen.transpose()*a;

                int idx_vel=para_.IndexPos();
                Vector3d vs(0,0,0);
                for(int j=idx_vel;j<idx_vel+3;j++) vs[j-idx_vel]=sat_info->brd_vel[j-idx_vel]-x[j];
                double rate=vs.dot(e)+OMGE_GPS/CLIGHT*(sat_info->brd_vel[1]*rover_xyz[0]+sat_info->brd_pos[1]*x[idx_vel]
                                                       -sat_info->brd_vel[0]*rover_xyz[1]-sat_info->brd_pos[0]*x[idx_vel+1]);
#endif

                int idx_clk_drift=para_.IndexClk(SYS_INDEX_GPS);
                double a=CLIGHT*sat_info->brd_clk[1];
                omc=lam*sat_info->cor_D[f]-(rate+x[idx_clk_drift]-CLIGHT*sat_info->brd_clk[1]);
                meas_var=GnssMeasVar(C,GNSS_OBS_DOPPLER,*sat_info);

                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var);

                for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-e[j]:(j==idx_clk_drift)?1.0:0.0);
                int idx;
                if(sys==SYS_BDS) {idx=idx_clk_drift+SYS_INDEX_BDS;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_BDS]++;}
                else if(sys==SYS_GAL) {idx=idx_clk_drift+SYS_INDEX_GAL;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_GAL]++;}
                else if(sys==SYS_GLO) {idx=idx_clk_drift+SYS_INDEX_GLO;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_GLO]++;}
                else if(sys==SYS_QZS) {idx=idx_clk_drift+SYS_INDEX_QZS;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_QZS]++;}
                else sys_mask[SYS_INDEX_GPS]++;

                num_doppler++;

                char buff[MAX_BUFF]={'\0'};
                sprintf(buff,"%s omc=%12.4f, var=%7.3f el=%3.1f clk_drift=%14.3f sat_clk_drift=%14.3f doppler=%10.3f rate=%10.3f",
                        sat_info->sat.sat_.id.c_str(),omc,meas_var,sat_info->el_az[0]*R2D,sys==SYS_GPS?x[idx_clk_drift]:x[idx],sat_info->brd_clk[1]*CLIGHT,-lam*doppler,rate);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;
            }
        }

        int idx_clk=para_.IndexClk(SYS_INDEX_GPS);
        for(int i=0;i<NSYS;i++){
            if(sys_mask[i]) continue;
            omcs.push_back(0.0);
            for(int j=0;j<num_full_x_;j++) H.push_back(j==i+idx_clk?1.0:0.0);
            meas_var_vec.push_back(0.01);
            num_doppler++;
        }

        L=Map<VectorXd>(omcs.data(),num_doppler);
        H_mat=Map<MatrixXd>(H.data(),num_full_x_,num_doppler);
        Map<VectorXd> var_vec(meas_var_vec.data(),num_doppler);
        R_mat=var_vec.asDiagonal();
        H.clear();omcs.clear();meas_var_vec.clear();

        return num_doppler;
    }

    /*
    *  Function : Estimate doppler velocity
    */
    void cSppSolver::EstDopVel(Vector3d rover_xyz) {
        VectorXd L,x(num_full_x_);
        MatrixXd H,R,Q;
        int nv,stat=0;

        for(int i=0;i<num_full_x_;i++) x[i]=0.0;
        for(int i=0;i<iter_;i++){
            if((nv=DopplerRes(spp_conf_,H,R,L,x,rover_xyz))<4){
                break;
            }

//            cout<<H.transpose()<<endl;
//            cout<<R<<endl;
//            cout<<L<<endl;
            stat=lsq_.Adjustment(L,H,R,x,Q,nv,num_full_x_);//最小二乘求解

            if(stat){
                for(int j=0;j<3;j++) igcoupled_sol_.vel[j]=x[j];//求解速度更新进最终结果
                break;
            }
        }
    }

    /*
    * Function : Calculate the dops
    */
    double cSppSolver::Dops() {
        tSatInfoUnit *sat_info=nullptr;
        MatrixXd H_mat(4,num_valid_sat_),Q,P;
        vector<double>H(4*num_valid_sat_,0.0);
        int n=0;

        for(int i=0;i<epoch_sat_info_collect_.size();i++){//逐卫星遍历
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED) continue;
            double cos_el=cos(sat_info->el_az[0]);
            double sin_el=sin(sat_info->el_az[0]);
            H[4*n]=cos_el*sin(sat_info->el_az[1]);
            H[1+4*n]=cos_el*cos(sat_info->el_az[1]);
            H[2+4*n]=sin_el;
            H[3+4*n++]=1.0;
        }
        if(n<4) return 0.0;
        H_mat=Map<MatrixXd>(H.data(),4,n);
        Q=H_mat*H_mat.transpose();
        P=Q.inverse(); 
        igcoupled_sol_.dops[0]=SQRT(P.trace());
        igcoupled_sol_.dops[1]=SQRT(P(0,0)+P(1,1)+P(2,2));
        igcoupled_sol_.dops[2]=SQRT(P(0,0)+P(1,1));
        igcoupled_sol_.dops[3]=SQRT(P(2,2));

        return igcoupled_sol_.dops[1];
    }

    /*
    * Function : Validate the solution using dops and unit weight STD
    */
    bool cSppSolver::ValidateSol(tIGCOUPLEDConf C) {
        tSatInfoUnit* sat_info= nullptr;

        for(int i=0;i<epoch_sat_info_collect_.size();i++){//逐卫星遍历,舍弃低高度角的卫星
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat==SAT_USED&&sat_info->el_az[0]<C.gnssC.ele_min*D2R){
                sat_info->stat=SAT_LOW_EL;
                continue;
            }
        }

        string buff;
        int num_no_used=0;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED){
                num_no_used++;
                if(num_no_used==1) buff=epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)+" SATELLITE NO USED: "+sat_info->sat.sat_.id+"("+kGnssSatStatStr[sat_info->stat+1]+") ";
                else{
                    buff+=sat_info->sat.sat_.id+"("+kGnssSatStatStr[sat_info->stat+1]+") ";
                }
            }
        }
        CLOG_IF(num_no_used>0,DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;

        double vv=lsq_.unit_weight_STD_*((num_L_>num_full_x_)?(num_L_-num_full_x_):num_L_);
        if(num_L_>num_full_x_&&vv>kChiSqr[num_L_-num_full_x_-1]){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" SPP SOLVE FAILED, CHI-SQR OVERRUN vv: "<<vv<<" THRESHOLD: "<<kChiSqr[num_L_-num_full_x_-1];
            return false;
        }

        double pdop=Dops();

        if(pdop>C.gnssC.max_pdop){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" SPP SOLVE FAILED, PDOP OVERRUN pdop: "<<pdop<<" THRESHOLD: "<<C.gnssC.max_pdop;
            return false;
        }

        return true;
    }
/*
* Function : RAIM-FDE algorithm to detect faulty satellite 
*/
    bool cSppSolver::RaimFde() {
        int i,j,k,exc_sat_idx=0;
        double rms_e=0.0;
        int stat=0;

        for(i=0;i<epoch_sat_info_collect_.size();i++){//同历元内逐卫星遍历

            exc_sat_idx=i;
            epoch_sat_info_collect_[exc_sat_idx].stat=SAT_NO_USE;

            if(!Estimator(spp_conf_)){//spp失败
                epoch_sat_info_collect_[exc_sat_idx].stat=SAT_USED;
                continue;
            }
            if(omc_L_.size()<5){
                epoch_sat_info_collect_[exc_sat_idx].stat=SAT_USED;
            }
            for(k=0;k<omc_L_.size();k++){
                rms_e+=SQR(omc_L_[k]);
            }
            rms_e=sqrt(rms_e/omc_L_.size());

            if(rms_e>100.0){
                epoch_sat_info_collect_[exc_sat_idx].stat=SAT_USED;
            }
            stat=1;
            break;
        }
        return stat;
    }

    /*
    * Function : 
    */
    bool cSppSolver::PostResidualQc(vector<double>omcs,vector<double>R) {

        bool flag=false;
        int i,sat_idx,f;
        vector<double>v,avg_v,fabs_v,norm_v;
        double s0=0.0,el,k0=ReNorm(0.95),k1=ReNorm(0.99);

        for(i=0;i<vflag_.size();i++){
            v.push_back(omcs[i]);
            fabs_v.push_back(fabs(omcs[i]));
            s0+=v[i];
        }

#if 0
        if(lsq_.unit_weight_STD_>0.3){
            auto max_v=max_element(begin(fabs_v),end(fabs_v));
            int idx_max_v=distance(begin(fabs_v),max_v);
            sat_idx=vflag_[idx_max_v]>>4&0xF;
            epoch_sat_info_collect_[sat_idx].stat=SAT_NO_USE;
            flag=true;
        }
        else{
            for(i=0;i<v.size();i++) avg_v.push_back(v[i]-s0/v.size());
            MatMul("NT",1,1,v.size(),1.0/(v.size()-1),avg_v.data(),avg_v.data(),0.0,&s0);
            for(i=0;i<v.size();i++) {
                sat_idx=vflag_[i]>>4&0xF;
                el=epoch_sat_info_collect_[sat_idx].el_az[0];
                if(epoch_sat_info_collect_[sat_idx].c_var_factor[0]!=1.0) continue;
                norm_v.push_back(fabs_v[i]/sin(el));
            }

            auto max_norm_v=max_element(begin(norm_v),end(norm_v));
            int idx_max_norm_v=distance(begin(norm_v),max_norm_v);
            sat_idx=vflag_[idx_max_norm_v]>>4&0xF;
            if(*max_norm_v>k1){
                epoch_sat_info_collect_[sat_idx].stat=SAT_NO_USE;
                flag=true;
            }
            else if(*max_norm_v>k0){
                double fact=(k0/(*max_norm_v))*SQR((k1-(*max_norm_v))/(k1-k0));
                epoch_sat_info_collect_[sat_idx].c_var_factor[0]=1.0/fact;
                flag=true;
            }
        }
#endif
        for(i=0;i<v.size();i++) avg_v.push_back(v[i]-s0/v.size());
        MatMul("NT",1,1,v.size(),1.0/(v.size()-1),avg_v.data(),avg_v.data(),0.0,&s0);
        for(i=0;i<v.size();i++) {
            sat_idx=vflag_[i]>>4&0xF;
            el=epoch_sat_info_collect_[sat_idx].el_az[0];
            if(epoch_sat_info_collect_[sat_idx].c_var_factor[0]!=1.0) continue;
            norm_v.push_back(fabs_v[i]/sin(el));
        }

        auto max_norm_v=max_element(begin(norm_v),end(norm_v));
        int idx_max_norm_v=distance(begin(norm_v),max_norm_v);
        sat_idx=vflag_[idx_max_norm_v]>>4&0xF;
        if(*max_norm_v>k1){
            epoch_sat_info_collect_[sat_idx].stat=SAT_NO_USE;
            flag=true;
        }
        else if(*max_norm_v>k0){
            double fact=(k0/(*max_norm_v))*SQR((k1-(*max_norm_v))/(k1-k0));
            epoch_sat_info_collect_[sat_idx].c_var_factor[0]=1.0/fact;
            flag=true;
        }

        v.clear();avg_v.clear();fabs_v.clear(),norm_v.clear();
        return flag;
    }

    /*
    * Function : GNSS obeservation residuals
    */
    int cSppSolver::GnssObsRes(int post, tIGCOUPLEDConf C,double* x) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        for(int i=0;i<NSYS;i++) sys_mask_[i]=0;

        vector<double>H,omcs,meas_var_vec;
        int sys;
        double omc,r,meas_var;
        bool have_large_res=false;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_xyz,ve;

        rover_xyz<<x[0],x[1],x[2];

        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;

        int num_used_frq=para_.GetGnssUsedFrqs();//获取使用的频率个数
        if(num_used_frq<=0) return 0;

        vector<int>sat_idx;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);

            CLOG_IF(i==0,DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SINGLE POINT POSITION RESIDUAL"<<"("<<epoch_idx_<<"-"<<post<<")"<<": "<<sat_info->t_tag.GetTimeStr(1)
                                       <<setprecision(12)<<" "<<"X="<<rover_xyz[0]<<" "<<"Y="<<rover_xyz[1]<<" "<<"Z="<<rover_xyz[2];

            sys=sat_info->sat.sat_.sys;//卫星系统
            if(sat_info->stat!=SAT_USED) continue;

            for(int f=0;f<num_used_frq;f++){//逐频率遍历
                double cor_P=0.0;
                if(spp_conf_.gnssC.ion_opt==ION_IF){//单频无电离层组合
                    cor_P=sat_info->cor_if_P[0];
                }
                else{
                    cor_P=sat_info->cor_P[f];
                }

                if(cor_P==0.0){
                    sat_info->stat=SAT_NO_CODE;
                    continue;
                }

                r=GeoDist(sat_info->brd_pos,rover_xyz,sig_vec);//卫地距
                if(SatElAz(rover_blh,sig_vec,sat_info->el_az)<=0.0){
                    continue;
                }
                double el=sat_info->el_az[0]*R2D;
                if(sat_info->el_az[0]*R2D<C.gnssC.ele_min){//低高度角
                    continue;
                }

#if 0
                if(sat_info->sat.sat_.sys==SYS_BDS&&sat_info->sat.sat_.prn<19){
                    gnss_err_corr_.BD2MultipathModel(sat_info);//北2去除多径效应
                }
#endif
                gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
                gnss_err_corr_.ion_model_.GetIonError();
                gnss_err_corr_.ion_model_.UpdateSatInfo();
                gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
                gnss_err_corr_.trp_model_.GetSaasTrp(0.7, nullptr, nullptr);
                gnss_err_corr_.trp_model_.UpdateSatInfo();

                double sag_err=gnss_err_corr_.SagnacCorr(sat_info->brd_pos,rover_xyz);//萨格纳克速度改正
                double rec_clk=x[para_.IndexClk(SYS_INDEX_GPS)];//接收机钟差（待估参数）
                double sat_clk=sat_info->brd_clk[0]*CLIGHT;//卫星钟差（广播星历）
                double trp_err=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];//对流层误差萨寺塔摩宁模型
                double ion_err=sat_info->ion_delay[0];//无电离层组合为0
                double cod_bia=sat_info->code_bias[f];//码偏差

                omc=cor_P-(r+sag_err+(rec_clk-CLIGHT*sat_info->brd_clk[0])+sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0]+sat_info->ion_delay[0]);
                if(epoch_idx_>20&&fabs(omc)>5000){ //500
//                    sat_info->stat=SAT_NO_USE;
//                    continue;
                }
                double a=GnssMeasVar(C,GNSS_OBS_CODE,*sat_info);//随机模型，观测值精度
                meas_var=GnssMeasVar(C,GNSS_OBS_CODE,*sat_info)+sat_info->brd_eph_var+sat_info->trp_var+sat_info->ion_var;
                meas_var*=sat_info->c_var_factor[f];
                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var); 

                int idx_clk=para_.IndexClk(SYS_INDEX_GPS),idx=0;
                for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-sig_vec[j]:(j==idx_clk)?1.0:0.0);
                if(spp_conf_.gnssC.est_bd3_isb){//如果估计北斗3的isb
                    if(sys==SYS_BDS&&sat_info->sat.sat_.prn<19) {idx=idx_clk+SYS_INDEX_BDS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_BDS]++;}
                    else if(sys==SYS_GAL) {idx=idx_clk+SYS_INDEX_GAL;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GAL]++;}
                    else if(sys==SYS_GLO) {idx=idx_clk+SYS_INDEX_GLO;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GLO]++;}
                    else if(sys==SYS_QZS) {idx=idx_clk+SYS_INDEX_QZS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_QZS]++;}
                    else if(sys==SYS_BDS&&sat_info->sat.sat_.prn>18) {
                        idx=idx_clk+NSYS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[NSYS]++;
                    }
                    else sys_mask_[SYS_INDEX_GPS]++;
                }
                else{
                    if(sys==SYS_BDS) {idx=idx_clk+SYS_INDEX_BDS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_BDS]++;}
                    else if(sys==SYS_GAL) {idx=idx_clk+SYS_INDEX_GAL;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GAL]++;}
                    else if(sys==SYS_GLO) {idx=idx_clk+SYS_INDEX_GLO;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GLO]++;}
                    else if(sys==SYS_QZS) {idx=idx_clk+SYS_INDEX_QZS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_QZS]++;}
                    else sys_mask_[SYS_INDEX_GPS]++;
                }
                num_L_++;

                if(post==0) sat_info->prior_res_P[f]=omc;

                char buff[MAX_BUFF]={'\0'};
                sprintf(buff,"%s omc=%12.4f, var=%7.3f el=%3.1f dtr=%14.3f dts=%14.9f trp=%6.3f ion=%6.3f cbias=%6.3f obs=%14.3f range=%14.3f",
                        sat_info->sat.sat_.id.c_str(),omc,meas_var,sat_info->el_az[0]*R2D,sys==SYS_GPS?rec_clk:x[idx],sat_clk/CLIGHT*1E6,trp_err,ion_err,cod_bia,cor_P,r+sag_err);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;

                if(f==0){
                    sat_idx.push_back(i);
                    num_valid_sat_++;
                }

                vflag_.push_back((sat_info->sat.sat_.no<<8)|(i<<4)|(f));
            }
        }

        if(epoch_idx_>10&&post&&spp_conf_.gnssC.res_qc){
            have_large_res=PostResidualQc(omcs,meas_var_vec);
        }

        int idx_clk=para_.IndexClk(SYS_INDEX_GPS);
        for(int i=0;i<para_.NumClk();i++){
            if(sys_mask_[i]) continue;
            omcs.push_back(0.0);
            for(int j=0;j<num_full_x_;j++) H.push_back(j==i+idx_clk?1.0:0.0);
            meas_var_vec.push_back(0.01);
            num_L_++;
        }

        omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
        Map<VectorXd> var_vec(meas_var_vec.data(),num_L_);
        R_=var_vec.asDiagonal();
        H.clear();omcs.clear();meas_var_vec.clear();

        return post==1?have_large_res:num_valid_sat_<=0;//true : 有大残差或者无有效卫星
    }

    /*
    * Function : SPP estimator using lsq
    */
    bool cSppSolver::Estimator(tIGCOUPLEDConf C) {
        bool stat;
        bool valid_flag=false;
        igcoupled_sol_.stat=SOL_NONE;
        tSatInfoUnit* sat_info= nullptr;

        VectorXd x=full_x_;
        MatrixXd Px=full_Px_;
        Vector3d re,ve;
        for(int i=0;i<iter_;i++){

            if(GnssObsRes(0,spp_conf_,x.data())) continue;
            if(num_valid_sat_<4){
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<"SPP NO ENOUGH VALID SATELLITE "<<num_valid_sat_;
                return false;
            }
            igcoupled_sol_.valid_sat_num=num_valid_sat_;
            stat=lsq_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_);//最小二乘估计

            if(GnssObsRes(i+1,spp_conf_,x.data())){
                x=full_x_;
                continue;
            }

            if(stat){
                full_x_=x;
                full_Px_=Px;
                valid_flag = ValidateSol(C);
                if(valid_flag){
                    igcoupled_sol_.stat=SOL_SPP;
                    igcoupled_sol_.sigma=lsq_.unit_weight_STD_;
                }
                else{
                    igcoupled_sol_.stat=SOL_NONE;
                }
                break;
            }
            if(i>=iter_){
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<"SPP ITER OVERRUN";
            }
        }

        return valid_flag;
    }

    /*
    * Function : Initialize spp solover
    */
    void cSppSolver::InitSolver(tIGCOUPLEDConf C) {
        spp_conf_=C;
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_=new cOutSol(spp_conf_);
        out_->InitOutSol(spp_conf_,spp_conf_.fileC.sol);
        out_->WriteHead();
        out_->ref_sols_=ref_sols_;

        spp_conf_.mode=MODE_SPP;
        spp_conf_.mode_opt=MODE_OPT_KINEMATIC;
        para_=cParSetting(spp_conf_);
        num_full_x_=para_.GetIGCOUPLEDPar(spp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    /*
    * Function : The process of spp solver
    */
    bool cSppSolver::SolverProcess(tIGCOUPLEDConf C,int idx) {
        char buff[MAX_BUFF]={'\0'};
        if(idx==-1) InitSolver(C);

        int i=0,num_epochs=rover_obs_.epoch_num;

		double outage_time = 0;
		if(spp_conf_.gnssC.use_outage){ // gnss_outage
			outage_time = spp_conf_.gnssC.outage_time;
		}

        if(idx!=-1){
            num_epochs=idx+1;
        }
        for(i=idx==-1?0:idx;i<num_epochs;i++){

            epoch_idx_+=1;
            epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
			if(outage_time != 0){
				double sow = epoch_sat_obs_.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS);
				if(sow > (outage_time + spp_conf_.gnssC.outage_len) && spp_conf_.gnssC.outage_period > spp_conf_.gnssC.outage_len){ // reinit outage time
					outage_time += static_cast<int>((sow - outage_time) / spp_conf_.gnssC.outage_len) * spp_conf_.gnssC.outage_len;
				}
				if(sow >= outage_time && sow <= outage_time + spp_conf_.gnssC.outage_len){ // outage
					epoch_fail_++;
					epoch_sat_info_collect_.clear();
					CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID) << " GNSS OUTAGE AT " << sow << "(s)" <<endl;
					continue;
				}
			}
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"START SOLVING: "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER: "<<epoch_sat_obs_.sat_num;

            UpdateGnssObs(spp_conf_,epoch_sat_obs_,REC_ROVER);
            if(SolverEpoch()){
                epoch_ok_++;
                SolutionUpdate();
                sprintf(buff,"%s SPP SOLVE SUCCESS POS: %14.3f %14.3f %14.3f VEL: %6.3f %6.3f %6.3f PDOP: %3.1f TOTAL SAT: %3d USED SAT: %3d",
                        igcoupled_sol_.t_tag.GetTimeStr(1).c_str(),igcoupled_sol_.pos[0],igcoupled_sol_.pos[1],igcoupled_sol_.pos[2],igcoupled_sol_.vel[0],igcoupled_sol_.vel[1],igcoupled_sol_.vel[2],
                       igcoupled_sol_.dops[1],epoch_sat_obs_.sat_num,num_valid_sat_);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;

                out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
                epoch_sat_info_collect_.clear();
            }
            else{
                epoch_fail_++;
                epoch_sat_info_collect_.clear();
            }

        }
        CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"TOTAL EPOCH: "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
    }

    /*
    * Function : Solver processing by epoch
    */
    bool cSppSolver::SolverEpoch() {
        gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);//星历误差改正
        Vector3d rr,ve;
        if(tc_mode_){//如果INS-GNSS紧组合
            RemoveLever(cur_imu_info_,spp_conf_.insC.lever,rr,ve);//移除杆壁效应
            for(int i=0;i<3;i++) full_x_[i]=rr[i];//用改正后的imu位置代替初值
        }
        else{
            rr<<full_x_[0],full_x_[1],full_x_[2];
        }
        CorrGnssObs(spp_conf_,rr);//改正GNSS观测值

        if(Estimator(spp_conf_)){//如果估计成功
            igcoupled_sol_.stat=SOL_SPP;
            igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;

            if(spp_conf_.gnssC.use_doppler){//如果使用多普勒观测值
                EstDopVel(rr);
            }
            return true;
        }
        else{//没有估计成功
            igcoupled_sol_.stat=SOL_NONE;
            igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;

            // failure detection and exclution
            if(num_valid_sat_>=6){
                if(RaimFde()){
                    igcoupled_sol_.stat=SOL_SPP;
                    if(spp_conf_.gnssC.use_doppler){
                        EstDopVel(rr);
                    }
                    return true;
                }
            }

            if(epoch_sat_info_collect_.size()) igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            return false;
        }
    }

    /*
    * Function : Update the solution using estimated x(position and clk)
    */
    bool cSppSolver::SolutionUpdate() {
        int m=0;
        int num_pos=para_.NumPos();
        for(int i=0;i<num_pos;i++) igcoupled_sol_.pos[i]=full_x_[i];

        int ic=para_.IndexClk(SYS_INDEX_GPS);
        if(sys_mask_[0]==0){
            for(m=1;m<para_.NumClk();m++){
                if(sys_mask_[m]!=0) break;
            }
            for(int i=1;i<para_.NumClk();i++){
                if(i==m) igcoupled_sol_.clk_error[i]=full_x_[ic+i];
                else{
                    if(sys_mask_[i]==0){
                        igcoupled_sol_.clk_error[i]=0.0;
                        continue;
                    }
                    igcoupled_sol_.clk_error[i]=full_x_[ic+i]-full_x_[ic+m]; //ISB
                }
            }
        }
        else{
            for(int i=0;i<para_.NumClk();i++){
                if(sys_mask_[i]){
                    igcoupled_sol_.clk_error[i]=full_x_[ic+i];
                }
                else igcoupled_sol_.clk_error[i]=0.0;
            }
        }
    }

    /*
    * Function : The constructor and destructor of class cTdcpSolver
    */
    cTdcpSolver::cTdcpSolver() {}

    cTdcpSolver::cTdcpSolver(tIGCOUPLEDConf conf) {
        tdcp_conf_=conf;
    }

    cTdcpSolver::~cTdcpSolver() {}

    /*
    * Function : The Constructor and destructor of class cPppSolver 
    */
    cPppSolver::cPppSolver() {}

    cPppSolver::cPppSolver(tIGCOUPLEDConf C) {
        ppp_conf_=C;
        para_=cParSetting(ppp_conf_);
        num_full_x_=para_.GetIGCOUPLEDPar(ppp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);

        num_real_x_fix_=para_.GetRealFixParNum(ppp_conf_);
        real_x_fix_=VectorXd::Zero(num_real_x_fix_);
        real_Px_fix_=MatrixXd::Zero(num_real_x_fix_,num_real_x_fix_);
    }

    cPppSolver::~cPppSolver() {}

    /*
    * Function : Initialize the ppp solver
    */
    void cPppSolver::InitSolver(tIGCOUPLEDConf C) {
        ppp_conf_=C;

        if(C.site_idx==1){
            cReadGnssPreEph clk_reader(C.fileC.clk,nav_);//读取精密钟差
            clk_reader.Reading(1);

            if(C.gnssC.ar_mode==AR_PPP_AR&&C.gnssC.ar_prod==AR_PROD_IRC_CNES&&!C.fileC.pr_clk.empty()){// IRC 固定模糊度
                cReadGnssPreEph code_clk_reader(C.fileC.pr_clk,nav_);//读取整数钟
                code_clk_reader.irc_pr_clk_=true;
                code_clk_reader.Reading(1);
            }

            cReadGnssPreEph orb_reader(C.fileC.sp3[1], nav_);//读取精密星历
            for(int i=0;i<3;i++){
                if(C.fileC.sp3[i].empty()) continue;
                orb_reader.file_=C.fileC.sp3[i];
                orb_reader.Reading(0);
            }

            if(C.gnssC.tid_opt){//潮汐模型
                cReadGnssErp erp_reader(C.fileC.erp,nav_);//读取地球自转参数erp
                erp_reader.Reading();

                cReadGnssOcean blq_reader(C.fileC.blq, nav_,C.site_name,REC_ROVER);
                blq_reader.Reading();
            }

            if((C.mode==MODE_PPP||C.mode_opt==MODE_OPT_PPP)&&C.gnssC.ar_mode==AR_PPP_AR&&C.gnssC.ar_prod==AR_PROD_FCB_SGG){//FCB固定模糊的
                cReadFcb fcb_reader(C.fileC.fcb,nav_);//读取FCB
                fcb_reader.Reading();
            }

            para_=cParSetting(C);
            gnss_err_corr_.InitGnssErrCorr(C,&nav_);//初始化大气误差
        }

        cReadGnssAntex atx_reader(C.fileC.atx,nav_);//读取天线信息文件
        atx_reader.Reading();
        atx_reader.AlignAntPar2Sat(C,*rover_obs_.GetStartTime(),nav_.sta_paras,nav_.sat_ant,nav_.rec_ant);

        if(out_) delete out_;
        out_=new cOutSol(ppp_conf_);
        out_->InitOutSol(ppp_conf_,ppp_conf_.fileC.sol);
        out_->WriteHead();
        out_->ref_sols_=ref_sols_;
        if(ppp_conf_.solC.out_err_fmt&&!ppp_conf_.fileC.snx.empty()){
            if(!GetRefPosFrmSnx(ppp_conf_,ppp_conf_.fileC.snx,out_->ref_sol_.pos)){
                ppp_conf_.solC.out_err_fmt=false;
            }
        }

        tIGCOUPLEDConf spp_conf=C;
        spp_conf.mode=MODE_SPP;
        spp_conf.mode_opt=MODE_OPT_KINEMATIC;
        spp_conf.gnssC.ion_opt=ION_KLB;
        spp_conf.gnssC.trp_opt=TRP_SAAS;
        spp_conf.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf.gnssC.eph_opt=EPH_BRD;
        if(spp_solver_) delete spp_solver_;
        spp_solver_=new cSppSolver(spp_conf);
        spp_solver_->spp_conf_=spp_conf;
        spp_solver_->nav_=nav_;
        spp_solver_->ref_sols_=ref_sols_;
        spp_solver_->InitSolver(spp_conf);
        spp_solver_->spp_conf_.gnssC.res_qc=false;
    }

    /*
    * Function : The process of spp solver
    */
    bool cPppSolver::SolverProcess(tIGCOUPLEDConf C,int idx) {
        bool stat=false;
        double rate=0.0;
        if(idx==-1) InitSolver(C);

        int i=0,num_epochs=rover_obs_.epoch_num;
        if(idx!=-1){
            num_epochs=idx+1;
        }

        if(C.site_idx!=1){
            epoch_ok_=0;
            epoch_idx_=0;
            epoch_fail_=0;
            ReInitPppSolver(C);
        }
		double outage_time = 0.0;
		if(ppp_conf_.gnssC.use_outage) outage_time = ppp_conf_.gnssC.outage_time;

        if(C.filter_type==FILTER_FORWARD){//前向滤波
            for(i=idx==-1?0:idx;i<num_epochs;i++){
				if(outage_time != 0){ //gnss outage
					double sow = rover_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS);
					if(sow > (outage_time + ppp_conf_.gnssC.outage_len) && 
						ppp_conf_.gnssC.outage_period > ppp_conf_.gnssC.outage_len){ // reinit outage time
							outage_time += static_cast<int>((sow - outage_time) / ppp_conf_.gnssC.outage_period) * 
											ppp_conf_.gnssC.outage_period;
						}
						if(sow >= outage_time && sow <= outage_time + ppp_conf_.gnssC.outage_len){// outage
							CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS OUTAGE AT" << sow << "(s) " << endl;
							continue;
						}
				}
                stat=SolverStart(i,idx);
            }
            if(idx!=-1) return stat;
            else CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (FORWARD): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
        }
        else if(C.filter_type==FILTER_BACKWARD){//后向滤波
            for(i=idx==-1?num_epochs-1:idx;i>=0;i--){
                stat=SolverStart(i,idx);
            }
            if(idx!=-1) return stat;
            else CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (BACKWARD): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
        }
        else if(C.filter_type==FILTER_COMBINED){//混合滤波
            for(i=idx==-1?0:idx;i<num_epochs;i++){
                SolverStart(i,idx);
                solf_.push_back(igcoupled_sol_);
            }

            epoch_ok_=0;
            epoch_idx_=0;
            epoch_fail_=0;
            ReinitSolver(C);
            for(i=idx==-1?num_epochs-1:idx;i>=0;i--){
                SolverStart(i,idx);
                solb_.push_back(igcoupled_sol_);
            }

            CombFbSol(C);
            solf_.clear();
            solb_.clear();
            CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (COMBINED): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
        }
    }

    /*
    * Function : The start of filter solver
    */
    bool cPppSolver::SolverStart(int i,int idx) {
        char buff[MAX_BUFF]={'\0'};

        epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"START "<<(ppp_conf_.filter_type==FILTER_FORWARD?"FORWARD ":"BACKWARD ") <<" PPP SOLVING "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER "<<epoch_sat_obs_.sat_num;

        epoch_idx_+=1;
        if(epoch_idx_==1) filter_start_=epoch_sat_obs_.obs_time;
        UpdateGnssObs(ppp_conf_,epoch_sat_obs_,REC_ROVER);
        igcoupled_sol_.observed_sat_num=epoch_sat_obs_.sat_num;
        InitEpochSatInfo(epoch_sat_info_collect_);

        if(SolverEpoch()){
            SolutionUpdate();
            sprintf(buff,"%s PPP SOLVE SUCCESS POS: %14.3f %14.3f %14.3f VEL: %6.3f %6.3f %6.3f PDOP: %3.1f TOTAL SAT: %3d USED SAT: %3d",
                    igcoupled_sol_.t_tag.GetTimeStr(1).c_str(),igcoupled_sol_.pos[0],igcoupled_sol_.pos[1],igcoupled_sol_.pos[2],igcoupled_sol_.vel[0],igcoupled_sol_.vel[1],igcoupled_sol_.vel[2],
                    igcoupled_sol_.dops[1],epoch_sat_obs_.sat_num,num_valid_sat_);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;
            if(idx!=-1){
                epoch_sat_info_collect_.clear();
                return true;
            }
        }
        if(idx!=-1) {
            epoch_sat_info_collect_.clear();
            return false;
        }
        epoch_sat_info_collect_.clear();
    }

    /*
    * Function : Solver processing by epoch
    */
    bool cPppSolver:: SolverEpoch() {
        InitSppSolver();
        if(spp_solver_->SolverEpoch()){

            Spp2Ppp();//先做SPP
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PPP-SPP SOLVE SUCCESS ROVER, SPP POS "<<spp_solver_->igcoupled_sol_.pos.transpose()<<" DOPPLER  VEL "<<spp_solver_->igcoupled_sol_.vel.transpose();

            if(ppp_conf_.gnssC.eph_opt==EPH_PRE){
                gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);//星历误差改正
            }

            ClockJumpRepair();//钟跳探测和修复
            CorrGnssObs(ppp_conf_,spp_solver_->igcoupled_sol_.pos);//改正GNSS观测值
            PppCycSlip(ppp_conf_);//PPP周跳检测
            StateTimeUpdate(ppp_conf_);//状态时刻更新

            if(Estimator(ppp_conf_)){//估计成功
                epoch_ok_++;
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP SOLVE SUCCESS ";
                return true;
            }
            else{
                epoch_fail_++;
                igcoupled_sol_=spp_solver_->igcoupled_sol_;
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP SOLVE FAILED ";
                return false;
            }
        }
        else{//spp初始化失败
            epoch_fail_++;
            igcoupled_sol_.stat=SOL_NONE;
            if(epoch_sat_info_collect_.size()) igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP-SPP SOLVE FAILED ";
            spp_solver_->epoch_sat_info_collect_.clear();
            return false;
        }
    }

    /*
    * Function : The estimator of PPP (using KF filter)
    */
    bool cPppSolver::Estimator(tIGCOUPLEDConf C) {

        // residuals
        VectorXd x;
        MatrixXd Px;
        vector<int>par_idx;
        vector<double>back_values;
        igcoupled_sol_.epoch_idx=epoch_idx_;

        Vector3d re,ve;
        if(tc_mode_){//INS-GNSS紧组合
            RemoveLever(cur_imu_info_,C.insC.lever,re,ve);//杆壁改正
        }
        else{
            re<<full_x_[0],full_x_[1],full_x_[2];
        }

        tImuInfoUnit cor_imu_info;
        int iter;
        for(iter=0;iter<max_iter_;iter++){//迭代（这里为什么要迭代？）
            cor_imu_info=cur_imu_info_;
            x=full_x_;
            Px=full_Px_;

            if(!GnssObsRes(0,C,x.data(),re)){
                par_idx.clear(),back_values.clear();
                continue;
            }

            kf_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_);// 卡尔曼滤波器

            if(tc_mode_){//紧组合
                CloseLoopState(x,&cor_imu_info);//用卡尔曼滤波闭环修正imu位置信息
                RemoveLever(cor_imu_info,C.insC.lever,re,ve);//杆壁改正，作为下一迭代周期的re
            }
            else{
                re<<x[0],x[1],x[2];
            }

            if(GnssObsRes(iter+1,C,x.data(),re)) {
                continue;
            }
            else{//当下一迭代周期的残差过大或有无效卫星，退出迭代
                if(tc_mode_){
                    cur_imu_info_=cor_imu_info;//更新imu信息
                    igcoupled_sol_.ins_stat=SOL_IG_TC;
                }
                full_x_=x;
                full_Px_=Px;
                igcoupled_sol_.stat=SOL_PPP;
                igcoupled_sol_.sigma=kf_.unit_weight_STD_;
                UpdateSatInfo(epoch_sat_info_collect_);//更新卫星信息
                break;
            }

            if(iter>max_iter_){//如果迭代次数过多，仍没推出迭代，则返回错误信息
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP ITERATION OVERRUN";
                return false;
            }
        }

        if(igcoupled_sol_.stat==SOL_PPP&&C.gnssC.ar_mode==AR_PPP_AR){//模糊度固定
            VectorXd xa=x;
#if 0
            // 轮换参考星策略,某些历元存在lamada搜索失败
            if(ResolvePppAmbRotRef(2,x, Px)){
                if(!GnssObsRes(5,C,x.data(),re)){
                    if(tc_mode_){
                        CloseLoopState(x,&cur_imu_info_);
                        RemoveLever(cur_imu_info_,C.insC.lever,re,ve);

                    }
                    else{
                        re<<full_x_[0],full_x_[1],full_x_[2];
                    }
                    igcoupled_sol_.stat=SOL_FIX;
                }
            }
            else{
                int a=1;
            }

#else
            // 固定参考星策略
            if(ResolvePppAmbFixRef(2,x, Px)){
                if(!GnssObsRes(5,C,x.data(),re)){
                    if(tc_mode_){
                        CloseLoopState(x,&cur_imu_info_);
                        RemoveLever(cur_imu_info_,C.insC.lever,re,ve);

                    }
                    else{
                        re<<full_x_[0],full_x_[1],full_x_[2];
                    }
                    igcoupled_sol_.stat=SOL_FIX;
                }
            }
            else{
                int a=1;
            }

#endif
        }

        int sat_no;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){//每颗卫星
            for(int f=0;f<MAX_GNSS_USED_FRQ_NUM;f++){//每个频率
                if(!epoch_sat_info_collect_[i].vsat[f]) continue;
                sat_no=epoch_sat_info_collect_[i].sat.sat_.no;
                previous_sat_info_[sat_no-1].outc[f]=0;
                previous_sat_info_[sat_no-1].lock[f]++;
                if(previous_sat_info_[sat_no-1].lock[f]<0||previous_sat_info_[sat_no-1].fix[f]>=2){
                    previous_sat_info_[sat_no-1].lock[f]++;
                }
            }
        }

        vflag_.clear();
        return igcoupled_sol_.stat;
    }

    /*
    * Function : Update the PPP solution after solver processing
    */
    bool cPppSolver::SolutionUpdate() {
//        igcoupled_sol_.stat=SOL_PPP;
        igcoupled_sol_.valid_sat_num=num_valid_sat_;
        igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;

        if(igcoupled_sol_.stat==SOL_PPP){
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) igcoupled_sol_.pos[i]=full_x_[i];//位置
            int num_clk=para_.NumClk();
            for(int i=para_.IndexClk(SYS_INDEX_GPS);i<para_.IndexClk(SYS_INDEX_GPS)+num_clk;i++) igcoupled_sol_.clk_error[i-para_.IndexClk(SYS_INDEX_GPS)]=full_x_[i];//钟差
            int num_ifb=para_.NumIfb();
            for(int i=para_.IndexIfb(SYS_INDEX_GPS);i<para_.IndexIfb(SYS_INDEX_GPS)+num_ifb;i++) igcoupled_sol_.rec_ifcb[i-para_.IndexIfb(SYS_INDEX_GPS)]=full_x_[i];//ifb
            igcoupled_sol_.zenith_trp_delay[0]=gnss_err_corr_.trp_model_.zenith_trp_[0]; //ZTD
            int num_trp=para_.NumTrp();
            for(int i=para_.IndexTrp();i<para_.IndexTrp()+num_trp;i++) igcoupled_sol_.zenith_trp_delay[i-para_.IndexTrp()+1]=full_x_[i];//ZWD

            int ip=para_.IndexPos();
            for(int i=0;i<3;i++) igcoupled_sol_.q_pos[ip+i]=full_Px_(ip+i,ip+i);//位置精度
        }
        else{
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) igcoupled_sol_.pos[i]=real_x_fix_[i];
            int num_clk=para_.NumClk();
            for(int i=para_.IndexClk(SYS_INDEX_GPS);i<para_.IndexClk(SYS_INDEX_GPS)+num_clk;i++) igcoupled_sol_.clk_error[i-para_.IndexClk(SYS_INDEX_GPS)]=full_x_[i];
            int num_ifb=para_.NumIfb();
            for(int i=para_.IndexIfb(SYS_INDEX_GPS);i<para_.IndexIfb(SYS_INDEX_GPS)+num_ifb;i++) igcoupled_sol_.rec_ifcb[i-para_.IndexIfb(SYS_INDEX_GPS)]=full_x_[i];
            igcoupled_sol_.zenith_trp_delay[0]=gnss_err_corr_.trp_model_.zenith_trp_[0]; //ZTD
            int num_trp=para_.NumTrp();
            for(int i=para_.IndexTrp();i<para_.IndexTrp()+num_trp;i++) igcoupled_sol_.zenith_trp_delay[i-para_.IndexTrp()+1]=full_x_[i];

            int ip=para_.IndexPos();
            for(int i=0;i<3;i++) igcoupled_sol_.q_pos[ip+i]=real_Px_fix_(ip+i,ip+i);
        }

        if(ppp_conf_.filter_type!=FILTER_COMBINED){
            if(ppp_conf_.solC.out_sol) out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
            if(ppp_conf_.solC.out_stat&&ppp_conf_.solC.sol_fmt) out_->WriteSatStat(&igcoupled_sol_,previous_sat_info_);
        }
    }

    /*
    * Function : Update the ambiguity
    */
    void cPppSolver::AmbUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumAmb()<=0) return;

        int ia;
        double dt=0.0;
        double bias[64]={0},ion;
        int slip[64]={0};
        int idx_flag[64]={0};
        tSatInfoUnit* sat_info= nullptr;
        int num_use_frqs=para_.GetGnssUsedFrqs();
        if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF) num_use_frqs=1;          //这里为什么是1

        for(int f=0;f<num_use_frqs;f++){
            for(int i=0;i<MAX_SAT_NUM;i++){
                if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                    if(previous_sat_info_[i].outc[0]>5||previous_sat_info_[i].rejc[f]>2){
                        ia=para_.IndexAmb(0,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                        ia=para_.IndexAmb(1,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                        ia=para_.IndexAmb(2,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    }
                }
                else{
                    if(previous_sat_info_[i].outc[f]>5||previous_sat_info_[i].rejc[f]>2){
                        ia=para_.IndexAmb(f,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    }
                };

            }

            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                sat_info=&epoch_sat_info_collect_.at(i);
                ia=para_.IndexAmb(f,sat_info->sat.sat_.no);
                bias[i]=0.0;
                idx_flag[i]=0.0;

                if(sat_info->stat!=SAT_USED) continue;

                // ionosphere-free combination
                if(C.gnssC.ion_opt==ION_IF){
                    if(C.gnssC.frq_opt==FRQ_SINGLE){
                        // single frequency
                        if(sat_info->raw_P[f]!=0.0||sat_info->raw_L[f]!=0.0){
                            bias[i]=sat_info->raw_L[f]*sat_info->lam[f]-sat_info->raw_P[f];
                            slip[i]=sat_info->slip[f];
                            idx_flag[i]=0;
                        }
                    }
                    else if(C.gnssC.frq_opt==FRQ_DUAL){
                        // dual frequency
                        if(sat_info->cor_if_L[f]==0.0||sat_info->cor_if_P[f]==0.0) continue;
                        bias[i]=sat_info->cor_if_L[f]-sat_info->cor_if_P[f];
                        slip[i]=sat_info->slip[f];
                        idx_flag[i]=0;
                    }
                    else if(C.gnssC.frq_opt==FRQ_TRIPLE){
                        // triple frequency with one IF-combination, while missing triple observables,use dual-frequency(L1_L2/L1_L5) instead
                        if(sat_info->cor_if_P[0]!=0.0&&sat_info->cor_if_L[0]!=0.0){
                            // L1_L2_L5
                            bias[i]=sat_info->cor_if_L[0]-sat_info->cor_if_P[0];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=0;
                            sat_info->tf_if_idx[0]=1;
                        }
                        else if(sat_info->cor_if_P[1]!=0.0&&sat_info->cor_if_L[1]!=0.0){
                            // L1_L2
                            bias[i]=sat_info->cor_if_L[1]-sat_info->cor_if_P[1];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=1;
                            sat_info->tf_if_idx[1]=1;
                        }
                        else if(sat_info->cor_if_P[2]!=0.0&&sat_info->cor_if_L[2]!=0.0){
                            // L1_L5
                            bias[i]=sat_info->cor_if_L[2]-sat_info->cor_if_P[2];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=2;
                            sat_info->tf_if_idx[2]=1;
                        }
                        else{
                            sat_info->stat=SAT_NO_USE;
                            continue;
                        }
                    }
                }
                else if(C.gnssC.ion_opt==ION_IF_DUAL&&C.gnssC.frq_opt==FRQ_TRIPLE){
                    // L1_L2
                    if(sat_info->cor_if_P[f]!=0.0&&sat_info->cor_if_L[f]!=0.0){
                        bias[i]=sat_info->cor_if_L[f]-sat_info->cor_if_P[f];
                        slip[i]=sat_info->slip[f];
                    }
                    else{
//                        sat_info->stat=SAT_NO_USE;
                        continue;
                    }
                }
                else if(sat_info->raw_L[f]!=0.0&&sat_info->raw_P[f]!=0.0){
                    slip[i]=sat_info->slip[f];
                    if(sat_info->raw_P[0]==0.0||sat_info->raw_P[f]==0.0) ion=0.0;
                    else{
                        ion=(sat_info->raw_P[0]-sat_info->raw_P[1])/(1.0-SQR(sat_info->lam[1]/sat_info->lam[0]));
                    }
                    bias[i]=sat_info->raw_L[f]*sat_info->lam[f]-sat_info->raw_P[f]+2.0*ion*SQR(sat_info->lam[f]/sat_info->lam[0]);
                }

                if(full_x_[ia]==0.0||slip[i]||bias[i]==0.0) continue;
            }

            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                sat_info=&epoch_sat_info_collect_.at(i);
                if(previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_.long_time==0.0){
                    dt=ppp_conf_.gnssC.sample_rate;
                }
                else dt=spp_solver_->igcoupled_sol_.t_tag.TimeDiff(previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);

                // for UC: f==0 L1, f==1 L2, f==2 L3
                // for DF-IF: f==0 L1_L2/L1_L5
                // for TF-IF: f==0 L1_L2_L5, f==1 L1_L2, f==2 L1_L5
                // for TF-IF-DUAL: f==0 L1_L2, f==1 L1_L5
                if(idx_flag[i]==0){
                    ia=para_.IndexAmb(f,sat_info->sat.sat_.no);
                }
                else{
                    ia=para_.IndexAmb(idx_flag[i],sat_info->sat.sat_.no);
                }

                                // UC-L5                                                                    // TF-IF       L1_L2_L5       L1_L5
                if((C.gnssC.frq_opt==FRQ_TRIPLE&&f==FRQ_TRIPLE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF&&(idx_flag[i]==0||idx_flag[i]==2))){
                    //  L1_L2_L5/L1_L5/L5
                    full_Px_(ia,ia)+=3E-7*fabs(dt);
                }
                else if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF_DUAL&&f==1){
                    full_Px_(ia,ia)+=3E-5*fabs(dt);
                }
                else {
                    // L1/L2/L1_L2 psd=0.0
                    full_Px_(ia,ia)+=C.gnssC.ait_psd[0]*fabs(dt);
                }

                if(reset_flag_){
                    InitX(bias[i],SQR(60),ia,full_x_.data(),full_Px_.data());
                }
                else{
                    if(bias[i]==0.0||(full_x_[ia]!=0.0&&!slip[i]&&full_x_[ia]!=DIS_FLAG)) continue;
                    InitX(bias[i],SQR(60),ia,full_x_.data(),full_Px_.data());
                }

                string s="UC-L";
                if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                    if(sat_info->tf_if_idx[0]==1){
                        s="TF-L_123";
                    }
                    else if(sat_info->tf_if_idx[1]==1){
                        s="TF-L_12";
                    }
                    else if(sat_info->tf_if_idx[2]==1){
                        s="TF-L_15";
                    }
                }
                else if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF_DUAL){
                    if(f==0){
                        s="TF-IF_DUAL-L12";
                    }
                    else if(f==1){
                        s="TF-IF_DUAL-L13";
                    }
                }
                else if(C.gnssC.frq_opt==FRQ_DUAL&&C.gnssC.ion_opt==ION_IF){
                    s="DF-L_12";
                }
                else if(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF){
                    s="SF-L";
                }
                else{
                    s+=to_string(f+1);
                }
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_.at(i).t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_.at(i).sat.sat_.id<<" "<<s<<" AMBIGUITY INITIALIZED "<<bias[i];
            }
        }
    }

    /*
    * Function : Update ion-parameter
    */
    void cPppSolver::IonUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumIon()<=0) return;

        int ii;
        tSatInfoUnit* sat_info= nullptr;
        double ion;
        for(int i=0;i<MAX_SAT_NUM;i++){
            ii=para_.IndexIon(i+1);
            if((full_x_[ii]!=0.0&&previous_sat_info_[i].outc[0]>120)||full_x_[ii]==DIS_FLAG){
                full_x_[ii]=0.0;
            }
        }

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED) continue;
            ii=para_.IndexIon(sat_info->sat.sat_.no);
            if(full_x_[ii]==0.0||reset_flag_){
                ii=para_.IndexIon(sat_info->sat.sat_.no);
                if(sat_info->cor_P[0]==0.0||sat_info->cor_P[1]==0.0){
                    sat_info->stat=SAT_NO_USE;
                    continue;
                }
                ion=(sat_info->cor_P[0]-sat_info->cor_P[1])/(1.0-SQR(sat_info->lam[1]-sat_info->lam[0]));
                InitX(ion,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            }
            else{
                tt=spp_solver_->igcoupled_sol_.t_tag.TimeDiff(previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
                double sinel=sin(MAX(sat_info->el_az[0],10.0*D2R));
                full_Px_(ii,ii)+=C.gnssC.ait_psd[1]/sinel*fabs(tt);
//                InitX(full_x_[ii],SQR(60.0),ii,full_x_.data(),full_Px_.data());

            }
        }
    }

    void cPppSolver::TrpUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumTrp()<=0) return;

        int it=para_.IndexTrp();
        if(full_x_[it]==0.0||reset_flag_){
            full_x_[it]=0.175;
            full_Px_(it,it)=SQR(0.3);
            if(C.gnssC.trp_opt==TRP_EST_GRAD){
                full_x_[it+1]=full_x_[it+2]=1E-6;
                full_Px_(it+1,it+1)=full_Px_(it+2,it+2)=SQR(0.01);
            }
            return;
        }
        else{
            full_Px_(it,it)+=C.gnssC.ait_psd[2]*tt;
            if(C.gnssC.trp_opt==TRP_EST_GRAD){
                full_Px_(it+1,it+1)+=10E-14*tt;
                full_Px_(it+2,it+2)+=10E-14*tt;
            }
        }
    }

    void cPppSolver::GloIfcbUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumGloIfcb()<=0) return;

        int ii,j;
        if(C.gnssC.glo_ifcb_opt==GLO_IFCB_LNF){
            ii=para_.IndexGloIfcb(1);
            if(full_x_[ii]==0.0||reset_flag_) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_QUAD){
            ii=para_.IndexGloIfcb(1);
            if(full_x_[ii]==0.0||reset_flag_) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            ii=para_.IndexGloIfcb(2);
            if(full_x_[ii]==0.0||reset_flag_) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT){
            double dtr=spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GLO];
            int ic=para_.IndexClk(SYS_INDEX_GLO);
            InitX(0.0,0.0,ic,full_x_.data(),full_Px_.data());
            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                if(epoch_sat_info_collect_[i].stat!=SAT_USED) continue;
                ii=para_.IndexGloIfcb(epoch_sat_info_collect_[i].sat.sat_.prn);
                if(full_x_[ii]==0.0||reset_flag_){
                    InitX(dtr,SQR(60.0),ii,full_x_.data(),full_Px_.data());
                }
                else full_x_(ii,ii)+=SQR(0.001)*tt;
            }
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1FRQ){
            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                if(epoch_sat_info_collect_[i].stat!=SAT_USED) continue;
                for(j=1;j<14;j++){
                    if(GLO_FRQ_NUM[j-1]==nav_.brd_glo_eph[epoch_sat_info_collect_[i].brd_eph_index].frq) break;
                }
                ii=para_.IndexGloIfcb(j);
                if(full_x_[ii]==0.0||reset_flag_) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            }
        }
    }

    void cPppSolver::ClkUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumClk()<=0) return;
        int ic=para_.IndexClk(SYS_INDEX_GPS);
        bool init=false;
        double psd=0.0;
        int num_sys=NSYS;
        if(ppp_conf_.gnssC.est_bd3_isb) num_sys+=1;

        if(SQR(full_x_[ic])+SQR(full_x_[ic+1])+SQR(full_x_[ic+2])+SQR(full_x_[ic+3])+SQR(full_x_[ic+4])<=0.0||reset_flag_){
            init=true;
        }

        if(init){
            for(int i=0;i<num_sys;i++){
                // GPS_clock, GC_ISB, GE_ISB, GR_ISB, GJ_ISB
                if(spp_solver_->igcoupled_sol_.clk_error[i]==0.0) continue;
                InitX(spp_solver_->igcoupled_sol_.clk_error[i],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
            }
        }
        else{
            if(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GPS]!=0.0){
                // GPS clock error reinitialize
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GPS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // BDS,GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_BDS;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i]==0.0){
                        InitX(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i-SYS_INDEX_GPS,ic+i-SYS_INDEX_GPS)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_BDS]!=0.0){
                // BDS clock error reinitialize
                ic+=SYS_INDEX_BDS;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_BDS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GAL;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i-SYS_INDEX_BDS]==0.0){
                        InitX(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i-SYS_INDEX_BDS,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i-SYS_INDEX_BDS,ic+i-SYS_INDEX_BDS)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_GAL]!=0.0){
                // GAL clock error reinitialize
                ic+=SYS_INDEX_GAL;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GAL],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GLO;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i-SYS_INDEX_GAL,ic+i-SYS_INDEX_GAL)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_GLO]!=0.0){
                // GLO clock error reinitialize
                ic+=SYS_INDEX_GLO;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GLO],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // QZS ISB modeling
                for(int i=SYS_INDEX_QZS;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i-SYS_INDEX_GLO,ic+i-SYS_INDEX_GLO)+=SQR(psd)*tt;
                }
                return ;
            }
        }
    }

    void cPppSolver::PhaseClkUpdate(tIGCOUPLEDConf C, double tt) {
        if(para_.NumPhaseClk()<=0) return;
        int ic=para_.IndexPhaseClk(SYS_INDEX_GPS);
        bool init=false;
        double psd=0.0;
        int num_sys=NSYS;
        if(ppp_conf_.gnssC.est_bd3_isb) num_sys+=1;

        if(SQR(full_x_[ic])+SQR(full_x_[ic+1])+SQR(full_x_[ic+2])+SQR(full_x_[ic+3])+SQR(full_x_[ic+4])<=0.0||reset_flag_){
            init=true;
        }

        if(init){
            for(int i=0;i<num_sys;i++){
                // GPS_clock, GC_ISB, GE_ISB, GR_ISB, GJ_ISB
                if(spp_solver_->igcoupled_sol_.clk_error[i]==0.0) continue;
                InitX(spp_solver_->igcoupled_sol_.clk_error[i],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
            }
        }
        else{
            if(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GPS]!=0.0){
                // GPS clock error reinitialize
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GPS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // BDS,GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_BDS;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i]==0.0){
                        InitX(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i-SYS_INDEX_GPS,ic+i-SYS_INDEX_GPS)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_BDS]!=0.0){
                // BDS clock error reinitialize
                ic+=SYS_INDEX_BDS;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_BDS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GAL;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i-SYS_INDEX_BDS]==0.0){
                        InitX(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i-SYS_INDEX_BDS,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i-SYS_INDEX_BDS,ic+i-SYS_INDEX_BDS)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_GAL]!=0.0){
                // GAL clock error reinitialize
                ic+=SYS_INDEX_GAL;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GAL],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GLO;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i-SYS_INDEX_GAL,ic+i-SYS_INDEX_GAL)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(igcoupled_sol_.clk_error[SYS_INDEX_GLO]!=0.0){
                // GLO clock error reinitialize
                ic+=SYS_INDEX_GLO;
                InitX(spp_solver_->igcoupled_sol_.clk_error[SYS_INDEX_GLO],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // QZS ISB modeling
                for(int i=SYS_INDEX_QZS;i<num_sys;i++){
                    if(spp_solver_->igcoupled_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i-SYS_INDEX_GLO,ic+i-SYS_INDEX_GLO)+=SQR(psd)*tt;
                }
                return ;
            }
        }
    }

    void cPppSolver::IfbUpdate(tIGCOUPLEDConf C, double tt) {
        if(para_.NumIfb()<=0) return;

        int i,ifcb=para_.IndexIfb(SYS_INDEX_GPS);

        for(i=0;i<NSYS;i++){
            ifcb+=i;
            if(sys_mask_[i]){
                if(full_x_[ifcb]==0.0||reset_flag_){
                    InitX(0.01,SQR(60),ifcb,full_x_.data(),full_Px_.data());
                }
                else{
                    full_Px_(ifcb,ifcb)+=10E-8*tt;
                }
            }
        }
    }

    void cPppSolver::PosUpdate(tIGCOUPLEDConf C) {
        if(tc_mode_) {
            return;
        }
        int ip=para_.IndexPos();
        Eigen::Vector3d var(SQR(60.0),SQR(60.0),SQR(60.0));

        if((SQR(full_x_[0])+SQR(full_x_[1])+SQR(full_x_[2]))==0.0||reset_flag_){
            for(int i=0;i<3;i++) full_x_[i]=spp_solver_->igcoupled_sol_.pos[i];
            full_Px_.block<3,3>(0,0)=var.asDiagonal();
            return;
        }

        if(C.mode_opt==MODE_OPT_STATIC){
            // nothing
        }
        else if(C.mode_opt==MODE_OPT_KINE_SIM||C.mode_opt==MODE_OPT_KINEMATIC||C.mode_opt==MODE_OPT_PPP){
            for(int i=0;i<3;i++) full_x_[i]=spp_solver_->igcoupled_sol_.pos[i];
            full_Px_.block<3,3>(0,0)=var.asDiagonal();
        }
    }


    void cPppSolver::StateTimeUpdate(tIGCOUPLEDConf C) {
        double tt=fabs(spp_solver_->igcoupled_sol_.t_tag.TimeDiff(igcoupled_sol_.t_tag.t_));
        if(C.gnssC.restart_gap>0&&epoch_sat_info_collect_[0].t_tag.TimeDiff(filter_start_.t_)>C.gnssC.restart_gap*3600.0){
            reset_flag_=true;
            filter_start_=epoch_sat_info_collect_[0].t_tag;
        }

        //position
        PosUpdate(C);

        //clock
        ClkUpdate(C,tt);

        //phase clock
        PhaseClkUpdate(C,tt);

        //dcb

        //ifb
        IfbUpdate(C,tt);

        //glo ifcb
        GloIfcbUpdate(C,tt);

        //trp
        TrpUpdate(C,tt);

        //ion
        IonUpdate(C,tt);

        //amb
        AmbUpdate(C,tt);

        if(reset_flag_) reset_flag_=false;
    }

    void cPppSolver::InitSppSolver() {
        spp_solver_->epoch_idx_=epoch_idx_;
        spp_solver_->epoch_sat_info_collect_=epoch_sat_info_collect_;
        spp_solver_->cur_imu_info_=cur_imu_info_;
        spp_solver_->tc_mode_=tc_mode_;
    }

    void cPppSolver::Spp2Ppp() {
        spp_solver_->SolutionUpdate();
        epoch_sat_info_collect_.clear();
        epoch_sat_info_collect_=spp_solver_->epoch_sat_info_collect_;
        spp_solver_->epoch_sat_info_collect_.clear();
        for(int i=0;i<para_.NumClk();i++) sys_mask_[i]=spp_solver_->sys_mask_[i];
    }

    /*
    * Function : PPP cycle slip detect 
    */
    void cPppSolver::PppCycSlip(tIGCOUPLEDConf C) {
        tSatInfoUnit* sat_info= nullptr;
        cTime t=igcoupled_sol_.t_tag;
        double dt=C.gnssC.sample_rate;

        if(t.t_.long_time!=0.0) dt=epoch_sat_info_collect_[0].t_tag.TimeDiff(t.t_);

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);

            // LLI 探测周跳
            //gnss_obs_operator_.LliCycleSlip(C,*sat_info,previous_sat_info_[sat_info->sat.sat_.no-1],para_.GetGnssUsedFrqs(),dt,REC_ROVER);

            if(ppp_conf_.gnssC.frq_opt==FRQ_SINGLE){
                // 单频使用伪距－相位
                gnss_obs_operator_.CodeMinuPhase(C,*sat_info,previous_sat_info_[sat_info->sat.sat_.no-1]);
            }
            else{//多频使用MW-GF
                gnss_obs_operator_.MwCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, nullptr,previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
                gnss_obs_operator_.GfCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, nullptr);
                gnss_obs_operator_.SmoothMw(C,sat_info, nullptr);
            }
        }
    }

    /*
    * Function : PPP clock jump detect and repair
    */
    void cPppSolver::ClockJumpRepair() {
        int i;
        tSatInfoUnit sat_info,pre_sat_info;
        double d1,d2,d3,d4,delta0=0,delta1=0;
        int cj_num=0,valid_sat_num=0;

        for(i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=epoch_sat_info_collect_[i];
            if(sat_info.sat.sat_.sys!=SYS_GPS) continue;

            if(sat_info.raw_P[0]*sat_info.raw_P[1]*sat_info.raw_L[0]*sat_info.raw_L[1]==0.0) continue;
            pre_sat_info=previous_sat_info_[sat_info.sat.sat_.no-1];

            if(pre_sat_info.raw_P[0]*pre_sat_info.raw_P[1]*pre_sat_info.raw_L[0]*pre_sat_info.raw_L[1]==0.0) continue;

            valid_sat_num++;
            d1=sat_info.raw_P[0]-pre_sat_info.raw_P[0];
            d2=sat_info.raw_P[1]-pre_sat_info.raw_P[1];
            d3=(sat_info.raw_L[0]-pre_sat_info.raw_L[0])*sat_info.lam[0];
            d4=(sat_info.raw_L[1]-pre_sat_info.raw_L[1])*sat_info.lam[1];

            if(fabs(d1-d3)>290000){
                delta0+=d1-d3;
                delta1+=d2-d4;
                cj_num++;
            }
        }

        double cj_f1,cj_f2,clk_jump=0.0;
        if(cj_num!=0&&cj_num==valid_sat_num){
            d1=delta0/cj_num;
            d2=delta1/cj_num;
            cj_f1=d1/CLIGHT*1000.0;
            cj_f2=Round(cj_f1);
            if(fabs(cj_f1-cj_f2)<2.5E-2){
                clk_jump=(int)cj_f2;
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" RECEIVER CLOCK JUMP";
            }
        }

        for(i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=epoch_sat_info_collect_[i];
            if(sat_info.sat.sat_.sys!=SYS_GPS) continue;

            if(sat_info.raw_L[0]!=0.0){
                epoch_sat_info_collect_[i].raw_L[0]+=clk_jump*CLIGHT/1000.0/sat_info.lam[0];
            }
            if(sat_info.raw_L[1]!=0.0){
                epoch_sat_info_collect_[i].raw_L[1]+=clk_jump*CLIGHT/1000.0/sat_info.lam[1];
            }

        }
    }

    /*
    * Function : GNSS observation residuals for PPP
    * -Args :
    *        int       post                     I                    index
    *        tIGCOUPLEDConf C                        I                    configure file
    *        double    x                        I                    parameters to be estimated
    *        Vector3d  re                       I                    position of imu
    * -Return :
    *        true : 
    *        false : 
    */
    int cPppSolver::GnssObsRes(int post, tIGCOUPLEDConf C, double *x,Vector3d re) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        for(int i=0;i<NSYS;i++) sys_mask_[i]=0;

        bool have_larger_res=false;
        char buff[MAX_BUFF]={'\0'};
        vector<double>H,omcs,meas_var_vec,idxs,types,larger_omcs;
        vector<int>frqs;
        double omc,r,meas,meas_var,ion=0.0,amb=0.0,fact;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_xyz;
        if(tc_mode_) rover_xyz=re;
        else{
            rover_xyz<<x[0],x[1],x[2];
        }

        cTime t=epoch_sat_info_collect_[0].t_tag;
        Vector3d dr;
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<(post==0?" PRIOR(":(" POST("))<<post<<")"<<" RESIDUAL, PRIOR ROVER COORDINATE "<<rover_xyz.transpose();

        if(C.gnssC.tid_opt>TID_OFF){//潮汐改正
            gnss_err_corr_.tid_model_.TidCorr(t.Gpst2Utc(),rover_xyz,dr);
            rover_xyz+=dr;
        }

        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return 0;
        if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF) num_used_frq=1;//三频无电离层组合

        double sat_clk,cbias,trp_del,alpha,beta;
        int m,idx_isb,idx_gps_clk,idx_ifcb,ii=0,ia;
        double base_rec_clk=0.0,isb=0.0,rec_clk=0.0,glo_ifcb=0.0;

        bool irc_flag=(C.gnssC.ar_mode==AR_PPP_AR&&C.gnssC.ar_prod==AR_PROD_IRC_CNES)?true:false;//模糊度固定是否使用IRC

        for(int i=0;i<epoch_sat_info_collect_.size();i++){//每颗卫星

            sat_info=&epoch_sat_info_collect_.at(i);

            if(sat_info->stat!=SAT_USED) continue;

            if((r=GeoDist(sat_info->pre_pos,rover_xyz,sig_vec))<=0.0){
                sat_info->stat=SAT_NO_USE;
                continue;
            }

            if(SatElAz(rover_blh,sig_vec,sat_info->el_az)<C.gnssC.ele_min*D2R){
                sat_info->stat=SAT_LOW_EL;
                continue;
            }

            //误差改正
            sat_info->sagnac=gnss_err_corr_.SagnacCorr(sat_info->pre_pos,rover_xyz);
            sat_info->shapiro=gnss_err_corr_.ShapiroCorr(sat_info->sat.sat_.sys,sat_info->pre_pos,rover_xyz);
            gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.trp_model_.GetTrpError(0.0,x,para_.IndexTrp());
            gnss_err_corr_.trp_model_.UpdateSatInfo();
            gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.ion_model_.GetIonError();
            gnss_err_corr_.ion_model_.UpdateSatInfo();

            int obs_type,frq;
            for(int f=0;f<num_used_frq*num_used_obs_type;f++){//每种观测值、每个频率
                frq=f/num_used_obs_type;
                obs_type=f%num_used_obs_type;

                if(irc_flag&&obs_type==GNSS_OBS_PHASE){
                    idx_gps_clk=para_.IndexPhaseClk(SYS_INDEX_GPS);
                }
                else{
                    idx_gps_clk=para_.IndexClk(SYS_INDEX_GPS);
                }

                //ionosphere-free
                if(C.gnssC.ion_opt==ION_IF){
                    if(C.gnssC.frq_opt==FRQ_TRIPLE){
                        if(sat_info->tf_if_idx[0]==1){
                            // L1_L2_L5
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[0]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[0]:sat_info->cor_D[0]);
                        }
                        else if(sat_info->tf_if_idx[1]==1){
                            // L1_L2
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[1]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[1]:sat_info->cor_D[1]);
                        }
                        else if(sat_info->tf_if_idx[2]==1){
                            // L1_L5
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[2]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[2]:sat_info->cor_D[2]);
                        }
                    }
                    else meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                }
                else if(C.gnssC.ion_opt==ION_IF_DUAL&&C.gnssC.frq_opt==FRQ_TRIPLE){
                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                }
                else{
                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_L[frq]:sat_info->cor_D[frq]);
                }

                if(meas==0.0){
#if 0
                    sat_info->stat=GNSS_OBS_CODE?sat_info->stat=SAT_NO_CODE:sat_info->stat=SAT_NO_PHASE;
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" MISSING  OBS"<<(obs_type==GNSS_OBS_CODE?"P":"L")<<frq+1;
#endif
                    continue;
                }
                ion=0.0;
                amb=0.0;

                //clock and inter-system bias
                int bd3_flag=(sat_info->sat.sat_.sys==SYS_BDS&&sat_info->sat.sat_.prn>18)?1:0;
                if(tc_mode_){
                    if(!post) for(int j=0;j<num_full_x_;j++) H.push_back(j<3?sig_vec[j]:0.0);
                }
                else{
                    if(!post) for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-sig_vec[j]:0.0);
                }

                if(ppp_conf_.gnssC.est_bd3_isb){
                    for(m=0;m<para_.NumClk();m++){
                        if(x[m+idx_gps_clk]!=0.0) break;
                    }
                    if(sat_info->sat.sat_.sys_idx==m&&!bd3_flag){
                        //clock
                        if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) break;
                        if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                        base_rec_clk=x[m+idx_gps_clk];
                    }
                    else{
                        for(idx_isb=0;idx_isb<para_.NumClk();idx_isb++){
                            if(idx_isb==m){
                                if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                                if(idx_isb==SYS_INDEX_BDS&&bd3_flag){
                                    if(!post) H[NSYS+idx_gps_clk+num_L_*num_full_x_]=1.0;
                                    isb=x[NSYS+idx_gps_clk];
                                }
                            }
                            else{
                                if(idx_isb!=sat_info->sat.sat_.sys_idx) continue;
                                if(x[idx_gps_clk+idx_isb]==0.0) continue;
                                if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) continue;
                                if(bd3_flag){
                                    if(!post) H[NSYS+idx_gps_clk+num_L_*num_full_x_]=1.0;
                                    isb=x[NSYS+idx_gps_clk];
                                }
                                else{
                                    if(!post) H[idx_isb+idx_gps_clk+num_L_*num_full_x_]=1.0;
                                    isb=x[idx_isb+idx_gps_clk];
                                }
                            }
                        }
                    }
                }
                else{
                    for(m=0;m<para_.NumClk();m++){
                        if(x[m+idx_gps_clk]!=0.0) break;
                    }
                    if(sat_info->sat.sat_.sys_idx==m){
                        if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) break;
                        if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                        base_rec_clk=x[m+idx_gps_clk];
                    }
                    else{
                        for(idx_isb=0;idx_isb<para_.NumClk();idx_isb++){
                            if(idx_isb==m){
                                if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                            }
                            else{
                                if(idx_isb!=sat_info->sat.sat_.sys_idx) continue;
                                if(x[idx_gps_clk+idx_isb]==0.0) continue;
                                if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) continue;
                                if(!post) H[idx_isb+idx_gps_clk+num_L_*num_full_x_]=1.0;
                                isb=x[idx_isb+idx_gps_clk];
                            }
                        }
                    }
                }

                rec_clk=base_rec_clk+isb;
#if 1
                // inter-frequency clock bias
                     // TF-UC-L5                                                                       //  TF-IF-DUAL L1_L5
                if((C.gnssC.frq_opt>=FRQ_TRIPLE&&frq==FRQ_TRIPLE&&obs_type==GNSS_OBS_CODE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF_DUAL&&frq==1&&obs_type==GNSS_OBS_CODE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF&&(sat_info->tf_if_idx[1]||sat_info->tf_if_idx[2])&&obs_type==GNSS_OBS_CODE)){
                    int ifcb=para_.IndexIfb(sat_info->sat.sat_.sys_idx);
                    rec_clk+=x[ifcb];
                    if(!post){
                        H[ifcb+num_L_*num_full_x_]=1.0;
                    }
                }
#endif
                // glonass inter-frequency code bias
                if(sat_info->sat.sat_.sys==SYS_GLO&&obs_type==GNSS_OBS_CODE){
                    int jj;
                    double frq=nav_.brd_glo_eph[sat_info->brd_eph_index].frq;
                    if(C.gnssC.glo_ifcb_opt==GLO_IFCB_LNF){
                        idx_ifcb=para_.IndexGloIfcb(1);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=frq;
                        glo_ifcb=frq*full_x_[idx_ifcb];
                        rec_clk+=glo_ifcb;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_QUAD){
                        idx_ifcb=para_.IndexGloIfcb(1);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=frq;
                        glo_ifcb=frq*full_x_[idx_ifcb];
                        idx_ifcb=para_.IndexGloIfcb(2);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=SQR(frq);
                        glo_ifcb+=SQR(frq)*full_x_[idx_ifcb];
                        rec_clk+=glo_ifcb;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT){
                        idx_ifcb=para_.IndexGloIfcb(sat_info->sat.sat_.prn);
                        rec_clk=full_x_[idx_ifcb];
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=1.0;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1FRQ){
                        for(jj=1;jj<14;jj++){
                            if(frq==GLO_FRQ_NUM[jj]) break;
                        }
                        idx_ifcb=para_.IndexGloIfcb(jj);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=1.0;
                    }
                }

                // tropospheric delay
                int it=para_.IndexTrp();
                if(!post) H[it+num_L_*num_full_x_]=sat_info->trp_wet_delay[1];
                if(C.gnssC.trp_opt==TRP_EST_GRAD&&!post){
                    H[it+1+num_L_*num_full_x_]=sat_info->trp_wet_delay[2];
                    H[it+2+num_L_*num_full_x_]=sat_info->trp_wet_delay[3];
                }

                // ionospheric delay
                if(C.gnssC.ion_opt>ION_IF_DUAL){
                    ii=para_.IndexIon(sat_info->sat.sat_.no);
                    fact=SQR(sat_info->frq[0]/sat_info->frq[frq])*(obs_type==GNSS_OBS_PHASE?(-1.0):1.0);
                    ion=x[ii]*fact;

                    if(x[ii]==0.0){
                        CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"NO VALUE FOR IONOSPHERIC PARAMETER";
                        continue;
                    }
                    if(!post) H[ii+num_L_*num_full_x_]=fact;
                }

                //Ambiguity
                if(obs_type==GNSS_OBS_PHASE||(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF&&obs_type==GNSS_OBS_CODE)){
                    if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                        if(sat_info->tf_if_idx[0]==1){
                            ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                        }
                        else if(sat_info->tf_if_idx[0]==0&&sat_info->tf_if_idx[1]!=0){
                            ia=para_.IndexAmb(1,sat_info->sat.sat_.no);
                        }
                        else if(sat_info->tf_if_idx[0]==0&&sat_info->tf_if_idx[2]!=0){
                            ia=para_.IndexAmb(2,sat_info->sat.sat_.no);
                        }
                    }
                    else{
                        ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                    }

                    amb=x[ia];

                    if(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF&&obs_type==GNSS_OBS_CODE) amb*=0.5;
                    if(amb==0.0) {
                        CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"NO VALUE FOR AMBIGUITY PARAMETER";
                        continue;
                    }
                    if(!post) H[ia+num_L_*num_full_x_]=1.0;
                }

                sat_clk=sat_info->pre_clk[0]*CLIGHT;
                if(irc_flag&&obs_type==GNSS_OBS_CODE) {
                    sat_clk=sat_info->pre_pr_clk[0]*CLIGHT;
                }

                trp_del=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];
                cbias=0.0;
                alpha=SQR(sat_info->frq[0])/(SQR(sat_info->frq[0])-SQR(sat_info->frq[1]));
                beta=-SQR(sat_info->frq[1])/(SQR(sat_info->frq[0])-SQR(sat_info->frq[1]));

                if(obs_type==GNSS_OBS_CODE){
                    if(C.gnssC.ion_opt==ION_IF){
                        cbias=alpha*sat_info->code_bias[0]+beta*sat_info->code_bias[1];
                    }
                    else cbias=sat_info->code_bias[f];
                }

                omc=meas-(r+sat_info->sagnac+rec_clk-sat_clk+trp_del+ion+amb-sat_info->shapiro);
                meas_var=GnssMeasVar(C,obs_type==0?GNSS_OBS_CODE:(obs_type==1?GNSS_OBS_PHASE:GNSS_OBS_DOPPLER),*sat_info)+sat_info->trp_var;
                if(obs_type==GNSS_OBS_CODE) meas_var*=previous_sat_info_[sat_info->sat.sat_.no-1].c_var_factor[frq];
                else if(obs_type==GNSS_OBS_PHASE) meas_var*=previous_sat_info_[sat_info->sat.sat_.no-1].p_var_factor[frq];

                if(obs_type==GNSS_OBS_CODE){
                    sprintf(buff,"%4d el=%3.1f var=%9.5f omc=%8.4f cor_obs=%14.3f range=%14.3f rec_clk=%12.3f sat_clk=%12.3f trp=%7.3f ion=%7.3f shaprio=%7.3f cbias=%7.3f",
                            epoch_idx_,sat_info->el_az[0]*R2D,meas_var,omc,meas,r+sat_info->sagnac,rec_clk,sat_clk,trp_del,ion,sat_info->shapiro,cbias);
                }
                else if(obs_type==GNSS_OBS_PHASE){
                    sprintf(buff,"%4d el=%3.1f var=%9.5f omc=%8.4f cor_obs=%14.3f range=%14.3f rec_clk=%12.3f sat_clk=%12.3f trp=%7.3f ion=%7.3f shaprio=%7.3f amb  =%7.3f phw=%7.3f",
                            epoch_idx_,sat_info->el_az[0]*R2D,meas_var,omc,meas,r+sat_info->sagnac,rec_clk,sat_clk,trp_del,ion,sat_info->shapiro,amb,sat_info->phase_wp);
                }

                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->sat.sat_.id<<" "<<(C.gnssC.ion_opt==ION_IF?((obs_type==GNSS_OBS_CODE?"IF_P":"IF_L")):((obs_type==GNSS_OBS_CODE?"P":"L")))<<frq+1<<" "<<buff;

//                if(!post&&fabs(omc)>((sat_info->sat.sat_.sys==SYS_GAL||sat_info->sat.sat_.sys==SYS_BDS)?100.0:C.gnssC.max_prior)){
//                    if(obs_type==GNSS_OBS_CODE) sat_info->stat=SAT_PRI_RES_C;
//                    else sat_info->stat=SAT_PRI_RES_P;
//                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" DETECTING OUTLIER IN PRIOR RESIDUALS "<<(obs_type==GNSS_OBS_CODE?"P":"L")<<frq+1<<"= "<<omc;
//                    return 0;
//                }

                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var);
                double el=sat_info->el_az[0];
                if(post&&fabs(omc)>sqrt(meas_var)/sin(el)*3.0){
                    have_larger_res=true;
                    larger_omcs.push_back(omc);
                    idxs.push_back(i);
                    frqs.push_back(frq);
                    types.push_back(obs_type);
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<"("<<epoch_idx_<<") "<<sat_info->sat.sat_.id<<" "<<(obs_type==GNSS_OBS_CODE?"P":"L")
                                <<frq+1<<" LARGER POST RESIDUAL res="<<omc<<" thres="<<sqrt(meas_var)*5.0;
                }

                if(!post){
                    if(obs_type==GNSS_OBS_CODE)  sat_info->prior_res_P[frq]=omc;
                    if(obs_type==GNSS_OBS_PHASE){
                        sat_info->prior_res_L[frq]=omc;
                    }

                }
                else{
                    if(obs_type==GNSS_OBS_CODE){
                        if(ppp_conf_.gnssC.ion_opt==ION_IF&&ppp_conf_.gnssC.frq_opt==FRQ_TRIPLE){
                            if(sat_info->tf_if_idx[0]==1){
                                sat_info->post_res_P[frq]=omc;
                            }
                            else if(sat_info->tf_if_idx[1]==1){
                                sat_info->post_res_P[1]=omc;
                            }
                        }
                        else sat_info->post_res_P[frq]=omc;
                    }
                    if(obs_type==GNSS_OBS_PHASE){
                        if(ppp_conf_.gnssC.ion_opt==ION_IF&&ppp_conf_.gnssC.frq_opt==FRQ_TRIPLE){
                            if(sat_info->tf_if_idx[0]==1){
                                sat_info->post_res_L[frq]=omc;
                                sat_info->float_amb[frq]=amb;
                            }
                            else if(sat_info->tf_if_idx[1]==1){
                                sat_info->post_res_L[1]=omc;
                                sat_info->float_amb[1]=amb;
                            }
                        }
                        else{
                            sat_info->post_res_L[frq]=omc;
                            sat_info->float_amb[frq]=amb;
                        }

                    }
                }

                sat_info->vsat[frq]=1;
                previous_sat_info_[sat_info->sat.sat_.no-1].vsat[frq]=1;

                vflag_.push_back((i<<8)|(obs_type<<4)|(frq));

                num_L_++;
                if(f==0) num_valid_sat_++;
            }
        }

        if(!post){
            omc_L_=Map<VectorXd>(omcs.data(),num_L_);
            H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
            Map<VectorXd> var_vec(meas_var_vec.data(),num_L_);
            R_=var_vec.asDiagonal();
        }

        if(post&&have_larger_res&&larger_omcs.size()>0){
            for(int i=0;i<larger_omcs.size();i++) larger_omcs[i]=fabs(larger_omcs[i]);
            auto max_omc=max_element(larger_omcs.begin(),larger_omcs.end());
            int idx_max_omc=distance(begin(larger_omcs),max_omc);
            double el=epoch_sat_info_collect_[idxs[idx_max_omc]].el_az[0]*R2D;
            epoch_sat_info_collect_[idxs[idx_max_omc]].stat=types[idx_max_omc]==GNSS_OBS_CODE?SAT_POS_RES_C:SAT_POS_RES_P;
            epoch_sat_info_collect_[idxs[idx_max_omc]].rejc[frqs[idx_max_omc]]++;
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[idxs[idx_max_omc]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idxs[idx_max_omc]].sat.sat_.id
                        <<" "<<(types[idx_max_omc]==GNSS_OBS_CODE?"P":"L")<<frqs[idx_max_omc]+1<<" EXCLUDED BY POST RESIDUAL res="<<larger_omcs[idx_max_omc]<<" el="<<el;
            have_larger_res=true;
        }

        if(!have_larger_res&&post&&num_valid_sat_>4&&ppp_conf_.gnssC.res_qc&&!(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF)){
            if(PppResidualQc(post,omcs,meas_var_vec)){
                have_larger_res=true;
            }
        }

        H.clear();omcs.clear();meas_var_vec.clear();types.clear();frqs.clear();idxs.clear();larger_omcs.clear();

        return post?have_larger_res:num_valid_sat_;
    }

    /*
    * Function : Reinitialize PPP solver without ambiguity
    */
    void cPppSolver::ReInitPppSolver(tIGCOUPLEDConf C) {
        ppp_conf_=C;
        para_=cParSetting(ppp_conf_);
        num_full_x_=para_.GetIGCOUPLEDPar(ppp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);

        num_real_x_fix_=para_.GetRealFixParNum(ppp_conf_);//不包括模糊度的参数个数
        real_x_fix_=VectorXd::Zero(num_real_x_fix_);
        real_Px_fix_=MatrixXd::Zero(num_real_x_fix_,num_real_x_fix_);
    }

    /*
    * Function : 
    * -Args :
    *        int            iter                   I                Number of irerations
    *        VectorXd       x                      I                Parameters to be estimated
    *        vector<int>    par_idx                O                the index of disabled parameters
    *        vector<double> back_values            O                the value of disabled parameters 
    * -Return :
    *         None
    */
    void cPppSolver::DisableX(int iter,VectorXd &x, vector<int> &par_idx, vector<double> &back_values) {

        int ia=0,ii=0,num_used_frq=para_.GetGnssUsedFrqs();//使用频率的个数
        int i,j;
        if(iter!=-1){
            par_idx.clear();
            back_values.clear();
        }

        if(0<iter&&iter<10){//迭代次数
            for(i=0;i<epoch_sat_info_collect_.size();i++){//每颗卫星
                for(j=0;j<num_used_frq;j++){//每个频率
                    if(!epoch_sat_info_collect_[i].vsat[j]){
                        if(ppp_conf_.gnssC.ion_opt>=ION_EST){//估计电离层参数
                            ii=para_.IndexIon(epoch_sat_info_collect_[i].sat.sat_.no);
                            back_values.push_back(x[ii]);
                            par_idx.push_back(ii);
                            x[ii]=DIS_FLAG;
                        }
                        ia=para_.IndexAmb(j,epoch_sat_info_collect_[i].sat.sat_.no);
                        back_values.push_back(x[ia]);
                        par_idx.push_back(ia);
                        x[ia]=DIS_FLAG;
                    }
                }
            }
        }

        if(iter==-1){
            for(i=0;i<par_idx.size();i++){
                x[par_idx[i]]=back_values[i];
            }
        }
    }

    /*
    * Function : PPP residual QC check
    * -Args :
    *        int            iter                  I              the number of iterations 
    *        vector<double> omcs                  I              observations minus calculations 
    *        vector<double> var                   I              variance
    * -Return :
    *         bool         
    */
    bool cPppSolver::PppResidualQc(int iter, vector<double> omcs, vector<double> var) {
        bool flag=false;
        int i,obs_type,f,idx_sat,sat_no,ia;
        vector<double>code_v,phase_v,norm_code_v,norm_phase_v;
        vector<int>idx_code_v,idx_phase_v;
        double el;
        double thres_code=3.0,thres_phase=0.03;
        double k0=1.5,k1=2.0;
        bool reduce_weight=false;

        if(iter>3) return false;

        for(i=0;i<vflag_.size();i++){
            obs_type=vflag_[i]>>4&0xF;
            idx_sat=vflag_[i]>>8&0xFF;
            if(obs_type==GNSS_OBS_CODE){//码伪距
                idx_code_v.push_back(i);
                code_v.push_back(fabs(omcs[i]));
                norm_code_v.push_back(fabs(omcs[i])/SQRT(var[i]));
            }
            else{//相位
                idx_phase_v.push_back(i);
                phase_v.push_back(fabs(omcs[i]));
                norm_phase_v.push_back(fabs(omcs[i])/SQRT(var[i]));
            }
        }

        bool code_flag=false;
        auto max_code_v=max_element(begin(code_v),end(code_v));//残差中的最大值
        int idx_max_code_v=distance(begin(code_v),max_code_v);//最大值的索引
        idx_sat=vflag_[idx_code_v[idx_max_code_v]]>>8&0xFF;//最大值残差的卫星索引
        el=epoch_sat_info_collect_[idx_sat].el_az[0];//对应的高度角

        if(epoch_sat_info_collect_[idx_sat].sat.sat_.sys==SYS_BDS){//如果是北斗卫星
            if(epoch_sat_info_collect_[idx_sat].sat.sat_.prn<=5){//GEO
                thres_code=10.0;
            }else{
                thres_code=6.0;
            }
        }

        if(*max_code_v>thres_code/sin(el)){
            code_flag=true;
        }

        bool norm_code_flag=false;
        auto max_norm_code_v=max_element(begin(norm_code_v),end(norm_code_v));
        int idx_max_norm_code_v=distance(begin(norm_code_v),max_norm_code_v);
        if(*max_norm_code_v>k1){
            norm_code_flag=true;
        }

        /* step0: exclude bad satellite due to bad pseudorange*/

        // step_01: check post-fit pseudorange residual,highest priority,exclude satellite
        if(code_flag){
            flag=true;
            idx_sat=vflag_[idx_code_v[idx_max_code_v]]>>8&0xFF;
            epoch_sat_info_collect_[idx_sat].stat=SAT_NO_USE;
            el=epoch_sat_info_collect_[idx_sat].el_az[0]*R2D;
            f=vflag_[idx_code_v[idx_max_code_v]]&0xF;
            if(iter==1) epoch_sat_info_collect_[idx_sat].rejc[f]++;
            sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
            ia=para_.IndexAmb(f,sat_no);
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" EXCULDED BY CODE RESIDUAL, el="<<el;
            return flag;
        }

        // step_02: check post-fit norm pseudorange residual, if true, exclude satellite
        if(norm_code_flag){
            flag=true;
            idx_sat=vflag_[idx_code_v[idx_max_norm_code_v]]>>8&0xFF;
            epoch_sat_info_collect_[idx_sat].stat=SAT_NO_USE;
            el=epoch_sat_info_collect_[idx_sat].el_az[0]*R2D;
            f=vflag_[idx_code_v[idx_max_norm_code_v]]&0xF;
            if(iter==1) epoch_sat_info_collect_[idx_sat].rejc[f]++;
            sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
            ia=para_.IndexAmb(f,sat_no);
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" EXCULDED BY NORM CODE RESIDUAL, el="<<el;
            return flag;
        }

        // step_03: if pass step_01 and step_02, judge for reducing the weight according norm pseudorange residual
        double fact=k0/(*max_norm_code_v)*SQR((k1-*max_norm_code_v)/(k1-k0));
        if(k0<=*max_norm_code_v&&*max_norm_code_v<=k1){
            obs_type=vflag_[idx_max_norm_code_v]>>4&0xF;
            idx_sat=vflag_[idx_max_norm_code_v]>>8&0xFF;
            f=vflag_[idx_code_v[idx_max_norm_code_v]]&0xF;
            previous_sat_info_[epoch_sat_info_collect_[idx_sat].sat.sat_.no-1].c_var_factor[f]*=1.0/fact;
            reduce_weight=true;
        }

        /* step1: exclude bad satellite due to bad phase*/
        // if check pseudorange OK
        if(phase_v.size()){

            // step_11: check post-fit phase residual
            auto max_phase_v=max_element(begin(phase_v),end(phase_v));
            int idx_max_phase_v=distance(begin(phase_v),max_phase_v);
            idx_sat=vflag_[idx_phase_v[idx_max_phase_v]]>>8&0xFF;
            f=vflag_[idx_phase_v[idx_max_phase_v]]&0xF;
            el=epoch_sat_info_collect_[idx_sat].el_az[0];
            if(epoch_sat_info_collect_[idx_sat].sat.sat_.sys==SYS_BDS){
                thres_phase=0.06;
            }

            if(*max_phase_v>thres_phase/sin(el)){
                int sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
                if(iter==1) epoch_sat_info_collect_[idx_sat].rejc_phase[f]++;
                if(iter==1) epoch_sat_info_collect_[idx_sat].rejc[f]++;
                int ia=para_.IndexAmb(f,sat_no);
                el=epoch_sat_info_collect_[idx_sat].el_az[0]*R2D;
                // over two consecutive epoch to reset amb,otherwise set the weight to zero
                if(epoch_sat_info_collect_[idx_sat].rejc_phase[f]>0){
                    InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
//                    previous_sat_info_[epoch_sat_info_collect_[idx_sat].sat.sat_.no-1].p_var_factor[f]*=10000;

                    epoch_sat_info_collect_[idx_sat].rejc_phase[f]=0;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" RESET AMB BY PHASE RESIDUAL, el="<<el;

                }
                else{
                    previous_sat_info_[epoch_sat_info_collect_[idx_sat].sat.sat_.no-1].p_var_factor[f]*=10000;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" REDUCE WEIGHT BY NORM PHASE RESIDUAL, el="<<el;
                }
                return true;
            }

            // step_12: check post-fit norm phase residual
            auto max_norm_phase_v=max_element(begin(norm_phase_v),end(norm_phase_v));
            int idx_max_norm_phase_v=distance(begin(norm_phase_v),max_norm_phase_v);
            idx_sat=vflag_[idx_phase_v[idx_max_norm_phase_v]]>>8&0xFF;
            el=epoch_sat_info_collect_[idx_sat].el_az[0]*R2D;
            f=vflag_[idx_phase_v[idx_max_norm_phase_v]]&0xF;
            if(*max_norm_phase_v>k1){
                int sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
                if(iter==1) epoch_sat_info_collect_[idx_sat].rejc_phase[f]++;
                if(iter==1) epoch_sat_info_collect_[idx_sat].rejc[f]++;
                int ia=para_.IndexAmb(f,sat_no);
                // over tow consecutive epoch to reset amb, otherwise set the weight to zero
                if(epoch_sat_info_collect_[idx_sat].rejc_phase[f]>1){
                    InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                    epoch_sat_info_collect_[idx_sat].rejc_phase[f]=0;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" RESET AMB BY NORM PHASE RESIDUAL, el="<<el;
                    return true;
                }
                else{
                    previous_sat_info_[epoch_sat_info_collect_[idx_sat].sat.sat_.no-1].p_var_factor[f]*=10000;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idx_sat].sat.sat_.id<<" REDUCE WEIGHT BY NORM PHASE RESIDUAL, el="<<el;
                    return true;
                }

            }

            // step_13: if pass step_11 and step_12, judge for reducing the weight according norm phase residual
            fact=k0/(*max_norm_phase_v)*SQR((k1-*max_norm_phase_v)/(k1-k0));
            if(k0<=*max_norm_phase_v&&*max_norm_phase_v<=k1){
                obs_type=vflag_[idx_max_norm_phase_v]>>4&0xF;
                idx_sat=vflag_[idx_max_norm_phase_v]>>8&0xFF;
                f=vflag_[idx_code_v[idx_max_norm_phase_v]]&0xF;
                previous_sat_info_[epoch_sat_info_collect_[idx_sat].sat.sat_.no-1].p_var_factor[f]*=1.0/fact;
                reduce_weight=true;
            }
        }
        return reduce_weight;
    }
    /*
    * Function : Resolve PPP ambiguity
    */
    bool cPppSolver::ResolvePppAmbRotRef(int nf, VectorXd& x,MatrixXd& P) {
        double elmask,var=0.0;
        int i,j,m=0,stat;
        if(num_valid_sat_<=0||ppp_conf_.gnssC.ion_opt!=ION_IF||nf<2) return false;
        elmask=ppp_conf_.gnssC.ar_el_mask>0.0?ppp_conf_.gnssC.ar_el_mask:ppp_conf_.gnssC.ele_min;
        vector<int>sat_no1,sat_no2,fix_wls;
        vector<double> res_wls;
        double res_wl;

        AverageLcAmb();

        for(i=0;i<3;i++) var+=SQRT(full_Px_(i,i));
        var=var/3.0;
        if(var>0.15){ //0.25
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"POSITIONING VARIANCE TOO LARGE FOR PPK AR var="<<var;
            igcoupled_sol_.num_ar_sat=0;
            return false;
        }

        int sat1,sat2,fix_wl_amb;
        tSatInfoUnit *sat_info_i=nullptr,*sat_info_j= nullptr;
        for(i=0;i<epoch_sat_info_collect_.size()-1;i++){
            sat_info_i=&epoch_sat_info_collect_.at(i);
            if(sat_info_i->sat.sat_.sys!=SYS_GPS) continue;
            for(j=i+1;j<epoch_sat_info_collect_.size();j++){
                sat_info_j=&epoch_sat_info_collect_.at(j);
                if(sat_info_j->sat.sat_.sys!=SYS_GPS) continue;
                if(!sat_info_i->vsat[0]||!sat_info_j->vsat[0]||sat_info_i->el_az[0]*R2D<elmask||sat_info_j->el_az[0]*R2D<elmask) continue;

                sat1=sat_info_i->sat.sat_.no;
                sat2=sat_info_j->sat.sat_.no;
                if(FixWlAmb(sat1,sat2,&fix_wl_amb,&res_wl)){
                    m++;
                    sat_no1.push_back(sat1);
                    sat_no2.push_back(sat2);
                    fix_wls.push_back(fix_wl_amb);
                    res_wls.push_back(res_wl);
                }
            }
        }

        if(ppp_conf_.gnssC.ar_mode==AR_PPP_AR&&m>4){
            if(ppp_conf_.gnssC.ar_prod==AR_PROD_FCB_SGG){
                stat=FixNlAmbILS_FCB(sat_no1.data(),sat_no2.data(),fix_wls.data(),res_wls.data(),m,x,P);
            }
            else if(ppp_conf_.gnssC.ar_prod==AR_PROD_IRC_CNES){
                stat=FixNlAmbILS_IRC(sat_no1.data(),sat_no2.data(),fix_wls.data(),res_wls.data(),m,x,P);
            }
            else
                return false;
        }
        else return false;

        return stat;
    }

    /*
    * Function : Fix PPP solution after fix ambiguity
    * -Args :
    *       int      sat1                      I                one satellite index
    *       int      sat2                      I                another satellite index
    *       double   NC                        I
    *       int      n                         I                number of parameters
    *       VectorXd x                         IO               parameters
    *       MatrixXd P                         IO               precision informations
    * -Return :
    *       bool
    */
    bool cPppSolver::FixPppSol(int *sat1, int *sat2, double *NC, int n,VectorXd& x,MatrixXd& P) {
        VectorXd v(n);
        MatrixXd H(num_full_x_,n),R(n,n);
        int i,j,k;
        if(n<=0) return false;

        H=MatrixXd::Zero(num_full_x_,n);
        R=MatrixXd::Zero(n,n);
        // constraints to fixed ambiguities
        for(i=0;i<n;i++){
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);
            v[i]=(NC[i]-(x[j]-x[k]));
            H(j,i)=1.0;
            H(k,i)=-1.0;
            R.data()[i+i*n]=SQR(0.0001);
        }

        kf_.Adjustment(v,H,R,x,P,n,num_full_x_);

        for(i=0;i<num_real_x_fix_;i++){
            real_x_fix_[i]=x[i];
            for(j=0;j<num_real_x_fix_;j++){
                real_Px_fix_.data()[i+j*num_real_x_fix_]=real_Px_fix_.data()[j+i*num_real_x_fix_]=P.data()[i+j*num_full_x_];
            }
        }

        for(i=0;i<n;i++){
            previous_sat_info_[sat1[i]-1].lc_amb.flags[sat2[i]-1]=1;
            previous_sat_info_[sat2[i]-1].lc_amb.flags[sat1[i]-1]=1;

        }
        return true;
    }

    /*
    * Function : Independent satellite pairs check
    */
    bool cPppSolver::AmbLinearDependCheck(int sat1, int sat2, int *flag, int *max_flg) {
        int i;

        if(flag[sat1-1]==0&&flag[sat2-1]==0){
            flag[sat1-1]=flag[sat2-1]=++(*max_flg);
        }
        else if(flag[sat1-1]==0&&flag[sat2-1]!=0){
            flag[sat1-1]=flag[sat2-1];
        }
        else if(flag[sat1-1]!=0&&flag[sat2-1]==0){
            flag[sat2-1]=flag[sat1-1];
        }
        else if(flag[sat1-1]>flag[sat2-1]){
            for(i=0;i<MAX_SAT_NUM;i++) if(flag[i]==flag[sat2-1]) flag[i]=flag[sat1-1];
        }
        else if(flag[sat1-1]<flag[sat2-1]){
            for(i=0;i<MAX_SAT_NUM;i++) if(flag[i]==flag[sat1-1]) flag[i]=flag[sat2-1];
        }
        else{
            return false;
        }

        return true;
    }

    /*
    * Function : find independent satellite pairs for ambiguity fixing
    */
    int cPppSolver::SelectAmb(int *sat1, int *sat2, double *N, double *var, int n) {
        int i,j=0,flag[MAX_SAT_NUM]={0},max_flg=0;

        for(i=0;i<n;i++) for(j=1;j<n-1;j++){
                if(var[j]>=var[j-1]) continue;
                SWAP_I(sat1[j],sat1[j-1]);
                SWAP_I(sat2[j],sat2[j-1]);
                SWAP_D(N[j],N[j-1]);
                SWAP_D(var[j],var[j-1]);
            }

        for(i=j=0;i<n;i++){
            if(!AmbLinearDependCheck(sat1[i],sat2[i],flag,&max_flg)) continue;
            sat1[j]=sat1[i];
            sat2[j]=sat2[i];
            N[j]=N[i];
            var[j++]=var[i];
        }
        return j;
    }

    /*
    * Function : Calculate the mean of linear combination ambiguity
    */
    void cPppSolver::AverageLcAmb() {
        tSatInfoUnit* sat_info= nullptr;
        tLcAmb* amb= nullptr;
        double LC1=0,LC2=0,LC3=0,var1,var2,var3;

        int i,j;
        for(i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            amb=&previous_sat_info_[sat_info->sat.sat_.no-1].lc_amb;

            if(sat_info->el_az[0]*R2D<ppp_conf_.gnssC.ele_min) continue;

            double bias=sat_info->cp_bias[0];
            double P1=sat_info->raw_P[0],P2=sat_info->raw_P[1];
            double L1=sat_info->raw_L[0],L2=sat_info->raw_L[1];
            double lam1=sat_info->lam[0],lam2=sat_info->lam[1];
            double MW=gnss_obs_operator_.GnssObsMwComb(P1,P2,L1,L2,lam1,lam2);
            LC1=gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,-1,0,*sat_info,GNSS_OBS_PHASE,nullptr)-gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,1,0,*sat_info,GNSS_OBS_CODE,&var1);
            //LC2
            //LC3
            if(sat_info->slip[0]||sat_info->slip[1]||sat_info->slip[2]||amb->n[0]==0.0||
               fabs(amb->t_tag[0].TimeDiff(epoch_sat_info_collect_[0].t_tag.t_))>ppp_conf_.gnssC.sample_rate){
                for(j=0;j<3;j++){
                    amb->n[j]=0.0;
                    amb->lc_amb[j]=0.0;
                    amb->var_amb[j]=0.0;
                }
                amb->fix_cnt=0;
                for(j=0;j<MAX_SAT_NUM;j++) amb->flags[j]=0;
            }
            if(LC1){
                amb->n[0]+=1;
                amb->lc_amb[0]+=(LC1-amb->lc_amb[0])/amb->n[0];
                amb->var_amb[0]+=(var1-amb->var_amb[0])/amb->n[0];
//                if(sat_info->sat.sat_.no==22) cout<<amb->lc_amb[0]<<endl;
            }
            if(LC2){
                amb->n[1]+=1;
                amb->lc_amb[1]+=(LC1-amb->lc_amb[1])/amb->n[1];
                amb->var_amb[1]+=(var1-amb->var_amb[1])/amb->n[1];
            }
            if(LC3){
                amb->n[2]+=1;
                amb->lc_amb[2]+=(LC1-amb->lc_amb[2])/amb->n[2];
                amb->var_amb[2]+=(var1-amb->var_amb[2])/amb->n[2];
            }
            amb->t_tag[0]=epoch_sat_info_collect_[0].t_tag;
        }
    }

    /*
    * Function : Fix the ambiguity of widelane
    */
    bool cPppSolver::FixWlAmb(int sat1,int sat2,int *sd_fix_wl,double *res_wl) {
        bool fix_flag=false;
        double sd_wl_amb,var_wl,lam_wl=gnss_obs_operator_.LinearCombLam(1,-1,0,previous_sat_info_[sat1-1]);
        tLcAmb *lc_amb1=&previous_sat_info_[sat1-1].lc_amb;
        tLcAmb *lc_amb2=&previous_sat_info_[sat2-1].lc_amb;

        if(!lc_amb1->n[0]||!lc_amb2->n[0]) return false;

        double C=1.0;
        if(ppp_conf_.gnssC.ar_prod==AR_PROD_IRC_CNES){
            C=-1.0;
        }else if(ppp_conf_.gnssC.ar_prod == AR_PROD_OSB_WUH || ppp_conf_.gnssC.ar_prod==AR_PROD_OSB_CNES){
            C=0.0;
        }

        double wl_bias1=nav_.wide_line_bias[sat1-1];
        double wl_bias2=nav_.wide_line_bias[sat2-1];
        double wl1=(lc_amb1->lc_amb[0])/lam_wl-wl_bias1*C;  //cycle
        double wl2=(lc_amb2->lc_amb[0])/lam_wl-wl_bias2*C;
        sd_wl_amb=wl1-wl2;

        *sd_fix_wl=(int)floor(sd_wl_amb+0.5);
        *res_wl=*sd_fix_wl-sd_wl_amb;
        var_wl=(lc_amb1->var_amb[0]/lc_amb1->n[0]+lc_amb2->var_amb[0]/lc_amb2->n[0])/SQR(lam_wl);

        bool fix_flag1=fabs(*sd_fix_wl-sd_wl_amb)<=0.25;
        bool fix_flag2=IntAmbConfFunc(*sd_fix_wl,sd_wl_amb,sqrt(var_wl))>=0.9999;
        fix_flag=fix_flag1&&fix_flag2;

        char buff[1024]={'\0'};
        sprintf(buff,"%s TO FIX SD NW AMB sd_wl_amb=%s-%s=(%6.3f-%6.3f)+(%4.3f-%4.3f)=%6.3f  round(sd_wl_amb)=%5d var=%4.2f %s",
                previous_sat_info_[sat1-1].t_tag.GetTimeStr(1).c_str(),previous_sat_info_[sat1-1].sat.sat_.id.c_str(),previous_sat_info_[sat2-1].sat.sat_.id.c_str(),
                wl1-wl_bias1,wl2-wl_bias2,wl_bias1,wl_bias2,sd_wl_amb,*sd_fix_wl,var_wl,fix_flag?"YES":"NO");
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;

        return fix_flag;
    }

    /*
    * Function : Fix narrowlane ambiguity using IRC
    */
    bool cPppSolver::FixNlAmbILS_IRC(int *sat1, int *sat2, int *fix_wls,double *res_wls, int n,VectorXd& x,MatrixXd& P) {
        double lam_nl=gnss_obs_operator_.LinearCombLam(1,1,0,previous_sat_info_[sat1[0]-1]);
        double lam1=0,lam2=0,alpha=0,beta=0,res_nl;
        vector<double>res_nls;
        int i,j,k,m=0,flag[MAX_SAT_NUM]={0},max_flag=0,info;
        char out_buff[1024]={'\0'};
        bool fix_flag=false;

        VectorXd B11(n),NC(n),s(2);
        MatrixXd N1(n,2),D(num_full_x_,n);
        B11=VectorXd::Zero(n);NC=VectorXd::Zero(n);
        N1=MatrixXd::Zero(n,2);D=MatrixXd::Zero(num_full_x_,n);

        for(i=0;i<n;i++){
            if(!AmbLinearDependCheck(sat1[i],sat2[i],flag,&max_flag)) continue;
            out_buff[0]='\0';
            lam1=previous_sat_info_[sat1[i]-1].lam[0];
            lam2=previous_sat_info_[sat1[i]-1].lam[1];
            alpha=SQR(lam2)/(SQR(lam2)-SQR(lam1));
            beta=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);

            B11[m]=(x[j]-x[k]+beta*lam2*fix_wls[i])/lam_nl;
            N1(m,0)=(int)floor(B11[m]+0.5);
            res_nl=N1(m,0)-B11[m];

            fix_flag=fabs(res_nl) <= 0.10;
            sprintf(out_buff,"%s TO FIX NL AMB sd_nl_amb(%s-%s)=sd_if_amb-sd_wl_amb=%6.3f-%6.3f=%6.3f round(sd_nl_amb)=%5d %s",
                    epoch_sat_info_collect_[0].t_tag.GetTimeStr(1).c_str(),previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),
                    (full_x_[j]-full_x_[k])/lam_nl,beta*lam2*fix_wls[i]/lam_nl,B11[m],(int)floor(B11[m]+0.5),fix_flag?"YES":"NO");
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<out_buff;

            if(!fix_flag) continue;

            D.data()[j+m*num_full_x_]=1.0/lam_nl;
            D.data()[k+m*num_full_x_]=-1.0/lam_nl;
            sat1[m]=sat1[i];
            sat2[m]=sat2[i];
            res_wls[m]=res_wls[i];
            fix_wls[m++]=fix_wls[i];
            res_nls.push_back(res_nl);
        }
        if(m<2) return false;

        MatrixXd E(m,num_full_x_),Q(m,m);
        MatMul("TN",m,num_full_x_,num_full_x_,1.0,D.data(),full_Px_.data(),0.0,E.data());
        MatMul("NN",m,m,num_full_x_,1.0,E.data(),D.data(),0.0,Q.data());

        VectorXd B1;
        B1=Map<VectorXd>(B11.data(),m); //float_amb

        char ss[20];
        string buff;
        for(i=0;i<B1.size();i++){
            sprintf(ss,"%5.2f  ",B1[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        for(i=0;i<B1.size();i++){
            sprintf(ss,"%5.2f  ",Q(i,i));
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   Q_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        N1=MatrixXd::Zero(n,2);
        if((info=lambda_.IntegerAmb(B1,Q,N1,m,2,s))){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP LAMBDA ERROR";
            return false;
        }

        for(i=0;i<m;i++){
            sprintf(ss,"%5.2f  ",N1.col(0)[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(1): "<<buff;
        buff.clear();
        ss[0]='\0';

        for(i=0;i<m;i++){
            sprintf(ss,"%5.2f  ",N1.col(1)[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(2): "<<buff;
        buff.clear();
        ss[0]='\0';


        if(s[0]<=0.0||s[1]/s[0]<ppp_conf_.gnssC.ar_thres[0]){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<< "PPP-AR VALIDATION FAILED, ratio="<<s[0]/s[1];
            return false;
        }


        igcoupled_sol_.ratio=MIN(s[1]/s[0],999.9);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP-AR VALIDATION OK, ratio="<<igcoupled_sol_.ratio;

        // nl to iono-free ambiguity
        for(i=0;i<m;i++){
            out_buff[0]='\0';
            NC[i]=alpha*lam1*N1(i,0)+beta*lam2*(N1(i,0)-fix_wls[i]);
//            NC[i]=lam_nl*(N1(i,0))-lam2*beta*(fix_wls[i]);
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);
            previous_sat_info_[sat2[i]-1].res_wl=res_wls[i];
            previous_sat_info_[sat2[i]-1].res_nl=res_nls[i];
            previous_sat_info_[sat2[i]-1].fix[0]=2;
            sprintf(out_buff,"REFACTOR SD_IF_AMB(%s-%s): float_if_amb=%6.3f-%6.3f=%6.3f fixed_if_amb=%6.3f",
                    previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),x[j],x[k],x[j]-x[k],NC[i]);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<out_buff;
        }

        return FixPppSol(sat1,sat2,NC.data(),m,x,P);
    }

    bool cPppSolver::MatchNlFcb(int sat1, int sat2, double *nl_fcb1,double *nl_fcb2) {
        cTime obs_time=epoch_sat_info_collect_[0].t_tag;
        int i;
        double fcb1=0.0,fcb2=0.0;
        bool stat=false;

        for(i=0;i<nav_.nl_fcbs.size();i++){
            if(obs_time.TimeDiff(nav_.nl_fcbs[i].ts.t_)>=0&&obs_time.TimeDiff(nav_.nl_fcbs[i].te.t_)<0){
                fcb1=nav_.nl_fcbs[i].nl_fcb[sat1-1];
                fcb2=nav_.nl_fcbs[i].nl_fcb[sat2-1];
                stat=true;
                break;
            }
        }
        if(nl_fcb1) *nl_fcb1=fcb1;
        if(nl_fcb2) *nl_fcb2=fcb2;
        return stat;
    }

    /*
    * Function : Fix narrowlane ambiguity using FCB
    */
    bool cPppSolver::FixNlAmbILS_FCB(int *sat1, int *sat2, int *fix_wls, double *res_wls, int n, VectorXd &x,MatrixXd &P) {
        double lam_nl=gnss_obs_operator_.LinearCombLam(1,1,0,previous_sat_info_[sat1[0]-1]);
        double lam1=0,lam2=0,alpha=0,beta=0,res_nl;
        vector<double>res_nls,sd_nl_fcbs;
        int i,j,k,m=0,flag[MAX_SAT_NUM]={0},max_flag=0,info,num_sys=1,round_nl;
        char buff[1024]={'\0'};
        double float_nl,nl_fcb1=0.0,nl_fcb2=0.0;
        bool fix_flag=false;

        vector<double> float_nl_ambs;
        MatrixXd D(num_full_x_,n);
        D=MatrixXd::Zero(num_full_x_,n);
        for(int s=0;s<num_sys;s++){
            for(i=0;i<n;i++){
                if(!AmbLinearDependCheck(sat1[i],sat2[i],flag,&max_flag)) continue;
                lam1=previous_sat_info_[sat1[i]-1].lam[0];
                lam2=previous_sat_info_[sat1[i]-1].lam[1];
                alpha=SQR(lam2)/(SQR(lam2)-SQR(lam1));
                beta=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
                j=para_.IndexAmb(0,sat1[i]);
                k=para_.IndexAmb(0,sat2[i]);
                if(!MatchNlFcb(sat1[i],sat2[i],&nl_fcb1,&nl_fcb2)) continue;

                float_nl =(x[j]-x[k]+beta*lam2*fix_wls[i])/lam_nl-(nl_fcb1-nl_fcb2);
                round_nl=(int)floor(float_nl+0.5);
                res_nl=float_nl-round_nl;
                fix_flag=fabs(res_nl)<=0.15;
                sprintf(buff,"%s TO FIX NL AMB sd_nl_amb(%s+%s)=sd_if_amb+sd_wl_amb=(%6.3f-%6.3f)+(%4.3f-%4.3f)=%6.3f round(sd_nl_amb)=%5d %s",
                        epoch_sat_info_collect_[0].t_tag.GetTimeStr(1).c_str(),previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),
                        (x[j]-x[k])/lam_nl,beta*lam2*fix_wls[i]/lam_nl,nl_fcb1,nl_fcb2,float_nl,round_nl,fix_flag?"YES":"NO");
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;

                if(!fix_flag) continue;

                D.data()[j+m*num_full_x_]=1.0/lam_nl;
                D.data()[k+m*num_full_x_]=-1.0/lam_nl;

                float_nl_ambs.push_back(float_nl);
                sat1[m]=sat1[i];
                sat2[m]=sat2[i];
                res_wls[m]=res_wls[i];
                fix_wls[m++]=fix_wls[i];
                res_nls.push_back(res_nl);
                sd_nl_fcbs.push_back((nl_fcb1-nl_fcb2));
            }
        }
        if(m<2) return false;

        MatrixXd E(m,num_full_x_),Qnl(m,m);
        MatMul("TN",m,num_full_x_,num_full_x_,1.0,D.data(),P.data(),0.0,E.data());
        MatMul("NN",m,m,num_full_x_,1.0,E.data(),D.data(),0.0,Qnl.data()); // 浮点窄项模糊度方差

        string str;
        for(i=0;i<float_nl_ambs.size();i++){
            sprintf(buff,"%5.2f  ",float_nl_ambs[i]);
            str+=buff;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(0): "<<str;
        str.clear();
        buff[0]='\0';
        for(i=0;i<float_nl_ambs.size();i++){
            sprintf(buff,"%5.2f  ",Qnl(i,i));
            str+=buff;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"    Q_AMB(0): "<<str;
        str.clear();
        buff[0]='\0';

        VectorXd s(2),B;
        MatrixXd fix_nl_ambs(m,2);
        B=Map<VectorXd>(float_nl_ambs.data(),m);
        if(!(info=lambda_.IntegerAmb(B,Qnl,fix_nl_ambs,m,2,s))){
            igcoupled_sol_.ratio=MIN(s[1]/s[0],999.9);

            if(s[0]>0.0&&ppp_conf_.gnssC.ar_thres[0]<igcoupled_sol_.ratio){
                for(i=0;i<m;i++){
                    sprintf(buff,"%5.2f  ",fix_nl_ambs.col(0)[i]);
                    str+=buff;
                }
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(1): "<<str;
                str.clear();
                buff[0]='\0';

                for(i=0;i<m;i++){
                    sprintf(buff,"%5.2f  ",fix_nl_ambs.col(1)[i]);
                    str+=buff;
                }
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(2): "<<str;
                str.clear();
                buff[0]='\0';
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP-AR VALIDATION OK, ratio="<<igcoupled_sol_.ratio;
            }
            else{
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<< "PPP-AR VALIDATION FAILED, ratio="<<igcoupled_sol_.ratio;
                return false;
            }
        }
        else{
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP LAMBDA ERROR";
            return false;
        }

        // nl to iono-free ambiguity
        VectorXd fix_if_ambs(m);
        for(i=0;i<m;i++){
            buff[0]='\0';
            fix_if_ambs[i]=0.0;
            fix_if_ambs[i]=lam_nl*(fix_nl_ambs(i,0)+sd_nl_fcbs[i])-lam2*beta*(fix_wls[i]);
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);
            sprintf(buff,"REFACTOR SD_IF_AMB(%s-%s): float_if_amb=%6.3f-%6.3f=%6.3f fixed_if_amb=%6.3f",
                    previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),x[j],x[k],x[j]-x[k],fix_if_ambs[i]);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<buff;
        }

        return FixPppSol(sat1,sat2,fix_if_ambs.data(),m,x,P);
    }

    /*
    * Function : Fix PPP ambiguity 
    */
    bool cPppSolver::ResolvePppAmbFixRef(int nf,VectorXd& x,MatrixXd& P){
        double elmask,var=0.0,res_wl;
        int i,j,m=0,stat,nb=0;
        if(num_valid_sat_<=0||ppp_conf_.gnssC.ion_opt!=ION_IF) return false;
        elmask=ppp_conf_.gnssC.ar_el_mask>0.0?ppp_conf_.gnssC.ar_el_mask:ppp_conf_.gnssC.ele_min;
        vector<int>sat_no1,sat_no2,fix_wls;
        vector<double>res_wls;

        AverageLcAmb();

        igcoupled_sol_.ratio=0.0;
        if(ppp_conf_.gnssC.ar_mode==AR_OFF||ppp_conf_.gnssC.ar_thres[0]<1.0){
            igcoupled_sol_.num_ar_sat=0;
            return 0;
        }

        if((nb=SelectFixSat(nullptr,1,1))<(ppp_conf_.gnssC.min_sat_num2fix-1)){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"NO ENOUGH VALID SATELLITE MAKING SINGLE DIFFERENCE0";
            return -1;
        }

        igcoupled_sol_.num_ar_sat=nb;

#if 1
        for(i=0;i<3;i++) var+=SQRT(full_Px_(i,i));
        var=var/3.0;
        if(var>0.05){ //0.25
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"POSITIONING VARIANCE TOO LARGE FOR PPP AR var="<<var;
            igcoupled_sol_.num_ar_sat=0;
            return false;
        }
#endif

        int ref_sat_no,sat_no,fix_wl_amb;
        for(i=0;i<sdambs_.size();i++){
            if(i==sdambs_[i].ref_sat_idx) continue;
            ref_sat_no=epoch_sat_info_collect_[sdambs_[i].ref_sat_idx].sat.sat_.no;
            sat_no=epoch_sat_info_collect_[sdambs_[i].sat_idx].sat.sat_.no;
            if(FixWlAmb(ref_sat_no,sat_no,&fix_wl_amb,&res_wl)){
                m++;
                sat_no1.push_back(ref_sat_no);
                sat_no2.push_back(sat_no);
                fix_wls.push_back(fix_wl_amb);
                res_wls.push_back(res_wl);
            }
        }

        if(ppp_conf_.gnssC.ar_mode==AR_PPP_AR&&m>3){
            if(ppp_conf_.gnssC.ar_prod==AR_PROD_FCB_SGG){
                stat=FixNlAmbILS_FCB(sat_no1.data(),sat_no2.data(),fix_wls.data(),res_wls.data(),m,x,P);
            }
            else if(ppp_conf_.gnssC.ar_prod==AR_PROD_IRC_CNES){
                stat=FixNlAmbILS_IRC(sat_no1.data(),sat_no2.data(),fix_wls.data(),res_wls.data(),m,x,P);
            }else if(ppp_conf_.gnssC.ar_prod==AR_PROD_OSB_WUH || ppp_conf_.gnssC.ar_prod==AR_PROD_OSB_CNES){//如果使用OSB zwl 
                stat= FixNlAmbILS_IRC(sat_no1.data(),sat_no2.data(),fix_wls.data(),res_wls.data(),m,x,P);
            }
            else
                return false;
        }
        else return false;

        return stat;
    }

    // -------------------------------------

    void cPppSolver::ResetAmb(double *bias, double *xa, int nb) {
        int i,n,m,f,index[MAX_SAT_NUM]={0},nv=0,nf=para_.GetGnssUsedFrqs();
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<num_real_x_fix_;i++) xa[i]=real_x_fix_[i];

#if 0
        for(m=0;m<NSYS;m++){
            for(f=0;f<nf;f++){
                for(n=i=0;i<MAX_SAT_NUM;i++){
                    sat_info=&previous_sat_info_[i];
                    if(sat_info->sat.sat_.sys_idx!=m||sat_info->fix[f]!=2) continue;
                    index[n++]=para_.IndexAmb(f,sat_info->sat.sat_.no);
                }

                if(n<2) continue;
                xa[index[0]]=full_x_[index[0]];

                for(i=1;i<n;i++){
                    xa[index[i]]=xa[index[0]]-bias[nv++];
                }
            }
        }
#endif
        int sat1,sat2,ib1,ib2;
        for(i=0;i<nb;i++){
            sat1=epoch_sat_info_collect_[sdambs_[i].ref_sat_idx].sat.sat_.no;
            sat2=epoch_sat_info_collect_[sdambs_[i].sat_idx].sat.sat_.no;
            f=0;

            ib1=para_.IndexAmb(f,sat1);
            ib2=para_.IndexAmb(f,sat2);
            // 基准卫星仍是浮点解
            xa[ib2]=xa[ib1]-bias[i];

            previous_sat_info_[sat1-1].fix_amb[f]=xa[ib1];
            previous_sat_info_[sat2-1].fix_amb[f]=-bias[i];
        }
    }

    bool cPppSolver::FixNlAmbILS1(int *sat1, int *sat2, int *fix_wls, int n,double *xa) {
        double lam_nl=gnss_obs_operator_.LinearCombLam(1,1,0,previous_sat_info_[sat1[0]-1]);
        double lam1=0,lam2=0,alpha=0,beta=0;
        int i,j,k,m=0,flag[MAX_SAT_NUM]={0},max_flag=0,info,fix_nl_amb;
        char out_buff[1024]={'\0'};
        int na=para_.GetRealFixParNum(ppp_conf_);

        VectorXd float_nl(n),s(2);
        float_nl=VectorXd::Zero(n);

        MatrixXd D(num_full_x_,num_full_x_);
        D=MatrixXd::Zero(num_full_x_,num_full_x_);

        for(i=0;i<na;i++){
            D(i,i)=1.0;
        }

        for(i=0;i<n;i++){
            out_buff[0]='\0';
            lam1=previous_sat_info_[sat1[i]-1].lam[0];
            lam2=previous_sat_info_[sat1[i]-1].lam[1];
            alpha=SQR(lam2)/(SQR(lam2)-SQR(lam1));
            beta=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);

            float_nl[m]=(full_x_[j]-full_x_[k]+beta*lam2*fix_wls[i])/lam_nl;
            fix_nl_amb=(int)floor(float_nl[m]+0.5);

            if(fabs(fix_nl_amb-float_nl[m])>0.25) continue;

            sprintf(out_buff,"%s TO FIX NL AMB sd_nl_amb(%s-%s)=sd_if_amb-sd_wl_amb=%6.3f-%6.3f=%6.3f round(sd_nl_amb)=%5d",
                    epoch_sat_info_collect_[0].t_tag.GetTimeStr(1).c_str(),previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),
                    (full_x_[j]-full_x_[k])/lam_nl,beta*lam2*fix_wls[i]/lam_nl,float_nl[m],fix_nl_amb);
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<out_buff;

            D.data()[j+(na+m)*num_full_x_]=1.0;
            D.data()[k+(na+m)*num_full_x_]=-1.0;
            sat1[m]=sat1[i];
            sat2[m]=sat2[i];
            fix_wls[m++]=fix_wls[i];
        }
        if(m<3) return false;

        MatrixXd fix_nl(n,2);
        fix_nl=MatrixXd::Zero(n,2);
        VectorXd B(m);
        B=Map<VectorXd>(float_nl.data(),m);

        char ss[20];
        string buff;
        for(i=0;i<B.size();i++){
            sprintf(ss,"%5.2f  ",B[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

#if 0
        if(!(info=lambda_.IntegerAmb(B1,Q,N1,m,2,s))){
            igcoupled_sol_.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
            if(igcoupled_sol_.ratio>999.9) igcoupled_sol_.ratio=999.9f;

            for(i=0;i<m;i++){
                sprintf(ss,"%5.2f  ",N1.col(0)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_AMB(1): "<<buff;
            buff.clear();
            ss[0]='\0';

            for(i=0;i<m;i++){
                sprintf(ss,"%5.2f  ",N1.col(1)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_AMB(2): "<<buff;
            buff.clear();
            ss[0]='\0';

            for(i=0;i<m;i++){
                out_buff[0]='\0';
                NC[i]=alpha*lam1*N1(i,0)+beta*lam2*(N1(i,0)-fix_wls[i]);
                j=para_.IndexAmb(0,sat1[i]);
                k=para_.IndexAmb(0,sat2[i]);
                sprintf(out_buff,"REFACTOR SD_IF_AMB(%s-%s): float_if_amb=%6.3f-%6.3f fixed_if_amb=%6.3f",
                        previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),full_x_[j],full_x_[k],NC[i]);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<out_buff;
            }

            VectorXd fix_float_diff(m);
            // validation by popular ratio-test of residuals
            if(s[0]<=0.0||s[1]/s[0]>=ppp_conf_.gnssC.ar_thres[0]){
                // init non phase-bias states and covariance with float solution values
                for(i=0;i<na;i++){
                    real_x_fix_[i]=full_x_[i];
                    for(j=0;j<na;j++) real_Px_fix_.data()[i+j*na]=full_Px_.data()[i+j*num_full_x_];
                }

                for(i=0;i<m;i++){
                    j=para_.IndexAmb(0,sat1[i]);
                    k=para_.IndexAmb(0,sat2[i]);
                    fix_float_diff[i]=(x[j]-x[k])-NC[i];
                }

                MatrixXd Qb_=Q.inverse();
                db=Qb_*fix_float_diff; //db=Qb^-1*(b0-b)
                real_x_fix_=real_x_fix_-Qab*db;   //xa=x-Qab-db
                QQ=Qab*Qb_;            // QQ=Qab*Qb^-1
                real_Px_fix_=real_Px_fix_-QQ*Qab.transpose();  //Pa=P-QQ*Qab^T

//                ReSetAmb(b.col(0).data(),xa,nb);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<" AR VALIDATION OK, ratio="<<s[1]/s[0];
            }
            else{
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<" RATIO VALIDATION FAILED ratio="<<s[1]/s[0];
                nb=0;
            }
        }
        else{
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"PPK LAMBDA ERROR";
            nb=0;
        }

#else
        int ny=na+m;
        VectorXd y(ny);
        MatrixXd Qyy(ny,ny),DP(ny,num_full_x_);
        y=D.transpose()*full_x_;
        MatMul("TN",ny,num_full_x_,num_full_x_,1.0,D.data(),full_Px_.data(),0.0,DP.data());
        MatMul("NN",ny,ny,num_full_x_,1.0,DP.data(),D.data(),0.0,Qyy.data());

        cout<<y.transpose()<<endl;

        MatrixXd Qbb(na,na),Qaa(m,m),Qba(na,m),Qnl(m,m);
        for(i=0;i<m;i++){
            for(j=0;j<m;j++){
                Qaa.data()[i+j*m]=Qyy.data()[(na+i)+(na+j)*ny];  //浮点解
                Qnl.data()[i+j*m]=Qaa.data()[i+j*m]/SQR(lam_nl); //窄项模糊度
            }
        }
        cout<<Qaa<<endl;
        for(i=0;i<na;i++) for(j=0;j<m;j++) Qba.data()[i+j*na]=Qyy.data()[i+(na+j)*ny]; //浮点解协方差

        for(i=0;i<B.size();i++){
            sprintf(ss,"%5.2f  ",Qnl(i,i));
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   Q_NL_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        if(!(info=lambda_.IntegerAmb(B,Qnl,fix_nl,m,2,s))){
            igcoupled_sol_.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
            if(igcoupled_sol_.ratio>999.9) igcoupled_sol_.ratio=999.9f;

            for(i=0;i<m;i++){
                sprintf(ss,"%5.2f  ",fix_nl.col(0)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(1): "<<buff;
            buff.clear();
            ss[0]='\0';

            for(i=0;i<m;i++){
                sprintf(ss,"%5.2f  ",fix_nl.col(1)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"SD_NL_AMB(2): "<<buff;
            buff.clear();
            ss[0]='\0';

            if(ppp_conf_.gnssC.ar_thres[0]>0.0&&3.0>igcoupled_sol_.ratio){
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP-AR VALIDATION FAILED, ratio="<<igcoupled_sol_.ratio;
                return false;
            }
            else{
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PP-AR VALIDATION OK, ratio="<<igcoupled_sol_.ratio;

                for(i=0;i<na;i++){
                    real_x_fix_[i]=full_x_[i];
                    for(j=0;j<na;j++) real_Px_fix_.data()[i+j*na]=full_Px_.data()[i+j*num_full_x_];
                }

                // nl to iono-free ambiguity
                VectorXd fix_if_amb(m),float_if_amb(m),fix_float_diff(m);
                fix_if_amb=VectorXd::Zero(m);
                float_if_amb=VectorXd::Zero(m);
                fix_float_diff=VectorXd::Zero(m);
                for(i=0;i<m;i++){
                    out_buff[0]='\0';
                    int nl_amb=fix_nl(i,0);
                    fix_if_amb[i]=alpha*lam1*fix_nl(i,0)+beta*lam2*(fix_nl(i,0)-fix_wls[i]);
                    j=para_.IndexAmb(0,sat1[i]);
                    k=para_.IndexAmb(0,sat2[i]);
                    float_if_amb[i]=y.data()[na+i];
                    fix_float_diff[i]=float_if_amb[i]-fix_if_amb[i];
                    sprintf(out_buff,"REFACTOR SD_IF_AMB(%s-%s): float_if_amb=%6.3f-%6.3f=%6.3f fixed_if_amb=%6.3f",
                            previous_sat_info_[sat1[i]-1].sat.sat_.id.c_str(),previous_sat_info_[sat2[i]-1].sat.sat_.id.c_str(),full_x_[j],full_x_[k],float_if_amb[i],fix_if_amb[i]);
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<out_buff;
                }

                MatrixXd Qaa_=Qaa.inverse();
                VectorXd db(m,1);
                cout<<Qaa_<<endl<<endl;
                cout<<fix_float_diff<<endl;
                db=Qaa_*fix_float_diff;
                cout<<db<<endl<<endl;
                cout<<Qba<<endl<<endl;
                cout<<Qba*db<<endl;
                real_x_fix_=real_x_fix_-Qba*db;
                MatrixXd QQ(na,m);
                QQ=Qba*Qaa_;
                real_Px_fix_=real_Px_fix_-QQ*Qba.transpose();
                ResetAmb(fix_if_amb.data(),xa,m);
            }
        }
        else{
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"PPP LAMBDA ERRPR";
            m=0;
        }
        return m;
//        return FixPppSol(sat1,sat2,fix_if_amb.data(),m,x,P);
#endif
    }

    int cPppSolver::SelectFixSat(double *D, int gps, int glo) {
        int i,j,k,m,f=0,nb=0,na=num_real_x_fix_,nf=para_.GetGnssUsedFrqs(),no_fix;
        double fix[MAX_SAT_NUM],ref[MAX_SAT_NUM],s[2];
        tSatInfoUnit *ref_sat_info= nullptr;
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<MAX_SAT_NUM;i++) for(j=0;j<MAX_GNSS_USED_FRQ_NUM;j++){
                previous_sat_info_[i].fix[j]=0;
        }

        int sat1,sat2;
        if(sdambs_.size()) sdambs_.clear();
        tSdAmb sd_amb={0};

        int sat_no;
        // 最高高度角
        int ref_sat_idx;
        for(i=0,ref_sat_idx=-1;i<epoch_sat_info_collect_.size();i++){
            if(epoch_sat_info_collect_[i].stat!=SAT_USED) continue;
            if(!previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1].vsat[0]) continue;
            if(ref_sat_idx<0||epoch_sat_info_collect_[ref_sat_idx].el_az[0]<epoch_sat_info_collect_[i].el_az[0]) ref_sat_idx=i;
        }
        if(ref_sat_idx<0) return 0;

        ref_sat_info=&previous_sat_info_[epoch_sat_info_collect_[ref_sat_idx].sat.sat_.no-1];
        for(i=0;i<epoch_sat_info_collect_.size();i++){
            if(i==ref_sat_idx) continue;
            sat_info=&previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1];
            k=para_.IndexAmb(f,epoch_sat_info_collect_[ref_sat_idx].sat.sat_.no);
            j=para_.IndexAmb(f,epoch_sat_info_collect_[i].sat.sat_.no);

            ref_sat_info->fix[f]=2;
            if(full_x_[j]==0.0||full_x_[k]==0.0) continue;
            if(sat_info->slip[f]) continue;
            no_fix=(m==0&&gps==0)||(m==1&&ppp_conf_.gnssC.bds_ar_mode==0)||(m==3&&glo==0);
            if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&sat_info->vsat[f]&&sat_info->el_az[0]*R2D>=ppp_conf_.gnssC.ar_el_mask&&!no_fix){
                sd_amb.ref_sat_idx=ref_sat_idx;
                sd_amb.sat_idx=i;
                sd_amb.f=0;
                sdambs_.push_back(sd_amb);
                ref[nb]=ref_sat_info->sat.sat_.no;
                fix[nb++]=sat_info->sat.sat_.no;
                sat_info->fix[f]=2;
            }
            else sat_info->fix[f]=1;
        }

        if(nb>0){
            VectorXd ref_sats,fix_sats;
            ref_sats=Map<VectorXd>(ref,nb,1);
            fix_sats=Map<VectorXd>(fix,nb,1);

            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"REF SATELLITES: "<<ref_sats.transpose();
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"FIX SATELLITES: "<<fix_sats.transpose();
        }
        return nb;
    }

    /*
    * Function : Fix ambiguity using LAMBDA
    */
    int cPppSolver::ResolveAmbLambda(double *xa,int gps,int glo) {
        double elmask,var=0.0,res_wl;
        int i,j,m=0,stat,nb=0;
        if(num_valid_sat_<=0||ppp_conf_.gnssC.ion_opt!=ION_IF) return false;
        elmask=ppp_conf_.gnssC.ar_el_mask>0.0?ppp_conf_.gnssC.ar_el_mask:ppp_conf_.gnssC.ele_min;
        vector<int>sat_no1,sat_no2,fix_wls;
        vector<double>res_wls;

        igcoupled_sol_.ratio=0.0;
        if(ppp_conf_.gnssC.ar_mode==AR_OFF||ppp_conf_.gnssC.ar_thres[0]<1.0){
            igcoupled_sol_.num_ar_sat=0;
            return 0;
        }

#if 1
        for(i=0;i<3;i++) var+=SQRT(full_Px_(i,i));
        var=var/3.0;
        if(var>0.05){ //0.25
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"POSITIONING VARIANCE TOO LARGE FOR PPP AR var="<<var;
            igcoupled_sol_.num_ar_sat=0;
            return false;
        }
#endif

        if((nb=SelectFixSat(nullptr,gps,glo))<(ppp_conf_.gnssC.min_sat_num2fix-1)){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"NO ENOUGH VALID SATELLITE MAKING SINGLE DIFFERENCE0";
            return -1;
        }

        igcoupled_sol_.num_ar_sat=nb;

        AverageLcAmb();

        int ref_sat_no,sat_no,fix_wl_amb;
        for(i=0;i<sdambs_.size();i++){
            if(i==sdambs_[i].ref_sat_idx) continue;
            ref_sat_no=epoch_sat_info_collect_[sdambs_[i].ref_sat_idx].sat.sat_.no;
            sat_no=epoch_sat_info_collect_[sdambs_[i].sat_idx].sat.sat_.no;
            if(FixWlAmb(ref_sat_no,sat_no,&fix_wl_amb,&res_wl)){
                m++;
                sat_no1.push_back(ref_sat_no);
                sat_no2.push_back(sat_no);
                fix_wls.push_back(fix_wl_amb);
                res_wls.push_back(res_wl);
            }
        }

        if(ppp_conf_.gnssC.ar_mode==AR_PPP_AR&&m>0){
            stat=FixNlAmbILS1(sat_no1.data(),sat_no2.data(),fix_wls.data(),m,xa);
        }
        else return false;

        return stat;
    }

    bool cPppSolver::ResolverPppAmb1(int nf, double *xa) {
        int nb=0,f,i,ar=0,lockc[MAX_GNSS_USED_FRQ_NUM],arsats[64]={0};
        bool exc_flag=false;
        tSatInfoUnit *sat_info= nullptr;

        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PREVIOUS EPOCH AR: ratio1="<<pre_epoch_ar_ratio1_<<" ratio2="<<pre_epoch_ar_ratio2_<<" ar_sat_num="<<igcoupled_sol_.num_ar_sat;

        // 若前一历元没有ＡＲ成功
        nf=1;
        if(pre_epoch_ar_ratio2_<ppp_conf_.gnssC.ar_thres[0]&&igcoupled_sol_.num_ar_sat>=ppp_conf_.gnssC.min_sat_num2drop){
            // previous epoch ar failed and satellite enough to drop
            for(f=0;f<nf;f++){
                for(i=0;i<epoch_sat_info_collect_.size();i++){
                    sat_info=&previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1];
                    if(sat_info->vsat[f]&&sat_info->lock[f]>=0&&sat_info->el_az[0]*R2D>=ppp_conf_.gnssC.ele_min){
                        arsats[ar++]=i;
                    }
                }
            }

            if(exc_sat_index_<ar){
                i=epoch_sat_info_collect_[arsats[exc_sat_index_]].sat.sat_.no;
                for(f=0;f<nf;f++){
                    lockc[f]=previous_sat_info_[i-1].lock[f];
                    previous_sat_info_[i-1].lock[f]=-igcoupled_sol_.num_ar_sat;
                }
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[i-1].sat.sat_.id<<" EXCLUDE BY AR";
                exc_flag=true;
            }
            else exc_sat_index_=0;
        }

        double ratio1;
        int dly;
        bool rerun=false;
        if(ppp_conf_.gnssC.glo_ar_mode!=GLO_AR_FIXHOLD){

            nb=ResolveAmbLambda(xa,1,0);

            ratio1=igcoupled_sol_.ratio;
            if(ppp_conf_.gnssC.ar_filter){
                if(nb>=0&&pre_epoch_ar_ratio2_>=ppp_conf_.gnssC.ar_thres[0]&&(igcoupled_sol_.ratio<ppp_conf_.gnssC.ar_thres[0])||
                   (igcoupled_sol_.ratio<ppp_conf_.gnssC.ar_thres[0]*1.1&&igcoupled_sol_.ratio<pre_epoch_ar_ratio1_/2.0)){
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1);
                    dly=2;
                    for(i=0;i<epoch_sat_info_collect_.size();i++){
                        for(f=0;f<nf;f++){
                            if(previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1].fix[f]!=2) continue;
                            if(previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1].lock[f]==0){
                                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1].sat.sat_.id<<" L"<<f+1<<" IS NEW OBSERVED";
                                previous_sat_info_[epoch_sat_info_collect_[i].sat.sat_.no-1].lock[f]=-ppp_conf_.gnssC.min_lock2fix-dly;
                                dly+=2;
                                rerun=true;
                            }
                        }
                    }
                }

                if(rerun){
                    if((nb=ResolveAmbLambda(xa,1,0))<2){
                        pre_epoch_ar_ratio1_=pre_epoch_ar_ratio2_=ratio1;
                        return false;
                    }
                }
            }
        }
        else{
            ratio1=0.0;
            nb=0;
        }

        if(exc_flag&&(igcoupled_sol_.ratio<ppp_conf_.gnssC.ar_thres[0])&&igcoupled_sol_.ratio<1.5*pre_epoch_ar_ratio2_){
            i=epoch_sat_info_collect_[arsats[exc_sat_index_++]].sat.sat_.no;
            for(f=0;f<nf;f++) previous_sat_info_[i-1].lock[f]=lockc[f];
        }

        pre_epoch_ar_ratio1_=ratio1>0?ratio1:igcoupled_sol_.ratio;
        pre_epoch_ar_ratio2_=igcoupled_sol_.ratio;

        if(nb>0) return true;
        else return false;
    }

    // -------------------------------------
    /*
    * Function : The constructor and detructor function of class cPpkSolver
    */
    cPpkSolver::cPpkSolver() {}

    cPpkSolver::cPpkSolver(tIGCOUPLEDConf C) {
        ppk_conf_=C;
//        ppk_conf_.mode=MODE_PPK;
        para_=cParSetting(ppk_conf_);
        num_full_x_=para_.GetIGCOUPLEDPar(ppk_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);

        num_real_x_fix_=para_.GetRealFixParNum(ppk_conf_);
        real_x_fix_=VectorXd::Zero(num_real_x_fix_);
        real_Px_fix_=MatrixXd::Zero(num_real_x_fix_,num_real_x_fix_);
    }

    cPpkSolver::~cPpkSolver() {}

    void cPpkSolver::InitSolver(tIGCOUPLEDConf C) {
        ppk_conf_=C;

        cReadGnssObs base_reader(C.fileC.base,nav_,base_obs_,REC_BASE);
        base_reader.SetGnssSysMask(C.gnssC.nav_sys);
        base_reader.Reading();

        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_=new cOutSol(ppk_conf_);
        out_->InitOutSol(ppk_conf_,ppk_conf_.fileC.sol);
        out_->WriteHead();

        tIGCOUPLEDConf spp_conf=C;
        spp_conf.mode=MODE_SPP;
        spp_conf.mode_opt=MODE_OPT_KINEMATIC;
        spp_conf.gnssC.ion_opt=ION_KLB;
        spp_conf.gnssC.trp_opt=TRP_SAAS;
        spp_conf.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf.gnssC.res_qc=0;
        spp_solver_=new cSppSolver(spp_conf);
        spp_solver_->spp_conf_=spp_conf;
        spp_solver_->nav_=nav_;
        spp_solver_->InitSolver(spp_conf);
    }

    bool cPpkSolver::SolverProcess(tIGCOUPLEDConf C,int idx) {
        bool stat=false;
        if(idx==-1) InitSolver(C);

        int i=idx,num_epochs=rover_obs_.epoch_num;
        if(idx!=-1){
            num_epochs=idx+1;
        }

		double outage_time = 0;
		if(ppk_conf_.gnssC.use_outage) outage_time = ppk_conf_.gnssC.outage_time;

        if(C.filter_type==FILTER_FORWARD){
            for(i=idx==-1?0:idx;i<num_epochs;i++){
				if(outage_time != 0){ //gnss outage
                     double sow = rover_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS);
                     if(sow > (outage_time + ppk_conf_.gnssC.outage_len) && 
                         ppk_conf_.gnssC.outage_period > ppk_conf_.gnssC.outage_len){ // reinit outage time
                             outage_time += static_cast<int>((sow - outage_time) / ppk_conf_.gnssC.outage_period) * 
                                             ppk_conf_.gnssC.outage_period;
                         }
                         if(sow >= outage_time && sow <= outage_time + ppk_conf_.gnssC.outage_len){// outage
                             CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS OUTAGE AT" << sow << "(s) " << endl;
                             continue;
                         }
                 }

                stat=SolverStart(i,idx);
            }
            if(idx!=-1) return stat;
            else CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (FORWARD): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
        }
        else if(C.filter_type==FILTER_BACKWARD){
            for(i=idx==-1?num_epochs-1:idx;i>=0;i--){
                stat=SolverStart(i,idx);
            }
            if(idx!=-1) return stat;
            else CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (BACKWARD): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
        }
        else if(C.filter_type==FILTER_COMBINED){
            for(i=idx==-1?0:idx;i<num_epochs;i++){
                SolverStart(i,idx);
                solf_.push_back(igcoupled_sol_);
            }

            epoch_ok_=0;
            epoch_fail_=0;
            epoch_idx_=0;
            ReinitSolver(C);
            for(i=idx==-1?num_epochs-1:idx;i>=0;i--){
                SolverStart(i,idx);
                solb_.push_back(igcoupled_sol_);
            }

            CombFbSol(C);
            CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<" TOTAL EPOCH (COMBINED): "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
            solf_.clear();
            solb_.clear();
        }


        CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"TOTAL EPOCH: "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
    }

    bool cPpkSolver:: SolverStart(int i, int idx) {
        epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"START PPK SOLVING "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER "<<epoch_sat_obs_.sat_num;

        if(MatchBaseObs(epoch_sat_obs_.obs_time)){

            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"MATCH BASE STATION OBSERVATIONS, BASE SATELLITE NUMBER "<<base_epoch_sat_obs_.sat_num;

            if(SolverEpoch()){
                SolutionUpdate();
                if(idx!=-1){
                    epoch_sat_info_collect_.clear();
                    base_sat_info_collect_.clear();
                    rover_res.clear();base_res.clear();
                    return true;
                }
            }

            if(idx!=-1){
                epoch_sat_info_collect_.clear();
                base_sat_info_collect_.clear();
                rover_res.clear();base_res.clear();
                return false;
            }
        }
        else{
            igcoupled_sol_.stat=SOL_NONE;
            igcoupled_sol_.t_tag=epoch_sat_obs_.obs_time;
            epoch_fail_++;
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"MATCH BASE STATION OBSERVATIONS FAILED";
        }

        epoch_sat_info_collect_.clear();
        base_sat_info_collect_.clear();
        rover_res.clear();base_res.clear();
    }

    bool cPpkSolver::SolverEpoch() {
        epoch_idx_+=1;
        base_xyz_=ppk_conf_.gnssC.rb;

        UpdateGnssObs(ppk_conf_,epoch_sat_obs_,REC_ROVER);
        UpdateGnssObs(ppk_conf_,base_epoch_sat_obs_,REC_BASE);
        InitEpochSatInfo(epoch_sat_info_collect_);

        InitSppSolver();
        if(spp_solver_->SolverEpoch()){
            Spp2Ppk();

            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PPK-SPP SOLVE SUCCESS, SPP POS "<<spp_solver_->igcoupled_sol_.pos.transpose()<<" DOPPLER  VEL "<<spp_solver_->igcoupled_sol_.vel.transpose();

            gnss_err_corr_.eph_model_.EphCorr(base_sat_info_collect_);
            if(ppk_conf_.gnssC.eph_opt==EPH_PRE){
                gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);
            }

            if(Estimator(ppk_conf_)){
                epoch_ok_++;
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK SOLVE SUCCESS";
                return true;
            }else{
                epoch_fail_++;
                igcoupled_sol_=spp_solver_->igcoupled_sol_;
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK SOLVE FAILED";
                return false;
            }
        }
        else{
            epoch_fail_++;
            igcoupled_sol_.stat=SOL_NONE;
            igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            spp_solver_->epoch_sat_info_collect_.clear();
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK-SPP SOLVE FAILED";
            return false;
        }

    }

    bool cPpkSolver::Estimator(tIGCOUPLEDConf C) {

        // select common satellites
        vector<int>ir,ib,cmn_sat_no;

        SelectCmnSat(C,ir,ib,cmn_sat_no);
        // cycle slip
        PpkCycleSlip(C,ir,ib,cmn_sat_no);
        // state time_update
        StateTimeUpdate(C,ir,ib,cmn_sat_no);

        // zero residual for base and rover
        if(!GnssZeroRes(C,REC_BASE,ib,full_x_.data(),C.gnssC.rb)){
            return false;
        }

        Vector3d rover_xyz,ve;

        static int ref_sat[NSYS][2*MAX_GNSS_USED_FRQ_NUM]={0};
        VectorXd x;
        tImuInfoUnit cor_imu_info;
        qc_iter_=0;
        for(int iter=0;iter<5;iter++){
            cor_imu_info=cur_imu_info_;
            x=full_x_;
            MatrixXd Px=full_Px_;

            if(tc_mode_){
                RemoveLever(cur_imu_info_,C.insC.lever,rover_xyz,ve);
            }
            else{
                rover_xyz<<full_x_[0],full_x_[1],full_x_[2];
            }

            if(!GnssZeroRes(C,REC_ROVER,ir,x.data(),rover_xyz)){
                break;
            }

            if(!GnssDdRes(0,C,ir,ib,cmn_sat_no,x.data(),ref_sat)){
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"MAKE DOUBLE DIFFERENCE RESIDUAL FAILED";
                return false;
            }

            if(!kf_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_)){
                igcoupled_sol_.stat=SOL_SPP;
                return false;
            }

            if(tc_mode_){
                CloseLoopState(x,&cor_imu_info);
                RemoveLever(cor_imu_info,C.insC.lever,rover_xyz,ve);
            }
            else{
                rover_xyz<<x[0],x[1],x[2];
            }

            if(GnssZeroRes(C,REC_ROVER,ir,x.data(),rover_xyz)){
                if(GnssDdRes(1,C,ir,ib,cmn_sat_no,x.data(),ref_sat)){

                    if(tc_mode_){
                        igcoupled_sol_.stat=SOL_FLOAT;
                        igcoupled_sol_.ins_stat=SOL_IG_TC;
                        cur_imu_info_=cor_imu_info;
                    }
                    else{
                        igcoupled_sol_.stat=SOL_FLOAT;
                    }
                    full_x_=x;
                    full_Px_=Px;
                    UpdateSatInfo(epoch_sat_info_collect_);
                    break;
                }

            }
            else{
                igcoupled_sol_.stat=SOL_SPP;
                return false;
            }
        }

        //PPK-AR
#if 1
        VectorXd xa=x;
        if(igcoupled_sol_.stat==SOL_FLOAT){
            if(ResolvePpkAmb(cmn_sat_no,para_.GetGnssUsedFrqs(),xa.data())){

                if(GnssZeroRes(C,REC_ROVER,ir,xa.data(),rover_xyz)){
                    GnssDdRes(5,C,ir,ib,cmn_sat_no,xa.data(),ref_sat);
                    igcoupled_sol_.stat=SOL_FIX;
                }

                if(++num_continuous_fix_>=C.gnssC.min_fix2hold){
                    if(C.gnssC.ar_mode==AR_FIX_HOLD||C.gnssC.glo_ar_mode==GLO_AR_FIXHOLD){
                        HoldAmb(cmn_sat_no,xa.data());
                    }
                }
            }
            else{
                num_continuous_fix_=0;
            }
        }

#endif
        int sat_no,ia;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            for(int f=0;f<MAX_GNSS_USED_FRQ_NUM;f++){
                if(!epoch_sat_info_collect_[i].vsat[f]) continue;
                sat_no=epoch_sat_info_collect_[i].sat.sat_.no;
                ia=para_.IndexAmb(f,sat_no);
                previous_sat_info_[sat_no-1].outc[f]=0;
                if(previous_sat_info_[sat_no-1].lock[f]<0||previous_sat_info_[sat_no-1].fix[f]>=2){
                    previous_sat_info_[sat_no-1].lock[f]++;
                }
            }
        }

        ir.clear();ib.clear();cmn_sat_no.clear();ddambs_.clear();vflag_.clear();
        return igcoupled_sol_.stat;
    }

    bool cPpkSolver::SolutionUpdate() {
        if(igcoupled_sol_.stat==SOL_SPP) return false;
        else if(igcoupled_sol_.stat==SOL_FLOAT){
            igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) igcoupled_sol_.pos[i]=full_x_[i];
            for(int i=0;i<num_pos;i++) igcoupled_sol_.q_pos[i]=full_Px_(i,i);
        }
        else if(igcoupled_sol_.stat==SOL_FIX){
            igcoupled_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) igcoupled_sol_.pos[i]=real_x_fix_[i];
            for(int i=0;i<num_pos;i++) igcoupled_sol_.q_pos[i]=full_Px_(i,i);
        }
        for(int i=0;i<3;i++) spp_solver_->full_x_[i]=full_x_[i];
        igcoupled_sol_.valid_sat_num=num_valid_sat_;

        if(ppk_conf_.solC.out_sol) out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
        if(ppk_conf_.solC.out_stat) out_->WriteSatStat(&igcoupled_sol_,previous_sat_info_);
    }

    void cPpkSolver::InitSppSolver() {
        spp_solver_->epoch_idx_=epoch_idx_;
        spp_solver_->epoch_sat_info_collect_=epoch_sat_info_collect_;
        spp_solver_->cur_imu_info_=cur_imu_info_;
        spp_solver_->tc_mode_=tc_mode_;
    }

    void cPpkSolver::Spp2Ppk() {
        spp_solver_->SolutionUpdate();
        epoch_sat_info_collect_.clear();
        epoch_sat_info_collect_=spp_solver_->epoch_sat_info_collect_;
        spp_solver_->epoch_sat_info_collect_.clear();
    }

    bool cPpkSolver::GnssZeroRes(tIGCOUPLEDConf C, RECEIVER_INDEX rec,vector<int>sat_idx,double* x,Vector3d rr) {

        vector<tSatInfoUnit> *sat_collect;
        vector<double> *res;

        sat_collect=rec==REC_ROVER?&epoch_sat_info_collect_:&base_sat_info_collect_;
        res=rec==REC_ROVER?&rover_res:&base_res;
        double omc,r,meas;
        tSatInfoUnit* sat_info= nullptr;

        Vector3d rec_xyz;
        if(rec==REC_ROVER){
            rec_xyz=rr;
        }
        else{
            rec_xyz=base_xyz_;
        }

        Vector3d rover_blh=Xyz2Blh(rec_xyz);

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return false;
        int nf=num_used_frq*num_used_obs_type;

        rover_res.resize(sat_idx.size()*num_used_frq*num_used_obs_type,0.0);
        base_res.resize(sat_idx.size()*num_used_frq*num_used_obs_type,0.0);

        char buff[1024]={'\0'};
//        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<(rec==REC_ROVER?"ROVER ":"BASE ")<<"STATION ZERO RESIDUAL(omc = obs - (r - dtr^s + trp)): "<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1);
        for(int i=0;i<sat_idx.size();i++){
            buff[0]='\0';
            sat_info=&sat_collect->at(sat_idx[i]);
            if(sat_info->stat!=SAT_USED){
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" NO USED("<<kGnssSatStatStr[sat_info->stat+1]<<")";
                continue;
            }

            if((r=GeoDist(sat_info->brd_pos,rec_xyz,sat_info->sig_vec))<=0.0){
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" NO USED";
                continue;
            }

            if(SatElAz(rover_blh,sat_info->sig_vec,sat_info->el_az)<C.gnssC.ele_min*D2R){
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" LOW ELE";
                continue;
            }

            r+=gnss_err_corr_.SagnacCorr(sat_info->brd_pos,rec_xyz);
            gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.trp_model_.GetTrpError(0.0,x,para_.IndexTrp());
            gnss_err_corr_.trp_model_.UpdateSatInfo();
            gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.ion_model_.GetIonError();
            gnss_err_corr_.ion_model_.UpdateSatInfo();
//            gnss_err_corr_.ant_model_.SatPcvCorr(sat_info,rover_blh, nullptr);
//            gnss_err_corr_.ant_model_.RecAntCorr(sat_info, nullptr,REC_ROVER);

            int obs_code,frq;
            for(int f=0;f<num_used_frq*num_used_obs_type;f++){
                //f%num_used_obs_type==0, L, ==1 P, ==2 D,
                //f/num_used_obs_type==1, f1,==2 f2,==3 f3
                frq=f>=num_used_frq?f-num_used_frq:f;
                obs_code=f<num_used_frq?0:1;

                //ionosphere-free
                if(C.gnssC.ion_opt==ION_IF){
//                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                    meas=obs_code?sat_info->cor_if_P[frq]:sat_info->cor_if_L[frq];
                }
                    //uncombined and undifference
                else{
                    meas=obs_code?sat_info->raw_P[frq]:sat_info->raw_L[frq]*sat_info->lam[frq];
//                    meas=obs_type==GNSS_OBS_CODE?sat_info->raw_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->raw_L[frq]*sat_info->lam[frq]:sat_info->cor_D[frq]);
                }

                if(meas==0.0) {
//                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" MISSING "<<(obs_code?"P":"L")<<frq+1;
                    continue;
                }

                double sat_clk=C.gnssC.eph_opt==EPH_BRD?sat_info->brd_clk[0]*CLIGHT:sat_info->pre_clk[0]*CLIGHT;
                double trp_del=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];
                double ion_del=sat_info->ion_delay[0];
                omc=meas-(r-sat_clk+trp_del);
                res->at(i*num_used_frq*num_used_obs_type+f)=omc;

                sprintf(buff,"%s %s%d omc=%6.3f obs=%12.3f r=%12.3f trp=%5.3f dts=%12.3f",sat_info->sat.sat_.id.c_str(),(obs_code?"P":"L"),frq+1,omc,meas,r,trp_del,sat_clk);
//                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->sat.sat_.id<<" "<<buff;
            }
        }

        res= nullptr;
        return true;
    }

    int cPpkSolver::GnssDdRes(int post, tIGCOUPLEDConf C,vector<int>ir,vector<int>ib,vector<int>cmn_sat_no, double *x,int refsat[NSYS][2*MAX_GNSS_USED_FRQ_NUM]) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        char buff[MAX_BUFF]={'\0'};

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return false;

        int nf=num_used_frq*num_used_obs_type,ref_sat_idx,j,ia,ref_ia,it,ref_it,ii,ref_ii,need_iter=1;
        int frq,obs_code,nobs[NSYS][nf],iobs=0;
        double omc,ambi,ambj,threshadj;
        vector<double> omcs,H,Ri,Rj,R;
        tSatInfoUnit* sat_info= nullptr, *ref_sat_info= nullptr;

        for(int i=0;i<NSYS;i++){
            for(int k=0;k<nf;k++){
                nobs[i][k]=0;
            }
        }

        for(int isys=0;isys<NSYS;isys++){

            for(int ifrq=0;ifrq<nf;ifrq++){
                ambi=0.0;ambj=0.0;
                frq=ifrq>=num_used_frq?ifrq-num_used_frq:ifrq;
                obs_code=ifrq<num_used_frq?0:1;

                for(j=0,ref_sat_idx=-1;j<ir.size();j++){
                    int sys=epoch_sat_info_collect_.at(ir[j]).sat.sat_.sys_idx;
                    if(sys!=isys) continue;
                    if(!ValidObs(j,num_used_frq,ifrq)) continue;
                    if(ref_sat_idx<0||epoch_sat_info_collect_.at(ir[ref_sat_idx]).el_az[0]<epoch_sat_info_collect_.at(ir[j]).el_az[0]) ref_sat_idx=j;
                }
                if(ref_sat_idx<0) continue;
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<kGnssSysStr[isys+1]<<(post?" POST":" PRIOR")<<" DOUBLE DIFFERENCE, "<<" REFERENCE SATELLITE "<<epoch_sat_info_collect_[ir[ref_sat_idx]].sat.sat_.id<<" "<<epoch_sat_info_collect_[ir[ref_sat_idx]].el_az[0]*R2D;

                if(refsat){
                    refsat[isys][ifrq]=cmn_sat_no[ref_sat_idx];
                }

                for(int isat=0;isat<cmn_sat_no.size();isat++){
                    sat_info=&epoch_sat_info_collect_[ir[isat]];
                    ref_sat_info=&epoch_sat_info_collect_[ir[ref_sat_idx]];

                    if(isat==ref_sat_idx) {
                        if(post==1&&ifrq<num_used_frq){
                            ref_sat_info->outc[frq]=0;
                            ref_sat_info->lock[frq]++;
                            ref_sat_info->vsat[frq]=1;
                        }
                        continue;
                    }

                    if(sat_info->stat!=SAT_USED) continue;
                    int sys=epoch_sat_info_collect_.at(ir[isat]).sat.sat_.sys_idx;
                    if(sys!=isys) continue;
                    if(!ValidObs(isat,num_used_frq,ifrq)){
                        continue;
                    }

                    if(C.gnssC.check_dual_phase&&C.gnssC.frq_opt>=FRQ_DUAL&&!obs_code){
                        if(!CheckDualFrq(*sat_info)){
                            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" NO DUAL PHASE OBSERVATIONS";
                            continue;
                        }
                    }

                    omc=(rover_res[ref_sat_idx*nf+ifrq]-base_res[ref_sat_idx*nf+ifrq])-(rover_res[isat*nf+ifrq]-base_res[isat*nf+ifrq]);

                    double df=0.0,glo_bia=0.0;
                    if(sat_info->sat.sat_.sys==SYS_GLO&&ref_sat_info->sat.sat_.sys==SYS_GLO){
                        if(C.gnssC.glo_ar_mode==GLO_AR_AUTO&&frq<2){
                            df=(CLIGHT/ref_sat_info->lam[frq]-CLIGHT/sat_info->lam[frq])/(ifrq==0?FREQ_GLO_D1:FREQ_GLO_D2);
                            glo_bia=df*x[para_.IndexGloIfpb()+frq];
                            omc-=glo_bia;
                        }
                        else if(C.gnssC.glo_ar_mode==GLO_AR_FIXHOLD&&frq<2){
                            glo_bia=previous_sat_info_[ref_sat_info->sat.sat_.no-1].glo_bias[frq]*ref_sat_info->lam[frq]-previous_sat_info_[sat_info->sat.sat_.no-1].glo_bias[frq]*sat_info->lam[frq];
                            omc-=glo_bia;
                        }
                    }

                    if(C.gnssC.trp_opt>=TRP_EST_WET){
                        // ???
                    }

                    if(C.gnssC.ion_opt==ION_EST){
                        // ???
                    }

                    if(!obs_code){
                        ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                        ref_ia=para_.IndexAmb(frq,ref_sat_info->sat.sat_.no);
                        if(C.gnssC.ion_opt!=ION_IF){
                            ambi=ref_sat_info->lam[frq]*x[ref_ia];ambj=sat_info->lam[frq]*x[ia];
                        }
                        else{
                            // ????
                        }
                        omc-=(ambi-ambj);
                    }

                    threshadj=obs_code||(full_Px_(ia,ia)>SQR(30/2))||(full_Px_(ref_ia,ref_ia)>SQR(30/2))?300:1.0;
                    if(C.gnssC.max_inno>0.0&&fabs(omc)>C.gnssC.max_inno*threshadj){
                        sat_info->rejc[frq]++;
                        CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<(obs_code?"P":"L")<<frq+1<<" INNOVATION OVERRUN inno="<<omc<<" threshold="<<threshadj;
                        continue;
                    }

                    omcs.push_back(omc);
                    Ri.push_back(GnssMeasVar(C,(obs_code?GNSS_OBS_CODE:GNSS_OBS_PHASE),*ref_sat_info));
                    double var=GnssMeasVar(C,(obs_code?GNSS_OBS_CODE:GNSS_OBS_PHASE),*sat_info);
                    if(obs_code){
                        var*=previous_sat_info_[sat_info->sat.sat_.no-1].c_var_factor[frq];
                    }
                    else{
                        var*=previous_sat_info_[sat_info->sat.sat_.no-1].p_var_factor[frq];
                    }
                    Rj.push_back(var);
                    if(!post){
                        //position
                        if(tc_mode_){
                            for(int i=0;i<num_full_x_;i++) H.push_back(i<3?(ref_sat_info->sig_vec[i]-sat_info->sig_vec[i]):0.0);
                        }
                        else{
                            for(int i=0;i<num_full_x_;i++) H.push_back(i<3?(-ref_sat_info->sig_vec[i]+sat_info->sig_vec[i]):0.0);
                        }

                        // glo ifpb
                        if(sat_info->sat.sat_.sys==SYS_GLO&&ref_sat_info->sat.sat_.sys==SYS_GLO){
                            if(C.gnssC.glo_ar_mode==GLO_AR_AUTO&&frq<2){
                                int glo_bias_idx=para_.IndexGloIfpb()+frq;
                                H[glo_bias_idx+num_full_x_*num_L_]=df;
                            }
                        }

                        //tropospheric delay
                        if(C.gnssC.ion_opt>=ION_EST){
                            //???
                        }
                        //ionospheric delay
                        if(C.gnssC.trp_opt>=TRP_EST_WET){
                            //????
                        }
                        //ambiguity
                        if(!obs_code){
                            if(C.gnssC.ion_opt!=ION_IF){
                                H[ref_ia+num_full_x_*num_L_]=ref_sat_info->lam[frq];
                                H[ia+num_full_x_*num_L_]=-sat_info->lam[frq];
                            }
                            else{
                                H[ref_ia+num_full_x_*num_L_]=1.0;
                                H[ia+num_full_x_*num_L_]=-1.0;
                            }
                        }
                        if(obs_code){
                            sat_info->prior_res_P[frq]=omcs.back();
                        }
                        else if(!obs_code){
                            sat_info->prior_res_L[frq]=omcs.back();
                        }
                    }
                    else{
                        if(obs_code){
                            sat_info->post_res_P[frq]=omcs.back();
                        }
                        else if(!obs_code&&post!=5){
                            sat_info->post_res_L[frq]=omcs.back();
                            if(C.gnssC.ion_opt!=ION_IF){
                                sat_info->float_amb[frq]=ambj/sat_info->lam[frq];
                                ref_sat_info->float_amb[frq]=ambi/ref_sat_info->lam[frq]; //cycle;
                            }
                            else{
                                sat_info->float_amb[frq]=ambj;    //m
                                ref_sat_info->float_amb[frq]=ambi;
                            }
                        }
                    }

                    sat_info->vsat[frq]=1;
                    sat_info->res_idx[frq]=num_L_;
                    previous_sat_info_[sat_info->sat.sat_.no-1].vsat[frq]=1;

                    sprintf(buff,"(%9.5f - %9.5f)-(%9.5f - %9.5f)-(%9.5f - %9.5f)=%9.5f el=%3.1f Ri=%10.6f Rj=%10.6f glo_bias=%5.3f",
                            rover_res[ref_sat_idx*nf+ifrq],base_res[ref_sat_idx*nf+ifrq],rover_res[isat*nf+ifrq],base_res[isat*nf+ifrq],ambi,ambj,omcs.back(),sat_info->el_az[0]*R2D,Ri.back(),Rj.back(),glo_bia);
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<(obs_code?"P":"L")<<frq+1<<": (ROVER_"<<ref_sat_info->sat.sat_.id<<"-BASE_"<<ref_sat_info->sat.sat_.id<<")-(ROVER_"<<sat_info->sat.sat_.id<<"-BASE_"<<sat_info->sat.sat_.id<<")-amb"
                              <<" = "<<buff;

                    vflag_.push_back((ref_sat_info->sat.sat_.no<<16|(sat_info->sat.sat_.no<<8)|(obs_code<<4)|(frq)));
                    num_L_++;
                    if(ifrq==0) num_valid_sat_++;
                    nobs[isys][iobs]++;
                }
                iobs++;
            }
        }

        R.resize(num_L_*num_L_,0.0);
        for(int isys=0,nt=0;isys<NSYS;isys++){
            for(int f=0;f<nf;f++){
                for(int j=0;j<nobs[isys][f];j++){
                    for(int k=0;k<nobs[isys][f];k++){
                        R[nt+k+(nt+j)*num_L_]=j==k?Ri[nt+j]+Rj[nt+j]:Ri[nt+j];
                    }
                }
                nt+=nobs[isys][f];
            }
        }

        if(!post){
            H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
            R_=Map<MatrixXd>(R.data(),num_L_,num_L_);
            omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        }

        if(post&&post!=5&&C.gnssC.res_qc){
            qc_iter_++;
           if(PpkResidualQc(post,ir,cmn_sat_no,omcs,Rj)){
               need_iter=0;
           }
        }

        H.clear();Ri.clear();Rj.clear();R.clear();omcs.clear();
        return post?need_iter:num_L_;
    }

    bool cPpkSolver::ValidObs(int i, int nf, int f) {
        return (rover_res[i*nf*2+f]!=0.0&&base_res[i*nf*2+f]!=0.0&&(f<nf||rover_res[f-nf+i*nf*2]!=0.0&&base_res[f-nf+i*nf*2]!=0.0));
    }

    bool cPpkSolver::MatchBaseObs(cTime t) {
        double sow1,sow2;
        int i,week=0,wod=0,info=false;

        for(i=base_idx_-100<0?0:base_idx_-10;base_idx_<base_obs_.GetGnssObs().size();i++){
            sow1=base_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(&week,&wod,SYS_GPS);
            sow2=t.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<DTTOL){
                base_idx_=i;info=true;
                base_epoch_sat_obs_=base_obs_.GetGnssObs().at(base_idx_);
                break;
            }
            else if((sow1-sow2)>2.0*DTTOL){
                info=false;break;
            }
        }

        return info;
    }

    int cPpkSolver::SelectCmnSat(tIGCOUPLEDConf C, vector<int> &ir, vector<int> &ib, vector<int> &cmn_sat_no) {
        string buff;
        tSatInfoUnit rover,base;
        for(int i=0,j=0;i<epoch_sat_info_collect_.size()&&j<base_sat_info_collect_.size();i++,j++){
            rover=epoch_sat_info_collect_.at(i);
            base=base_sat_info_collect_.at(j);
            if(rover.sat.sat_.no<base.sat.sat_.no) j--;
            else if(rover.sat.sat_.no>base.sat.sat_.no) i--;
            else if(rover.el_az[0]>=C.gnssC.ele_min*D2R){
                if(rover.sat.sat_.no!=base.sat.sat_.no) continue;
                ir.push_back(i);
                ib.push_back(j);
                cmn_sat_no.push_back(rover.sat.sat_.no);
            }
        }

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            buff+=epoch_sat_info_collect_.at(i).sat.sat_.id+" ";
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PPK ROVER STATION OBSERVED SATELLITES: "<<epoch_sat_info_collect_.size()<<" "<<buff;
        buff.clear();

        for(int i=0;i<base_sat_info_collect_.size();i++){
            buff+=base_sat_info_collect_.at(i).sat.sat_.id+" ";
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PPK BASE STATION OBSERVED SATELLITES : "<<base_sat_info_collect_.size()<<" "<<buff;
        buff.clear();

        for(int i=0;i<ir.size();i++){
            buff+=epoch_sat_info_collect_.at(ir[i]).sat.sat_.id+" ";
        }

        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PPK ROVER AND BASE COMMON SATELLITE  : "<<ir.size()<<" "<<buff;

        return ir.size();
    }

    void cPpkSolver::PpkCycleSlip(tIGCOUPLEDConf C,vector<int>& iu,vector<int>& ib,vector<int>& cmn_sat_no) {

        tSatInfoUnit* sat_info= nullptr;
        tSatInfoUnit* base_sat= nullptr;
        cTime t=igcoupled_sol_.t_tag;
        double dt=C.gnssC.sample_rate*(C.filter_type?-1.0:1.0);
        int f;

        if(t.t_.long_time!=0.0) dt=spp_solver_->igcoupled_sol_.t_tag.TimeDiff(igcoupled_sol_.t_tag.t_);
//        if(t.t_.long_time!=0.0) dt=epoch_sat_info_collect_[0].t_tag.TimeDiff(t.t_);

        for(int i=0;i<cmn_sat_no.size();i++){
            sat_info=&epoch_sat_info_collect_.at(iu[i]);
            base_sat=&base_sat_info_collect_.at(ib[i]);

            for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++) previous_sat_info_[sat_info->sat.sat_.no-1].slip[j]&=0xFC;
            gnss_obs_operator_.LliCycleSlip(C,*sat_info,previous_sat_info_[sat_info->sat.sat_.no-1],2,dt,REC_ROVER);

            if(gnss_obs_operator_.LliCycleSlip(C,*base_sat,previous_sat_info_[sat_info->sat.sat_.no-1],2,dt,REC_BASE)){
                for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++) sat_info->slip[j]|=base_sat->slip[j];
            }
            gnss_obs_operator_.MwCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, base_sat,previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
            gnss_obs_operator_.GfCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, base_sat);

            gnss_obs_operator_.SmoothMw(C,sat_info, base_sat);

//            for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++) previous_sat_info_[sat_info->sat.sat_.no-1].slip[j]=sat_info->slip[j];
        }


    }

    void cPpkSolver::StateTimeUpdate(tIGCOUPLEDConf C,vector<int>& ir,vector<int>& ib,vector<int>& cmn_sat_no) {
        double tt=spp_solver_->igcoupled_sol_.t_tag.TimeDiff(igcoupled_sol_.t_tag.t_);
        //position
        PosUpdate(C);

        // glo ifb
        GloIfpbUpdate(C,tt);

        //tropospheric delay
        TrpUpdate(C,tt);
        //ionospheric delay
        IonUpdate(C,tt);
        //ambiguity
        AmbUpdate(C,tt,ir,ib,cmn_sat_no);
    }

    void cPpkSolver::PosUpdate(tIGCOUPLEDConf  C) {
        if(para_.NumPos()<=0) return;
        if(tc_mode_) return;

        Vector3d q(SQR(30),SQR(30),SQR(30));
        int ip=para_.IndexPos();
        if((SQR(full_x_[ip])+SQR(full_x_[ip+1])+SQR(full_x_[ip+1]))==0.0){
            for(int i=0;i<3;i++) InitX(spp_solver_->igcoupled_sol_.pos[i],SQR(30.0),ip+i,full_x_.data(),full_Px_.data());
        }else{
            if(C.mode_opt==MODE_OPT_STATIC) return;
            if(C.mode_opt==MODE_OPT_KINE_SIM||C.mode_opt==MODE_OPT_KINEMATIC||C.mode_opt==MODE_OPT_PPK){
                for(int i=0;i<3;i++) InitX(spp_solver_->igcoupled_sol_.pos[i],SQR(30.0),ip+i,full_x_.data(),full_Px_.data());
            }
        }
    }

    void cPpkSolver::GloIfpbUpdate(tIGCOUPLEDConf C,double tt) {
        if(para_.NumGloIfpb()<=0) return;

        int i,j;
        for(i=0;i<2;i++){
            j=para_.IndexGloIfpb()+i;

            if(full_x_[j]==0.0){
                InitX(C.gnssC.ar_thres[2]+1E-6,C.gnssC.ar_thres[3],j,full_x_.data(),full_Px_.data());
            }
            else if(num_continuous_fix_>C.gnssC.min_fix2hold){
                return;
            }
            else{
                full_Px_(j,j)+=SQR(C.gnssC.ar_thres[4])*fabs(tt);
            }
        }
    }

    void cPpkSolver::TrpUpdate(tIGCOUPLEDConf C, double tt) {
        if(para_.NumTrp()<=0) return;
    }

    void cPpkSolver::IonUpdate(tIGCOUPLEDConf C, double tt) {
        if(para_.NumIon()<=0) return;
    }

    void cPpkSolver::AmbUpdate(tIGCOUPLEDConf C, double tt,vector<int>& ir,vector<int>& ib,vector<int>& cmn_sat_no) {
        if(para_.NumAmb()<=0) return;

        int i,j,reset,ia,slip,rejc;
        double offset;
        for(int f=0;f<para_.GetGnssUsedFrqs();f++){

            for(i=0;i<MAX_SAT_NUM;i++){
                ia=para_.IndexAmb(f,i+1);
                reset=previous_sat_info_[i].outc[f]>C.gnssC.max_out;
                if(C.gnssC.ar_mode==AR_INST&&full_x_[ia]!=0.0){
                    InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                }
                else if(reset&&full_x_[ia]!=0.0){
                    InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    previous_sat_info_[i].outc[f]=0;
                }
                if(C.gnssC.ar_mode!=AR_INST&&reset){
                    previous_sat_info_[i].lock[f]=-C.gnssC.min_lock2fix;
                }
            }

            for(i=0;i<cmn_sat_no.size();i++){
                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                full_Px_(ia,ia)+=SQR(C.gnssC.ait_psd[0])*fabs(tt);
                slip=epoch_sat_info_collect_[ir[i]].slip[f];
                rejc=epoch_sat_info_collect_[ir[i]].rejc[f];
                if(C.gnssC.ar_mode==AR_INST||(!(slip&1)&&rejc<2)) continue;
                full_x_[ia]=0.0;
                epoch_sat_info_collect_[ir[i]].rejc[f]=0;
                epoch_sat_info_collect_[ir[i]].lock[f]=-C.gnssC.min_lock2fix;
                if(previous_sat_info_[ir[i]-1].sat.sat_.sys!=SYS_GLO) previous_sat_info_[ir[i]-1].glo_bias[f]=0.0;
            }

            double cp,pr,amb,com_offset;
            vector<double>bias;
            for(i=0,j=0,offset=0.0;i<cmn_sat_no.size();i++){
                if(C.gnssC.ion_opt!=ION_IF){
                    if(C.gnssC.frq_opt==FRQ_DUAL&&C.gnssC.check_dual_phase){
                        if(!CheckDualFrq(epoch_sat_info_collect_[ir[i]])){
                            epoch_sat_info_collect_[ir[i]].stat=SAT_NO_USE;
                            bias.push_back(0.0);
                            continue;
                        }
                    }
                    cp=gnss_obs_operator_.GnssSdObs(epoch_sat_info_collect_[ir[i]],base_sat_info_collect_[ib[i]],f,GNSS_OBS_PHASE);
                    pr=gnss_obs_operator_.GnssSdObs(epoch_sat_info_collect_[ir[i]],base_sat_info_collect_[ib[i]],f,GNSS_OBS_CODE);
                    if(cp==0.0||pr==0.0){
                        if(cp==0.0){
                            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[ir[i]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[i]].sat.sat_.id<<" MISSING L"<<f+1<<" FOR INITIALIZING AMBIGUITY";
                        }
                        if(pr==0.0){
                            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[ir[i]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[i]].sat.sat_.id<<" MISSING P"<<f+1<<" FOR INITIALIZING AMBIGUITY";
//                            epoch_sat_info_collect_[ir[i]].stat=SAT_NO_PR;
                        }
                        bias.push_back(0.0);
                        continue;
                    }
                    amb=cp-pr/epoch_sat_info_collect_[ir[i]].lam[f];
                    bias.push_back(amb);
                }
                else{
                    // ????
                }

                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                if(full_x_[ia]!=0.0){
                    offset+=amb-full_x_[ia]*epoch_sat_info_collect_[ir[i]].lam[f];
                    j++;
                }
            }

            com_offset=j>0?offset/j:0.0;

            for(i=0;i<cmn_sat_no.size();i++){
                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                if(bias[i]==0.0||full_x_[ia]!=0.0) continue;
                InitX((bias[i]-com_offset),SQR(30),ia,full_x_.data(),full_Px_.data());
                epoch_sat_info_collect_[ir[i]].lock[f]=-C.gnssC.min_lock2fix;
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_.at(ir[i]).t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_.at(ir[i]).sat.sat_.id<<" L"<<f+1<<" AMBIGUITY INITIALIZED "<<(bias[i]-com_offset);
            }
            bias.clear();
        }
    }

    // can improved AR rate
    bool cPpkSolver::CheckDualFrq(tSatInfoUnit &sat_info) {
        for(int i=0;i<2;i++){
            if(sat_info.raw_L[i]==0.0) return false;
        }
        return true;
    }

    bool cPpkSolver::PpkResidualQc(int iter,vector<int>ir, vector<int> cmn_sat_no,vector<double>& omcs, vector<double>R) {
        bool flag=false;
        vector<double>v_code,v_phase,norm_v_code,norm_v_phase;
        vector<int> code_idx,phase_idx;
        int i,sat1,sat2,f,f1,f2,obs_type,sat_idx;
        double el;
        double thres_code=3.0,thres_phase=0.03;
        double k0=2.0,k1=3.0;
        bool reduce_weight=false;

        if(!kf_.qc_flag) return false;
        if(qc_iter_>1) return false;

        for(i=0;i<omcs.size();i++){
            obs_type=(vflag_[i]>>4)&0xF;
            sat1=(vflag_[i]>>16)&0xFF;
            sat2=(vflag_[i]>>8)&0xFF;
            f=(vflag_[i])&0xF;

            if(obs_type==1){
                code_idx.push_back(i);
                v_code.push_back(fabs(omcs[i]));
                norm_v_code.push_back(fabs(omcs[i])/SQRT(fabs(R[i])));
            }
            else if(obs_type==0){
                phase_idx.push_back(i);
                v_phase.push_back(fabs(omcs[i]));
                norm_v_phase.push_back(fabs(omcs[i])/SQRT(fabs(R[i])));
            }
        }

        // step-1: test pseudorange residual
        bool code_flag=false;
        auto max_v_code=max_element(begin(v_code),end(v_code));
        int idx_max_v_code=distance(begin(v_code),max_v_code);

        sat1=(vflag_[code_idx[idx_max_v_code]]>>16)&0xFF;
        sat2=(vflag_[code_idx[idx_max_v_code]]>>8)&0xFF;
        f=(vflag_[code_idx[idx_max_v_code]])&0xF;
        int j1=0;
        for(j1=0;j1<cmn_sat_no.size();j1++){
            if(cmn_sat_no[j1]==sat2) break;
        }
        el=epoch_sat_info_collect_[ir[j1]].el_az[0];

        if(epoch_sat_info_collect_[ir[j1]].sat.sat_.sys==SYS_BDS){
            if(epoch_sat_info_collect_[ir[j1]].sat.sat_.prn<=5){
                thres_code=10.0;
            }
            else{
                thres_code=6.0;
            }
        }

        if(*max_v_code>thres_code/sin(el)){
            code_flag=true;
        }

        // step-2: test standardized pseudorange residual
        bool norm_code_flag=false;
        auto max_norm_v_code=max_element(begin(norm_v_code),end(norm_v_code));
        int idx_max_norm_v_code=distance(begin(norm_v_code),max_norm_v_code);
        sat1=(vflag_[code_idx[idx_max_norm_v_code]]>>16)&0xFF;
        sat2=(vflag_[code_idx[idx_max_norm_v_code]]>>8)&0xFF;
        f=(vflag_[code_idx[idx_max_norm_v_code]])&0xF;

        int j2=0;
        for(j2=0;j2<cmn_sat_no.size();j2++){
            if(cmn_sat_no[j2]==sat2) break;
        }
        if(*max_norm_v_code>k1){
            norm_code_flag=true;
        }

        if(code_flag){
            flag=true;
            el=epoch_sat_info_collect_[ir[j1]].el_az[0];
            epoch_sat_info_collect_[ir[j1]].stat=SAT_NO_USE;
            f=(vflag_[code_idx[idx_max_v_code]])&0xF;
            if(iter==1) epoch_sat_info_collect_[ir[j1]].rejc[f]++;
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" P"<<f+1<<" EXCULDED BY CODE RESIDUAL, res="<<*max_v_code<<" el="<<el*R2D;
            return flag;
        }

        if(norm_code_flag){
            flag=true;
            el=epoch_sat_info_collect_[ir[j2]].el_az[0];
            epoch_sat_info_collect_[ir[j2]].stat=SAT_NO_USE;
            f=(vflag_[code_idx[idx_max_norm_v_code]])&0xF;
            if(iter==1) epoch_sat_info_collect_[ir[j2]].rejc[f]++;
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j2]].sat.sat_.id<<" P"<<f+1<<" EXCULDED BY NORM CODE RESIDUAL, norm_res="<<*max_norm_v_code<<" el="<<el*R2D;
            return flag;
        }

        int sat_no;
        double fact=k0/(*max_norm_v_code)*SQR((k1-*max_norm_v_code)/(k1-k0));
        if(k0<=*max_norm_v_code&&*max_norm_v_code<=k1){
            sat_no=epoch_sat_info_collect_[ir[j2]].sat.sat_.no;
            f=(vflag_[code_idx[idx_max_norm_v_code]])&0xF;
            previous_sat_info_[sat_no].c_var_factor[f]*=1.0/fact;
            reduce_weight=true;
        }

#if 0
        if(v_phase.size()){
            auto max_v_phase=max_element(begin(v_phase),end(v_phase));
            int idx_max_v_phase=distance(begin(v_phase),max_v_phase);

            sat1=(vflag_[phase_idx[idx_max_v_phase]]>>16)&0xFF;
            sat2=(vflag_[phase_idx[idx_max_v_phase]]>>8)&0xFF;
            f1=(vflag_[phase_idx[idx_max_v_phase]])&0xF;
            for(j1=0;j1<cmn_sat_no.size();j1++){
                if(cmn_sat_no[j1]==sat2) break;
            }
            el=epoch_sat_info_collect_[ir[j1]].el_az[0];

            if(epoch_sat_info_collect_[ir[j1]].sat.sat_.sys==SYS_BDS){
                thres_phase=0.06;
            }

            int ia;
            if(*max_v_phase>thres_phase/sin(el)){
                sat_no=epoch_sat_info_collect_[ir[j1]].sat.sat_.no;
                if(iter==1) epoch_sat_info_collect_[ir[j1]].rejc[f]++;
                if(iter==1) epoch_sat_info_collect_[ir[j1]].rejc_phase[f]++;
                ia=para_.IndexAmb(f,sat_no);
                if(epoch_sat_info_collect_[ir[j1]].rejc_phase[f]>1){
                    InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                    epoch_sat_info_collect_[ir[j1]].rejc_phase[f]=0;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" L"<<f+1<<" RESET AMB BY PHASE RESIDUAL,  res="<<*max_v_phase<<" el="<<el*R2D;
                }
                else{
                    previous_sat_info_[epoch_sat_info_collect_[ir[j1]].sat.sat_.no-1].p_var_factor[f]*=10000;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" L"<<f+1<<" REDUCE WEIGHT BY NORM PHASE RESIDUAL, res="<<*max_v_phase<<" el="<<el*R2D;
                }
                return true;
            }

            // step-4: test standardized phase residual
            auto max_norm_v_phase=max_element(begin(norm_v_phase),end(norm_v_phase));
            int idx_max_norm_v_phase=distance(begin(norm_v_phase),max_norm_v_phase);
            sat1=(vflag_[phase_idx[idx_max_norm_v_phase]]>>16)&0xFF;
            sat2=(vflag_[phase_idx[idx_max_norm_v_phase]]>>8)&0xFF;
            f2=(vflag_[phase_idx[idx_max_norm_v_phase]])&0xF;

            for(j2=0;j2<cmn_sat_no.size();j2++){
                if(cmn_sat_no[j2]==sat2) break;
            }

            if(*max_norm_v_phase>k1){
                sat_no=epoch_sat_info_collect_[ir[j2]].sat.sat_.no;
                if(iter==1) epoch_sat_info_collect_[ir[j2]].rejc_phase[f]++;
                if(iter==1) epoch_sat_info_collect_[ir[j2]].rejc[f]++;
                ia=para_.IndexAmb(f2,sat_no);
                el=epoch_sat_info_collect_[ir[j2]].el_az[0];

                if(epoch_sat_info_collect_[ir[j2]].rejc_phase[f]>1){
                    InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                    epoch_sat_info_collect_[ir[j2]].rejc_phase[f]=0;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j2]].sat.sat_.id<<" L"<<f+1
                    <<" RESET AMB BY NORM PHASE RESIDUAL,norm_res="<<*max_norm_v_phase<<" el="<<el*R2D<<" res="<<v_phase[idx_max_norm_v_phase];
                    return true;
                }
                else{
//                    epoch_sat_info_collect_[ir[j2]].stat=SAT_NO_USE;
                    previous_sat_info_[epoch_sat_info_collect_[ir[j2]].sat.sat_.no-1].p_var_factor[f]*=10000;
                    CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j2]].sat.sat_.id<<" L"<<f+1
                    <<" REDUCE WEIGHT BY NORM PHASE RESIDUAL,norm_res="<<*max_norm_v_phase<<" el="<<el*R2D<<" res="<<v_phase[idx_max_norm_v_phase];
                    return true;
                }
            }

//            fact=k0/(*max_norm_v_phase)*SQR((k1-*max_norm_v_phase)/(k1-k0));
//            if(k0<=*max_norm_v_phase&&*max_norm_v_phase<=k1){
//                sat_no=epoch_sat_info_collect_[ir[j2]].sat.sat_.no;
//                f=(vflag_[code_idx[idx_max_norm_v_code]])&0xF;
//                previous_sat_info_[sat_no].p_var_factor[f]*=1.0/fact;
//                reduce_weight=true;
//            }
        }
#endif

        v_code.clear();v_phase.clear();norm_v_code.clear();norm_v_phase.clear();code_idx.clear();phase_idx.clear();
        return reduce_weight;
    }

    void cPpkSolver::ReSetAmb(double *bias, double *xa,int nb) {
        int i,n,m,f,index[MAX_SAT_NUM]={0},nv=0,nf=para_.GetGnssUsedFrqs();
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<num_real_x_fix_;i++) xa[i]=real_x_fix_[i];

#if 0
        for(m=0;m<NSYS;m++){
            for(f=0;f<nf;f++){
                for(n=i=0;i<MAX_SAT_NUM;i++){
                    sat_info=&previous_sat_info_[i];
                    if(sat_info->sat.sat_.sys_idx!=m||sat_info->fix[f]!=2) continue;
                    index[n++]=para_.IndexAmb(f,sat_info->sat.sat_.no);
                }

                if(n<2) continue;
                xa[index[0]]=full_x_[index[0]];

                for(i=1;i<n;i++){
                    xa[index[i]]=xa[index[0]]-bias[nv++];
                }
            }
        }
#endif
        int sat1,sat2,ib1,ib2;
        for(i=0;i<nb;i++){
            sat1=ddambs_[i].ref_sat;
            sat2=ddambs_[i].sat;
            f=ddambs_[i].f;

            ib1=para_.IndexAmb(f,sat1);
            ib2=para_.IndexAmb(f,sat2);
            // 基准卫星仍是浮点解
            xa[ib2]=xa[ib1]-bias[i];

            previous_sat_info_[sat1-1].fix_amb[f]=xa[ib1];
            previous_sat_info_[sat2-1].fix_amb[f]=-bias[i];
        }
    }

    int cPpkSolver::DdMat(double *D, int gps, int bds, int glo,int gal,double el_mask) {
        int i,j,k,f,nb=0,na=num_real_x_fix_,nf=para_.GetGnssUsedFrqs(),no_fix;
        double fix[MAX_SAT_NUM],ref[MAX_SAT_NUM],s[2];
        tSatInfoUnit *ref_sat_info= nullptr;
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<MAX_SAT_NUM;i++) for(j=0;j<MAX_GNSS_USED_FRQ_NUM;j++){
            previous_sat_info_[i].fix[j]=0;
        }

        for(i=0;i<na;i++) D[i+i*num_full_x_]=1.0;

#if 0
        for(m=0;m<NSYS;m++){

            no_fix=(m==0&&gps==0)||(m==1&&ppk_conf_.gnssC.bds_ar_mode==0)||(m==3&&glo==0);

            for(f=0,k=na;f<nf;f++,k+=MAX_SAT_NUM){

                for(i=k;i<k+MAX_SAT_NUM;i++){
                    sat_info=&previous_sat_info_[i-k];
                    if(full_x_[i]==0.0||sat_info->sat.sat_.sys_idx!=m||!sat_info->vsat[f]){
                        continue;
                    }
                    if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&sat_info->el_az[0]*R2D>ppk_conf_.gnssC.ar_el_mask&&!no_fix){
                        sat_info->fix[f]=2;
                        break;
                    }
                    else sat_info->fix[f]=1;
                }

                if(sat_info->fix[f]!=2) continue;

                for(j=k;j<k+MAX_SAT_NUM;j++){
                    sat_info=&previous_sat_info_[j-k];
                    if(i==j||full_x_[j]==0.0||sat_info->sat.sat_.sys_idx!=m||!sat_info->vsat[f]) continue;
                    if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&previous_sat_info_[i-k].vsat[f]&&sat_info->el_az[0]*R2D>=ppk_conf_.gnssC.ar_el_mask&&!no_fix){
                        D[i+(na+nb)*num_full_x_]=1.0;
                        D[j+(na+nb)*num_full_x_]=-1.0;
                        ref[nb]=i-k+1;
                        fix[nb++]=j-k+1;
                        sat_info->fix[f]=2;
                    }
                    else sat_info->fix[f]=1;
                }

            }
        }
#endif
        int sat1,sat2,sys=SYS_GPS;
        if(ddambs_.size()) ddambs_.clear();
        tDdAmb dd_amb={0};
        for(i=0;i<vflag_.size();i++){
            if(((vflag_[i]>>4)&0xF)==1) continue;
            sat1=(vflag_[i]>>16)&0xFF;
            sat2=(vflag_[i]>> 8)&0xFF;
            f=vflag_[i]&0xF;
            if(previous_sat_info_[sat1-1].sat.sat_.sys==SYS_BDS){
                continue;
            }
            sys=previous_sat_info_[sat1-1].sat.sat_.sys;

            no_fix=(sys==SYS_GPS&&gps==0)||(sys==SYS_BDS&&bds==0)||(sys==SYS_GLO&&glo==0)||(sys==SYS_GAL&&gal==0);
            ref_sat_info=&previous_sat_info_[sat1-1];
            sat_info=&previous_sat_info_[sat2-1];
            k=para_.IndexAmb(f,sat1);
            j=para_.IndexAmb(f,sat2);

            ref_sat_info->fix[f]=2;
            if(full_x_[j]==0.0||full_x_[k]==0.0) continue;

            if(sat_info->slip[f]) continue;

            if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&sat_info->vsat[f]&&sat_info->el_az[0]*R2D>=el_mask&&!no_fix){
                D[k+(na+nb)*num_full_x_]= 1.0;
                D[j+(na+nb)*num_full_x_]=-1.0;

                dd_amb.ref_sat=ref_sat_info->sat.sat_.no;
                dd_amb.sat=sat_info->sat.sat_.no;
                dd_amb.f=f;
                ddambs_.push_back(dd_amb);
                ref[nb]=ref_sat_info->sat.sat_.no;
                fix[nb++]=sat_info->sat.sat_.no;
                sat_info->fix[f]=2;
            }
            else sat_info->fix[f]=1;
        }

        if(nb>0){
            VectorXd ref_sats,fix_sats;
            ref_sats=Map<VectorXd>(ref,nb,1);
            fix_sats=Map<VectorXd>(fix,nb,1);

            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"REF SATELLITES: "<<ref_sats.transpose();
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"FIX SATELLITES: "<<fix_sats.transpose();
        }
        return nb;
    }

    int cPpkSolver::ResolveAmbLambda(double *xa, int gps, int bds,int glo,int gal, int qzs,double el_mask) {
        int i,j,ny,nb,info,na=num_real_x_fix_;
        VectorXd QQb(MAX_SAT_NUM);
        double var=0.0;
        VectorXd s(2);
        string buff;

        igcoupled_sol_.ratio=0.0;
        if(ppk_conf_.gnssC.ar_mode==AR_OFF||ppk_conf_.gnssC.ar_thres[0]<1.0){
            igcoupled_sol_.num_ar_sat=0;
            return 0;
        }

        // skip AR if position variance too lager to avoid false fix
#if 1
        for(i=0;i<3;i++) var+=full_Px_(i,i);
        var=var/3.0;
        if(var>ppk_conf_.gnssC.ar_thres[1]){ //0.25
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"POSITIONING VARIANCE TOO LARGE FOR PPK AR var="<<var;
            igcoupled_sol_.num_ar_sat=0;
            return 0;
        }
#endif
        // create single to DD transformation matrix,used to translate SD phase biases to DD
        MatrixXd D(num_full_x_,num_full_x_);
        D=Eigen::MatrixXd::Zero(num_full_x_,num_full_x_);

        if((nb=DdMat(D.data(),gps,bds,glo,gal,el_mask))<(ppk_conf_.gnssC.min_sat_num2fix-1)){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"NO ENOUGH VALID SATELLITE MAKING DOUBLE DIFFERENCE MATRIX";
            return -1;
        }
        igcoupled_sol_.num_ar_sat=nb;

        ny=na+nb;
        VectorXd y(ny);
        MatrixXd Qy(ny,ny),DP(ny,num_full_x_);
        VectorXd db(nb,1);
        MatrixXd b(nb,2),Qb(nb,nb),Qab(na,nb),QQ(na,nb);

        y=D.transpose()*full_x_;
        MatMul("TN",ny,num_full_x_,num_full_x_,1.0,D.data() ,full_Px_.data(),0.0,DP.data());   /* DP=D'*P */
        MatMul("NN",ny,ny,num_full_x_,1.0,DP.data(),D.data(),0.0,Qy.data());

        //phase-bias covariance(Qb) and real-parameters to bias covariance(Qab)
        for(i=0;i<nb;i++){
            QQb[i]=Qy.data()[(na+i)+(na+i)*ny];
            for(j=0;j<nb;j++){
                Qb.data()[i+j*nb]=Qy.data()[(na+i)+(na+j)*ny];
            }
        }
        for(i=0;i<na;i++) for(j=0;j<nb;j++) Qab.data()[i+j*na]=Qy.data()[i+(na+j)*ny];

        VectorXd float_amb(nb),fix_float_diff(nb);
        float_amb=Map<VectorXd>(y.data()+na,nb,1); //DD ambiguity

        char ss[20];
        for(i=0;i<float_amb.size();i++){
            sprintf(ss,"%5.2f  ",float_amb[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"DD_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        for(i=0;i<float_amb.size();i++){
            sprintf(ss,"%7.5f  ",QQb[i]);
            buff+=ss;
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<" Q_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        if(!(info=lambda_.IntegerAmb(float_amb,Qb,b,nb,2,s))){
            igcoupled_sol_.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
            if(igcoupled_sol_.ratio>999.9) igcoupled_sol_.ratio=999.9f;

            for(i=0;i<nb;i++){
                sprintf(ss,"%5.2f  ",b.col(0)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"DD_AMB(1): "<<buff;
            buff.clear();
            ss[0]='\0';

            for(i=0;i<nb;i++){
                sprintf(ss,"%5.2f  ",b.col(1)[i]);
                buff+=ss;
            }
            CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"DD_AMB(2): "<<buff;
            buff.clear();
            ss[0]='\0';

            // validation by popular ratio-test of residuals
            if(s[0]<=0.0||s[1]/s[0]>=ppk_conf_.gnssC.ar_thres[0]){
                // init non phase-bias states and covariance with float solution values
                for(i=0;i<na;i++){
                    real_x_fix_[i]=full_x_[i];
                    for(j=0;j<na;j++) real_Px_fix_.data()[i+j*na]=full_Px_.data()[i+j*num_full_x_];
                }

                for(i=0;i<nb;i++){
                    fix_float_diff[i]=y.data()[na+i]-b.col(0).data()[i];
                }

                MatrixXd Qb_=Qb.inverse();
                db=Qb_*fix_float_diff; //db=Qb^-1*(b0-b)
                real_x_fix_=real_x_fix_-Qab*db;   //xa=x-Qab-db
                QQ=Qab*Qb_;            // QQ=Qab*Qb^-1
                real_Px_fix_=real_Px_fix_-QQ*Qab.transpose();  //Pa=P-QQ*Qab^T

                ReSetAmb(b.col(0).data(),xa,nb);
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<" AR VALIDATION OK, ratio="<<s[1]/s[0];
            }
            else{
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<" RATIO VALIDATION FAILED ratio="<<s[1]/s[0];
                nb=0;
            }
        }
        else{
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"PPK LAMBDA ERROR";
            nb=0;
        }

        return nb;
    }

    // function from rtklib
    bool cPpkSolver::ResolvePpkAmb(vector<int> cmn_sat, int nf, double *xa) {
        int nb=0,f,i,ar=0,lockc[MAX_GNSS_USED_FRQ_NUM],arsats[64]={0};
        bool exc_flag=false;
        tSatInfoUnit *sat_info= nullptr;

        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"PREVIOUS EPOCH AR: ratio1="<<pre_epoch_ar_ratio1<<" ratio2="<<pre_epoch_ar_ratio2<<" ar_sat_num="<<igcoupled_sol_.num_ar_sat;

        // 若前一历元没有AR成功
        if(pre_epoch_ar_ratio2<ppk_conf_.gnssC.ar_thres[0]&&igcoupled_sol_.num_ar_sat>=ppk_conf_.gnssC.min_sat_num2drop){
            // previous epoch ar failed and satellite enough to drop
            for(f=0;f<nf;f++){
                for(i=0;i<cmn_sat.size();i++){
                    sat_info=&previous_sat_info_[cmn_sat[i]-1];
                    if(sat_info->vsat[f]&&sat_info->lock[f]>=0&&sat_info->el_az[0]*R2D>=ppk_conf_.gnssC.ele_min){
                        arsats[ar++]=i;
                    }
                }
            }

            if(exc_sat_index<ar){
                i=cmn_sat[arsats[exc_sat_index]];
                for(f=0;f<nf;f++){
                    lockc[f]=previous_sat_info_[i-1].lock[f];
                    previous_sat_info_[i-1].lock[f]=-igcoupled_sol_.num_ar_sat;
                }
                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[i-1].sat.sat_.id<<" EXCLUDE BY AR";
                exc_flag=true;
            }
            else exc_sat_index=0;
        }

        double ratio1;
        int dly;
        bool rerun=false;
        if(ppk_conf_.gnssC.glo_ar_mode!=GLO_AR_FIXHOLD||amb_hold_flag){

            nb=ResolveAmbLambda(xa,ppk_conf_.gnssC.ar_mode,ppk_conf_.gnssC.bds_ar_mode,ppk_conf_.gnssC.glo_ar_mode,ppk_conf_.gnssC.gal_ar_mode,0,ppk_conf_.gnssC.ar_el_mask);

            ratio1=igcoupled_sol_.ratio;
            if(ppk_conf_.gnssC.ar_filter){
                if((nb>=0&&pre_epoch_ar_ratio2>=ppk_conf_.gnssC.ar_thres[0]&&(igcoupled_sol_.ratio<ppk_conf_.gnssC.ar_thres[0]))||
                        (igcoupled_sol_.ratio<ppk_conf_.gnssC.ar_thres[0]*1.1&&igcoupled_sol_.ratio<pre_epoch_ar_ratio1/2.0)){
                    dly=2;
                    for(i=0;i<cmn_sat.size();i++){
                        for(f=0;f<nf;f++){
                            if(previous_sat_info_[cmn_sat[i]-1].fix[f]!=2) continue;
                            if(previous_sat_info_[cmn_sat[i]-1].lock[f]==0){
                                CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" AR FILTER "<<previous_sat_info_[cmn_sat[i]-1].sat.sat_.id<<" L"<<f+1<<" IS NEW OBSERVED";
                                previous_sat_info_[cmn_sat[i]-1].lock[f]=-ppk_conf_.gnssC.min_lock2fix-dly;
                                dly+=2;
                                rerun=true;
                            }
                        }
                    }
                }

                if(rerun){
                    if((nb=ResolveAmbLambda(xa,ppk_conf_.gnssC.ar_mode,ppk_conf_.gnssC.bds_ar_mode,ppk_conf_.gnssC.glo_ar_mode,ppk_conf_.gnssC.gal_ar_mode,0,ppk_conf_.gnssC.ar_el_mask))<2){
                        pre_epoch_ar_ratio1=pre_epoch_ar_ratio2=ratio1;

                        if(ppk_conf_.gnssC.par_ar){
                            double el=ppk_conf_.gnssC.ar_el_mask;
                            if(igcoupled_sol_.num_ar_sat>3){
                                for(int j=0;j<igcoupled_sol_.num_ar_sat;j++){
                                    el+=5.0;
                                    if(igcoupled_sol_.num_ar_sat-j<3) return false;
                                    ResolveAmbLambda(xa,ppk_conf_.gnssC.ar_mode,ppk_conf_.gnssC.bds_ar_mode,ppk_conf_.gnssC.glo_ar_mode,ppk_conf_.gnssC.gal_ar_mode,0,el);
                                    if(igcoupled_sol_.num_ar_sat<3){
                                        return false;
                                    }
                                }
                            }
                        }
                        else
                            return false;
                    }
                }
            }
            pre_epoch_ar_ratio1=ratio1;
        }
        else{
            ratio1=0.0;
            nb=0;
        }

        if(ppk_conf_.gnssC.glo_ar_mode==GLO_AR_FIXHOLD){

        }

        if(exc_flag&&(igcoupled_sol_.ratio<ppk_conf_.gnssC.ar_thres[0])&&igcoupled_sol_.ratio<1.5*pre_epoch_ar_ratio2){
            i=cmn_sat[arsats[exc_sat_index++]];
            for(f=0;f<nf;f++) previous_sat_info_[i-1].lock[f]=lockc[f];
        }

        pre_epoch_ar_ratio1=ratio1>0?ratio1:igcoupled_sol_.ratio;
        pre_epoch_ar_ratio2=igcoupled_sol_.ratio;

        if(nb>0) return true;
        else return false;
    }

    void cPpkSolver::HoldAmb(vector<int> cmn_sat, double *xa) {
        MatrixXd H,R;
        int i,j,n,nb=num_full_x_-num_real_x_fix_,nf=para_.GetGnssUsedFrqs();
        int index[MAX_SAT_NUM]={0};

        H=MatrixXd::Zero(nb,num_full_x_);
        vector<double>v;
        int m,f,nv=0;
        for(m=0;m<NSYS;m++){
            for(f=0;f<nf;f++){
                for(n=i=0;i<MAX_SAT_NUM;i++){
                    if(previous_sat_info_[i].sat.sat_.sys_idx!=m||previous_sat_info_[i].fix[f]!=2||
                    previous_sat_info_[i].el_az[0]*R2D<ppk_conf_.gnssC.hold_er_mask){
                        continue;
                    }
                    index[n++]=para_.IndexAmb(f,i+1);
                    previous_sat_info_[i].fix[f]=3;
                }

                for(i=1;i<n;i++){
                    v.push_back((xa[index[0]]-xa[index[i]])-(full_x_[index[0]]-full_x_[index[i]]));
                    H.data()[index[0]+nv*num_full_x_]=1.0;
                    H.data()[index[i]+nv*num_full_x_]=-1.0;
                    nv++;
                }
            }
        }

        if(nv<ppk_conf_.gnssC.min_sat_num2hold){
            return;
        }
        amb_hold_flag=true;

        R=MatrixXd::Zero(nv,nv);
        for(i=0;i<nv;i++) R(i,i)=0.1;

        VectorXd l;
        l=Map<VectorXd>(v.data(),nv);

        kf_.Adjustment(l,H,R,full_x_,full_Px_,nv,num_full_x_);

        if(ppk_conf_.gnssC.glo_ar_mode!=GLO_AR_FIXHOLD) return;

        double sum=0.0,dd=0.0;
        for(f=0;f<nf;f++){
            i=-1;

            for(j=nv=0;j<MAX_SAT_NUM;j++){
                if(previous_sat_info_[j].sat.sat_.sys_idx==SYS_GLO&&previous_sat_info_[j].vsat[f]&&previous_sat_info_[j].lock[f]>=0){
                    if(i<0){
                        i=j;
                        index[nv++]=j;
                    }
                    else{
                        int ia1=para_.IndexAmb(f,j+1);
                        int ia2=para_.IndexAmb(f,i+1);
                        dd=full_x_[ia1]-full_x_[ia2];
                        dd=0.01*(dd-Round(dd));
                        full_x_[ia1]-=dd;
                        previous_sat_info_[j].glo_bias[f]+=dd;
                        sum+=dd;
                        index[nv++]=j;
                    }
                }
            }
        }


    }

    cFusionSolver::cFusionSolver() {}

    cFusionSolver::cFusionSolver(tIGCOUPLEDConf C) {
        fs_conf_=C;
        num_full_x_=para_.GetIGCOUPLEDPar(C);
        full_x_=VectorXd::Zero(num_full_x_);
    }

    cFusionSolver::~cFusionSolver() {}

    bool cFusionSolver::InputImuData(int ws) {
        int i,week;

        if(imu_index_++<0||imu_index_>=imu_data_.data_.size()) return false;

        // prepare imu data for static detect
//        for(i=imu_index_;i>=0&&i<ws&&i<imu_data_.data_.size();i++){
//            imu_data_zd_.push_back(imu_data_.data_.at(i));
//        }

        cur_imu_info_.t_tag=imu_data_.data_[imu_index_].t_tag;

        if(fs_conf_.insC.imu_type==IMU_M39){ // 应该是：如果使用的IMU是M39类型的需要加上GPS周的时间
            if(fs_conf_.mode_opt==MODE_OPT_GSOF){
                gnss_sols_[0].t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
            }
            else if(fs_conf_.mode_opt>MODE_OPT_SOL){
                rover_obs_.GetGnssObs()[0].obs_time.Time2Gpst(&week, nullptr,SYS_GPS);
            }

            cur_imu_info_.t_tag+=week*604800.0;
        }

        cur_imu_info_.raw_gyro=imu_data_.data_[imu_index_].gyro;
        cur_imu_info_.raw_acce=imu_data_.data_[imu_index_].acce;
        return true;
    }
    //匹配基准站和移动站的观测值
    bool cFusionSolver::MatchGnssObs(bool isBack) {
        double sow1,sow2;
        int i,week=0,wod=0,info=false;

        for(i=rover_idx_-100<0?0:rover_idx_-10;rover_idx_<rover_obs_.GetGnssObs().size();i++){
            sow1=rover_obs_.GetGnssObs()[i].obs_time.Time2Gpst(&week,&wod,SYS_GPS);
            sow2=cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<0.5/fs_conf_.insC.sample_rate){
                rover_idx_=i;info=true;

                break;
            }
            else if((!isBack && (sow1-sow2)>1.0/fs_conf_.insC.sample_rate) ||
					(isBack && (sow2 - sow1) > 1.0 / fs_conf_.insC.sample_rate)){
                info=false;break;
            }
        }

        if(info&&fs_conf_.mode_opt==MODE_OPT_PPK){
            for(i=base_idx_-100<0?0:base_idx_-10;base_idx_<gnss_solver_->base_obs_.GetGnssObs().size();i++){
                sow1=gnss_solver_->base_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(&week,&wod,SYS_GPS);
                sow2=rover_obs_.GetGnssObs()[rover_idx_].obs_time.Time2Gpst(nullptr, nullptr,SYS_GPS);
                if(fabs(sow1-sow2)<DTTOL){
                    base_idx_=i;info=true;
                    break;
                }
                else if((!isBack && (sow1-sow2)>2.0*DTTOL) || 
						(isBack && (sow2 - sow1) > 2.0 * DTTOL)){
                    info=false;break;
                }
            }
        }

        return info;
    }
    //匹配GNSS位置信息的值
    bool cFusionSolver::MatchGnssSol(bool isBack) {
        double sow1,sow2;
        bool info=false;
        int i;

        for(i=gnss_sol_idx>3?gnss_sol_idx-3:0;i<gnss_sols_.size()&&i>=0;i++){
            sow1=gnss_sols_[i].t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            sow2=cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<0.5/fs_conf_.insC.sample_rate){
                info=true;gnss_sol_idx=i;
                break;
            }
            else if((!isBack && (sow1-sow2)>1.0/fs_conf_.insC.sample_rate)|| 
					(isBack && (sow2 - sow1) > 1.0 / fs_conf_.insC.sample_rate)){
                info=false;break;
            }
        }
        return info;
    }

    //将IMU解转换成最终解
    void cFusionSolver::InsSol2IgcoupledSol(tImuInfoUnit &imu_sol, tSolInfoUnit &igcoupled_sol) {
        igcoupled_sol_.t_tag=imu_sol.t_tag;
        igcoupled_sol.pos=imu_sol.re;
        igcoupled_sol.vel=imu_sol.ve;
        igcoupled_sol.att=imu_sol.rpy;

        igcoupled_sol.accl_bias=imu_sol.ba;
        igcoupled_sol.gyro_bias=imu_sol.bg;

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        Matrix3d Pp,Pv,Pa;
        //将权矩阵从Px中取出
        if(fs_conf_.mode == MODE_IGTC){
            Pp=gnss_solver_->full_Px_.block<3,3>(ip,ip);
            Pv=gnss_solver_->full_Px_.block<3,3>(iv,iv);
            Pa=gnss_solver_->full_Px_.block<3,3>(ia,ia);
        }else if(fs_conf_.mode == MODE_IGLC){
			Pp = full_Px_.block<3,3>(ip,ip);
			Pv = full_Px_.block<3,3>(iv,iv);
			Pa = full_Px_.block<3,3>(ia,ia);
		}

        //取出位置、速度和姿态的精度信息
        for(int i=0;i<3;i++) igcoupled_sol.q_pos[i]=Pp(i,i);
        igcoupled_sol.q_pos[3]=Pp(0,1);
        igcoupled_sol.q_pos[4]=Pp(1,2);
        igcoupled_sol.q_pos[5]=Pp(2,0);

        for(int i=0;i<3;i++) igcoupled_sol.q_vel[i]=Pv(i,i);
        igcoupled_sol.q_vel[3]=Pv(0,1);
        igcoupled_sol.q_vel[4]=Pv(1,2);
        igcoupled_sol.q_vel[5]=Pv(2,0);

        for(int i=0;i<3;i++) igcoupled_sol.q_att[i]=Pa(i,i);
        igcoupled_sol.q_att[3]=Pa(0,1);
        igcoupled_sol.q_att[4]=Pa(1,2);
        igcoupled_sol.q_att[5]=Pa(2,0);
    }
    //由速度计算航向角yaw（速度矢量粗对准）
    double cFusionSolver::Vel2Yaw(Vector3d vn) {
        return atan2(vn[1],fabs(vn[0])<1E-4?1E-4:vn[0]);
    }
    //位置之差除以时间间隔计算平均速度
    Vector3d cFusionSolver::Pos2Vel(tSolInfoUnit &sol1, tSolInfoUnit &sol2) {
        Vector3d vel={0,0,0};
        double t=sol1.t_tag.TimeDiff(sol2.t_tag.t_);

        vel=(sol1.pos-sol2.pos)/t;
        return vel;
    }

    //将GNSS位置和速度信息更新为前一时刻的INS的位置和速度以及RPY
    bool cFusionSolver::GnssSol2Ins(Vector3d re,Vector3d ve) {
        if(re.norm()==0.0||ve.norm()==0.0) return false;
        Vector3d blh=Xyz2Blh(re),rn;
        Vector3d wiee(0,0,OMGE_GPS);
        Matrix3d Cne;

		if(cur_imu_info_.t_tag.TimeDiff(pre_imu_info_.t_tag.t_) < 0) wiee = {0,0,-OMGE_GPS};
        Cne=CalcCen(blh,COORD_NED).transpose();
        pre_imu_info_.rn=blh;
        pre_imu_info_.vn=Cne.transpose()*ve;

        pre_imu_info_.rpy[2]=Vel2Yaw(pre_imu_info_.vn);
        pre_imu_info_.Cbn=Euler2RotationMatrix(pre_imu_info_.rpy);
        pre_imu_info_.Cbe=Cne*pre_imu_info_.Cbn;

        // lever correct
        pre_imu_info_.re=re-pre_imu_info_.Cbe*fs_conf_.insC.lever;
        Matrix3d T=VectorSkew(cur_imu_info_.raw_gyro*cur_imu_info_.dt);//当前时刻初始角增量的旋转对称
        Matrix3d Omge=VectorSkew(wiee); //地球自传参数的旋转对称
        pre_imu_info_.ve=ve-pre_imu_info_.Cbe*T*fs_conf_.insC.lever+Omge*pre_imu_info_.Cbe*fs_conf_.insC.lever;//对速度进行杆壁改正
        pre_imu_info_.rn=Xyz2Blh(pre_imu_info_.re);//更新改正后的导航系下的位置和速度
        pre_imu_info_.vn=Cne.transpose()*pre_imu_info_.ve;

        return true;
    }

    //INS粗对准
    bool cFusionSolver::InsAlign() {
        int num_sols=5;
        static vector<tSolInfoUnit>sols;
        SOL_STAT  align_level=SOL_FIX;
        if(fs_conf_.mode_opt==MODE_OPT_PPP){
            align_level=SOL_PPP;
        }
        else if(fs_conf_.mode_opt==MODE_OPT_SPP){
            align_level=SOL_SPP;
        }

        if(fs_conf_.insC.ins_align==ALIGN_GNSS_OBS){

            if(gnss_alignor_->SolverProcess(gnss_conf_,rover_idx_)) {
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"USING GNSS TO ALIGN INS("<<sols.size()<<"): "<<gnss_alignor_->igcoupled_sol_.pos.transpose();
                if(gnss_alignor_->igcoupled_sol_.stat==align_level){
                    sols.push_back(gnss_alignor_->igcoupled_sol_);
                }
            }
            else{
                sols.clear();
                return false;
            }

        }
        else if(fs_conf_.insC.ins_align==ALIGN_GNSS_SOL){
            sols.push_back(gnss_sols_.at(gnss_sol_idx));
        }

        if(sols.size()>=num_sols){
            if(fs_conf_.insC.ins_align==ALIGN_GNSS_OBS){
                for(int i=0;i<sols.size()-1;i++){
                    if(sols[i+1].t_tag.TimeDiff(sols[i].t_tag.t_)>fs_conf_.gnssC.sample_rate){
                        sols.clear();
                        return false;
                    }
                }
            }
            else{
                for(int i=0;i<sols.size();i++){
                    //if(sols[i].stat>SOL_FIX||sols[i].stat==SOL_NONE){
                      //  sols.clear();
                        //return false;
                    //}
                }
            }

            Vector3d ve;
            if(sols.back().vel.norm()==0.0){
                ve=Pos2Vel(sols.back(),*(sols.end()-2));
            }
            else ve=sols.back().vel;

            if(ve.norm()<5.0||cur_imu_info_.raw_gyro.norm()>30.0*D2R){
                sols.clear();
                return false;
            }


            if(GnssSol2Ins(sols.back().pos,ve)){
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" INS INITIALIZATION OK";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pre_imu_info_.re.transpose()<<" m";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT VELOCITY(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.vn.transpose()<<" m/s";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT ATTITUTE(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.rpy.transpose()*R2D<<" deg";
                igcoupled_sol_.ins_stat=SOL_INS_INIT;
                InsSol2IgcoupledSol(pre_imu_info_,igcoupled_sol_);
                return true;
            }
            else sols.clear();
        }
        return false;
    }
    //初始化求解模式、读取IMU数据
    void cFusionSolver::InitSolver(tIGCOUPLEDConf C) {
        fs_conf_=C;
        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
		if(C.filter_type != FILTER_COMBINED){
			out_=new cOutSol(C);
            out_->InitOutSol(C,C.fileC.sol);
            out_->WriteHead();
		}
        //用初始unc初始化Px
        InitInsPx(C,num_full_x_,full_Px_);

        if(fs_conf_.mode==MODE_IGTC) tc_mode_=true;

        if(fs_conf_.mode_opt==MODE_OPT_SPP){
            tIGCOUPLEDConf spp_conf=C;
            if(C.mode==MODE_IGLC){
                spp_conf.mode=MODE_SPP;
            }
            spp_conf.solC.out_sol=false;
			spp_conf.filter_type = FILTER_FORWARD;
            gnss_solver_=new cSppSolver(spp_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(spp_conf);
            gnss_conf_=spp_conf;

            if(C.insC.ins_align==ALIGN_GNSS_OBS){
                spp_conf.mode=MODE_SPP;
                spp_conf.mode_opt=MODE_OPT_KINEMATIC;
                gnss_alignor_=new cSppSolver(spp_conf);
                gnss_alignor_->nav_=nav_;
                gnss_alignor_->InitSolver(spp_conf);
                gnss_alignor_->rover_obs_=rover_obs_;
            }
        }
        else if(fs_conf_.mode_opt==MODE_OPT_PPP){
            tIGCOUPLEDConf ppp_conf=C;
            if(C.mode==MODE_IGLC) ppp_conf.mode=MODE_PPP;
            ppp_conf.solC.out_sol=false;
            ppp_conf.solC.out_head=false;
			ppp_conf.filter_type = FILTER_FORWARD;
            gnss_solver_=new cPppSolver(ppp_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(ppp_conf);
            gnss_conf_=ppp_conf;

            if(C.insC.ins_align==ALIGN_GNSS_OBS){
                ppp_conf.mode=MODE_PPP;
                ppp_conf.mode_opt=MODE_OPT_KINEMATIC;
                gnss_alignor_=new cPppSolver(ppp_conf);
                gnss_alignor_->nav_=nav_;
                gnss_alignor_->InitSolver(ppp_conf);
                gnss_alignor_->rover_obs_=rover_obs_;

            }
        }
        else if(fs_conf_.mode_opt==MODE_OPT_PPK){
            tIGCOUPLEDConf ppk_conf=C;
            if(C.mode==MODE_IGLC) ppk_conf.mode=MODE_PPK;
            ppk_conf.solC.out_sol=false;
            ppk_conf.solC.out_head=false;
			ppk_conf.filter_type = FILTER_FORWARD;
            gnss_solver_=new cPpkSolver(ppk_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(ppk_conf);
            gnss_conf_=ppk_conf;

            if(C.insC.ins_align==ALIGN_GNSS_OBS){
                ppk_conf.mode=MODE_PPK;
                ppk_conf.mode_opt=MODE_OPT_KINEMATIC;
                gnss_alignor_=new cPpkSolver(ppk_conf);
                gnss_alignor_->nav_=nav_;
                gnss_alignor_->InitSolver(ppk_conf);
                gnss_alignor_->rover_obs_=rover_obs_;
            }
        }
        else if(fs_conf_.mode_opt==MODE_OPT_GSOF){
            cDecodeGsof gsof_decoder;
            gsof_decoder.DecodeGsof(C.fileC.gsof,gnss_sols_);
			fs_conf_.gnssC.sample_rate = fabs(gnss_sols_[0].t_tag.TimeDiff(gnss_sols_[1].t_tag.t_));
        }
        else if(fs_conf_.mode_opt==MODE_OPT_SOL){
            ReadSol(C,C.fileC.gsof,gnss_sols_);
			fs_conf_.gnssC.sample_rate = fabs(gnss_sols_[0].t_tag.TimeDiff(gnss_sols_[1].t_tag.t_));
        }

        if(C.insC.imu_type==IMU_M39){
            cDecodeImuM39 m39_decoder;
            m39_decoder.DecodeM39(C.fileC.imu,C.insC,imu_data_.data_);
        }
        else if(C.insC.imu_type==IMU_NOVTEL_A1||C.insC.imu_type==IMU_NOVTEL_CPT||C.insC.imu_type==IMU_POS){
            int week=0;
            if(C.insC.imu_type==IMU_POS && !rover_obs_.GetGnssObs().empty()) rover_obs_.GetGnssObs()[0].obs_time.Time2Gpst(&week, nullptr,SYS_GPS);
			else if(C.insC.imu_type >= IMU_POS && !gnss_sols_.empty()) gnss_sols_[0].t_tag.Time2Gpst(&week, nullptr, SYS_GPS);
            cReadImu imu_reader(C.fileC.imu,week);
            imu_reader.SetImu(C);
            imu_reader.Reading();
            imu_data_=*imu_reader.GetImus();
        }

		if(C.filter_type == FILTER_BACKWARD){
			for(int i = imu_data_.data_.size() - 1; i >= 0; i--){
				imu_data_.data_[i].gyro *= -1;
				imu_data_.data_[i].t_tag = i == 0? (imu_data_.data_[i].t_tag + -1 / fs_conf_.insC.sample_rate) : imu_data_.data_[i-1].t_tag;
			}
			std::reverse(imu_data_.data_.begin(), imu_data_.data_.end());
			std::reverse(gnss_sols_.begin(), gnss_sols_.end());
			if(C.mode_opt == MODE_OPT_SOL || C.mode_opt == MODE_OPT_GSOF){}
			else{
				rover_obs_.ReverseObs();
				gnss_alignor_->rover_obs_ = rover_obs_;
			}
		}

        if(fs_conf_.mode_opt>MODE_OPT_SOL) gnss_solver_->rover_obs_=rover_obs_;
        if(tc_mode_) gnss_solver_->full_Px_=full_Px_;
    }

	bool cFusionSolver::SolverProcess(tIGCOUPLEDConf C, int idx){
		if(C.filter_type == FILTER_FORWARD || C.filter_type == FILTER_BACKWARD){
			CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<(C.filter_type == FILTER_FORWARD?"FORWARD":("BACKWARD"))<<"FILTER BEGIN"<<endl;
			process(C, 0);
		}else if(C.filter_type == FILTER_COMBINED){
			solb_.clear();
			solf_.clear();
			CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"COMBINED / FORWARD FILTER BEGIN"<<endl;
			process(C, 0);
			imu_index_ = 0;
			gnss_sol_idx = 0;
			rover_idx_ = 0;
			ins_mech_idx = 0;
			full_x_.setZero();
			full_Px_.setZero();
			imu_data_.data_.clear();
			gnss_sols_.clear();
			CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID) << "COMBINED / BACKWARD FILTER BEGIN"<<endl;
			process(C, 1);
			CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID) << "COMBINED F/B SOLUTIONS"<<endl;
			CombFbSol(C);
			solf_.clear();
			solb_.clear();

		}
	}

    bool cFusionSolver::process(tIGCOUPLEDConf C,int idx) {

		if(C.filter_type == FILTER_COMBINED && idx == 1){
			tIGCOUPLEDConf tmpC = C;
			tmpC.filter_type = FILTER_BACKWARD;
			InitSolver(tmpC);
			fs_conf_ = C;
		}else{
			InitSolver(C);
		}
        tSatInfoUnit sat_info;

        int gnss_obs_flag=false,ins_align=false,gnss_sol_flag=false;

		double outage_time = 0;
		if(fs_conf_.gnssC.use_outage) outage_time = fs_conf_.gnssC.outage_time;

		//反向需调整outage time
		if((C.filter_type == FILTER_BACKWARD || idx == 1) && fs_conf_.gnssC.use_outage){
			//找到IMU和GNSS中sow比较小的那个
			double sow1 = imu_data_.data_[0].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
			double sow2;
			if(fs_conf_.mode_opt == MODE_OPT_SOL || fs_conf_.mode_opt == MODE_OPT_GSOF){
				sow2 = gnss_sols_[0].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
			}else{
				sow2 = rover_obs_.GetGnssObs()[0].obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS);
			}
			sow2 = MIN(sow2, sow1);
			while(outage_time + fs_conf_.gnssC.outage_len + fs_conf_.gnssC.outage_period < sow2){
				outage_time += fs_conf_.gnssC.outage_period;
				if(fs_conf_.gnssC.outage_period <= fs_conf_.gnssC.outage_len) break;
			}
			outage_time += fs_conf_.gnssC.outage_len;
			fs_conf_.gnssC.outage_time = outage_time;
		}
        while(InputImuData(5)){
            if(fs_conf_.mode>MODE_INS){
                if(fs_conf_.mode_opt==MODE_OPT_GSOF||fs_conf_.mode_opt==MODE_OPT_SOL){
				 	if(fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) gnss_sol_flag = MatchGnssSol(true);
					else gnss_sol_flag = MatchGnssSol(false);
					if(gnss_sol_flag && outage_time != 0){
						double sow_sol = gnss_sols_[gnss_sol_idx].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
						if(((fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) && sow_sol < (outage_time - fs_conf_.gnssC.outage_len) && 
								fs_conf_.gnssC.outage_period > fabs(fs_conf_.gnssC.outage_len)) || 
								((fs_conf_.filter_type == FILTER_FORWARD || (fs_conf_.filter_type == FILTER_COMBINED && idx == 0)) && 
								sow_sol >(outage_time + fs_conf_.gnssC.outage_len) && fs_conf_.gnssC.outage_period > fs_conf_.gnssC.outage_len)){
									outage_time = fs_conf_.gnssC.outage_time + static_cast<int>((sow_sol - fs_conf_.gnssC.outage_time)/fs_conf_.gnssC.outage_period) * fs_conf_.gnssC.outage_period;
								}

						if(((fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) && sow_sol <= outage_time && sow_sol >= outage_time - fs_conf_.gnssC.outage_len)||
							((fs_conf_.filter_type == FILTER_FORWARD || (fs_conf_.filter_type == FILTER_COMBINED && idx == 0)) && 
							sow_sol >= outage_time + fs_conf_.gnssC.outage_len)){
								gnss_sol_flag = false;
								CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS OUTAGE AT"<<sow_sol<<"(s)";
							}
					}
                }
                else{
					if(fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) gnss_obs_flag = MatchGnssObs(true);
                    else gnss_sol_flag = MatchGnssObs(false);
                    if(gnss_sol_flag && outage_time != 0){
                        double sow_sol = gnss_sols_[gnss_sol_idx].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
                        if(((fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) && sow_sol < (outage_time - fs_conf_.gnssC.outage_len) && 
                                fs_conf_.gnssC.outage_period > fabs(fs_conf_.gnssC.outage_len)) || 
                                ((fs_conf_.filter_type == FILTER_FORWARD || (fs_conf_.filter_type == FILTER_COMBINED && idx == 0)) && 
                                sow_sol >(outage_time + fs_conf_.gnssC.outage_len) && fs_conf_.gnssC.outage_period > fs_conf_.gnssC.outage_len)){
                                    outage_time = fs_conf_.gnssC.outage_time + static_cast<int>((sow_sol - fs_conf_.gnssC.outage_time)/fs_conf_.gnssC.outage_period) * fs_conf_.gnssC.outage_period;
                                  
                        if(((fs_conf_.filter_type == FILTER_BACKWARD || idx == 1) && sow_sol <= outage_time && sow_sol >= outage_time - fs_conf_.gnssC.outage_len)||
                            ((fs_conf_.filter_type == FILTER_FORWARD || (fs_conf_.filter_type == FILTER_COMBINED && idx == 0)) && 
                            sow_sol >= outage_time + fs_conf_.gnssC.outage_len)){
                                gnss_obs_flag = false;
                                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS OUTAGE AT"<<sow_sol<<"(s)";
                            }
                       }
                   }
                }
                //如果观测值或者GNSS位置解匹配成功
                if(gnss_sol_flag||gnss_obs_flag){
                    //如果没有对准或上次对准失败
                    if(!ins_align){
                        //如果使用GNSS观测值
                        if(gnss_obs_flag) ins_align=InsAlign();
                        else {
                            //如果位置的模为0，即位置不存在
                            if(gnss_sols_[gnss_sol_idx].pos.norm()==0.0) continue;

                            // initial alignment
                            ins_align=InsAlign();
                        }
                        pre_imu_info_.t_tag=cur_imu_info_.t_tag;
                        pre_imu_info_.cor_gyro=pre_imu_info_.raw_gyro=cur_imu_info_.raw_gyro;
                        pre_imu_info_.cor_acce=pre_imu_info_.raw_acce=cur_imu_info_.raw_acce;
                        epoch_sat_obs_.epoch_data.clear();
                        base_epoch_sat_obs_.epoch_data.clear();
                        continue;
                    }
                    //如果对准成功
                    // measurement update
                    if(C.mode==MODE_IGLC){
                        LooseCouple(C);
                    }
                    else if(C.mode==MODE_IGTC) TightCouple(C);
                }
                else{ // 位置或观测值没有匹配成功
                    if(!ins_align){
                        pre_imu_info_=cur_imu_info_;
                        continue;
                    }
                    // time update
                    igcoupled_sol_.stat=SOL_NONE;
                    if(ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
                        igcoupled_sol_.ins_stat=SOL_INS_MECH;
                        InsSol2IgcoupledSol(cur_imu_info_,igcoupled_sol_);//用INS解代替最终解
                    }
                    StateTimeUpdate();
                }
            }
            else{
                // pure ins mech
                igcoupled_sol_.stat=SOL_NONE;
                if(ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
                    igcoupled_sol_.ins_stat=SOL_INS_MECH;
                    InsSol2IgcoupledSol(cur_imu_info_,igcoupled_sol_);
                }
            }

            pre_imu_info_=cur_imu_info_;

			if(C.filter_type == FILTER_COMBINED){
				if(C.solC.out_ins_mech_frq != 0 && ins_mech_idx%C.solC.out_ins_mech_frq==0&&igcoupled_sol_.ins_stat==SOL_INS_MECH){
					igcoupled_sol_.valid_sat_num = 0;
					if(idx == 0){
						solf_.push_back(igcoupled_sol_);
					}
					else if(idx == 1){
						solb_.push_back(igcoupled_sol_);
					}
				}
				else if(igcoupled_sol_.ins_stat != SOL_INS_MECH){
					if(idx == 0){
						solf_.push_back(igcoupled_sol_);
					}
					else if(idx == 1){
						solb_.push_back(igcoupled_sol_);
					}
				}
			}else{
            	if(C.solC.out_ins_mech_frq!=0&&ins_mech_idx%C.solC.out_ins_mech_frq==0&&igcoupled_sol_.ins_stat==SOL_INS_MECH){
            	    igcoupled_sol_.valid_sat_num=0;
                	out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
           	 	}
            	else if(igcoupled_sol_.ins_stat!=SOL_INS_MECH) out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
				}
            if(C.solC.out_ins_mech_frq!=0&&ins_mech_idx%C.solC.out_ins_mech_frq==0&&igcoupled_sol_.ins_stat==SOL_INS_MECH){
                igcoupled_sol_.valid_sat_num=0;
                out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
            }
            else if(igcoupled_sol_.ins_stat!=SOL_INS_MECH) out_->WriteSol(igcoupled_sol_,epoch_sat_info_collect_);
        }
    }

    bool cFusionSolver::SolverEpoch() {
        gnss_solver_->epoch_idx_=epoch_idx_;
        gnss_solver_->epoch_sat_obs_=epoch_sat_obs_;
        if(gnss_solver_->SolverEpoch()){
            gnss_solver_->out_->WriteSol(gnss_solver_->igcoupled_sol_,epoch_sat_info_collect_);
        }
    }

    bool cFusionSolver::SolutionUpdate() {
        // ????
    }

    void cFusionSolver::DisableX(VectorXd &x, int idx) {
        if(idx>para_.NumClPar()) return;
        else{
            x[idx]=DIS_FLAG;
        }
    }

    void cFusionSolver::StateTimeUpdate() {
        double dt=fabs(cur_imu_info_.dt=cur_imu_info_.t_tag.TimeDiff(pre_imu_info_.t_tag.t_));
        int nx=para_.GetInsTransParNum(fs_conf_);//所有INS相关参数（位置速度姿态零偏比例因子非正交）
        MatrixXd F,Q;
        F=MatrixXd::Zero(nx,nx);

        //状态转移矩阵
        F=ins_mech_.StateTransferMat(fs_conf_,pre_imu_info_,cur_imu_info_,nx,dt);

//        cout<<F<<endl<<endl;

#if 0
        Q=InitQ(fs_conf_,dt,nx);
#else
        Q=InitPrecQ(fs_conf_,dt,nx,cur_imu_info_.Cbe);
#endif
        VectorXd x;
        MatrixXd Px=MatrixXd::Zero(nx,nx);
        int i,j;
        if(tc_mode_){//紧组合
            x=Map<VectorXd>(gnss_solver_->full_x_.data(),nx);

            for(i=0;i<nx;i++) for(j=0;j<nx;j++){
                Px.data()[i+j*nx]=gnss_solver_->full_Px_.data()[i+j*num_full_x_];
            }
//            Px=gnss_solver_->full_Px_.block<15,15>(0,0);
        }
        else{//松组合
            x=Map<VectorXd>(full_x_.data(),nx);
            for(i=0;i<nx;i++) for(j=0;j<nx;j++){
                Px.data()[i+j*nx]=full_Px_.data()[i+j*num_full_x_];
            }
//            Px=full_Px_.block<15,15>(0,0);
        }

        if(fabs(dt)>60.0){
            InitInsPx(fs_conf_,nx,Px);//用配置文件初始化权矩阵Px
        }
        else{
            PropVariance(F,Q,nx,Px);
        }

        for(int i=0;i<nx;i++) x[i]=1E-20;
        if(tc_mode_){
            gnss_solver_->full_x_.segment(0,nx)=Map<VectorXd>(x.data(),nx);

            for(i=0;i<nx;i++) for(j=0;j<nx;j++){
                gnss_solver_->full_Px_.data()[i+j*num_full_x_]=Px.data()[i+j*nx];
            }
//            gnss_solver_->full_Px_.block<15,15>(0,0)=Px;
        }
        else{
            full_x_.segment(0,nx)=Map<VectorXd>(x.data(),nx);
            for(i=0;i<nx;i++) for(j=0;j<nx;j++){
                full_Px_.data()[i+j*num_full_x_]=Px.data()[i+j*nx];
            }
//            full_Px_.block<15,15>(0,0)=Px;
        }
    }
    // PQ = Px + 1/2Q
    // FPF = F（PQ）F^T
    // Px = FPF + 1/2 Q
    void cFusionSolver::PropVariance(MatrixXd& F,MatrixXd& Q,int nx,MatrixXd& Px) {

        MatrixXd prior_P=Px;
        MatrixXd PQ(nx,nx),FPF(nx,nx);

        PQ=MatrixXd::Zero(nx,nx);FPF=MatrixXd::Zero(nx,nx);
        int i,j;
//        PQ.block<15,15>(0,0)=prior_P.block<15,15>(0,0)+0.5*Q.block<15,15>(0,0);

        for(i=0;i<nx;i++){
            for(j=0;j<nx;j++) PQ.data()[i+j*nx]=prior_P.data()[i+j*nx]+0.5*Q.data()[i+j*nx];
        }

        FPF=F*PQ*F.transpose();

#if 0
        cout<<std::fixed<<setprecision(12)<<prior_P<<endl<<endl;

        cout<<std::fixed<<setprecision(12)<<Q<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<PQ<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<F<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<FPF<<endl<<endl;
#endif

//        Px.block<15,15>(0,0)=FPF+0.5*Q;
        for(i=0;i<nx;i++){
            for(j=0;j<nx;j++) Px.data()[i+j*nx]=FPF.data()[i+j*nx]+0.5*Q.data()[i+j*nx];
        }
    }

    int cFusionSolver::BuildLcHVR(int post,tIGCOUPLEDConf C,tImuInfoUnit& imu_info,double *meas_pos,double *meas_vel,Vector3d& q_pos,Vector3d& q_vel) {
        Vector3d gnss_re,gnss_ve;
        Matrix3d Cbe=imu_info.Cbe;
        double omc=0.0;
        vector<double>omcs;
        num_L_=0;

        if(!meas_pos&&meas_vel) return 0;
        if(meas_pos) num_L_+=3;
        if(meas_vel) num_L_+=3;
        H_=MatrixXd::Zero(num_L_,num_full_x_);
        R_=MatrixXd::Zero(num_L_,num_L_);

        RemoveLever(imu_info,C.insC.lever,gnss_re,gnss_ve);

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        if(meas_pos){
            for(int i=0;i<3;i++){
                omc=meas_pos[i]-gnss_re[i];//测量值-观测值=残差
                omcs.push_back(omc);
            }

            if(!post){
                // position to position jacobian
                H_.block<3,3>(ip,ip)=-1.0*Matrix3d::Identity();

                // position to attitude jacobian
                H_.block<3,3>(ip,ia)=VectorSkew(Cbe*C.insC.lever);
            }
            R_.block<3,3>(ip,ip)=q_pos.asDiagonal();
        }

        if(meas_vel){
            for(int i=0;i<3;i++){
                omc=meas_vel[i]-gnss_ve[i];
                omcs.push_back(omc);
            }

            if(!post){
                // velocity to velocity jacobian
                H_.block<3,3>(iv,iv)=-1.0*Matrix3d::Identity();

                // velocity to attitude jacobian
                Vector3d wiee(0.0,0.0,OMGE_GPS);
                Vector3d T,W;
                T=Cbe*VectorSkew(imu_info.cor_gyro)*C.insC.lever;
                W=VectorSkew(wiee)*Cbe*C.insC.lever;
                H_.block<3,3>(iv,ia)=VectorSkew(T-W);

                // velocity to ba jacobian
                H_.block<3,3>(iv,iba)=Cbe*(VectorSkew(C.insC.lever));

            }
            R_.block<3,3>(3,3)=q_vel.asDiagonal();
        }

        omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        omcs.clear();

        return num_L_;
    }

    bool cFusionSolver::LcFilter(tIGCOUPLEDConf C) {
        double *meas_pos= nullptr,*meas_vel= nullptr;
        Vector3d q_pos,q_vel;

        if(C.mode_opt>MODE_OPT_SOL){//当用GNSS位置进行组合时，测量值为GNSS位置
            igcoupled_sol_.stat=gnss_solver_->igcoupled_sol_.stat;
            igcoupled_sol_.pos=gnss_solver_->igcoupled_sol_.pos;
            igcoupled_sol_.vel=gnss_solver_->igcoupled_sol_.vel;
            if(igcoupled_sol_.pos.norm()!=0.0){
                meas_pos=igcoupled_sol_.pos.data();
                q_pos<<igcoupled_sol_.q_pos[0],igcoupled_sol_.q_pos[1],igcoupled_sol_.q_pos[2];
            }
            if(igcoupled_sol_.vel.norm()!=0.0){
                meas_vel=igcoupled_sol_.vel.data();
                q_vel<<igcoupled_sol_.q_vel[0],igcoupled_sol_.q_vel[1],igcoupled_sol_.q_vel[2];
                if(q_vel.norm()==0.0){
                    q_vel<<0.1,0.1,0.1;
                }
            }
        }
        else{//没有GNSS位置组合时
            igcoupled_sol_.stat=gnss_sols_[gnss_sol_idx].stat;
            igcoupled_sol_.valid_sat_num=gnss_sols_[gnss_sol_idx].valid_sat_num;
            if(gnss_sols_[gnss_sol_idx].pos.norm()!=0.0){
                meas_pos=gnss_sols_[gnss_sol_idx].pos.data();
                q_pos<<(gnss_sols_[gnss_sol_idx].q_pos[0]),(gnss_sols_[gnss_sol_idx].q_pos[1]),(gnss_sols_[gnss_sol_idx].q_pos[2]);
            }
            if(gnss_sols_[gnss_sol_idx].vel.norm()!=0.0){
                meas_vel=gnss_sols_[gnss_sol_idx].vel.data();
                q_vel<<(gnss_sols_[gnss_sol_idx].q_vel[0]),(gnss_sols_[gnss_sol_idx].q_vel[1]),(gnss_sols_[gnss_sol_idx].q_vel[2]);
                if(q_vel.norm()==0.0){
                    q_vel<<0.1,0.1,0.1;
                }
            }
        }

        Vector3d ins_re=cur_imu_info_.re,ins_ve=cur_imu_info_.ve;

        Vector3d pos,vel;
        for(int i=0;i<3;i++){
            if(meas_pos) pos[i]=meas_pos[i];
            if(meas_vel) vel[i]=meas_vel[i];
        }
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS GNSS   FUSION(-): "<<cur_imu_info_.t_tag.GetTimeStr(3);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS MECH POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.re.transpose();
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS MECH VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.ve.transpose();
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS     POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pos.transpose()<<" "<<q_pos.transpose();
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS     VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<vel.transpose()<<" "<<q_vel.transpose();

        VectorXd x=full_x_;
        MatrixXd Px=full_Px_;
        tImuInfoUnit imu_corr=cur_imu_info_;
        bool stat=false;
        if(BuildLcHVR(0,C,cur_imu_info_,meas_pos,meas_vel,q_pos,q_vel)){

            if(R_.trace()>=10.0){
                for(int i=0;i<para_.IndexAtt()+para_.NumAtt();i++){
                    DisableX(x,i);
                }
                for(int i=0;i<para_.IndexBa()+para_.NumBa();i++){
                    DisableX(x,i);
                }
                for(int i=0;i<para_.IndexBg()+para_.NumBg();i++){
                    DisableX(x,i);
                }
            }

#if 0
            cout<<omc_L_.transpose()<<endl<<endl;
            cout<<R_<<endl<<endl;
            cout<<H_.transpose()<<endl;
            cout<<full_Px_<<endl;
#endif
            kf_.Adjustment(omc_L_,H_.transpose(),R_,x,Px,num_L_,num_full_x_);

            CloseLoopState(x,&imu_corr);

            if(BuildLcHVR(1,C,imu_corr,meas_pos,meas_vel,q_pos,q_vel)){

                stat=ValidSol(x, 10.0);
                if(stat){
                    full_x_=x;
                    full_Px_=Px;
                    cur_imu_info_=imu_corr;
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS GNSS   FUSION(+): "<<cur_imu_info_.t_tag.GetTimeStr(3);
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS MECH POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.re.transpose();
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS MECH VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.ve.transpose();
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS     POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pos.transpose()<<" "<<q_pos.transpose();
                    CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS     VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<vel.transpose()<<" "<<q_vel.transpose();
                }
            }
        }
        else{
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" MAKE OMC ERROR";
            stat=false;
        }

        return stat;
    }

    void cFusionSolver::ControlPx(int nx, MatrixXd &Px) {
        double var=0.0;
        for(int i=0;i<3;i++) {
            var+=SQRT(Px(i,i));
        }
        if((var/3)>SQR(100)) InitInsPx(fs_conf_,nx,Px);
    }

    bool cFusionSolver::ValidSol(VectorXd& x, double thres) {
        int flag=0;
        double fact=SQR(thres);
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();
        Vector3d att(x[ia],x[ia+1],x[ia+2]);
        Vector3d ba(x[iba],x[iba+1],x[iba+2]);
        Vector3d bg(x[ibg],x[ibg+1],x[ibg+2]);
        if(x[ia]==DIS_FLAG&&att.norm()>5.0*D2R){
            flag|=1;
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE ATTITUDE AND OVERRUN";
        }
        if(x[iba]==DIS_FLAG&&ba.norm()>1E5*MG2MS2){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE ACCE BIAS AND OVERRUN";
            flag|=1;
        }
        if(x[ibg]==DIS_FLAG&&bg.norm()>10.0*D2R){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE GYRO BIAS AND OVERRUN";
            flag|=1;
        }

        if(flag){
            return false;
        }

        int i,j;
        for(j=0,i=0;i<num_L_&&thres;i++){
            if(omc_L_[i]*omc_L_[i]<fact*R_.data()[i+i*num_L_]) continue;
            j++;
        }
        return j<num_L_;
    }

    bool cFusionSolver::LooseCouple(tIGCOUPLEDConf C) {
        igcoupled_sol_.stat=SOL_NONE;
        igcoupled_sol_.valid_sat_num=0;
        if(!ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
            return false;
        }
        StateTimeUpdate();

        igcoupled_sol_.ins_stat=SOL_INS_MECH;
        epoch_idx_++;

        if(fs_conf_.mode_opt>MODE_OPT_SOL){
//            cout<<gnss_solver_->full_x_.transpose()<<endl;
            if(gnss_solver_->SolverProcess(gnss_conf_,rover_idx_)){
                igcoupled_sol_=gnss_solver_->igcoupled_sol_;
                if(LcFilter(C)){
                    igcoupled_sol_.ins_stat=SOL_IG_LC;
                    InsSol2IgcoupledSol(cur_imu_info_,igcoupled_sol_);
                }
            }
        }
        else{
            if(LcFilter(C)){
                igcoupled_sol_.ins_stat=SOL_IG_LC;
                InsSol2IgcoupledSol(cur_imu_info_,igcoupled_sol_);
            }
        }
    }

    bool cFusionSolver::TightCouple(tIGCOUPLEDConf C) {
        StateTimeUpdate();

        if(!ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
            return false;
        }

        igcoupled_sol_.stat=SOL_NONE;
        igcoupled_sol_.ins_stat=SOL_INS_MECH;
        igcoupled_sol_.valid_sat_num=0;

//        cout<<gnss_solver_->full_Px_.block<15,15>(0,0)<<endl<<endl;
//        cout<<gnss_solver_->full_x_.transpose()<<endl<<endl;

        epoch_idx_++;
        gnss_solver_->cur_imu_info_=cur_imu_info_;
        gnss_solver_->tc_mode_= true;
        if(!gnss_solver_->SolverProcess(fs_conf_,rover_idx_)){
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<"WARNING";
        }
        cur_imu_info_=gnss_solver_->cur_imu_info_;
        InsSol2IgcoupledSol(cur_imu_info_,igcoupled_sol_);
        igcoupled_sol_.stat=gnss_solver_->igcoupled_sol_.stat;
        igcoupled_sol_.ins_stat=SOL_IG_TC;
    }

	cFGOSolver::cFGOSolver() {}
    cFGOSolver::cFGOSolver(tIGCOUPLEDConf C) {
        fgo_conf_=C;
        imu_idx_ = 0;
        gnss_sol_idx_ = 0;
        rover_idx_ = 0;
    }
    cFGOSolver::~cFGOSolver() {}


    Vector3d cFGOSolver::Pos2Vel(tSolInfoUnit &sol1, tSolInfoUnit &sol2) {
        Vector3d vel = {0,0,0};
        double t = fabs(sol1.t_tag.TimeDiff(sol2.t_tag.t_));

        vel = (sol1.pos - sol2.pos) / t;

        return vel;
    }

    double cFGOSolver::Vel2Yaw(Vector3d vn){
        return atan2(vn[1],fabs(vn[0]) < 1E-4?1E-4:vn[0]);
    }

    bool cFGOSolver::GnssSol2Ins(Vector3d re,Vector3d ve) {
        if(re.norm()==0.0||ve.norm()==0.0) return false;
        Vector3d blh=Xyz2Blh(re),rn;
        Vector3d wiee(0,0,OMGE_GPS);
        Matrix3d Cne;

        if(fgo_conf_.filter_type == FILTER_BACKWARD) wiee = {0,0,-OMGE_GPS};

        Cne=CalcCen(blh,COORD_NED).transpose();
        pre_imu_info_.rn=blh;
        pre_imu_info_.vn=Cne.transpose()*ve;

        pre_imu_info_.rpy[0] = atan2(-pre_imu_info_.cor_acce[1],-pre_imu_info_.cor_acce[2]);
        pre_imu_info_.rpy[1] = atan(pre_imu_info_.cor_acce[0] / Norm(pre_imu_info_.cor_acce.tail(2),2));

        pre_imu_info_.rpy[2]=Vel2Yaw(pre_imu_info_.vn);
        pre_imu_info_.Cbn=Euler2RotationMatrix(pre_imu_info_.rpy);
        pre_imu_info_.Cbe=Cne*pre_imu_info_.Cbn;

        // lever correct
        pre_imu_info_.re=re-pre_imu_info_.Cbe*fgo_conf_.insC.lever;
        Matrix3d T=VectorSkew(cur_imu_info_.raw_gyro*cur_imu_info_.dt);//当前时刻初始角增量的旋转对称
        Matrix3d Omge=VectorSkew(wiee); //地球自传参数的旋转对称
        pre_imu_info_.ve=ve-pre_imu_info_.Cbe*T*fgo_conf_.insC.lever+Omge*pre_imu_info_.Cbe*fgo_conf_.insC.lever;//对速度进行杆壁改正
        pre_imu_info_.rn=Xyz2Blh(pre_imu_info_.re);//更新改正后的导航系下的位置和速度
        pre_imu_info_.vn=Cne.transpose()*pre_imu_info_.ve;

        return true;
    }

    bool cFGOSolver::InsAlign() {
        int num_sols=5;
        static vector<tSolInfoUnit>sols;
        SOL_STAT  align_level=SOL_FIX;
        if(fgo_conf_.mode_opt==MODE_OPT_PPP){
            align_level=SOL_PPP;
        }
        else if(fgo_conf_.mode_opt==MODE_OPT_SPP){
            align_level=SOL_SPP;
        }

        if(fgo_conf_.insC.ins_align==ALIGN_GNSS_OBS){

            if(gnss_alignor_->SolverProcess(gnss_conf_,rover_idx_)) {
                CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"USING GNSS TO ALIGN INS("<<sols.size()<<"): "<<gnss_alignor_->igcoupled_sol_.pos.transpose();
                if(gnss_alignor_->igcoupled_sol_.stat==align_level){
                    sols.push_back(gnss_alignor_->igcoupled_sol_);
                }
            }
            else{
                sols.clear();
                return false;
            }

        }
        else if(fgo_conf_.insC.ins_align==ALIGN_GNSS_SOL){
            sols.push_back(gnss_sols_.at(gnss_sol_idx_));
        }

        if(sols.size()>=num_sols){
            if(fgo_conf_.insC.ins_align==ALIGN_GNSS_OBS){
                for(int i=0;i<sols.size()-1;i++){
                    if(fabs(sols[i+1].t_tag.TimeDiff(sols[i].t_tag.t_))>fgo_conf_.gnssC.sample_rate){
                        sols.clear();
                        return false;
                    }
                }
            }
            else{
                for(int i=0;i<sols.size();i++){
                    if(sols[i].stat>SOL_FIX||sols[i].stat==SOL_NONE){
                     //  sols.clear();
                      // return false; // zwl 2022.5.15
                    }
                }
            }

            Vector3d ve;
            if(sols.back().vel.norm()==0.0){
                ve=Pos2Vel(sols.back(),*(sols.end()-2));
            }
            else ve=sols.back().vel;

            if(ve.norm()<5.0||cur_imu_info_.raw_gyro.norm()>30.0*D2R){
                sols.clear();
                return false;
            }

            if(GnssSol2Ins(sols.back().pos,ve)){
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" INS INITIALIZATION OK";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pre_imu_info_.re.transpose()<<" m";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT VELOCITY(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.ve.transpose()<<" m/s";
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"INIT ATTITUTE(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.rpy.transpose()*R2D<<" deg";
                igcoupled_sol_.ins_stat=SOL_INS_INIT;
                InsSol2IgcoupledSol(pre_imu_info_,igcoupled_sol_);
                return true;
            }
            else sols.clear();
        }
        return false;
    }

    void cFGOSolver::InitSolver(tIGCOUPLEDConf C){
        fgo_conf_ = C;
        para_ = cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_ = new cOutSol(C);
        out_->InitOutSol(C,C.fileC.sol);
        out_->WriteHead();

        if(fgo_conf_.mode == MODE_IGTC) tc_mode_=true;
        if(fgo_conf_.use_image) image_mode_ = true;

        if(fgo_conf_.mode_opt == MODE_OPT_SPP){
            tIGCOUPLEDConf spp_conf = C;
         //  if(C.mode == MODE_IGLC){
                spp_conf.mode = MODE_SPP;
         //  }
            spp_conf.solC.out_sol = false;
            spp_conf.filter_type = FILTER_FORWARD;
            gnss_solver_ = new cSppSolver(spp_conf);
            gnss_solver_->nav_ = nav_;
            gnss_solver_->InitSolver(spp_conf);
            gnss_conf_ = spp_conf;
            if(C.insC.ins_align == ALIGN_GNSS_OBS){
                spp_conf.mode = MODE_SPP;
                spp_conf.mode_opt = MODE_OPT_KINEMATIC;
                gnss_alignor_ = new cSppSolver(spp_conf);
                gnss_alignor_->nav_ = nav_;
                gnss_alignor_->InitSolver(spp_conf);
                gnss_alignor_->rover_obs_=rover_obs_;
            }
        }
        else if(fgo_conf_.mode_opt == MODE_OPT_PPP){
            tIGCOUPLEDConf ppp_conf = C;
           // if(C.mode == MODE_IGLC) ppp_conf.mode = MODE_PPP;
            ppp_conf.mode = MODE_PPP;
            ppp_conf.solC.out_sol = false;
            ppp_conf.solC.out_head =false;
            ppp_conf.filter_type = FILTER_FORWARD;
            gnss_solver_ = new cPppSolver(ppp_conf);
            gnss_solver_->nav_ = nav_;
            gnss_solver_->InitSolver(ppp_conf);
            gnss_conf_ = ppp_conf;

            if(C.insC.ins_align == ALIGN_GNSS_OBS){
                ppp_conf.mode = MODE_PPP;
                ppp_conf.mode_opt = MODE_OPT_KINEMATIC;
                gnss_alignor_ = new cPppSolver(ppp_conf);
                gnss_alignor_->nav_ = nav_;
                gnss_alignor_->InitSolver(ppp_conf);
                gnss_alignor_->rover_obs_ = rover_obs_;
            }
        }
        else if(fgo_conf_.mode_opt == MODE_OPT_PPK){
            tIGCOUPLEDConf ppk_conf = C;
           // if(C.mode == MODE_IGLC) ppk_conf.mode = MODE_PPK;
            ppk_conf.mode = MODE_PPK;
            ppk_conf.solC.out_sol = false;
            ppk_conf.solC.out_head = false;
            ppk_conf.filter_type = FILTER_FORWARD;
            gnss_solver_ = new cPpkSolver(ppk_conf);
            gnss_solver_->nav_ = nav_;
            gnss_solver_->InitSolver(ppk_conf);
            gnss_conf_ = ppk_conf;

            if(C.insC.ins_align == ALIGN_GNSS_OBS){
                ppk_conf.mode = MODE_PPK;
                ppk_conf.mode_opt = MODE_OPT_KINEMATIC;
                gnss_alignor_ = new cPpkSolver(ppk_conf);
                gnss_alignor_->nav_ = nav_;
                gnss_alignor_->InitSolver(ppk_conf);
                gnss_alignor_->rover_obs_ = rover_obs_;
            }
        }
        else if(fgo_conf_.mode_opt == MODE_OPT_GSOF){
            cDecodeGsof gsof_decoder;
            gsof_decoder.DecodeGsof(C.fileC.gsof,gnss_sols_);
        }
        else if(fgo_conf_.mode_opt == MODE_OPT_SOL){
            ReadSol(C,C.fileC.gsof,gnss_sols_);
        }

        if(C.insC.imu_type == IMU_M39){
            cDecodeImuM39 m39_decoder;
            m39_decoder.DecodeM39(C.fileC.imu,C.insC,imu_data_.data_);
        }
        else if(C.insC.imu_type == IMU_NOVTEL_A1 || C.insC.imu_type == IMU_NOVTEL_CPT || C.insC.imu_type >= IMU_MTI_CSV){
            int week = 0;
            if(C.insC.imu_type >= IMU_POS && !rover_obs_.GetGnssObs().empty()) rover_obs_.GetGnssObs()[0].obs_time.Time2Gpst(&week, nullptr,SYS_GPS);
            else if(C.insC.imu_type >= IMU_POS && !gnss_sols_.empty()) gnss_sols_[0].t_tag.Time2Gpst(&week, nullptr, SYS_GPS);
            cReadImu imu_reader(C.fileC.imu,week);
            imu_reader.SetImu(C);
            imu_reader.Reading();
            imu_data_ = *imu_reader.GetImus();
        }
        if(fgo_conf_.mode_opt > MODE_OPT_SOL) {
            gnss_solver_->rover_obs_ = rover_obs_;
            gnss_solver_->tc_mode_ = false;
        }
       // if(tc_mode_) gnss_solver_->full_Px_ = full_Px_;
        GnssObsList_.clear();
        GnssSolList_.clear();
        gnssSolverList_.clear();
      //  numlList_.clear();
        EpochSatInfolist_.clear();
        preintegrationlist_.clear();
        timelist_.clear();
        statelist_.clear();
        statedatalist_.clear();
        if(fgo_conf_.fgo_window != 0){
            statelist_.resize(fgo_conf_.fgo_window + 1);
            statedatalist_.resize(fgo_conf_.fgo_window + 1);
        }else{
            statelist_.resize(2);
            statedatalist_.resize(2);
        }

        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INIT FGO SOLVER OK"<<endl;
    }

    bool cFGOSolver::MatchGnssSol() {
        double sow1, sow2;
        bool info = false;
        int i;

        for(i = gnss_sol_idx_ > 3 ? gnss_sol_idx_-3:0;i < gnss_sols_.size()&&i>=0;i++){
            sow1 = gnss_sols_[i].t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            sow2 = cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow2-sow1) < 0.5/fgo_conf_.insC.sample_rate){
                info = true;
                gnss_sol_idx_ = i;
                break;
            }
            else if((fgo_conf_.filter_type == FILTER_FORWARD || fgo_conf_.filter_type == FILTER_COMBINED) && (sow1 - sow2) > 1.0 / fgo_conf_.insC.sample_rate){
                info = false;
                break;
            }else if((fgo_conf_.filter_type == FILTER_BACKWARD)  && (sow2 - sow1) > 1.0 / fgo_conf_.insC.sample_rate ){
                info = false;
                break;
            }
        }
        return info;
    }

    bool cFGOSolver::MatchGnssObs() {
        double sow1,sow2;
        int i,week = 0, wod = 0, info = false;

        for(i = rover_idx_-100<0?0:rover_idx_-10; i < rover_obs_.GetGnssObs().size();i++){
            sow1 = rover_obs_.GetGnssObs()[i].obs_time.Time2Gpst(&week, &wod, SYS_GPS);
            sow2 = cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2) < 0.5 / fgo_conf_.insC.sample_rate){
                rover_idx_ = i;
                info = true;
                break;
            }
            else if( (fgo_conf_.filter_type == FILTER_FORWARD || fgo_conf_.filter_type == FILTER_COMBINED ) && (sow1-sow2) > 1.0/fgo_conf_.insC.sample_rate){
                info = false;break;
            }else if((fgo_conf_.filter_type == FILTER_BACKWARD)  && (sow2 - sow1) > 1.0 / fgo_conf_.insC.sample_rate ){
                info = false;
                break;
            }
        }
        return info;
    }

    void cFGOSolver::InsSol2IgcoupledSol(tImuInfoUnit &imu_sol, tSolInfoUnit &igcoupled_sol) {
        igcoupled_sol.t_tag=imu_sol.t_tag;
        igcoupled_sol.pos=imu_sol.re;
        igcoupled_sol.vel=imu_sol.ve;
        igcoupled_sol.att=imu_sol.rpy;

        igcoupled_sol.accl_bias=imu_sol.ba;
        igcoupled_sol.gyro_bias=imu_sol.bg;

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        Matrix3d Pp,Pv,Pa;
        //将权矩阵从Px中取出
        if(gnss_alignor_){
            Pp=gnss_alignor_->full_Px_.block<3,3>(ip,ip);
            Pv=gnss_alignor_->full_Px_.block<3,3>(iv,iv);
            Pa=gnss_alignor_->full_Px_.block<3,3>(ia,ia);
        }
        //取出位置、速度和姿态的精度信息
        for(int i=0;i<3;i++) igcoupled_sol.q_pos[i]=Pp(i,i);
        igcoupled_sol.q_pos[3]=Pp(0,1);
        igcoupled_sol.q_pos[4]=Pp(1,2);
        igcoupled_sol.q_pos[5]=Pp(2,0);

        for(int i=0;i<3;i++) igcoupled_sol.q_vel[i]=Pv(i,i);
        igcoupled_sol.q_vel[3]=Pv(0,1);
        igcoupled_sol.q_vel[4]=Pv(1,2);
        igcoupled_sol.q_vel[5]=Pv(2,0);

        for(int i=0;i<3;i++) igcoupled_sol.q_att[i]=Pa(i,i);
        igcoupled_sol.q_att[3]=Pa(0,1);
        igcoupled_sol.q_att[4]=Pa(1,2);
        igcoupled_sol.q_att[5]=Pa(2,0);
    }

    //将IMU解作为最终解
    void cFGOSolver::PreIntSol2IgcoupledSol(const cTime& time,const tFGOState& state, tSolInfoUnit &igcoupled_sol, Vector3d lever) {
        igcoupled_sol.t_tag = time;
        igcoupled_sol.pos = state.p;
        igcoupled_sol.vel = state.v;
        igcoupled_sol.att = Quaternion2Euler(state.q);


        igcoupled_sol.pos = igcoupled_sol.pos + state.q.toRotationMatrix() * lever;
        igcoupled_sol.vel =

        igcoupled_sol.accl_bias = state.ba;
        igcoupled_sol.gyro_bias = state.bg;

        int ip = para_.IndexPos();
        int iv = para_.IndexVel();
        int ia = para_.IndexAtt();

        Matrix3d Pp,Pv,Pa;
        //将权矩阵从Px中取出
        if(gnss_solver_){
            Pp = gnss_solver_->full_Px_.block<3,3>(ip,ip);
            Pv = gnss_solver_->full_Px_.block<3,3>(iv,iv);
            Pa = gnss_solver_->full_Px_.block<3,3>(ia,ia);
        }

        //取出位置、速度和姿态的精度信息
        for(int i = 0; i < 3; i++) igcoupled_sol.q_pos[i] = Pp(i,i);
        igcoupled_sol.q_pos[3] = Pp(0,1);
        igcoupled_sol.q_pos[4] = Pp(1,2);
        igcoupled_sol.q_pos[5] = Pp(2,0);

        for(int i = 0; i < 3; i++) igcoupled_sol.q_vel[i] = Pv(i,i);
        igcoupled_sol.q_vel[3] = Pv(0,1);
        igcoupled_sol.q_vel[4] = Pv(1,2);
        igcoupled_sol.q_vel[5] = Pv(2,0);

        for(int i = 0; i < 3; i++) igcoupled_sol.q_att[i] = Pa(i,i);
        igcoupled_sol.q_att[3] = Pa(0,1);
        igcoupled_sol.q_att[4] = Pa(1,2);
        igcoupled_sol.q_att[5] = Pa(2,0);
    };

    void cFGOSolver::ImuInterp(const tImuDataUnit &imu01, tImuDataUnit &imu00, tImuDataUnit &imu11, double mid) const {
        double time = mid;
        tImuDataUnit buff = imu01;

        double scale = fabs(buff.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) - time) / (1 / fgo_conf_.insC.sample_rate);

        imu00.t_tag += time - imu00.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
        imu00.gyro = buff.gyro * (1 - scale);
        imu00.acce   = buff.acce * (1 - scale);
        for(int i = 0; i < 3; i++) imu00.odo.vr[i] = buff.odo.vr[i] * (1 - scale);
        imu00.odo.odovel = buff.odo.odovel * (1 - scale);


        imu11.t_tag   = buff.t_tag;
        imu11.gyro = buff.gyro * scale;
        imu11.acce   = buff.acce * scale;
        imu11.odo.odovel = buff.odo.odovel * scale;
        for(int i = 0; i < 3; i++) imu11.odo.vr[i] = buff.odo.vr[i] * scale;
    }

    int cFGOSolver::isNeedInterpolation(tImuDataUnit &imu0, tImuDataUnit &imu1, double mid) {
        double time = mid;
        bool back = false;
        if(imu1.t_tag.TimeDiff(imu0.t_tag.t_) < 0) back = true;

        if(!back && imu0.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) < time && imu1.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) > time){
            double dt = time - imu0.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);

            // 前一个历元接近
            // close to the first epoch
            if (dt < 0.0001) {
                return -1;
            }

            // 后一个历元接近
            // close to the second epoch
            dt = imu1.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) - time;
            if (dt < 0.0001) {
                return 1;
            }

            // 需内插
            // need interpolation
            return 2;
        }else if(back && imu0.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) > time && imu1.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) < time){
            double dt = imu0.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) - time;

            // 前一个历元接近
            // close to the first epoch
            if (dt < 0.0001) {
                return -1;
            }

            // 后一个历元接近
            // close to the second epoch
            dt = time - imu1.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS);
            if (dt < 0.0001) {
                return 1;
            }

            // 需内插
            // need interpolation
            return 2;
        }
        return 0;
    }

    void cFGOSolver::StateTimeUpdate(tFGOState &state,VectorXd &x) {

    }

    double cFGOSolver::findOutage(double sow, vector<double> outageList) {
        if(!outageList.empty()){
            if(sow >= outageList[0] ||  sow <= outageList[0]){
                for(int i = 0; i < outageList.size();i++){
                    if(fgo_conf_.filter_type == FILTER_BACKWARD && sow <= outageList[i] && sow >= outageList[i] - fgo_conf_.gnssC.outage_len){
                        return fabs(sow - outageList[i]);
                    }else if((fgo_conf_.filter_type == FILTER_COMBINED || fgo_conf_.filter_type == FILTER_FORWARD) && sow >= outageList[i]
                    && sow <= outageList[i] + fgo_conf_.gnssC.outage_len){
                        return fabs(sow - outageList[i]);
                    }
                }
            }else return -1;
        }else return -1;
    }



    bool cFGOSolver::SolverProcess(tIGCOUPLEDConf C, int type) {

        InitSolver(C);

        tSatInfoUnit sat_info;
        vector<double> outageList;
        outageList.clear();

        int gnss_obs_flag = false, ins_align = false, gnss_sol_flag = false;

        int week = 0;
        //
        //使用原始观测值还是结果
        int gnss_obs = false, gnss_sol = false;
        double outage_time = 0;
        if(fgo_conf_.gnssC.use_outage) outage_time = fgo_conf_.gnssC.outage_time;

        //imu 数据初始化
            for (; imu_idx_ < imu_data_.data_.size(); imu_idx_++) {
                cur_imu_info_.t_tag = imu_data_.data_[imu_idx_].t_tag;
                cur_imu_info_.raw_gyro = imu_data_.data_[imu_idx_].gyro;
                cur_imu_info_.raw_acce = imu_data_.data_[imu_idx_].acce;
                if (fgo_conf_.mode > MODE_INS) {// 松紧组合
                    if (fgo_conf_.mode_opt == MODE_OPT_GSOF || fgo_conf_.mode_opt == MODE_OPT_SOL) {//用GNSS位置松组合
                        gnss_sol_flag = MatchGnssSol();
                        gnss_sol = true;
                    } else {
                        gnss_obs_flag = MatchGnssObs();
                        gnss_obs = true;

                    }
                    if (gnss_obs_flag || gnss_sol_flag) {
                        pre_imu_info_.t_tag = cur_imu_info_.t_tag;
                        pre_imu_info_.cor_gyro = pre_imu_info_.raw_gyro = cur_imu_info_.raw_gyro;
                        pre_imu_info_.cor_acce = pre_imu_info_.raw_acce = cur_imu_info_.raw_acce;
                        if (gnss_obs_flag) ins_align = InsAlign();
                        else {
                            if (gnss_sols_[gnss_sol_idx_].pos.norm() == 0.0) continue;
                            ins_align = InsAlign();
                        }
                        if (!ins_align) continue;
                    } else continue;
                    if (gnss_obs_flag) rover_idx_++;
                    if (gnss_sol_flag) gnss_sol_idx_++;
                    pre_imu_info_.ba = pre_imu_info_.bg = {0, 0, 0};
                    epoch_sat_obs_.epoch_data.clear();
                    base_epoch_sat_obs_.epoch_data.clear();
                    break;
                }
            }
        if(!ins_align){
            CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID) << "IMU DATA NOT ALIGN" << endl;
            return false;
        }
        tFGOState state_curr = {
                .t_tag = pre_imu_info_.t_tag,
                .p = pre_imu_info_.re,
                .q = RotationMatrix2Quaternion(pre_imu_info_.Cbe),
                .v = pre_imu_info_.ve,
                .bg = {0,0,0},
                .ba = {0,0,0},
                .sodo = 0.0,
                .abv = {fgo_conf_.insC.abv[1], fgo_conf_.insC.abv[2]},
        };

        tImuDataUnit imu_pre = imu_data_.data_[imu_idx_- 1];
        tImuDataUnit imu_cur = imu_data_.data_[imu_idx_++];

        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INITIALIZATION " << pre_imu_info_.t_tag.Time2Gpst(&week, nullptr, SYS_GPS) << " S";

        double sow = 0.0;
        if(gnss_obs) {//初始对准时刻的观测值
            GnssObsList_.push_back(rover_obs_.GetGnssObs()[rover_idx_ - 1]);
            if(gnss_solver_->SolverProcess(gnss_conf_, rover_idx_ - 1)){
                GnssSolList_.push_back(gnss_solver_->igcoupled_sol_);
                if(tc_mode_){// 有问题
                    //
				}

            }
            sow = round(GnssObsList_.back().obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS));
        }
        else if(gnss_sol){
            GnssSolList_.push_back(gnss_sols_[gnss_sol_idx_ - 1]);
            sow = round(GnssSolList_.back().t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS));
        }



        timelist_.push_back(sow);
        statelist_[0] = state_curr;
        statedatalist_[0] = cPreIntegration::stateToData(state_curr);


        //初始预积分
        preintegrationlist_.emplace_back(
                cPreIntegration::createPreintegration(fgo_conf_, imu_pre, state_curr)
                );

        //读取下一个整秒GNSS
        tSolInfoUnit gnssSol;
        tEpochSatUnit gnssObs;
        if(gnss_sol){
            gnssSol = gnss_sols_[gnss_sol_idx_++];
        }else if(gnss_obs){
            gnssObs = rover_obs_.GetGnssObs()[rover_idx_++];
        }

        //边缘化信息
        shared_ptr<cMarginalizationInfo> last_marginalization_info;
        vector<double *> last_marginalization_parameter_blocks;

        //下一个积分节点
        sow += fgo_conf_.integration_len;


        while(true) {
            if (imu_idx_ > imu_data_.data_.size()) break;


            // 加入IMU数据
            preintegrationlist_.back()->addNewImu(imu_cur);
            imu_pre = imu_cur;
            imu_cur = imu_data_.data_[imu_idx_++];

            if (((fgo_conf_.filter_type == FILTER_FORWARD || fgo_conf_.filter_type == FILTER_COMBINED) &&
                 imu_cur.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) > sow) ||
                (fgo_conf_.filter_type == FILTER_BACKWARD &&
                 imu_cur.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) < sow)) {
                //当前IMU数据时间等于GNSS数据时间, 读取新的GNSS
                if ((gnss_sol && fabs(gnssSol.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) - sow) <
                                 (1.0 / fgo_conf_.insC.sample_rate))
                    || gnss_obs && fabs(gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS) - sow) <
                                   (1.0 / fgo_conf_.insC.sample_rate)) {
                    if (gnss_sol) {
                        GnssSolList_.push_back(gnssSol);
                        gnssSol = gnss_sols_[gnss_sol_idx_++];
                        //固定阈值抗差
                        //  while (gnssSol.q_pos[0] > 0.04 || gnssSol.q_pos[1] > 0.04 || gnssSol.q_pos[2] > 0.04)
                        //   gnssSol = gnss_sols_[gnss_sol_idx_++];
                    } else if (gnss_obs) {
                        GnssObsList_.push_back(gnssObs);
                        if (gnss_solver_->SolverProcess(gnss_conf_, rover_idx_ - 1)) {
                            GnssSolList_.push_back(gnss_solver_->igcoupled_sol_);
                            if (tc_mode_) { 
                               // 
							}
                        }
                        gnssObs = rover_obs_.GetGnssObs()[rover_idx_++];
                    }


                    //中断配置
                    if ((fgo_conf_.filter_type == FILTER_FORWARD || fgo_conf_.filter_type == FILTER_COMBINED) &&
                        fgo_conf_.gnssC.use_outage) {
                        if (gnss_sol && lround(gnssSol.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS)) == outage_time) {
                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << " GNSS OUTAGE AT " << outage_time << " s" << endl;
                            while (gnssSol.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) <
                                   outage_time + fgo_conf_.gnssC.outage_len || lround(gnssSol.t_tag.Time2Gpst(
                                    nullptr, nullptr, SYS_GPS)) == outage_time + fgo_conf_.gnssC.outage_len) {
                                gnssSol = gnss_sols_[gnss_sol_idx_++];
                            }
                            if (fgo_conf_.gnssC.outage_period > fgo_conf_.gnssC.outage_len) {
                                outage_time += fgo_conf_.gnssC.outage_period;
                                outageList.push_back(outage_time);
                            }
                        } else if (gnss_obs &&
                                   lround(gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS)) == outage_time) {
                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << " GNSS OUTAGE AT " << outage_time << " s" << endl;

                            while (gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS) <
                                   outage_time + fgo_conf_.gnssC.outage_len ||
                                   lround(gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS)) ==
                                   outage_time + fgo_conf_.gnssC.outage_len) {
                                gnssObs = rover_obs_.GetGnssObs()[rover_idx_++];
                            }
                            if (fgo_conf_.gnssC.outage_period > fgo_conf_.gnssC.outage_len) {
                                outage_time += fgo_conf_.gnssC.outage_period;
                                outageList.push_back(outage_time);
                            }
                        }
                    } else if (fgo_conf_.filter_type == FILTER_BACKWARD && fgo_conf_.gnssC.use_outage) {
                        if (gnss_sol && lround(gnssSol.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS)) == outage_time) {
                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << " GNSS OUTAGE AT " << outage_time << " s" << endl;
                            while (gnssSol.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) >
                                   outage_time - fgo_conf_.gnssC.outage_len || lround(gnssSol.t_tag.Time2Gpst(
                                    nullptr, nullptr, SYS_GPS)) == outage_time - fgo_conf_.gnssC.outage_len) {
                                gnssSol = gnss_sols_[gnss_sol_idx_++];
                            }
                            if (fgo_conf_.gnssC.outage_period > fgo_conf_.gnssC.outage_len &&
                                (outage_time - fgo_conf_.gnssC.outage_period - fgo_conf_.gnssC.outage_len) >=
                                fgo_conf_.gnssC.outage_time) {
                                outage_time -= fgo_conf_.gnssC.outage_period;
                                outageList.push_back(outage_time);
                            }
                        } else if (gnss_obs &&
                                   lround(gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS)) == outage_time) {
                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << " GNSS OUTAGE AT " << outage_time << " s" << endl;
                            while (gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS) >
                                   outage_time - fgo_conf_.gnssC.outage_len ||
                                   lround(gnssObs.obs_time.Time2Gpst(nullptr, nullptr, SYS_GPS)) ==
                                   outage_time - fgo_conf_.gnssC.outage_len) {
                                gnssObs = rover_obs_.GetGnssObs()[rover_idx_++];
                            }
                            if (fgo_conf_.gnssC.outage_period > fgo_conf_.gnssC.outage_len) {
                                outage_time -= fgo_conf_.gnssC.outage_period;
                                outageList.push_back(outage_time);
                            }
                        }
                    }

                    if ((gnss_sol && gnss_sol_idx_ > gnss_sols_.size()) ||
                        (gnss_obs && rover_idx_ > rover_obs_.GetGnssObs().size())) {
                        if (gnss_sol) {
                            gnssSol.t_tag.t_.sec = 0.0;
                            gnssSol.t_tag.t_.long_time = time_t(0);
                        } else if (gnss_obs) {
                            gnssObs.obs_time.t_.sec = 0.0;
                            gnssSol.t_tag.t_.long_time = time_t(0);
                        }
                        CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << " GNSS OVER AT " << sow << " s" << endl;
                    }
                }

                // IMU 内插处理
                if (isNeedInterpolation(imu_pre, imu_cur, sow) == 1) {
                    preintegrationlist_.back()->addNewImu(imu_cur);

                    imu_pre = imu_cur;
                    imu_cur = imu_data_.data_[imu_idx_++];
                } else if (isNeedInterpolation(imu_pre, imu_cur, sow) == 2) {
                    ImuInterp(imu_cur, imu_pre, imu_cur, sow);
                    preintegrationlist_.back()->addNewImu(imu_pre);
                }

                //下一积分节点
                timelist_.push_back(sow);
                sow += fgo_conf_.integration_len;

                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "NEXT INTEGRATION SECOND " << sow << " S" << endl;

                //当前整秒状态加入到划窗中
                state_curr = preintegrationlist_.back()->currentState();

                // gnss 数据丢失过久重新初始化 ba 和 bg  zwl 2022.6.10
                if (gnss_sol && GnssSolList_.empty() || state_curr.ba.norm() > 1.0 ||
                    state_curr.bg.norm() > 1.0) {
                    state_curr.ba = state_curr.bg = {0, 0, 0};
                } else if (gnss_obs && GnssObsList_.empty() || state_curr.ba.norm() > 1.0 ||
                         state_curr.bg.norm() > 1.0)
                    state_curr.ba = state_curr.bg = {0, 0, 0};

                statelist_[preintegrationlist_.size()] = state_curr;
                statedatalist_[preintegrationlist_.size()] = cPreIntegration::stateToData(state_curr);
                //构建优化问题
                {
                    CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "CONSTRUCT FGO SOLVER " << endl;
                    ceres::Solver solver;
                    ceres::Problem problem;
                    ceres::Solver::Summary summary;
                    ceres::Solver::Options options;
                   options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
                    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
                    options.num_threads = 4;
                    options.max_num_iterations = 200;
                   // options.use_nonmonotonic_steps = true;
                    options.minimizer_progress_to_stdout = true;

                    //参数块
                    for (size_t k = 0; k <= preintegrationlist_.size(); k++) {
                        //位姿
                        ceres::LocalParameterization *parameterization = new (cPoseParameterization);
                        problem.AddParameterBlock(statedatalist_[k].pose, cPreIntegration::numPoseParameter(),
                                                  parameterization);
                        problem.AddParameterBlock(statedatalist_[k].mix, cPreIntegration::numMixParameter(fgo_conf_));
                        if (tc_mode_) {
						 //
                        }


                    }
                    CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "PARAMETERIZATION  OK" << endl;

                    //GNSS残差
                    if (!tc_mode_) {
                        int index = 0;
                        for (auto &data: GnssSolList_) {
                            auto factor = new cGnssFactor_SOL(data, fgo_conf_.insC.lever);
                            for (size_t i = index; i <= preintegrationlist_.size(); ++i) {
                                if (fabs(data.t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS) - timelist_[i]) <
                                    (1 / fgo_conf_.insC.sample_rate)) {
                                    problem.AddResidualBlock(factor, nullptr, statedatalist_[i].pose);
                                    index++;
                                    break;
                                }

                            }
                        }
                    }
                    else if (tc_mode_) {
					     //
                        }
                        CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "GNSS FGO FACTOR OK" << endl;

                        //预积分残差
                        for (size_t k = 0; k < preintegrationlist_.size(); k++) {
                            auto factor = new cPreIntegrationFactor(preintegrationlist_[k]);
                            problem.AddResidualBlock(factor, nullptr, statedatalist_[k].pose, statedatalist_[k].mix,
                                                     statedatalist_[k + 1].pose, statedatalist_[k + 1].mix);

                        }
                        {
                            //IMU 误差控制
                            auto factor = new ImuErrorFactor(*preintegrationlist_.rbegin());
                            problem.AddResidualBlock(factor, nullptr, statedatalist_[preintegrationlist_.size()].mix);

                        }
                        CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "IMU FGO FACTOR OK" << endl;

                        //边缘化残差
                        if (last_marginalization_info && last_marginalization_info->isValid()) {
                            auto factor = new cMarginalizationFactor(last_marginalization_info);
                            problem.AddResidualBlock(factor, nullptr, last_marginalization_parameter_blocks);
                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "MARGINAL FGO FACTOR OK" << endl;
                        }

                        //求解最小二乘
                        solver.Solve(options, &problem, &summary);

                        //输出优化过程
                        cout << summary.FullReport() << "\n";


                    }

                    if (preintegrationlist_.size() == static_cast<size_t>(fgo_conf_.fgo_window)) {
                        {
                            //边缘化
                            std::shared_ptr<cMarginalizationInfo> marginalization_info = std::make_shared<cMarginalizationInfo>();
                            if (last_marginalization_info && last_marginalization_info->isValid()) {
                                std::vector<int> marginilized_index;
                                for (size_t k = 0; k < last_marginalization_parameter_blocks.size(); k++) {
                                    if (last_marginalization_parameter_blocks[k] == statedatalist_[0].pose ||
                                        last_marginalization_parameter_blocks[k] == statedatalist_[0].mix) {
                                        marginilized_index.push_back(static_cast<int>(k));
                                    }
                                }


                                auto factor = std::make_shared<cMarginalizationFactor>(last_marginalization_info);
                                auto residual = std::make_shared<cResidualBlockInfo>(
                                        factor, nullptr, last_marginalization_parameter_blocks, marginilized_index);
                                marginalization_info->addResidualBlockInfo(residual);

                            }

                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "MARGINALIZATION PARAMETER OK" << endl;

                            //IMU残差
                            {
                                auto factor = std::make_shared<cPreIntegrationFactor>(preintegrationlist_[0]);
                                auto residual = std::make_shared<cResidualBlockInfo>(
                                        factor, nullptr,
                                        std::vector<double *>{statedatalist_[0].pose, statedatalist_[0].mix,
                                                              statedatalist_[1].pose,
                                                              statedatalist_[1].mix},
                                        std::vector<int>{0, 1});
                                marginalization_info->addResidualBlockInfo(residual);
                            }

                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "IMU MARGINALIZATION OK" << endl;

                            //GNSS残差
                            {
                                if (!tc_mode_) {
                                    if (fabs(
                                            timelist_[0] - GnssSolList_[0].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS)) <
                                        (1.0 / fgo_conf_.insC.sample_rate)) {
                                        auto factor = std::make_shared<cGnssFactor_SOL>(GnssSolList_[0],
                                                                                        fgo_conf_.insC.lever);
                                        auto residual = std::make_shared<cResidualBlockInfo>(
                                                factor, nullptr, std::vector<double *>{statedatalist_[0].pose},
                                                std::vector<int>{0});
                                        marginalization_info->addResidualBlockInfo(residual);
                                    }
                                }
                                else if (tc_mode_) {
								  //
                                }

                            }

                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "GNSS MARGINALIZATION OK" << endl;


                            // 边缘化处理
                            // do marginalization
                            marginalization_info->marginalization();

                            // 数据指针调整
                            // get new pointers
                            std::unordered_map<long, double *> address;
                            for (size_t k = 1; k <= preintegrationlist_.size(); k++) {
                                address[reinterpret_cast<long>(statedatalist_[k].pose)] = statedatalist_[k - 1].pose;
                                address[reinterpret_cast<long>(statedatalist_[k].mix)] = statedatalist_[k - 1].mix;
                                if (tc_mode_) {
                                //
								}
                            }
                            last_marginalization_parameter_blocks = marginalization_info->getParamterBlocks(address);
                            last_marginalization_info = std::move(marginalization_info);

                        }
                        CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "MARGINALIZATION OK" << endl;

                        // 滑窗处理
                        {
                            if (gnss_sol && lround(timelist_[0]) ==
                                            lround(GnssSolList_[0].t_tag.Time2Gpst(nullptr, nullptr, SYS_GPS))) {
                                GnssSolList_.erase(GnssSolList_.begin());
                            } else if (gnss_obs &&
                                       lround(timelist_[0]) == lround(GnssObsList_[0].obs_time.Time2Gpst(nullptr,
                                                                                                         nullptr,
                                                                                                         SYS_GPS))) {
                                GnssObsList_.erase(GnssObsList_.begin());
                                GnssSolList_.erase(GnssSolList_.begin());
                                if (tc_mode_) {
								     //
                                }
                            }
                            timelist_.erase(timelist_.begin());
                            preintegrationlist_.erase(preintegrationlist_.begin());

                            for (int k = 0; k < fgo_conf_.fgo_window; k++) {
                                statedatalist_[k] = statedatalist_[k + 1];
                                statelist_[k] = cPreIntegration::stateFromData(statedatalist_[k]);
                                if (tc_mode_) {
								//
								}
                            }
                            statelist_[fgo_conf_.fgo_window] = cPreIntegration::stateFromData(
                                    statedatalist_[fgo_conf_.fgo_window]);
                            if (tc_mode_){
							//
							}
                            state_curr = statelist_[fgo_conf_.fgo_window];

                            CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "WINDOW MOVE" << endl;
                        }
                    } else {
                        state_curr = cPreIntegration::stateFromData(statedatalist_[preintegrationlist_.size()]);
                        if (tc_mode_){
						//
						}
                    }

                    //输出结果
                    cTime time;
                    time.Gpst2Time(week, *timelist_.rbegin(), SYS_GPS);
                    PreIntSol2IgcoupledSol(time, state_curr, igcoupled_sol_, fgo_conf_.insC.lever);
                        out_->WriteSol(igcoupled_sol_, epoch_sat_info_collect_);


                    //建立新的预积分
                    preintegrationlist_.emplace_back(
                            cPreIntegration::createPreintegration(fgo_conf_, imu_pre, state_curr)
                    );
                }else{
                    auto integration = *preintegrationlist_.rbegin();
                    //输出结果
                    PreIntSol2IgcoupledSol(integration->endTime(), integration->currentState(), igcoupled_sol_,
                                      fgo_conf_.insC.lever);
                        out_->WriteSol(igcoupled_sol_, epoch_sat_info_collect_);
                    }

                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "STATE CURRENT" << endl;
                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "FGO POSITION :" << state_curr.p << "m" << endl;
                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "FGO VELOCITY :" << state_curr.v << "m/s" << endl;
                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "FGO ATTITUDE :" << state_curr.q << endl;
                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "FGO BA       :" << state_curr.ba << endl;
                CLOG(DEBUG, ELPP_CURR_FILE_LOGGER_ID) << "FGO BG       :" << state_curr.bg << endl;

            }




	}    
}
