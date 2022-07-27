#include "CmnFunc.h"
#include "ReadFiles.h"
#include "Solver.h"
INITIALIZE_EASYLOGGINGPP

using namespace IGCoupled;

static tIGCOUPLEDConf kConf;

static int AutoMatchFile(string rover_path)
{
    char f[1024]={'\0'};

    int dgps_flag=kConf.mode==MODE_DGNSS||kConf.mode==MODE_PPK||kConf.mode_opt==MODE_OPT_PPK;
    int ins=kConf.mode>=MODE_INS;

    if(dgps_flag){
        vector<string> splits=MultiSplitStr(rover_path,".");
        sprintf(f,"%s_base.%s",splits[0].c_str(),splits[1].c_str());
        if((access(f,0))==-1){
            CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<"BASE OBSERVATION FILE NO EXIST path="<<f;
            return 0;
        }
        kConf.fileC.base=f;
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"BASE OBSERVATION FILE: "<<f;
    }

    if(ins){
        vector<string> splits=MultiSplitStr(rover_path,".");
        sprintf(f,"%s.imu",splits[0].c_str());
        if((access(f,0))==-1){
            CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<"IMU FILE NO EXIST path="<<f;
            return 0;
        }
        kConf.fileC.imu=f;
    }

    if(kConf.mode_opt==MODE_OPT_KINEMATIC&&kConf.solC.out_err_fmt){
        vector<string> splits=MultiSplitStr(rover_path,".");
        sprintf(f,"%s.reff",splits[0].c_str());
        if((access(f,0))!=-1){
            kConf.fileC.kine_ref=f;
        }
    }

    cMatchFile match_file;
    match_file.InitMatchFile(kConf,FILEPATHSEP);
    if(!match_file.MatchFileAuto()){
        return 0;
    }

    return 1;
}

void PrintHelp(){
    fprintf(stdout,"IGCoupled Usage: \n"
                   "-C  configuration file path\n"
                   "-M  processing mode(SPP-KINE,PPP-KINE,PPP-STATIC,PPK-KINE,PPK-STATIC,IGLC-GSOF,IGLC-PPP,IGLC-PPK,IGTC-PPP,IGTC-PPK)\n"
                   "-S  selected GNSS(G,B,E,R,GBJ,GR,..)\n"
                   "-L  debug level(1:Debug 32:Warning 128:Info)\n");
}

static void CheckConf()
{
    if(kConf.mode==MODE_SPP||kConf.mode_opt==MODE_OPT_SPP){
        kConf.gnssC.ion_opt=ION_KLB;
        kConf.gnssC.trp_opt=TRP_SAAS;
        kConf.gnssC.eph_opt=EPH_BRD;
    }
    else if(kConf.mode==MODE_PPP||kConf.mode_opt==MODE_OPT_PPP){
        kConf.gnssC.eph_opt=EPH_PRE;
    }
    else if(kConf.mode==MODE_PPK||kConf.mode_opt==MODE_OPT_PPK){

    }
    if(kConf.gnssC.ion_opt==ION_IF_DUAL) kConf.gnssC.frq_opt=FRQ_TRIPLE;
}

static void LoadConf()
{
    int i;
    Config::Ptr_ config=Config::GetInstance();

    kConf.data_dir=config->Get<string>("data_dir");
    kConf.use_custom_dir=config->Get<int>("use_custom_dir");
    kConf.site_name=config->Get<string>("site_name");
    vector<int> epoch;
    epoch=config->GetArray<int>("prc_date");
    double ep[6]={0};
    for(int i=0;i<3;i++) ep[i]=epoch[i];
    kConf.prc_date.Epoch2Time(ep);
	kConf.estimator = static_cast<SOLVE_ESTIMATOR>(config->Get<int>("estimator"));
    kConf.filter_type= static_cast<FILTER_TYPE>(config->Get<int>("filter_type"));
	kConf.fgo_window = config->Get<int>("fgo_window");
	kConf.integration_len = config->Get<double>("integration_len");

    tGnssConf *gnssC=&kConf.gnssC;
    gnssC->restart_gap=config->Get<int>("restart_gap");
    gnssC->ele_min=config->Get<double>("ele_min");
    gnssC->frq_opt= static_cast<GNSS_FRQ_OPT>(config->Get<int>("frq_opt"));
    vector<int>frqs;
    frqs=config->GetArray<int>("gps_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_GPS][i]=frqs[i];
    frqs.clear();
    frqs=config->GetArray<int>("bds_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_BDS][i]=frqs[i];
    frqs.clear();
    frqs=config->GetArray<int>("gal_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_GAL][i]=frqs[i];
    frqs.clear();
    frqs=config->GetArray<int>("glo_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_GLO][i]=frqs[i];
    frqs.clear();
    frqs=config->GetArray<int>("qzs_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_QZS][i]=frqs[i];
    frqs.clear();
    frqs=config->GetArray<int>("bd3_frq");
    for(i=0;i<frqs.size();i++) gnssC->gnss_frq[SYS_INDEX_QZS+1][i]=frqs[i];
    frqs.clear();
    gnssC->use_bd3=config->Get<int>("use_bd3");
    gnssC->est_bd3_isb=config->Get<int>("est_bd3_isb");
    gnssC->adj_obs=config->Get<int>("adj_obs");
    gnssC->use_doppler=config->Get<int>("use_doppler");
    gnssC->code_phase_ratio=config->Get<double>("code_phase_ratio");
    vector<double>ratio;
    ratio=config->GetArray<double>("meas_err_factor");
    for(i=0;i<ratio.size();i++) gnssC->meas_err_factor[i]=ratio[i];
    gnssC->ac_opt= static_cast<GNSS_AC_OPT>(config->Get<int>("ac_opt"));
    gnssC->eph_opt= static_cast<GNSS_EPH_OPT>(config->Get<int>("eph_opt"));
    gnssC->ion_opt= static_cast<GNSS_ION_OPT>(config->Get<int>("ion_opt"));
    gnssC->trp_opt= static_cast<GNSS_TRP_OPT>(config->Get<int>("trp_opt"));
    gnssC->tid_opt= static_cast<GNSS_TID_OPT>(config->Get<int>("tid_opt"));
    gnssC->glo_ifcb_opt= static_cast<GLO_IFCB_OPT>(config->Get<int>("glo_ifcb_opt"));
    gnssC->sat_pcv=config->Get<int>("sat_pcv");
    gnssC->rec_ant=config->Get<int>("rec_ant");
    ratio.clear();
    ratio=config->GetArray<double>("cs_thres");
    for(i=0;i<ratio.size();i++) gnssC->cs_thres[i]=ratio[i];
    gnssC->max_pdop=config->Get<double>("max_pdop");
    gnssC->max_prior=config->Get<double>("max_prior");
    gnssC->max_inno=config->Get<double>("max_inno");
    gnssC->max_out=config->Get<int>("max_out");
    ratio.clear();
    ratio=config->GetArray<double>("ait_psd");
    for(i=0;i<ratio.size();i++) gnssC->ait_psd[i]=ratio[i];
    gnssC->check_dual_phase=config->Get<int>("check_dual_phase");
    gnssC->ar_mode= static_cast<GNSS_AR_MODE>(config->Get<int>("ar_mode"));
    gnssC->ar_prod= static_cast<GNSS_AR_PROD>(config->Get<int>("ar_prod"));
    gnssC->glo_ar_mode= static_cast<GLO_AR_MODE>(config->Get<int>("glo_ar_mode"));
    gnssC->bds_ar_mode= config->Get<int>("bds_ar_mode");
    gnssC->gal_ar_mode= config->Get<int>("gal_ar_mode");
    ratio.clear();
    ratio=config->GetArray<double>("ar_thres");
    for(i=0;i<ratio.size();i++) gnssC->ar_thres[i]=ratio[i];
    gnssC->ar_el_mask=config->Get<double>("ar_el_mask");
    gnssC->min_sat_num2fix=config->Get<int>("min_sat_num2fix");
    gnssC->min_sat_num2drop=config->Get<int>("min_sat_num2drop");
    gnssC->min_lock2fix=config->Get<int>("min_lock2fix");
    gnssC->min_sat_num2hold=config->Get<int>("min_sat_num2hold");
    gnssC->min_fix2hold=config->Get<int>("min_fix2hold");
    gnssC->hold_er_mask=config->Get<double>("hold_el_mask");
    gnssC->ar_filter=config->Get<int>("ar_filter");
    gnssC->par_ar=config->Get<int>("partial_ar");
    gnssC->res_qc= config->Get<int>("res_qc");
    ratio.clear();
    ratio=config->GetArray<double>("base_coord");
	gnssC->use_outage = false;
	gnssC->use_outage = config->Get<bool>("use_outage");
	gnssC->outage_time = config->Get<double>("outage_time");
	gnssC->outage_len = config->Get<double>("outage_len");
	gnssC->outage_period = config->Get<double>("outage_period");
    if(ratio[0]==1.0){
        // ecef
        for(i=0;i<3;i++) gnssC->rb[i]=ratio[i+1];
    }
    else if(ratio[0]==0.0){
        // blh
        Vector3d blh;
        blh[0]=ratio[1]*D2R;
        blh[1]=ratio[2]*D2R;
        blh[2]=ratio[3];
        gnssC->rb=Blh2Xyz(blh);
    }

    tInsConf *insC=&kConf.insC;
    insC->imu_type= static_cast<IMU_TYPE>(config->Get<int>("imu_type"));
    insC->coord_type= static_cast<IMU_COORD_TYPE>(config->Get<int>("coord_type"));
    insC->data_format= static_cast<IMU_DATA_FORMAT>(config->Get<int>("data_format"));
    insC->gyro_val_format= static_cast<GYRO_DATA_FORMAT>(config->Get<int>("gyro_val_format"));
    insC->sample_rate=config->Get<double>("sample_hz");
    insC->ins_align=static_cast<INS_ALIGN>(config->Get<int>("ins_align"));
    ratio.clear();
    ratio=config->GetArray<double>("lever");
    for(i=0;i<ratio.size();i++) insC->lever[i]=ratio[i];
    insC->correction_time_ba=config->Get<double>("correction_time_ba");
    insC->correction_time_bg=config->Get<double>("correction_time_bg");
    insC->init_pos_unc=config->Get<double>("init_pos_unc");
    insC->init_vel_unc=config->Get<double>("init_vel_unc");
    insC->init_att_unc=config->Get<double>("init_att_unc");
    insC->init_ba_unc=config->Get<double>("init_ba_unc");
    insC->init_bg_unc=config->Get<double>("init_bg_unc");
    insC->init_ra_unc = config->Get<double>("init_ra_unc");
    insC->init_rg_unc = config->Get<double>("init_rg_unc");
    insC->init_sa_unc = config->Get<double>("init_sa_unc");
    insC->init_sg_unc = config->Get<double>("init_sg_unc");
    insC->psd_acce=config->Get<double>("psd_acce");
    insC->psd_gyro=config->Get<double>("psd_gyro");
    insC->psd_ba=config->Get<double>("psd_ba");
    insC->psd_bg=config->Get<double>("psd_bg");
    int flag_err = config->Get<int>("err_model");
    if (flag_err == 1) insC->err_model = true;
    else insC->err_model = false;
    if (insC->err_model) {
        insC->est_ra = insC->est_rg = insC->est_sa = insC->est_sg = true;
    }
    
    insC->use_odo = config->Get<bool>("use_odo");
	ratio.clear();
	ratio = config->GetArray<double>("abv");
	for(i = 0; i < ratio.size(); i++) insC->abv[i] = ratio[i];
    ratio.clear();
    ratio=config->GetArray<double>("lodo");
	for(i = 0; i < ratio.size(); i++) insC->lodo[i] = ratio[i];
	insC->init_odo_unc = config->Get<double>("init_odo_unc");
	insC->odo_srw = config->Get<double>("odo_srw");

    tSolConf *solC=&kConf.solC;
    solC->out_sol=config->Get<int>("out_sol");
    solC->sol_coord= static_cast<COORDINATE_TYPE>(config->Get<int>("sol_coord"));
    solC->sol_fmt=config->Get<int>("sol_fmt");
    solC->out_head=config->Get<int>("out_head");
    solC->out_vel=config->Get<int>("out_vel");
    solC->out_trp=config->Get<int>("out_trp");
    solC->out_att=config->Get<int>("out_att");
    solC->out_ba=config->Get<int>("out_ba");
    solC->out_bg=config->Get<int>("out_bg");
    solC->out_stat=config->Get<int>("out_stat");
    solC->out_err_fmt=config->Get<int>("out_err_fmt");
    solC->out_ins_mech_frq=config->Get<int>("out_ins_mech_frq");
}

static int ParsePara(int arc,char *arv[], string& conf_file)
{
    string mode;
    const char *p;
    double ep[6]={0};
    int mask=SYS_NONE;
    int level=128;

    for(int i=0;i<arc;i++){
        if(!strcmp(arv[i],"-C")&&i+1<arc){
            conf_file=arv[++i];
            if(conf_file.empty()){
                fprintf(stderr,"CONFIGURATION FILE REQUIRED\n");
                return 0;
            }
        }
        else if(!strcmp(arv[i],"-M")&&i+1<arc){
            mode=arv[++i];
        }
        else if(!strcmp(arv[i],"-S")&&i+1<arc){
            p=arv[++i];
            for(;*p&&*p!=' ';p++){
                switch(*p){
                    case 'G': mask|=SYS_GPS;break;
                    case 'B': mask|=SYS_BDS;break;
                    case 'E': mask|=SYS_GAL;break;
                    case 'R': mask|=SYS_GLO;break;
                    case 'J': mask|=SYS_QZS;break;
                }
            }
            if(mask==SYS_NONE){
                fprintf(stderr,"SATELLITE SYSTEM SET ERROR: %s\n", p);
            }
        }
        else if(!strcmp(arv[i],"-L")&&i+1<arc){
            level=atoi(arv[++i]);
        }
    }
    kConf.gnssC.nav_sys=mask;

    vector<string>mode_opt;
    SplitString(mode,mode_opt,"-");
    if(!mode_opt.size()){
        fprintf(stderr,"PROCESS MODE ERROR\n");
        return 0;
    }
    for(int i=0;kIgcoupledModeStr;i++){
        if(mode_opt[0]==kIgcoupledModeStr[i]){
            kConf.mode= static_cast<IGCOUPLED_MODE>(i);break;
        }
        if(i>7) return 0;
    }
    for(int i=0;kIgcoupledModeOptStr;i++){
        if(mode_opt[1]==kIgcoupledModeOptStr[i]){
            kConf.mode_opt= static_cast<IGCOUPLED_MODE_OPT>(i);break;
        }
        if(i>7) return 0;
    }

    string logini_path = SetLogConfPath("");
    InitLog(arc,arv,logini_path,level);

    Config::Ptr_ config=Config::GetInstance();
    if(!config->Open(conf_file)){
        fprintf(stderr,"OPEN CONFIGURATION FILE ERROR file=%s\n", conf_file.c_str());
        return 0;
    }
    else{
        LoadConf();
    }

	return 1;
}

static int Processer()
{
    string data_dir;
    char f[1024]={'\0'};

    if(kConf.data_dir.empty()){
        CLOG(FATAL,ELPP_CURR_FILE_LOGGER_ID)<<"PROCESS DATA DIR IS EMPTY: dir="<<kConf.data_dir;
        return 0;
    }

    int year,week,doy;
    kConf.prc_date.Time2Epoch();
    year=kConf.prc_date.GetEpoch()[0];
    kConf.prc_date.Time2Gpst(&week,nullptr,SYS_GPS);
    doy=kConf.prc_date.Time2Doy();
    if(!kConf.use_custom_dir){
        sprintf(f,"%s%c%d%c%d%c%03d%cobs",kConf.data_dir.c_str(),FILEPATHSEP,year,FILEPATHSEP,week,FILEPATHSEP,doy,FILEPATHSEP);
        data_dir=f;
    }
    else{
        data_dir=kConf.data_dir;
    }

    CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"PROCESS DATA DIR: "<<data_dir;

    DIR *dir;
    if(!(dir=opendir(data_dir.c_str()))){
        CLOG(FATAL,ELPP_CURR_FILE_LOGGER_ID)<<"PROCESS DATA DIR OPEN ERROR";
        return 0;
    }

    cSolver *solver;
    switch(kConf.mode){
        case MODE_SPP:  solver=new cSppSolver(kConf);break;
        case MODE_PPP:  solver=new cPppSolver(kConf);break;
        case MODE_PPK:  solver=new cPpkSolver(kConf);break;
        //case MODE_IGLC: 
       // case MODE_IGTC:
                       // solver=new cFusionSolver(kConf);break;
    }

	if(kConf.mode >= IGCoupled::MODE_INS){
		switch(kConf.estimator){
			case IGCoupled::SOLVE_FGO : solver = new cFGOSolver(kConf);break;
			case IGCoupled::SOLVE_LSQ :
			case IGCoupled::SOLVE_KF  :
				solver = new cFusionSolver(kConf);break;
		}
	}
    // loose coupled with gsof / sol
    if(kConf.mode==MODE_IGLC&&(kConf.mode_opt==MODE_OPT_GSOF||kConf.mode_opt==MODE_OPT_SOL)){
        Config::Ptr_ config=Config::GetInstance();
        string file=config->Get<string>("imu");
		file = data_dir + file;
        if((access(file.c_str(),0))==-1){
            CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<"IMU FILE NO EXIST path="<<file;
            return 0;
        }
        kConf.fileC.imu=file;
        if(kConf.mode_opt==MODE_OPT_GSOF){
            file=config->Get<string>("gsof");
        }
        else if(kConf.mode_opt==MODE_OPT_SOL){
            file=config->Get<string>("gnss_sol");
        }
		file = data_dir + file;
        if((access(file.c_str(),0))==-1){
            CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<"GNSS SOLUTION FILE NO EXIST path="<<f;
            return 0;
        }
        kConf.fileC.gsof=file;
        kConf.fileC.sol=config->Get<string>("rslt");

        solver->SolverProcess(kConf,0);

        return true;
    }

    struct dirent *file;
    char *ext;
    kConf.site_idx=0;
    bool single_flag=kConf.site_name.empty()?false:true;
    while((file=readdir(dir))!= nullptr){
        if(strncmp(file->d_name,".",1)==0) continue;
        else if(strstr(file->d_name,"base")) continue;
        else if(strstr(file->d_name,"gsof")) continue;
        else if(!(ext=strrchr(file->d_name,'.'))) continue;
        else if(!strstr(ext+3,"o")) continue;
        else if(!kConf.site_name.empty()&&!strstr(file->d_name,kConf.site_name.c_str())) continue;

        f[0]='\0';
        sprintf(f,"%s%c%s",data_dir.c_str(),FILEPATHSEP,file->d_name);
        kConf.fileC.rover=f;
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"ROVER OBSERVATION FILE: "<<f;

        if(kConf.site_name.empty()){
            string name=file->d_name;
            kConf.site_name=name.substr(0,4);
        }
        kConf.site_idx++;

        CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"==== START PROCESS: "<<kConf.site_name;

        if(kConf.site_idx==1){
            if(!AutoMatchFile(kConf.fileC.rover)){
                kConf.site_name.clear();
                CLOG(INFO,ELPP_CURR_FILE_LOGGER_ID)<<"==== PROCESS FAIL, MISSING FILE: "<<kConf.site_name;
                continue;
            }
        }
        else{
            cMatchFile match_file;
            match_file.InitMatchFile(kConf,FILEPATHSEP);
            match_file.MatchOut();
        }

        long t1=clock();

        solver->InitReader(kConf);
        kConf.gnssC.sample_rate=solver->rover_obs_.GetGnssObs()[1].obs_time.TimeDiff(solver->rover_obs_.GetGnssObs()[0].obs_time.t_);
        solver->SolverProcess(kConf,-1);
        if(single_flag) break;


        kConf.site_name.clear();
        long t2=clock();
        double t=(double)(t2-t1)/CLOCKS_PER_SEC;
        cout<<"total(s): "<<t<<endl;
    }

    return 1;
}

int main(int argc, char** argv)
{
    string conf_file;

    if(!ParsePara(argc,argv,conf_file)){
        PrintHelp();
        return 0;
    }

    Processer();

    return 1;
}
