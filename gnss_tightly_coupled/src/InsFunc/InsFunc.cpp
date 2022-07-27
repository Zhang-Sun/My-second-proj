// realization of class and function in InsFunc.h
// Created by wlzhang on 7/4/22.
// 

#include "InsFunc.h"

extern const double Crf[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};

namespace IGCoupled {
    /*
    * Function : vector to skew symmetric matrix
    * -Args : 
    *        Vector3d vec          I        3-dimention vector
    * -Return:
    *        Matrix3d dcm                   skew symmetric matrix
    */
    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec){   // 向量求反对称矩阵
        Eigen::Matrix3d dcm=Matrix3d::Zero();

        dcm(0,1)=-vec(2);
        dcm(0,2)=vec(1);

        dcm(1,0)=vec(2);
        dcm(1,2)=-vec(0);

        dcm(2,0)=-vec(1);
        dcm(2,1)=vec(0);
        return dcm;
    }
    /*
    * Function : Quanternion to rotation matrix
    * -Args : 
    *        Quaterniond q           I      quaternion
    * -Return : 
    *        Matrix3d dcm                    rotation matrix
    */
    Eigen::Matrix3d Quaternion2RotationMatrix(const Eigen::Quaterniond& q){ // 四元数->旋转矩阵
        return q.toRotationMatrix();
    }
    /*
    * Function : Rotation matrix to quaternion 
    * -Args : 
    *       Matrix3d 
    */
    Eigen::Quaterniond RotationMatrix2Quaternion(const Eigen::Matrix3d& m){ // 旋转矩阵->四元数
        Eigen::Quaterniond q(m);
        return q;
    }

    Eigen::Quaterniond Euler2Quaternion(const Vector3d& rpy){  // 欧拉角->四元数
        Eigen::AngleAxisd roll_angle(Eigen::AngleAxisd(rpy(2),Eigen::Vector3d::UnitZ()));
        Eigen::AngleAxisd pitch_angle(Eigen::AngleAxisd(rpy(1),Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd yaw_angle(Eigen::AngleAxisd(rpy(0),Eigen::Vector3d::UnitX()));

        return roll_angle*pitch_angle*yaw_angle;
    }

    Eigen::Matrix3d Euler2RotationMatrix(const Vector3d& rpy){   // 欧拉角->旋转矩阵
        return Quaternion2RotationMatrix(Euler2Quaternion(rpy));
    }

    template <typename Derived>
    static Eigen::Matrix<Derived,3,1> EulerAngles(const Eigen::Matrix<Derived,3,3>& m){  //
        Eigen::Matrix<Derived, 3, 1> res;

        const size_t i = 2;
        const size_t j = 1;
        const size_t k = 2;
        typedef Eigen::Matrix<Derived, 2, 1> Vector2;

        res[2]=atan2(m(j,k),m(k,k));
        Derived c2=Vector2(m(i,i),m(i,j)).norm();
//        if(res[2]<Derived(0))
//        {
//            res[2]+=Derived(EIGEN_PI)*2;
//        }
        res[1]=atan2(-m(i,k),c2);
        Derived s1=sin(res[2]);
        Derived c1=cos(res[2]);
        res[0]=atan2(s1*m(k,i)-c1*m(j,i),c1*m(j,j)-s1*m(k,j));

        return res;
    }

    Eigen::Vector3d RotationMatrix2Euler(const Matrix3d &m){
//        return m.eulerAngles(2,1,0);
//        return EulerAngles(m);
        Vector3d rpy;
        rpy[0]=atan2(m(1,2),m(2,2));
        rpy[1]=-asin(m(0,2));
        rpy[2]=atan2(m(0,1),m(0,0));
        return rpy;
    }

    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q){
        return RotationMatrix2Euler(q.toRotationMatrix());
    }

    Eigen::Quaterniond RotationVector2Quaternion(const Vector3d& rv){
        Eigen::Quaterniond qfromrv;
        Vector3d rv_2=rv*0.5;
        double norm=rv_2.norm();
        qfromrv.w()=cos(norm);
        qfromrv.vec()=norm<1E-8?rv_2:(sin(norm)/norm)*rv_2;
        return qfromrv;
    }

	Eigen::Matrix4d quaternionleft(const Quaterniond &q) {
        Eigen::Matrix4d ans;
        ans(0, 0)             = q.w();
        ans.block<1, 3>(0, 1) = -q.vec().transpose();
        ans.block<3, 1>(1, 0) = q.vec();
        ans.block<3, 3>(1, 1) = q.w() * Eigen::Matrix3d::Identity() + VectorSkew(q.vec());
        return ans;
    }

     Eigen::Matrix4d quaternionright(const Quaterniond &p) {
        Eigen::Matrix4d ans;
        ans(0, 0)             = p.w();
        ans.block<1, 3>(0, 1) = -p.vec().transpose();
        ans.block<3, 1>(1, 0) = p.vec();
        ans.block<3, 3>(1, 1) = p.w() * Eigen::Matrix3d::Identity() - VectorSkew(p.vec());
        return ans;
    }


    cImuData::cImuData(){}

    cImuData::cImuData(IGCoupled::cTime *ts, IGCoupled::cTime *te){
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cImuData::~cImuData() {data_.clear();}

    void cImuData::SetImu(tInsConf C){
        imu_type_=C.imu_type;
        imu_coord_type_=C.coord_type;
        data_format_=C.data_format;
        gyro_format_=C.gyro_val_format;
        hz_=C.sample_rate;
    }

    void cImuData::SetImuType(IGCoupled::IMU_TYPE type) {imu_type_=type;}

    void cImuData::SetImuCoordType(IGCoupled::IMU_COORD_TYPE type) {imu_coord_type_=type;}

    void cImuData::SetTimeSpan(cTime *ts, cTime *te) {
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cInsMech::cInsMech() {}

    cInsMech::~cInsMech() {}

    void cInsMech::RotScullCorr(IGCoupled::tImuInfoUnit &pre_imu_info, IGCoupled::tImuInfoUnit &cur_imu_info, double dt,double *da,double *dv) {
        Vector3d pre_da,pre_dv,cur_da,cur_dv;
        int i;

        pre_da=pre_imu_info.cor_gyro*cur_imu_info.dt; //前一时刻的角增量
        pre_dv=pre_imu_info.cor_acce*cur_imu_info.dt; //前一时刻的速度增量

        cur_da=cur_imu_info.cor_gyro*dt; // 当前时刻的角增量
        cur_dv=cur_imu_info.cor_acce*dt; // 当前时刻的速度增量

        Vector3d t1,t2,t3,t4;
        CrossVec3(cur_da.data(),cur_dv.data(),t1.data()); // t1 = 当前时刻角增量 * 当前时刻速度增量
        CrossVec3(pre_da.data(),cur_dv.data(),t2.data()); // t2 = 前一时刻的角增量 *当前时刻速度增量
        CrossVec3(pre_dv.data(),cur_da.data(),t3.data()); // t3 = 前一时刻速度增量 * 当前时刻角增量
        CrossVec3(cur_da.data(),t1.data(),t4.data());     // t4 = 当前时刻角增量 * (当前时刻角增量 * 当前时刻速度增量)

        double a1,a2;
        double b=cur_da.norm();           // 当前角增量的模  
        if(fabs(b)>1E-6){
            a1=(1.0-cos(b))/SQR(b);
            a2=1.0/SQR(b)*(1.0-sin(b)/b);
        }
        else{
            a1=0.5-SQR(b)/24.0+SQR(SQR(b))/720.0;
            a2=1.0/6.0-SQR(b)/120.0+SQR(SQR(b))/5040.0;
        }

        for(i=0;i<3&&dv;i++){
            dv[i]=a1*t1[i]+a2*t4[i]+1.0/12.0*(t2[i]+t3[i]);  // 旋转效应补偿项 + 划桨效应补偿项
        }

        if(da){
            CrossVec3(pre_da,cur_da,da);
            for(i=0;i<3;i++) da[i]*=1.0/12.0;         //单子样的旋转矢量的前一周期部分
        }
    }

    Eigen::Quaterniond cInsMech::AttitudeUpdate(IGCoupled::tImuInfoUnit &pre_imu_info,
                                                IGCoupled::tImuInfoUnit &cur_imu_info,double dt,Vector3d da) {
        //当前时刻的角增量   前一时刻的角增量
        Vector3d theta_k(cur_imu_info.cor_gyro*dt),theta_k_1(pre_imu_info.cor_gyro*dt);

        //等效旋转矢量
        Vector3d cur_phi=theta_k+da; //单子样+前一周期
        //等效旋转矢量转换为四元数（即姿态更新的Rbb项）
        Quaterniond quat_bb=RotationVector2Quaternion(cur_phi);

        //地球自转角速度向量
        Vector3d wiee(0,0,-OMGE_GPS);
        //地球自转引起的角增量
        Vector3d zeta=wiee*dt;
        //姿态更新的第一项Ree项
        Quaterniond quat_ee=RotationVector2Quaternion(zeta);
        //前一时刻的姿态矩阵转换为四元数
        Quaterniond quat_k_1=RotationMatrix2Quaternion(pre_imu_info.Cbe.transpose()).conjugate();
        //姿态更新方程
        Quaterniond qbn_k=quat_ee*quat_k_1*quat_bb;

        return qbn_k.normalized(); // 归一化：多次姿态更新后四元数会失去规范化特性，因此每次进行姿态更新后都应该进行归一化操作
    }

    Eigen::Vector3d cInsMech::VelocityUpdate(IGCoupled::tImuInfoUnit &pre_imu_info,
                                             IGCoupled::tImuInfoUnit &cur_imu_info,double dt,Vector3d dv) {
        //地固系下前一时刻的位置和速度
        Vector3d pos=pre_imu_info.re, vel=pre_imu_info.ve;
        // 地球旋转矢量
        Vector3d wiee(0,0,OMGE_GPS);
        // 当前时刻角增量矢量                    前一时刻的角增量矢量
        Vector3d theta_k(cur_imu_info.cor_gyro*dt),theta_k_1(pre_imu_info.cor_gyro*dt);
        // 当前时刻的速度增量+旋转划桨效应       前一时刻的速度增量
        Vector3d vb_k(cur_imu_info.cor_acce*dt+dv),vb_k_1(pre_imu_info.cor_acce*dt);

        //位置地固系-》当地地理坐标系
        Vector3d coord_blh=Xyz2Blh(pos);
        //计算地固系-》导航系的坐标旋转矩阵
        Matrix3d Cen=CalcCen(coord_blh,COORD_NED);
        // 计算n系下的重力加速度并将其转到地固系
        Vector3d ge=Cen.transpose()*CalculateGravity(coord_blh,false);
        Vector3d omgea_n=wiee*2.0;
        // 由重力加速度引起的速度增量在e系下的
        Vector3d delta_gcor=(ge-omgea_n.cross(vel))*dt;

        /*
        // 计算卯酉圈曲率半径和子午圈曲率半径
        Vector3d R_E = RE_WGS84 / sqrt(1.0 - E_2 * SQR(sin(coord_blh[1])));
        Vector3d R_N = RE_WGS84 * (1.0 - E_2) / ((1.0 - E_2 * SQR(sin(coord_blh[1]))) * sqrt(1.0 - E_2 * SQR(sin(coord_blh[1]))));
        // 计算wenn
        Vector3d vn = Cen * vel;
        Vector3d wenn(vn[1] / (R_E + coord_blh[2]), -vn[0] / (R_N + coord_blh[2]), -vn[1] * tan(coord_blh[1]) / (R_E + coord_blh[2]));
        // 计算n系下由重力加速度引起的速度增量（考虑了科氏加速度）
        Vector3d delta_gcor_n = (CalculateGravity(coord_blh, false) - (Cen * omgea_n + wenn).cross(vn)) * dt;
       */
        


        Matrix3d Cee=Matrix3d::Identity()-VectorSkew(wiee*0.5*dt);

//        Vector3d vrot=theta_k.cross(vb_k)*0.5;
//        Vector3d vscul=(theta_k_1.cross(vb_k)+vb_k_1.cross(theta_k))/12.0;

        Quaterniond pre_quat=RotationMatrix2Quaternion(pre_imu_info.Cbe);
        Matrix3d Cbe=pre_quat.toRotationMatrix();
        Vector3d delta_ve=Cee*Cbe*(vb_k);

        // 暂时没加入修改后的delta_gcor
        return (vel+delta_gcor+delta_ve);
    }

    Eigen::Vector3d cInsMech::PositionUpdate(const tImuInfoUnit& pre_imu_info,const Vector3d& cur_vel, double dt) {
        Eigen::Vector3d pos;
        // 前一时刻和当前时刻的速度均值作为平均速度进行位置更新
        pos=(pre_imu_info.ve+cur_vel)*0.5*dt+pre_imu_info.re;
        return pos;
    }

    // INS机械编排
    bool cInsMech::InsMechanization(bool err_model,tImuInfoUnit &pre_imu_info, tImuInfoUnit &cur_imu_info,int idx) {

        for(int i=0;i<3;i++){
            if(isnan(pre_imu_info.re[i])||isnan(pre_imu_info.ve[i])||
               isinf(pre_imu_info.re[i])||isinf(pre_imu_info.ve[i])){
                CLOG(ERROR,ELPP_CURR_FILE_LOGGER_ID)<<cur_imu_info.t_tag.GetTimeStr(4)<<" "<<idx<<" NUMERIC ERROR";
                return false;
            }
        }
        TraceInsMechInfo(pre_imu_info,true,idx);

        double dt=fabs(cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_));
        if(dt>60.0||fabs(dt)<1E-6){
            cur_imu_info.dt=dt;
            cur_imu_info.pt=pre_imu_info.pt;
            CLOG(WARNING,ELPP_CURR_FILE_LOGGER_ID)<<"TIME DIFFERENCE TOO LARGER";
            return false;
        }

        //对imu输出数据进行误差补偿
        if(err_model){ 
            ErrModel(pre_imu_info, cur_imu_info);
        }else{
            cur_imu_info.cor_gyro=cur_imu_info.raw_gyro-pre_imu_info.bg;
            cur_imu_info.cor_acce=cur_imu_info.raw_acce-pre_imu_info.ba;
        }

        // 进行姿态、速度、位置更新
        Vector3d da,dv;
        RotScullCorr(pre_imu_info,cur_imu_info,dt,da.data(),dv.data());
        cur_imu_info.Cbe=Quaternion2RotationMatrix(AttitudeUpdate(pre_imu_info,cur_imu_info,dt,da));
        cur_imu_info.ve=VelocityUpdate(pre_imu_info,cur_imu_info,dt,dv);
        cur_imu_info.re=PositionUpdate(pre_imu_info,cur_imu_info.ve,cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_));

        Vector3d blh=Xyz2Blh(cur_imu_info.re);
        Matrix3d Cen=CalcCen(blh,COORD_NED);
        Matrix3d Cnb=cur_imu_info.Cbe.transpose()*Cen.transpose();
        cur_imu_info.rpy=RotationMatrix2Euler(Cnb);
        cur_imu_info.vn=Cen*cur_imu_info.ve;
        cur_imu_info.rn=Cen*cur_imu_info.re;

        TraceInsMechInfo(cur_imu_info,false,idx);
        cur_imu_info.dt=dt;
        cur_imu_info.pt=pre_imu_info.t_tag;

        return true;
    }

    // 考虑了比例因子和交叉耦合误差的误差模型
    void cInsMech::ErrModel(tImuInfoUnit& pre_imu_info, tImuInfoUnit& cur_imu_info) {
        Matrix3d Mai, Mgi;
        Vector3d T1,T2,Gf;

        Mai = MatrixXd::Identity(3, 3) + pre_imu_info.sa + pre_imu_info.ra;
        Mgi = MatrixXd::Identity(3, 3) + pre_imu_info.sg + pre_imu_info.rg;

        if (!MatInv(Mai.data(), 3) && !MatInv(Mgi.data(), 3)) {
                cur_imu_info.cor_acce = Mai * (cur_imu_info.raw_acce - pre_imu_info.ba);
                cur_imu_info.cor_gyro = Mgi * (cur_imu_info.raw_gyro - pre_imu_info.bg);
        }
    }
    //状态转移矩阵
    Eigen::MatrixXd cInsMech:: StateTransferMat(tIGCOUPLEDConf C,IGCoupled::tImuInfoUnit &pre_imu_info,
                                               IGCoupled::tImuInfoUnit &cur_imu_info,int nx,double dt) {
        using Eigen::Matrix3d;
        using Eigen::MatrixXd;
        using Eigen::Vector3d;

        auto &vel=cur_imu_info.ve;
        auto &fb=cur_imu_info.cor_acce;
        auto &wb=cur_imu_info.cor_gyro;
        Vector3d wiee(0,0,OMGE_GPS);
        auto &Cbe=cur_imu_info.Cbe;

		if(cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_) < 0) wiee = {0,0,-OMGE_GPS};

        MatrixXd F=MatrixXd::Zero(nx,nx);

        int ip=0;
        int iv=3;
        int ia=6;
        int iba=9;
        int ibg=12;
        int isa=0,isg=0,ira=0,irg=0,ilev=0;
        if(C.insC.est_sa) isa=ibg+3;
        else isa=ibg;
        if(C.insC.est_sg) isg=isa+3;
        else isg=isa;
        if(C.insC.est_ra) ira=isg+3;
        else ira=isg;
        if(C.insC.est_rg) irg=ira+6;
        else irg=ira;
        if(C.insC.est_level) ilev=irg+6;
        else ilev=irg;

        //position-velocity  Mpv
        F.block<3,3>(ip,iv)=Matrix3d::Identity();

        //velocity-velocity  Mvv
        F.block<3,3>(iv,iv)=(-2.0*VectorSkew(wiee));
        //velocity-attitude  Mva
        F.block<3,3>(iv,ia)=-VectorSkew(Cbe*fb);
        //velocity-ba
        F.block<3,3>(iv,iba)=-Cbe;

        //attitude-attitude  Maa
        F.block<3,3>(ia,ia)=-1.0*VectorSkew(wiee);
        //attitute-bg
        F.block<3,3>(ia,ibg)=Cbe;
//        cout<<F<<endl;

        //ba-ba
//        F.block<3,3>(iba,iba)=Matrix3d::Identity()*(-fabs(dt)/C.insC.correction_time_ba);
        F.block<3,3>(iba,iba)=Matrix3d::Zero();
        //bg-bg
//        F.block<3,3>(ibg,ibg)=Matrix3d::Identity()*(-fabs(dt)/C.insC.correction_time_bg);
        F.block<3,3>(ibg,ibg)=Matrix3d::Zero();

//        cout<<MatrixXd::Identity(nx,nx)+F*dt<<endl;

        return MatrixXd::Identity(nx,nx)+F*dt;
    }

    Eigen::Vector3d CalculateGravity(const Vector3d coord_blh,bool is_ecef){
        if(is_ecef){
            const double constant_J2=0.00108263;
            const double constant_J4=-2.37091222e-6;
            const double constant_J6=6.08347e-9;
            double p = sqrt(coord_blh(0) * coord_blh(0) + coord_blh(1) * coord_blh(1) + coord_blh(2) * coord_blh(2));
            double t = coord_blh(2) / p;
            double a_p = WGS84_EARTH_LONG_RADIUS/ p;
            double a1 = -WGS84_GM / p / p;
            double a2 = 1 + 1.5 * constant_J2 * a_p * a_p - (15.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (35.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a3 = -4.5 * constant_J2 * a_p * a_p + (75.0 / 4) * constant_J4 * pow(a_p,3) * a_p - (735.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a4 = -(175.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (2205.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a5 = -(1617.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double b1 = 3 * constant_J2 * a_p * a_p - (15.0 / 2) * constant_J4 * pow(a_p,3) * a_p + (105.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b2 = (35.0 / 2) * constant_J4 * pow(a_p,3) * a_p - (945.0 / 12) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b3 = (693.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double c1 = a2;
            double c2 = a3 - b1;
            double c3 = a4 - b2;
            double c4 = a5 - b3;
            double d1 = a2 + b1;
            double d2 = c2 + b2;
            double d3 = c3 + b3;
            double d4 = c4;
            Vector3d ge_vec;
            ge_vec(0) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(0) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(0);
            ge_vec(1) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(1) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(1);
            ge_vec(2) = (d1 + d2 * t * t + d3 * pow(t,3) * t + d4 * pow(t,3) * pow(t,3)) * coord_blh(2) * a1 / p;
            return ge_vec;
        }
        else{
            double gn = 9.7803267715 * (1 + 0.0052790414 * sin(coord_blh(0)) * sin(coord_blh(0)) + 0.0000232719 * pow(sin(coord_blh(0)),3) * sin(coord_blh(0)));
            gn += (-0.0000030876910891 + 0.0000000043977311 * sin(coord_blh(0)) * sin(coord_blh(0))) * coord_blh(2);
            gn += 0.0000000000007211 * coord_blh(2) * coord_blh(2);
            Vector3d gn_vec{0, 0, gn};
            return gn_vec;
        }
    }

    void cInsMech::TraceInsMechInfo(IGCoupled::tImuInfoUnit &imu_info,bool prior,int idx) {
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"INS MECHANIZATION"<<(prior?"- ":"+ ")<< "("<<idx<<"): "<<imu_info.t_tag.GetTimeStr(4);
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"GYRO VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_gyro.transpose()<<" rad/s";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"ACCE VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_acce.transpose()<<" m/s^2";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"ATTITUDE:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rpy.transpose()*R2D<<" deg";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"VELOCITY:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.vn.transpose()<<" m/s";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"POSITION:(e)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.re.transpose()<<" m";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"GYRO BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.bg.transpose()<<" rad/s";
        CLOG(DEBUG,ELPP_CURR_FILE_LOGGER_ID)<<"   "<<"ACCE BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.ba.transpose()<<" m/s^2";
    }

    cInsAlign::cInsAlign() {}

    cInsAlign::cInsAlign(cImuData imu_data,tIGCOUPLEDConf C) {
        imu_data_=&imu_data;
        C_=C;
    }

    cInsAlign::~cInsAlign() {}

    bool cInsAlign::CoarseAlign(tImuInfoUnit &imu_info) {
        Eigen::Matrix3d A = Matrix3d::Zero();
        Matrix3d B = Matrix3d::Zero();

        Vector3d pos = imu_info.re;
        Vector3d coord_blh = Xyz2Blh(pos);
        Matrix3d Cen = CalcCen(coord_blh, COORD_NED);
        Vector3d gn = CalculateGravity(coord_blh, false);
        Vector3d ge = Cen.transpose() * gn;
        double g = Norm(gn, 3);

        /* 误差改正 */
        imu_info.cor_gyro = imu_info.raw_gyro - imu_info.bg;
        imu_info.cor_acce = imu_info.raw_acce - imu_info.ba;

        /* 获取初始姿态*/
        A(0, 0) = -tan(coord_blh[0]) / g;
        A(0, 1) = 1.0 / (OMGE_GPS * cos(coord_blh[0]));

        A(1, 2) = -1.0 / (g * OMGE_GPS * cos(coord_blh[0]));

        A(2, 0) = -1.0 / g;

        CrossVec3(imu_info.cor_acce, imu_info.cor_gyro, gn);
        for (int i = 0; i < 3; i++) {
            B(0, i) = imu_info.cor_acce[i];
            B(1, i) = imu_info.cor_gyro[i];
            B(2, i) = gn[i];
        }

        imu_info.Cbn = A * B;
        imu_info.Cbe = Cen.transpose() * imu_info.Cbn;
        imu_info.Cbn.normalized();
        imu_info.Cbe.normalized();

    }

	cPreIntegrationBase::cPreIntegrationBase(tIGCOUPLEDConf C, tImuDataUnit &imu, tFGOState state)
    : C_(std::move(C)), current_state_(std::move(state)){
        start_time_ = imu.t_tag;
        end_time_ = imu.t_tag;

        imu_buffer_.clear();
        imu_buffer_.push_back(imu);


        gravity_ = CalculateGravity(state.p,true);
    }

    void cPreIntegrationBase::integration(tImuDataUnit &imu_pre, tImuDataUnit &imu_cur) {
        //区间时间累积、
        double dt = fabs(imu_cur.t_tag.TimeDiff(imu_pre.t_tag.t_));
        delta_time_ += dt;

        end_time_ = imu_cur.t_tag;
        current_state_.t_tag = imu_cur.t_tag;

        //连续状态积分。先位置速度，再姿态
        //位置速度
        Vector3d vb_k(imu_cur.acce * dt), theta_k(imu_cur.gyro *dt);
        Vector3d vb_k_1(imu_pre.acce * dt), theta_k_1(imu_pre.gyro * dt);

        Vector3d dvfb = vb_k + 0.5 * theta_k.cross(vb_k) +
                1.0 / 12.0 * (theta_k_1.cross(vb_k) + vb_k + vb_k_1.cross(theta_k));

        Vector3d dvel = current_state_.q.toRotationMatrix() * dvfb + gravity_ * dt; // 速度增量
        current_state_.p += dt * current_state_.v + 0.5 * dt * dvel;
        current_state_.v += dvel;

        //姿态

        Vector3d dtheta = theta_k + 1.0 / 12.0 * theta_k_1.cross(theta_k);// 等效旋转矢量
        current_state_.q *= RotationVector2Quaternion(dtheta);
        current_state_.q.normalize();

        //预积分
        dvel = delta_state_.q.toRotationMatrix() * dvfb;
        delta_state_.p += dt * delta_state_.v + 0.5 * dt * dvel;
        delta_state_.v += dvel;

        //姿态
        delta_state_.q *= RotationVector2Quaternion(dtheta);
        delta_state_.q.normalize();
    }

    void cPreIntegrationBase::addNewImu(const tImuDataUnit &imu) {
        imu_buffer_.push_back(imu);
        integrationProcess(imu_buffer_.size() - 1);
    }

    void cPreIntegrationBase::reintegration(tFGOState &state){
        current_state_ = std::move(state);
        resetState(current_state_);

        for(size_t k = 1; k < imu_buffer_.size(); k++){
            integrationProcess(k);
        }
    }

    tImuDataUnit cPreIntegrationBase::compensationBias(const tImuDataUnit &imu) const {
        tImuDataUnit imu_calib = imu;
        imu_calib.gyro -= delta_state_.bg;
        imu_calib.acce -= delta_state_.ba;

        return imu_calib;
    }

    tImuDataUnit cPreIntegrationBase::compensationScale(const tImuDataUnit &imu) const {
        tImuDataUnit imu_calib = imu;

        for(int k = 0; k < 3; k++){
            imu_calib.gyro[k] *= (1.0 - delta_state_.sg[k]);
            imu_calib.acce[k] *= (1.0 - delta_state_.sa[k]);
        }
        return imu_calib;
    }

    void cPreIntegrationBase::stateToData(const tFGOState &state, tIntegrationStateData &data){
        data.t_tag = state.t_tag;

        memcpy(data.pose, state.p.data(), sizeof(double) * 3);
        memcpy(data.pose + 3, state.q.coeffs().data(), sizeof(double) * 4);

        memcpy(data.mix, state.v.data(), sizeof(double) * 3);
        memcpy(data.mix + 3, state.bg.data(), sizeof(double) * 3);
        memcpy(data.mix + 6, state.ba.data(), sizeof(double) * 3);
    }

    void cPreIntegrationBase::stateFromData(const tIntegrationStateData &data, tFGOState &state) {
        state.t_tag = data.t_tag;

        memcpy(state.p.data(), data.pose, sizeof(double) * 3);
        memcpy(state.q.coeffs().data(), data.pose + 3, sizeof(double) * 4);
        state.q.normalize();

        memcpy(state.v.data(), data.mix, sizeof(double) * 3);
        memcpy(state.bg.data(), data.mix + 3, sizeof(double) * 3);
        memcpy(state.ba.data(), data.mix + 6, sizeof(double) * 3);
    }

    cPreIntegrationEarthOdo::cPreIntegrationEarthOdo(tIGCOUPLEDConf C, tImuDataUnit &imu, tFGOState state)
    : cPreIntegrationBase(C, imu, std::move(state)){

        //Reset state
        if(C_.insC.use_odo) resetState(state, 19);
        else resetState(state, 15);

        //set initial noise matrix
        setNoiseMatrix();

        //里程计参数
        if(C_.insC.use_odo)
        {
            cvb_ = Euler2RotationMatrix(C.insC.abv).transpose();
            lodo_ = C.insC.lodo;
        }

    }

    MatrixXd cPreIntegrationEarthOdo::evaluate(const tFGOState &state0, const tFGOState &state1, double *residuals) {

        sqrt_information_ = Eigen::LLT<Eigen::Matrix<double, NUM_STATE, NUM_STATE>>(covariance_.inverse()).matrixL().transpose();
        Eigen::Map<Eigen::Matrix<double, NUM_STATE, 1>> residual(residuals);
        if(C_.insC.use_odo){
            sqrt_information_ = Eigen::LLT<Eigen::Matrix<double, NUM_STATE_ODO, NUM_STATE_ODO>>(covariance_.inverse()).matrixL().transpose();
            Eigen::Map<Eigen::Matrix<double, NUM_STATE_ODO, 1>> residual(residuals);
        }

        Matrix3d dp_dbg   = jacobian_.block<3, 3>(0, 9);
        Matrix3d dp_dba   = jacobian_.block<3, 3>(0, 12);
        Matrix3d dv_dbg   = jacobian_.block<3, 3>(3, 9);
        Matrix3d dv_dba   = jacobian_.block<3, 3>(3, 12);
        Matrix3d dq_dbg   = jacobian_.block<3, 3>(6, 9);
        Vector3d ds_dsodo;
        Matrix3d ds_dbg;
        if(C_.insC.use_odo){
            ds_dsodo = jacobian_.block<3, 1>(15, 18);
            ds_dbg   = jacobian_.block<3, 3>(15, 9);
        }

        //零偏误差
        Vector3d dbg = state0.bg - delta_state_.bg;
        Vector3d dba = state0.ba - delta_state_.ba;
        double dsodo;
        if(C_.insC.use_odo) {
            dsodo = state0.sodo - delta_state_.sodo;
        }

        //位置补偿项
        Vector3d p_cor{0,0,0};
        for(const auto &pe : pe_){
            p_cor += (pe.second - state0.p) * pe.first;
        }
        p_cor = 2.0 * iewe_skew_ * p_cor;

        //速度补偿
        Vector3d v_cor;
        v_cor = 2.0 * iewe_skew_ * (state1.p - state0.p);

        // 姿态
        Vector3d dee    = -iewe_ * delta_time_;
        Quaterniond qee = RotationVector2Quaternion(dee);

        dpe_ = state1.p - state0.p - state0.v * delta_time_ - 0.5 * gravity_ * delta_time_ * delta_time_ + p_cor;
        dve_ = state1.v - state0.v - gravity_ * delta_time_ + v_cor;

        // 积分校正
        corrected_p_ = delta_state_.p + dp_dba * dba + dp_dbg * dbg;
        corrected_v_ = delta_state_.v + dv_dba * dba + dv_dbg * dbg;
        corrected_q_ = delta_state_.q * RotationVector2Quaternion(dq_dbg * dbg);

        if(C_.insC.use_odo) corrected_s_ = delta_state_.s + ds_dbg * dbg + ds_dsodo * dsodo;

        Quaterniond qeb0 = state0.q.inverse();
        Matrix3d ceb0    = qeb0.toRotationMatrix();
        qb0b1_           = state1.q.inverse() * qee * state0.q;

        // Residuals
        residual.block<3, 1>(0, 0)  = ceb0 * dpe_ - corrected_p_;
        residual.block<3, 1>(3, 0)  = ceb0 * dve_ - corrected_v_;
        residual.block<3, 1>(6, 0)  = 2 * (qb0b1_ * corrected_q_).vec();
        residual.block<3, 1>(9, 0)  = state1.bg - state0.bg;
        residual.block<3, 1>(12, 0) = state1.ba - state0.ba;

        if(C_.insC.use_odo){
            residual.block<3, 1>(15, 0) = ceb0 * (state1.p - state0.p) - corrected_s_;
            residual(18)                = state1.sodo - state0.sodo;
        }

        residual = sqrt_information_ * residual;

        return residual;
    }

    MatrixXd cPreIntegrationEarthOdo::residualJacobianPose0(const tFGOState &state0, const tFGOState &state1,
                                                            double *jacobian) {
        Eigen::Map<Eigen::Matrix<double, NUM_STATE, NUM_POSE, Eigen::RowMajor>> jaco(jacobian);

        if(C_.insC.use_odo){
            Eigen::Map<Eigen::Matrix<double, NUM_STATE_ODO, NUM_POSE, Eigen::RowMajor>> jaco(jacobian);
        }

        jaco.setZero();

        // 地球自转的位置补偿项
        double dt = 0;
        for (const auto &pn : pe_) {
            dt += pn.first;
        }

        Matrix3d ceb0 = state0.q.inverse().toRotationMatrix();

        jaco.block(0, 0, 3, 3) = -ceb0 - 2.0 * ceb0 * iewe_skew_ * dt;
        jaco.block(0, 3, 3, 3) = VectorSkew(ceb0 * dpe_);
        jaco.block(3, 0, 3, 3) = -2.0 * ceb0 * iewe_skew_;
        jaco.block(3, 3, 3, 3) = VectorSkew(ceb0 * dve_);
        jaco.block(6, 3, 3, 3) =
                (quaternionleft(qb0b1_) * quaternionright(corrected_q_)).bottomRightCorner<3, 3>();

        if(C_.insC.use_odo) {
            jaco.block(15, 0, 3, 3) = -ceb0;
            jaco.block(15, 3, 3, 3) = VectorSkew(ceb0 * (state1.p - state0.p));
        }
        jaco = sqrt_information_ * jaco;
        return jaco;
    }

    MatrixXd cPreIntegrationEarthOdo::residualJacobianPose1(const tFGOState &state0,
                                                           const tFGOState &state1, double *jacobian) {
        Eigen::Map<Eigen::Matrix<double, NUM_STATE, NUM_POSE, Eigen::RowMajor>> jaco(jacobian);
        if(C_.insC.use_odo) Eigen::Map<Eigen::Matrix<double, NUM_STATE_ODO, NUM_POSE, Eigen::RowMajor>> jaco(jacobian);
        jaco.setZero();

        Matrix3d ceb0 = state0.q.inverse().toRotationMatrix();

        jaco.block(0, 0, 3, 3)  = ceb0;
        jaco.block(3, 0, 3, 3)  = 2.0 * ceb0 * iewe_skew_;
        jaco.block(6, 3, 3, 3)  = -quaternionleft(qb0b1_ * corrected_q_).bottomRightCorner<3, 3>();

        if(C_.insC.use_odo) jaco.block(15, 0, 3, 3) = ceb0;

        jaco = sqrt_information_ * jaco;
        return jaco;
    }

    MatrixXd cPreIntegrationEarthOdo::residualJacobianMix0(const tFGOState &state0,
                                                       const tFGOState &state1, double *jacobian) {
        Eigen::Map<Eigen::Matrix<double, NUM_STATE, NUM_MIX, Eigen::RowMajor>> jaco(jacobian);
        if(C_.insC.use_odo) {
            Eigen::Map<Eigen::Matrix<double, NUM_STATE_ODO, NUM_MIX_ODO, Eigen::RowMajor>> jaco(jacobian);
        }
        jaco.setZero();

        Eigen::Matrix3d dp_dbg = jacobian_.block<3, 3>(0, 9);
        Eigen::Matrix3d dp_dba = jacobian_.block<3, 3>(0, 12);
        Eigen::Matrix3d dv_dbg = jacobian_.block<3, 3>(3, 9);
        Eigen::Matrix3d dv_dba = jacobian_.block<3, 3>(3, 12);
        Eigen::Matrix3d dq_dbg = jacobian_.block<3, 3>(6, 9);
        if(C_.insC.use_odo) {
            Vector3d ds_dsodo = jacobian_.block<3, 1>(15, 18);
            Matrix3d ds_dbg   = jacobian_.block<3, 3>(15, 9);

            jaco.block(15, 3, 3, 3) = -ds_dbg;
            jaco.block(15, 9, 3, 1) = -ds_dsodo;
            jaco(18, 9)             = -1.0;
        }


        Matrix3d ceb0 = state0.q.inverse().toRotationMatrix();

        jaco.block(0, 0, 3, 3)  = -ceb0 * delta_time_;
        jaco.block(0, 3, 3, 3)  = -dp_dbg;
        jaco.block(0, 6, 3, 3)  = -dp_dba;
        jaco.block(3, 0, 3, 3)  = -ceb0;
        jaco.block(3, 3, 3, 3)  = -dv_dbg;
        jaco.block(3, 6, 3, 3)  = -dv_dba;
        jaco.block(6, 3, 3, 3)  = quaternionleft(qb0b1_ * delta_state_.q).bottomRightCorner<3, 3>() * dq_dbg;
        jaco.block(9, 3, 3, 3)  = -Eigen::Matrix3d::Identity();
        jaco.block(12, 6, 3, 3) = -Eigen::Matrix3d::Identity();


        jaco = sqrt_information_ * jaco;
        return jaco;
    }


    MatrixXd cPreIntegrationEarthOdo::residualJacobianMix1(const tFGOState &state0,
                                                          const tFGOState &state1, double *jacobian) {
        Eigen::Map<Eigen::Matrix<double, NUM_STATE, NUM_MIX, Eigen::RowMajor>> jaco(jacobian);
        if(C_.insC.use_odo) Eigen::Map<Eigen::Matrix<double, NUM_STATE_ODO, NUM_MIX_ODO, Eigen::RowMajor>> jaco(jacobian);
        jaco.setZero();

        jaco.block(3, 0, 3, 3)  = state0.q.inverse().toRotationMatrix();
        jaco.block(9, 3, 3, 3)  = Eigen::Matrix3d::Identity();
        jaco.block(12, 6, 3, 3) = Eigen::Matrix3d::Identity();
        if(C_.insC.use_odo) jaco(18, 9)             = 1.0;

        jaco = sqrt_information_ * jaco;
        return jaco;
    }

    int cPreIntegrationEarthOdo::numResiduals() {
        if(!C_.insC.use_odo) return NUM_STATE;
        return NUM_STATE_ODO;
    }

    vector<int> cPreIntegrationEarthOdo::numBlocksParameters() {
        if(C_.insC.use_odo) return std::vector<int>{NUM_POSE, NUM_MIX_ODO, NUM_POSE, NUM_MIX_ODO};
        return std::vector<int>{NUM_POSE, NUM_MIX, NUM_POSE, NUM_MIX};
    }

    tIntegrationStateData cPreIntegrationEarthOdo::stateToData( const tFGOState &state) {
        tIntegrationStateData data;
        cPreIntegrationBase::stateToData(state, data);
        data.mix[9] = state.sodo;

        return data;
    }

    tFGOState cPreIntegrationEarthOdo::stateFromData( const tIntegrationStateData &data) {
        tFGOState state;
        cPreIntegrationBase::stateFromData(data, state);
        state.sodo = data.mix[9];

        return state;
    }

    void cPreIntegrationEarthOdo::constructState(const double *const *parameters, tFGOState &state0,
                                                 tFGOState &state1) {
        state0 = tFGOState {
                .p    = {parameters[0][0], parameters[0][1], parameters[0][2]},
                .q    = {parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]},
                .v    = {parameters[1][0], parameters[1][1], parameters[1][2]},
                .bg   = {parameters[1][3], parameters[1][4], parameters[1][5]},
                .ba   = {parameters[1][6], parameters[1][7], parameters[1][8]},
                .sodo = parameters[1][9],
        };

        state1 = tFGOState {
                .p    = {parameters[2][0], parameters[2][1], parameters[2][2]},
                .q    = {parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]},
                .v    = {parameters[3][0], parameters[3][1], parameters[3][2]},
                .bg   = {parameters[3][3], parameters[3][4], parameters[3][5]},
                .ba   = {parameters[3][6], parameters[3][7], parameters[3][8]},
                .sodo = parameters[3][9],
        };
    }

    void cPreIntegrationEarthOdo::integrationProcess(unsigned long index) {
        tImuDataUnit imu_pre = compensationBias(imu_buffer_[index - 1]);
        tImuDataUnit imu_cur = compensationBias(imu_buffer_[index]);

        // 区间时间累积
        double dt = fabs(imu_cur.t_tag.TimeDiff(imu_pre.t_tag.t_));
        delta_time_ += dt;

        end_time_           = imu_cur.t_tag;
        current_state_.t_tag = imu_cur.t_tag;

        // 连续状态积分, 先位置速度再姿态

        // 位置速度
        Vector3d vb_k(imu_cur.acce * dt), theta_k(imu_cur.gyro *dt);
        Vector3d vb_k_1(imu_pre.acce * dt), theta_k_1(imu_pre.gyro * dt);
        Vector3d dvfb = vb_k + 0.5 * theta_k.cross(vb_k) +
                        1.0 / 12.0 * (theta_k_1.cross(vb_k) + vb_k_1.cross(theta_k));
        // 哥氏项和重力项
        Vector3d coord_blh = Xyz2Blh(current_state_.p);
        Matrix3d Cen = CalcCen(coord_blh,COORD_NED);
        Vector3d gravity = Cen.transpose()* CalculateGravity(coord_blh,false);
        Vector3d dv_cor_g = (gravity - 2.0 * iewe_.cross(current_state_.v)) * dt;

        // 地球自转补偿项, 省去了enwn项
        Vector3d dee    = -iewe_ * dt;
        Quaterniond qee = RotationVector2Quaternion(dee);

        Vector3d dvel =
                0.5 * (Matrix3d::Identity() + qee.toRotationMatrix()) * current_state_.q.toRotationMatrix() * dvfb + dv_cor_g; //速度增量

        // 前后历元平均速度计算位置
        current_state_.p += dt * current_state_.v + 0.5 * dt * dvel;
        current_state_.v += dvel;

        // 缓存IMU时刻位置, 时间间隔为两个历元的间隔
        pe_.emplace_back(std::make_pair(dt, current_state_.p));

        // 姿态
        Vector3d dtheta = theta_k + 1.0 / 12.0 * theta_k_1.cross(theta_k);

        current_state_.q = qee * current_state_.q * RotationVector2Quaternion(dtheta);
        current_state_.q.normalize();

        // 预积分

        // 中间时刻的地球自转等效旋转矢量
        dee           = -(delta_time_ - 0.5 * dt) * iewe_;
        Matrix3d cbbe = (q0_.inverse() * RotationVector2Quaternion(dee) * q0_ * delta_state_.q).toRotationMatrix();

        // 里程增量

        if(C_.insC.use_odo){
            Vector3d dsodo = Vector3d(imu_cur.odo.vr[0], 0, 0); // TODO：需要理清这里的vr是什么
            delta_state_.s += cbbe * (cvb_ * dsodo * (1 + delta_state_.sodo) -
                                      RotationVector2Quaternion(theta_k).toRotationMatrix() * lodo_ + lodo_);
        }

        // 前后历元平均速度计算位置
        dvel = cbbe * dvfb;
        delta_state_.p += dt * delta_state_.v + 0.5 * dt * dvel;
        delta_state_.v += dvel;

        // 姿态
        delta_state_.q *= RotationVector2Quaternion(dtheta);
        delta_state_.q.normalize();

        // 更新系统状态雅克比和协方差矩阵
        updateJacobianAndCovariance(imu_pre, imu_cur);
    }


    void cPreIntegrationEarthOdo::resetState(const tFGOState &state) {
        resetState(state, NUM_STATE);
        if(C_.insC.use_odo) resetState(state, NUM_STATE_ODO);
    }

    void cPreIntegrationEarthOdo::updateJacobianAndCovariance(const tImuDataUnit &imu_pre, const tImuDataUnit &imu_cur) {
        // dp, dv, dq, dbg, dba

        Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(NUM_STATE, NUM_STATE);
        if(C_.insC.use_odo) phi = Eigen::MatrixXd::Zero(NUM_STATE_ODO, NUM_STATE_ODO);

        double dt = fabs(difftime(imu_cur.t_tag.t_.long_time,imu_pre.t_tag.t_.long_time) + imu_cur.t_tag.t_.sec - imu_pre.t_tag.t_.sec);


        Vector3d dee  = -iewe_ * delta_time_;
        Matrix3d cbb0 = -(q0_.inverse() * RotationVector2Quaternion(dee) * q0_ * delta_state_.q).toRotationMatrix();

        // jacobian

        Vector3d vb_k(imu_cur.acce * dt), theta_k(imu_cur.gyro *dt);
        Vector3d vb_k_1(imu_pre.acce * dt), theta_k_1(imu_pre.gyro * dt);
        // phi = I + F * dt
        phi.block<3, 3>(0, 0)   = Matrix3d::Identity();
        phi.block<3, 3>(0, 3)   = Matrix3d::Identity() * dt;
        phi.block<3, 3>(3, 3)   = Matrix3d::Identity();
        phi.block<3, 3>(3, 6)   = cbb0 * VectorSkew(vb_k);
        phi.block<3, 3>(3, 12)  = cbb0 * dt;
        phi.block<3, 3>(6, 6)   = Matrix3d::Identity() - VectorSkew(theta_k);
        phi.block<3, 3>(6, 9)   = -Matrix3d::Identity() * dt;
        phi.block<3, 3>(9, 9)   = Matrix3d::Identity() * (1 - dt / C_.insC.correction_time_bg);
        phi.block<3, 3>(12, 12) = Matrix3d::Identity() * (1 - dt / C_.insC.correction_time_ba);

        if(C_.insC.use_odo){
            Vector3d dsodo  = Vector3d(imu_cur.odo.vr[0], imu_cur.odo.vr[1], imu_cur.odo.vr[2]); // TODO: Vr 代表什么
            Vector3d stheta = cvb_ * dsodo * (1 + delta_state_.sodo) - theta_k.cross(lodo_);

            phi.block<3, 3>(15, 6)  = -cbb0 * VectorSkew(stheta);
            phi.block<3, 3>(15, 9)  = -cbb0 * VectorSkew(lodo_) * dt;
            phi.block<3, 3>(15, 15) = Matrix3d::Identity();
            phi.block<3, 1>(15, 18) = delta_state_.q.toRotationMatrix() * cvb_ * dsodo;
            phi(18, 18)             = 1.0;
        }


        jacobian_ = phi * jacobian_;

        // covariance

        Eigen::MatrixXd gt = Eigen::MatrixXd::Zero(NUM_STATE, NUM_NOISE);

        if(C_.insC.use_odo) gt = Eigen::MatrixXd::Zero(NUM_STATE_ODO, NUM_NOISE_ODO);

        gt.block<3, 3>(3, 3)   = cbb0;
        gt.block<3, 3>(6, 0)   = -Matrix3d::Identity();
        gt.block<3, 3>(9, 6)   = Matrix3d::Identity();
        gt.block<3, 3>(12, 9)  = Matrix3d::Identity();
        if(C_.insC.use_odo){
            gt.block<3, 3>(15, 0)  = cbb0 * VectorSkew(lodo_);
            gt.block<3, 3>(15, 12) = cbb0 * cvb_ * (1 + delta_state_.sodo);
            gt(18, 15)             = 1.0;
        }


        Eigen::MatrixXd Qk =
                0.5 * dt * (phi * gt * noise_ * gt.transpose() + gt * noise_ * gt.transpose() * phi.transpose());
        covariance_ = phi * covariance_ * phi.transpose() + Qk;
    }

    void cPreIntegrationEarthOdo::resetState(const tFGOState &state, int num) {
        delta_time_ = 0;
        delta_state_.p.setZero();
        delta_state_.q.setIdentity();
        delta_state_.v.setZero();
        delta_state_.s.setZero();
        delta_state_.bg   = state.bg;
        delta_state_.ba   = state.ba;
        delta_state_.sodo = state.sodo;

        jacobian_.setIdentity(num, num);
        covariance_.setZero(num, num);

        // 预积分起点的绝对姿态
        q0_ = current_state_.q;

        // 地球自转, 近似使用初始时刻位置
        iewe_ = {0,0,OMGE_GPS};
        if(C_.filter_type == FILTER_BACKWARD){
            iewe_ = {0,0,-OMGE_GPS};
        }
        iewe_skew_ = VectorSkew(iewe_);

        pe_.clear();
    }

    void cPreIntegrationEarthOdo::setNoiseMatrix() {
        noise_.setIdentity(NUM_NOISE, NUM_NOISE);
        if(C_.insC.use_odo) noise_.setIdentity(NUM_NOISE_ODO, NUM_NOISE_ODO);
        noise_.block<3, 3>(0, 0) *= C_.insC.psd_gyro; // nw
        noise_.block<3, 3>(3, 3) *= C_.insC.psd_acce; // na
        noise_.block<3, 3>(6, 6) *= C_.insC.psd_bg; // nbg
        noise_.block<3, 3>(9, 9) *= C_.insC.psd_ba; // nba

        if(C_.insC.use_odo){
            noise_(12, 12) *= C_.insC.odo_std[0] * C_.insC.odo_std[0];                    // nodo
            noise_(13, 13) *= C_.insC.odo_std[1] * C_.insC.odo_std[1];                    // nodo
            noise_(14, 14) *= C_.insC.odo_std[2] * C_.insC.odo_std[2];                    // nodo
            noise_(15, 15) *= C_.insC.odo_srw * C_.insC.odo_srw;  // nsodo
        }

    }

    int cPreIntegrationEarthOdo::imuErrorNumResiduals() {
        if(C_.insC.use_odo) return 7;
        return 6;
    }

    vector<int> cPreIntegrationEarthOdo::imuErrorNumBlocksParameters() {
        if(C_.insC.use_odo) return std::vector<int>{NUM_MIX_ODO};
        return std::vector<int>{NUM_MIX};
    }

    void cPreIntegrationEarthOdo::imuErrorEvaluate(const double *const *parameters, double *residuals) {
        // bg, ba
        residuals[0] = parameters[0][3] / IMU_GRY_BIAS_STD;
        residuals[1] = parameters[0][4] / IMU_GRY_BIAS_STD;
        residuals[2] = parameters[0][5] / IMU_GRY_BIAS_STD;
        residuals[3] = parameters[0][6] / IMU_ACC_BIAS_STD;
        residuals[4] = parameters[0][7] / IMU_ACC_BIAS_STD;
        residuals[5] = parameters[0][8] / IMU_ACC_BIAS_STD;
        if(C_.insC.use_odo){
            residuals[6] = parameters[0][9] / ODO_SCALE_STD;
        }
    }

    void cPreIntegrationEarthOdo::imuErrorJacobian(double *jacobian) {
        Eigen::Map<Eigen::Matrix<double, 6, NUM_MIX, Eigen::RowMajor>> jaco(jacobian);
        if(C_.insC.use_odo) Eigen::Map<Eigen::Matrix<double, 7, NUM_MIX_ODO, Eigen::RowMajor>> jaco(jacobian);

        jaco.setZero();

        jaco(0, 3) = 1.0 / IMU_GRY_BIAS_STD;
        jaco(1, 4) = 1.0 / IMU_GRY_BIAS_STD;
        jaco(2, 5) = 1.0 / IMU_GRY_BIAS_STD;
        jaco(3, 6) = 1.0 / IMU_ACC_BIAS_STD;
        jaco(4, 7) = 1.0 / IMU_ACC_BIAS_STD;
        jaco(5, 8) = 1.0 / IMU_ACC_BIAS_STD;
        if(C_.insC.use_odo){
            jaco(6, 9) = 1.0 / ODO_SCALE_STD;
        }
    }

    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt) {
        Vector3d gyro,acce;

        if(coord_type==IMU_COORD_RFU){
            gyro=imu_data.gyro;
            acce=imu_data.acce;
            MatMul("NN",3,1,3,1.0,Crf,gyro.data(),0.0,imu_data.gyro.data());
            MatMul("NN",3,1,3,1.0,Crf,acce.data(),0.0,imu_data.acce.data());
        }
        if(data_format==IMU_FORMAT_INCR){
            for(int j=0;j<3;j++) imu_data.acce[j]/=dt;
            for(int j=0;j<3;j++) imu_data.gyro[j]/=dt;
        }
        if(gyro_val_format==GYRO_FORMAT_DEG){
            for(int j=0;j<3;j++) imu_data.gyro[j]*=D2R;
        }
    }

}

