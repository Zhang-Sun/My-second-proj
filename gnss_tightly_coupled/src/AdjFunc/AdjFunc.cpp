//
// Created by wlzhang on 7/4/22.
// 

#include "AdjFunc.h"

#define QC

namespace IGCoupled {

    /*
    * Function : The constructor and destructor function of class cAdjuster
    */
    cAdjuster::cAdjuster() {}

    cAdjuster::~cAdjuster() {}

    int cAdjuster::Adjustment(VectorXd L, const MatrixXd H, const MatrixXd R, VectorXd &X, MatrixXd &Px, int nl, int nx) {}

    /*
    * Function : Remove the empty rows in matrix
    * -Args :
    *        MatrixXd H                      IO           matrix to be precessed
    * -Returns :
    *        None 
    */
    void cAdjuster::RemoveEmptyColumns(MatrixXd &H) {
        MatrixXd H_=H;
        vector<int>zip_idx;
        int i,j;
        for(i=0;i<H_.rows();i++){
            if(H_.row(i).sum()==0.0) continue;//
            zip_idx.push_back(i);//
        }

        MatrixXd HH;
        HH=MatrixXd::Zero(zip_idx.size(),H_.cols());

        for(i=0;i<zip_idx.size();i++){ //
            for(j=0;j<H_.cols();j++) HH(i,j)=H(zip_idx[i],j);
        }
        H=HH;
#if 0
        cout<<H.transpose()<<endl<<endl;
        cout<<HH.transpose()<<endl<<endl;
#endif
        zip_idx.clear();
    }
    /*
    * Function : The constructor and destructor function of class cLsqAdjuster
    */
    cLsqAdjuster::cLsqAdjuster() {}

    cLsqAdjuster::~cLsqAdjuster() {}
    /*
    * Function : Least square process
    * -Args :
    *        VectorXd L                     I             Observation vector 
    *        MatrixXd H                     I             Design matrix 
    *        MatrixXd R                     I             Cofactor matrix
    *        VectorXd X                     IO            Parameters to be estimated 
    *        MatrixXd Px                                  Weight matrix (not used)
    *        int      nl                    I             The number of observation
    *        int      nx                    I             The number of paremeters to be estimated 
    * -Returns :
    *        true (for success)
    */
    int cLsqAdjuster::Adjustment(VectorXd L, const MatrixXd H, const MatrixXd R, VectorXd &X, MatrixXd &Px, int nl, int nx) {
        bool flag=true;

        Eigen::MatrixXd W=R.inverse();//
        Eigen::MatrixXd HTWH=H*W*H.transpose();//
        Eigen::MatrixXd HTWL=H*W*L;// A^TPL

#if 0
        cout<<X.transpose()<<endl;
        cout<<H.transpose()<<endl;
#endif

        dx_=HTWH.inverse()*HTWL;//
        for(int i=0;i<nx;i++) X[i]+=dx_[i];//
        v_=H.transpose()*dx_-L;//
        unit_weight_STD_=(v_.transpose()*W*v_);//
        unit_weight_STD_=SQRT(unit_weight_STD_/((nl>nx)?nl-nx:nl));//

        flag=dx_.norm()<1E-4;

        return flag;
    }
    /*
    *  Function :°°The constructor and destructor function of class cKfAdjuster (Kalman filter)
    */
    cKfAdjuster::cKfAdjuster() {}

    cKfAdjuster::~cKfAdjuster() {}
    /*
    * Function : Kalman filter (observation update)
    * -Args :
    *        VectorXd L                   I         observation 
    *        MatrixXd H                   I         design matrix
    *        MatrixXd R                   I         Measurement noise covariance matrix 
    *        VectorXd X                   IO        parameter to be estimated
    *        MatrixXd Px                  IO        Error and covariance matrix 
    *        int      nl                  I         The number of observation
    *        int      nx                  I         The number of parameters
    * -Returns :
    *        true(success) false (failed)
    */
    int cKfAdjuster::Adjustment(VectorXd L, const MatrixXd H, const MatrixXd R, VectorXd &X, MatrixXd &Px, int nl, int nx) {

        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        vector<int>zip_idx;
        int i,j;
        qc_flag=false;
       
        /* init */
        if(Qvv_.data()) Qvv_.resize(0,0);
        if(v_.data()) v_.resize(0);
        unit_weight_STD_=0.0;

        /* remove disabled paremeter */
        for(i=0;i<nx;i++){
//            if(X[i]!=0.0&&Px(i,i)>0.0&&X[i]!=DIS_FLAG&&H.row(i).norm()) zip_idx.push_back(i);
            if(X[i]!=0.0&&fabs(Px(i,i))>0.0&&X[i]!=DIS_FLAG) zip_idx.push_back(i);

        }
        MatrixXd H_,Px_;
        VectorXd X_;
        H_=MatrixXd::Zero(zip_idx.size(),nl);Px_=MatrixXd::Zero(zip_idx.size(),zip_idx.size());X_=VectorXd::Zero(zip_idx.size());

        for(i=0;i<zip_idx.size();i++){
            X_[i]=X[zip_idx[i]];
            for(j=0;j<zip_idx.size();j++) Px_(i,j)=Px(zip_idx[i],zip_idx[j]);
            for(j=0;j<nl;j++) H_(i,j)=H(zip_idx[i],j);
        }

#if 0
        cout<<X_.transpose()<<endl;
        cout<<H_.transpose()<<endl;
        cout<<Px_.transpose()<<endl;
#endif
        //π€≤‚∏¸–¬
        MatrixXd Kk=(H_.transpose()*Px_*H_+R);
        if(MatInv(Kk.data(),Kk.cols())==-1) return 0;
        MatrixXd K=Px_*H_*Kk;//

        dx_=K*L; 
        for(int i=0;i<zip_idx.size();i++) X_[i]+=dx_[i];//
        MatrixXd I=MatrixXd::Identity(zip_idx.size(),zip_idx.size());
        Px_=(I-K*H_.transpose())*Px_;//

        v_=H_.transpose()*dx_-L;//
        MatrixXd Rr=R;
        if(MatInv(Rr.data(),R.cols())!=-1){
            unit_weight_STD_=v_.transpose()*Rr*v_;
            unit_weight_STD_=SQRT(unit_weight_STD_/((nl>nx)?nl-nx:nl));
        }

        for(i=0;i<zip_idx.size();i++){
            X[zip_idx[i]]=X_[i];
            for(int j=0;j<zip_idx.size();j++) Px(zip_idx[i],zip_idx[j])=Px_(i,j);
        }

#ifdef QC
        MatrixXd R_inv=R.inverse();
        RemoveEmptyColumns(H_);
//        cout<<"P:---------------"<<endl;
//        cout<<R_inv<<endl;
//        cout<<"(BtPB)^-1-------------"<<endl;
        MatrixXd M=H_*R_inv*H_.transpose();
        if(MatInv(M.data(),M.cols())!=-1){
            qc_flag= true;
            MatrixXd N=H_.transpose()*M*H_;
            Qvv_=(R-N);
//            cout<<"Qvv"<<endl;
//            cout<<Qvv_<<endl<<endl;
        }
#endif

        zip_idx.clear();

        return true;
    }
}
