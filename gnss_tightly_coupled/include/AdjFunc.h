// The class of LSQ and KF 
// Created by wlzhang  on 6/30/22.
// 

#ifndef FAST_LIO_IGCOUPLED_ADJFUNC_H
#define FAST_LIO_IGCOUPLED_ADJFUNC_H

#include "CmnFunc.h"

namespace IGCoupled{

    class cAdjuster {
    public:
        cAdjuster();
        virtual ~cAdjuster();

    public:
        virtual int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx);
        void RemoveEmptyColumns(MatrixXd& H);

    public:
        VectorXd dx_;
        VectorXd v_;
        MatrixXd Qvv_;
        double unit_weight_STD_=0.0;
        bool qc_flag=false;
    };

    // Lsq adjuster extend class cAdjuster 
    class cLsqAdjuster:public cAdjuster{
    public:
        cLsqAdjuster();
        ~cLsqAdjuster();

    public:
        int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx) override;
    };

    //Kalman fiter adjuster extend class cAdjuster
    class cKfAdjuster:public cAdjuster{
    public:
        cKfAdjuster();
        ~cKfAdjuster();

    public:
        int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx) override;
    };

}
#endif
