//
// Created by wlzhang on 7/1/22.
//

#ifndef FAST_LIO_IGCOUPLED_INSFUNC_H
#define FAST_LIO_IGCOUPLED_INSFUNC_H

#include <ceres/ceres.h>
#ifdef LOG
#undef LOG
#endif

#include "CmnFunc.h"
#include "GnssFunc.h"
#include <memory>

namespace IGCoupled {

    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec);
    Eigen::Matrix3d Quaternion2RotationMatrix(const Eigen::Quaterniond& q);
    Eigen::Quaterniond RotationMatrix2Quaternion(const Eigen::Matrix3d& m);
    Eigen::Quaterniond Euler2Quaternion(const Vector3d& rpy);
    Eigen::Matrix3d Euler2RotationMatrix(const Vector3d& rpy);
    Eigen::Vector3d RotationMatrix2Euler(const Matrix3d &m);
    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q);
    Eigen::Quaterniond RotationVector2Quaternion(const Vector3d& rv);
	Eigen::Matrix4d quaternionleft(const Quaterniond &q);
    Eigen::Matrix4d quaternionright(const Quaterniond &p);

    Eigen::Vector3d CalculateGravity(const Vector3d coord_blh,bool is_ecef);

    typedef struct{
        cTime t_tag;
        double dt;     // time difference related to increment distance
        double dr;     // increment of distance (m)
        double vr[3];  // wheel velocity in vehicle rear frame
		double odovel; //wheel velocity in w frame;
    }tOdoDataUnit; // struct of odometry measurement

    typedef struct {
        cTime t_tag;
        Vector3d gyro;
        Vector3d acce;

        unsigned int pps;
        unsigned int imu_cnt;

        short int odo_cnt;
        tOdoDataUnit odo;
    }tImuDataUnit; // struct of IMU data Unit 
  
    // class of IMU data and some function 
    class cImuData{
    public:
        cImuData();
        cImuData(cTime* ts,cTime* te);
        ~cImuData();

    public:
        void SetImu(tInsConf C);
        void SetImuType(IMU_TYPE type);
        void SetImuCoordType(IMU_COORD_TYPE type);
        void SetTimeSpan(cTime* ts, cTime* te);

    public:
        cTime ts_,te_;
        IMU_TYPE imu_type_;
        IMU_COORD_TYPE imu_coord_type_;
        IMU_DATA_FORMAT data_format_;
        GYRO_DATA_FORMAT gyro_format_;
        double hz_;
        vector<tImuDataUnit> data_;
    };


    typedef struct{
        cTime t_tag;
        Vector3d raw_gyro,raw_acce;
        Vector3d cor_gyro,cor_acce;

        Vector3d re,ve,ae;
        Vector3d rn,vn,an;
        Matrix3d Cbe,Cbn;
        Vector3d rpy;

        Vector3d ba,bg;
		Matrix3d ra, rg, sa, sg; // non-orthogonal between sensor axes and body frame for accl and gyro / scale factor error

        double dt;
        cTime pt;
    }tImuInfoUnit; // class of IMU information Unit
    
	//class of IMU Mechanization
    class cInsMech{
    public:
        cInsMech();
        ~cInsMech();

    public:
        bool InsMechanization(bool err_model,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int idx);
        Eigen::MatrixXd StateTransferMat(tIGCOUPLEDConf C,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int nx,double dt);

    private:
        void RotScullCorr(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,double *da,double *dv);
        Eigen::Quaterniond AttitudeUpdate(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,Vector3d da);
        Eigen::Vector3d VelocityUpdate(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,Vector3d dv);
        Eigen::Vector3d PositionUpdate(const tImuInfoUnit& pre_imu_info,const Vector3d& cur_vel,double dt);
        void TraceInsMechInfo(tImuInfoUnit &imu_info,bool prior,int idx);
		void ErrModel(tImuInfoUnit& pre_imu_info, tImuInfoUnit& cur_imu_info);

    };

    // class of IMU init alignment
    class cInsAlign{
    public:
        cInsAlign();
        cInsAlign(cImuData imu_data,tIGCOUPLEDConf C);
        ~cInsAlign();

    public:
        bool CoarseAlign(tImuInfoUnit& imu_info);
        bool GnssSolAlign();
        bool GnssObsAlign();

    private:
        cImuData *imu_data_;
        tIGCOUPLEDConf C_;

    public:
        int imu_idx=0;
    };

	// INS preintegration

	typedef struct tIntegrationStateData {
		cTime t_tag;
		
		double pose[7]; // 3 pose + 4 quaterniond
		

		//vel + bias : 3 + 6 = 9
		//vel + bias + sodo : 3 + 6 + 1 = 10
		//vel + bias + sodo + abv : 3 + 6 + 1 + 2 = 12
		//vel + bias + scale : 3 + 6 + 6 = 15
		//ve; + bias + sodo + scale : 3 + 6 + 1 + 6 = 16
		//vel + bias + sodo + scale + abv : 3 + 6 + 1 + 6 + 2 = 18
		double mix[18];

		double dt[NSYS + 1 + NSYS + NSYS + NUM_GLO_SAT + 2]; // clk + bd3_isb + phase_clk + ifb + glo_ifcb + glo_ifpb : 5 +1 + 5 + 5 + NUM_GLO_SAT + 2 
		double trp[6 + MAX_SAT_NUM]; //trp + ion : 6 + MAX_SAT_NUM
		double amb[MAX_SAT_NUM * 3]; //amb
		//f1 : MAX_SAT_NUM + f2 : MAX_SAT_NUM + f3 : MAX_SAT_NUM

	}tIntegrationStateData;

	class cPreIntegrationBase{
	
	public:
		cPreIntegrationBase(tIGCOUPLEDConf C, tImuDataUnit& imu, tFGOState state);
		virtual ~cPreIntegrationBase() = default;

		const tFGOState &currentState(){
			return current_state_;
		}

		const tFGOState &deltaState(){
			return delta_state_;
		}

		double deltaTime() const{
			return delta_time_;
		}

		cTime startTime() const {
			return start_time_;
		}

		cTime endTime() const {
			return end_time_;
		}

		const Vector3d &gravity(){
			return gravity_;
		}

		const vector<tImuDataUnit> &imuBuffer(){
			return imu_buffer_;
		}

		void addNewImu(const tImuDataUnit &imu);

		void reintegration(tFGOState &state);
	public:
        virtual Eigen::MatrixXd evaluate(const tFGOState &state0, const tFGOState &state1,
                                         double *residuals) = 0;
        virtual Eigen::MatrixXd residualJacobianPose0(const tFGOState &state0, const tFGOState &state1,
                                                      double *jacobian) = 0;
        virtual Eigen::MatrixXd residualJacobianPose1(const tFGOState &state0, const tFGOState &state1,
                                                      double *jacobian) = 0;
        virtual Eigen::MatrixXd residualJacobianMix0(const tFGOState &state0, const tFGOState &state1,
                                                     double *jacobian)  = 0;
        virtual Eigen::MatrixXd residualJacobianMix1(const tFGOState &state0, const tFGOState &state1,
                                                     double *jacobian)  = 0;

        virtual int numResiduals()                     = 0;
        virtual std::vector<int> numBlocksParameters() = 0;

        virtual void constructState(const double *const *parameters, tFGOState &state0,
                                    tFGOState &state1) = 0;

        virtual int imuErrorNumResiduals()                     = 0;
        virtual std::vector<int> imuErrorNumBlocksParameters() = 0;

        virtual void imuErrorEvaluate(const double *const *parameters, double *residuals) = 0;

        virtual void imuErrorJacobian(double *jacobian) = 0;

    protected:
        // need compensate bias
        virtual void updateJacobianAndCovariance(const tImuDataUnit &imu_pre, const tImuDataUnit &imu_cur) = 0;

        virtual void resetState(const tFGOState &state) = 0;
        virtual void integrationProcess(unsigned long index)   = 0;

        static void stateToData(const tFGOState &state, tIntegrationStateData &data);
        static void stateFromData(const tIntegrationStateData &data, tFGOState &state);

        void integration(tImuDataUnit &imu_pre, tImuDataUnit &imu_cur);

        tImuDataUnit compensationBias(const tImuDataUnit &imu) const;

        tImuDataUnit compensationScale(const tImuDataUnit &imu) const;

    public:
        static constexpr int NUM_POSE = 7;
    protected:
        static constexpr double IMU_GRY_BIAS_STD = 7200 / 3600.0 * M_PI / 180.0; // 7200 deg / hr
        static constexpr double IMU_ACC_BIAS_STD = 2.0e4 * 1.0e-5;               // 20000 mGal
        static constexpr double IMU_SCALE_STD    = 5.0e3 * 1.0e-6;               // 5000 PPM
        static constexpr double IMU_ACC_Z_SCALE  = 100;
        static constexpr double ODO_SCALE_STD    = 2.0e4 * 1.0e-6; // 0.02

        //const std::shared_ptr<IntegrationParameters> parameters_;
        tIGCOUPLEDConf C_;

        tFGOState current_state_; // 当前绝对状态
        tFGOState  delta_state_; // 当前预积分状态

        vector<tImuDataUnit> imu_buffer_;
        double delta_time_{0};
        cTime start_time_;
        cTime end_time_;

        Vector3d gravity_; // n系下的重力

        Eigen::MatrixXd jacobian_, covariance_;
        Eigen::MatrixXd noise_;
        Eigen::MatrixXd sqrt_information_;

        Quaterniond corrected_q_;
        Vector3d corrected_p_, corrected_v_;

	};

	class cPreIntegrationEarthOdo: public cPreIntegrationBase{

    public:
        cPreIntegrationEarthOdo(tIGCOUPLEDConf C, tImuDataUnit &imu, tFGOState state);

        MatrixXd evaluate(const tFGOState &state0, const tFGOState &state1, double *residuals) override;
        Eigen::MatrixXd residualJacobianPose0(const tFGOState &state0, const tFGOState &state1,
                                              double *jacobian) override;
        Eigen::MatrixXd residualJacobianPose1(const tFGOState &state0, const tFGOState &state1,
                                              double *jacobian) override;
        Eigen::MatrixXd residualJacobianMix0(const tFGOState &state0, const tFGOState &state1,
                                             double *jacobian) override;
        Eigen::MatrixXd residualJacobianMix1(const tFGOState &state0, const tFGOState &state1,
                                             double *jacobian) override;
        int numResiduals() override;
        vector<int> numBlocksParameters() override;

        static tIntegrationStateData stateToData( const tFGOState &state);
        static tFGOState stateFromData( const tIntegrationStateData &data);
        void constructState(const double *const *parameters, tFGOState &state0, tFGOState &state1) override;

        int imuErrorNumResiduals() override;
        vector<int> imuErrorNumBlocksParameters() override;
        void imuErrorEvaluate(const double *const *parameters, double *residuals) override;
        void imuErrorJacobian(double *jacobian) override;

    protected:
        void integrationProcess(unsigned long index) override;
        void resetState(const tFGOState &state) override;

        void updateJacobianAndCovariance(const tImuDataUnit &imu_pre, const tImuDataUnit &imu_cur) override;

    private:
        void resetState(const tFGOState &state, int num);
        void setNoiseMatrix();

    public:
        static constexpr int NUM_MIX = 9;
        static constexpr int NUM_MIX_ODO = 10;

    private:
        static constexpr int NUM_STATE = 15;
        static constexpr int NUM_STATE_ODO = 19;
        static constexpr int NUM_NOISE = 12;
        static constexpr int NUM_NOISE_ODO =16;

        static constexpr int NUM_ERROR_RESIDUAL = 7;

        Matrix3d cvb_;
        Vector3d lodo_;

        Vector3d corrected_s_;

        Quaterniond q0_;
        Vector3d iewe_;
        Matrix3d iewe_skew_;

        vector<std::pair<double, Vector3d>> pe_;
        Vector3d dpe_, dve_;
        Quaterniond qb0b1_;
    };

    class cPreIntegration{

    public:
        cPreIntegration() = default;

        static std::shared_ptr<cPreIntegrationBase> createPreintegration(const tIGCOUPLEDConf C, tImuDataUnit &imu0, const tFGOState &state) {
            std::shared_ptr<cPreIntegrationBase> preintegration;


            preintegration = std::make_shared<cPreIntegrationEarthOdo>(C, imu0, state);

            return preintegration;
        }

        static int numPoseParameter() {
            return cPreIntegrationBase::NUM_POSE;
        }

        static tIntegrationStateData stateToData(const tFGOState &state) {
            return cPreIntegrationEarthOdo::stateToData(state);
        }

        static tFGOState stateFromData(const tIntegrationStateData &data) {
                return cPreIntegrationEarthOdo::stateFromData(data);
            }

        static int numMixParameter(tIGCOUPLEDConf C) {
            int num = 0;
            if(C.insC.use_odo){
                num = cPreIntegrationEarthOdo::NUM_MIX_ODO;
            }else num = cPreIntegrationEarthOdo::NUM_MIX;
            return num;
        }

    };
    // Ins Preintegration Factor

    class cPreIntegrationFactor : public ceres::CostFunction{

    public:

        cPreIntegrationFactor() = delete;

        explicit cPreIntegrationFactor(std::shared_ptr<cPreIntegrationBase> preintegration)
        : preintegration_(std::move(preintegration)){
            // parameter
            *mutable_parameter_block_sizes() = preintegration_->numBlocksParameters();

            // residual
            set_num_residuals(preintegration_->numResiduals());
        }

        bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override{
            // construct state
            tFGOState state0, state1;
            preintegration_->constructState(parameters, state0, state1);

            // residual
            preintegration_->evaluate(state0, state1, residuals);

            if (jacobians) {
                if (jacobians[0]) {
                    preintegration_->residualJacobianPose0(state0, state1, jacobians[0]);
                }
                if (jacobians[1]) {
                    preintegration_->residualJacobianMix0(state0, state1, jacobians[1]);
                }
                if (jacobians[2]) {
                    preintegration_->residualJacobianPose1(state0, state1, jacobians[2]);
                }
                if (jacobians[3]) {
                    preintegration_->residualJacobianMix1(state0, state1, jacobians[3]);
                }
            }

            return true;
        }

    private:
        std::shared_ptr<cPreIntegrationBase> preintegration_;
    };

    class ImuErrorFactor : public ceres::CostFunction {

    public:
        explicit ImuErrorFactor(std::shared_ptr<cPreIntegrationBase> preintegration)
                : preintegration_(std::move(preintegration)) {

            *mutable_parameter_block_sizes() = preintegration_->imuErrorNumBlocksParameters();
            set_num_residuals(preintegration_->imuErrorNumResiduals());
        }

        bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {

            // parameters: vel[3], bg[3], ba[3]

            preintegration_->imuErrorEvaluate(parameters, residuals);

            if (jacobians) {
                if (jacobians[0]) {
                    preintegration_->imuErrorJacobian(jacobians[0]);
                }
            }

            return true;
        }

    private:
        std::shared_ptr<cPreIntegrationBase> preintegration_;
    };


    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt);
}
#endif //FAST_LIO_IGCOUPLED_INSFUNC_H
