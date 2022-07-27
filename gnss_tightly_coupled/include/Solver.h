// class of gnss solvers
// Created by wlzhang on 7/1/22.
//

#ifndef FAST_LIO_IGCOUPLED_SOLVER_H
#define FAST_LIO_IGCOUPLED_SOLVER_H

#include <utility>

#include "GnssFunc.h"
#include "GnssErrorModel.h"
#include "InsFunc.h"
#include "AdjFunc.h"
#include "GnssAR.h"
#include "OutSol.h"


#define POSE_LOCAL_SIZE 6
#define POSE_GLOBAL_SIZE 7

using namespace Eigen;
namespace IGCoupled{

//  GNSS Solver basic class
    class cSolver {
    public:
        cSolver();
        virtual ~cSolver();

    private:
        int GetSingalInd(tSatObsUnit& sat_obs,int *f);
        int AdjustObsInd(tSatObsUnit& sat_obs,int* i,int* j,int* k);

    public:
        void InitEpochSatInfo(vector<tSatInfoUnit>& sat_infos);
        void UpdateSatInfo(vector<tSatInfoUnit>& sat_infos);
        void UpdateGnssObs(tIGCOUPLEDConf C,tEpochSatUnit& epoch_sat,RECEIVER_INDEX rec);
        void CorrGnssObs(tIGCOUPLEDConf C,Vector3d& rr);
        virtual int GnssObsRes(int post,tIGCOUPLEDConf C,double* x);
        void CombFbSol(tIGCOUPLEDConf C);
        int Smoother(const VectorXd xf, const MatrixXd Qf, const VectorXd xb, const MatrixXd Qb, VectorXd& xs, MatrixXd& Qs, int n);

        virtual void InitSolver(tIGCOUPLEDConf C);
        virtual bool InitReader(tIGCOUPLEDConf C);
        virtual bool SolverProcess(tIGCOUPLEDConf C,int idx);
        virtual bool SolverStart(int i,int idx);
        virtual bool SolverEpoch();
        virtual bool Estimator(tIGCOUPLEDConf C);
        virtual bool SolutionUpdate();

        void ReinitSolver(tIGCOUPLEDConf C);
        void InitInsPx(tIGCOUPLEDConf C,int nx,MatrixXd& Px);
        void InitX(double xi,double var,int idx,double *x,double *xp);
        Eigen::MatrixXd InitQ(tIGCOUPLEDConf,double dt,int nx);
        Eigen::MatrixXd InitPrecQ(tIGCOUPLEDConf C,double dt,int nx,Matrix3d Cbe);

        void RemoveLever(const tImuInfoUnit& imu_info,Vector3d& lever,Vector3d& gnss_re,Vector3d& gnss_ve);
        void CloseLoopState(VectorXd& x,tImuInfoUnit* imu_info_corr);

    public:
        cGnssObsOperator gnss_obs_operator_;
        cGnssErrCorr gnss_err_corr_;
        cParSetting para_;
        cLsqAdjuster lsq_;
        cKfAdjuster kf_;
        cOutSol *out_;

    public:
        int epoch_idx_=0;
        int epoch_ok_=0;
        int epoch_fail_=0;
        tNav nav_;
        cGnssObs rover_obs_;
        cGnssObs base_obs_;
        vector<tSolInfoUnit> ref_sols_;

        tEpochSatUnit epoch_sat_obs_;
        tEpochSatUnit base_epoch_sat_obs_;

        vector<tSatInfoUnit> epoch_sat_info_collect_;
        vector<tSatInfoUnit> base_sat_info_collect_;

        vector<int>vflag_;
        tSolInfoUnit igcoupled_sol_;
        tSatInfoUnit previous_sat_info_[MAX_SAT_NUM];
        cLambda lambda_;

        vector<tSolInfoUnit> solf_,solb_;

    public:
        int num_full_x_;
        VectorXd full_x_;
        MatrixXd full_Px_;
        int num_valid_sat_;
        int num_L_;
        VectorXd omc_L_;
        MatrixXd H_;
        MatrixXd R_;
        MatrixXd F_; //for coupled
        int sys_mask_[NSYS+1]={0};
        int num_real_x_fix_;
        VectorXd real_x_fix_;
        MatrixXd real_Px_fix_;
        vector<double>code_sigma,phase_sigma;

        bool tc_mode_=false;
        tImuInfoUnit cur_imu_info_;
    };

// SPP solver
    class cSppSolver:public cSolver {
    public:
        cSppSolver();
        cSppSolver(tIGCOUPLEDConf conf);
        ~cSppSolver();

    private:
        void CorrDoppler(tSatInfoUnit& sat_info, Vector3d& rover_xyz,int f);
        Vector3d SatVelSagnacCorr(const Vector3d& sat_pos,const double tau);
        int DopplerRes(tIGCOUPLEDConf C,MatrixXd& H_mat, MatrixXd& R_mat,VectorXd& L,VectorXd& x,Vector3d rover_xyz);
        void EstDopVel(Vector3d rover_xyz);

        double Dops();
        bool ValidateSol(tIGCOUPLEDConf C);
        bool PostResidualQc(vector<double>omcs,vector<double>R);
        bool RaimFde();

    public:
        int GnssObsRes(int post,tIGCOUPLEDConf C,double* x) override;
        void InitSolver(tIGCOUPLEDConf C) override;
        bool SolverProcess(tIGCOUPLEDConf C,int idx) override;
        bool SolverEpoch() override;
        bool Estimator(tIGCOUPLEDConf C) override;
        bool SolutionUpdate() override;

    private:
        int iter_=10;
    public:
        tIGCOUPLEDConf spp_conf_;
    };

    class cTdcpSolver:public cSolver {
    public:
        cTdcpSolver();
        cTdcpSolver(tIGCOUPLEDConf conf);
        ~cTdcpSolver();

    public:
        tIGCOUPLEDConf  tdcp_conf_;
    };

// PPP Solver
    class cPppSolver:public cSolver {
    public:
        cPppSolver();
        cPppSolver(tIGCOUPLEDConf C);
        ~cPppSolver();

    public:
        void InitSolver(tIGCOUPLEDConf C) override;
        bool SolverProcess(tIGCOUPLEDConf C,int idx) override;
        bool SolverStart(int i,int idx) override ;
        bool SolverEpoch() override;
        bool Estimator(tIGCOUPLEDConf C) override;
        bool SolutionUpdate() override;

    private:
        void InitSppSolver();
        void Spp2Ppp();
        void PppCycSlip(tIGCOUPLEDConf C);
        void ClockJumpRepair();
        void PosUpdate(tIGCOUPLEDConf C);
        void ClkUpdate(tIGCOUPLEDConf C,double tt);
        void PhaseClkUpdate(tIGCOUPLEDConf C,double tt);
        void IfbUpdate(tIGCOUPLEDConf C,double tt);
        void GloIfcbUpdate(tIGCOUPLEDConf C,double tt);
        void TrpUpdate(tIGCOUPLEDConf C,double tt);
        void IonUpdate(tIGCOUPLEDConf C,double tt);
        void AmbUpdate(tIGCOUPLEDConf C,double tt);
        void StateTimeUpdate(tIGCOUPLEDConf C);
        int GnssObsRes(int post,tIGCOUPLEDConf C,double *x,Vector3d re);
        void ReInitPppSolver(tIGCOUPLEDConf C);

        void DisableX(int iter,VectorXd& x,vector<int>&par_idx,vector<double>& back_values);

        // quality control
        bool PppResidualQc(int iter,vector<double>omcs,vector<double>var);

        // PPP-AR
        bool ResolvePppAmbRotRef(int nf,VectorXd& x,MatrixXd& P);
        // PPP-AR select high-ele satellite
        bool ResolvePppAmbFixRef(int nf,VectorXd& x,MatrixXd& P);

        bool FixPppSol(int *sat1,int *sat2,double *NC,int n,VectorXd& x,MatrixXd& P);
        bool AmbLinearDependCheck(int sat1,int sat2,int *flag,int *max_flg);
        int SelectAmb(int *sat1,int *sat2,double *N,double *var,int n);

        // fix wide-lane ambiguity
        void AverageLcAmb();
        bool FixWlAmb(int sat1,int sat2,int *fix_wl,double *res_wl);

        // fix narrow-lane ambiguity with CNES IRC
        bool FixNlAmbILS_IRC(int *sat1,int *sat2,int *fix_wls,double *res_wls,int n,VectorXd& x,MatrixXd& P);

        // fix narrow-lane ambiguity with SGG FCB
        bool MatchNlFcb(int sat1,int sat2,double *nl_fcb1,double *nl_fcb2);
        bool FixNlAmbILS_FCB(int *sat1,int *sat2,int *fix_wls,double *res_wls,int n,VectorXd& x,MatrixXd& P);

        // ppp-ar use ppk-ar function
        void ResetAmb(double *bias,double *xa,int nb);
        bool FixNlAmbILS1(int *sat1,int *sat2,int *fix_wls,int n,double *xa);
        int SelectFixSat(double *D,int gps,int glo);
        int ResolveAmbLambda(double *xa,int gps,int glo);
        bool ResolverPppAmb1(int nf,double *xa);
    private:
        tIGCOUPLEDConf ppp_conf_;
        int max_iter_=8;
        vector<tSdAmb> sdambs_;
        cTime filter_start_;
        bool reset_flag_=false;

    public:
        cSppSolver *spp_solver_;
        int exc_sat_index_=0;
        double pre_epoch_ar_ratio1_=0.0;
        double pre_epoch_ar_ratio2_=0.0;

    };

// PPk solver
    class cPpkSolver:public cSolver {
    public:
        cPpkSolver();
        cPpkSolver(tIGCOUPLEDConf C);
        ~cPpkSolver();

    public:
        void InitSolver(tIGCOUPLEDConf C) override;
        bool SolverProcess(tIGCOUPLEDConf C,int idx) override;
        bool SolverStart(int i,int idx) override;
        bool SolverEpoch() override;
        bool Estimator(tIGCOUPLEDConf C) override;
        bool SolutionUpdate() override;

    private:
        void InitSppSolver();
        void Spp2Ppk();
        bool GnssZeroRes(tIGCOUPLEDConf C,RECEIVER_INDEX rec,vector<int>sat_idx,double* x,Vector3d rr);
        int GnssDdRes(int post,tIGCOUPLEDConf C,vector<int>ir,vector<int>ib,vector<int>cmn_sat_no,double* x,int refsat[NSYS][2*MAX_GNSS_USED_FRQ_NUM]);
        bool ValidObs(int i,int nf,int f);
        bool MatchBaseObs(cTime t);
        int SelectCmnSat(tIGCOUPLEDConf C,vector<int>& ir,vector<int>& iu,vector<int>& cmn_sat_no);
        void PpkCycleSlip(tIGCOUPLEDConf C,vector<int>& iu,vector<int>& ib,vector<int>& cmn_sat_no);
        void StateTimeUpdate(tIGCOUPLEDConf C,vector<int>& iu,vector<int>& ib,vector<int>& cmn_sat_no);
        void PosUpdate(tIGCOUPLEDConf C);
        void GloIfpbUpdate(tIGCOUPLEDConf C,double tt);
        void TrpUpdate(tIGCOUPLEDConf C,double tt);
        void IonUpdate(tIGCOUPLEDConf C,double tt);
        void AmbUpdate(tIGCOUPLEDConf C, double tt,vector<int>& ir,vector<int>& ib,vector<int>& cmn_sat_no);

        bool CheckDualFrq(tSatInfoUnit& sat_info);
        bool PpkResidualQc(int iter,vector<int>ir,vector<int> cmn_sat,vector<double>& omcs, vector<double>R);

        void ReSetAmb(double *bias,double *xa,int nb);
        int DdMat(double *D,int gps,int bds,int glo,int gal,double el_mask);
        int ResolveAmbLambda(double *xa,int gps, int bds,int glo,int gal, int qzs,double el_mask);
        bool ResolvePpkAmb(vector<int>cmn_sat,int nf,double *xa);
        void HoldAmb(vector<int>cmn_sat,double *xa);

    private:
        cSppSolver *spp_solver_;
        tIGCOUPLEDConf ppk_conf_;
        Vector3d base_xyz_;

    private:
        double tt_=0.0;
        int base_idx_=0;
        vector<double>base_res,rover_res;
        vector<tDdAmb> ddambs_;
        int exc_sat_index=0;
        int num_continuous_fix_=0;
        bool amb_hold_flag=false;
        double pre_epoch_ar_ratio1=0.0;
        double pre_epoch_ar_ratio2=0.0;
        int qc_iter_=0;
    };

// Fusion Solver
    class cFusionSolver:public cSolver {
    public:
        cFusionSolver();
        cFusionSolver(tIGCOUPLEDConf C);
        ~cFusionSolver();

    private:
        bool InputImuData(int ws);
        bool MatchGnssObs(bool isBack);
        bool MatchGnssSol(bool isBack);
        void InsSol2IgcoupledSol(tImuInfoUnit &imu_sol,tSolInfoUnit &igcoupled_sol);

        double Vel2Yaw(Vector3d vn);
        bool GnssSol2Ins(Vector3d re,Vector3d ve);
        Vector3d Pos2Vel(tSolInfoUnit& sol1,tSolInfoUnit& sol2);
        bool InsAlign();

    public:
        void InitSolver(tIGCOUPLEDConf C) override;
        bool SolverProcess(tIGCOUPLEDConf C,int idx) override;
        bool SolverEpoch() override;
        bool SolutionUpdate() override;

    private:
        void DisableX(VectorXd& x,int idx);
        void StateTimeUpdate();
        void PropVariance(MatrixXd& F,MatrixXd& Q,int nx,MatrixXd& Px);
        int BuildLcHVR(int post,tIGCOUPLEDConf C,tImuInfoUnit& imu_info,double *meas_pos,double *meas_vel,Vector3d& q_pos,Vector3d& q_vel);
        bool LcFilter(tIGCOUPLEDConf C);
        bool ValidSol(VectorXd& x, double thres);
        void ControlPx(int nx,MatrixXd& Px);
		bool process(tIGCOUPLEDConf C, int idx);

    public:
        bool LooseCouple(tIGCOUPLEDConf C);
        bool TightCouple(tIGCOUPLEDConf C);

    private:
        cInsMech ins_mech_;
        cSolver *gnss_alignor_;
        cSolver *gnss_solver_;
        tIGCOUPLEDConf fs_conf_;
        tIGCOUPLEDConf gnss_conf_;

    private:
        cImuData imu_data_;
        vector<tSolInfoUnit> gnss_sols_;

    private:
        int imu_index_=0;
        int rover_idx_=0;
        int base_idx_=0;
        int gnss_sol_idx=0;
        int ins_mech_idx=0;
        tImuInfoUnit cur_imu_info_={0};
        tImuInfoUnit pre_imu_info_={0};
        vector<tImuDataUnit> imu_data_zd_;
    };

	//FGO Solver
	class cFGOSolver: public cSolver{
    public:
        cFGOSolver();
        explicit cFGOSolver(tIGCOUPLEDConf C);
        ~cFGOSolver() override;

    public:
        void InitSolver(tIGCOUPLEDConf C) override;
        bool SolverProcess(tIGCOUPLEDConf C, int idx) override;
        bool process(tIGCOUPLEDConf C, int type);

    private:
		void State2Data(tFGOState& state, tIntegrationStateData &data);
        void StateFromData(tIntegrationStateData &data, tFGOState &state);
        void StateTimeUpdate(tFGOState &state, VectorXd& x);
        void x2state(cSolver* solver, VectorXd x, tFGOState& state);
        bool InsAlign();
        bool MatchGnssSol();
        bool MatchGnssObs();
        Eigen::Vector3d Pos2Vel(tSolInfoUnit &sol1, tSolInfoUnit &sol2);
        double Vel2Yaw(Vector3d vn);
        bool GnssSol2Ins(Vector3d re, Vector3d ve);
        void InsSol2IgcoupledSol(tImuInfoUnit &imu_sol,tSolInfoUnit &igcp_sol);
        void PreIntSol2IgcoupledSol(const cTime& time, const tFGOState& state, tSolInfoUnit &igcp_sol, Vector3d lever);
        static int isNeedInterpolation(tImuDataUnit &imu0, tImuDataUnit &imu1, double mid);
        void ImuInterp(const tImuDataUnit& imu01, tImuDataUnit &imu00,tImuDataUnit &imu11, double mid) const;
        void ComFBSol(tIGCOUPLEDConf C);
        double findOutage(double sow, vector<double> outageList);
    private:
        cSolver *gnss_alignor_;
        cSolver *gnss_solver_;
        tIGCOUPLEDConf fgo_conf_;
        tIGCOUPLEDConf gnss_conf_;

    private:
        cImuData imu_data_;
        vector<tSolInfoUnit> gnss_sols_; //

        vector<std::shared_ptr<cPreIntegrationBase>> preintegrationlist_;
        vector<std::shared_ptr<cSolver>> gnssSolverList_;
        vector<tFGOState> statelist_;
        vector<tIntegrationStateData> statedatalist_;
        vector<vector<tSatInfoUnit>> EpochSatInfolist_;




        vector<tSolInfoUnit> GnssSolList_; // 存放窗口的 sol结果
        vector<tEpochSatUnit> GnssObsList_; // 存放窗口的 观测值

        vector<double> timelist_; // 存放窗口内的gps秒
    private:
        bool image_mode_;
        int imu_idx_;
        int gnss_sol_idx_;
        int rover_idx_;
        tImuInfoUnit cur_imu_info_;
        tImuInfoUnit pre_imu_info_;


    };

	//
	//残差块
    class cResidualBlockInfo{

    public:
        cResidualBlockInfo(std::shared_ptr<ceres::CostFunction> cost_function,
                          std::shared_ptr<ceres::LossFunction> loss_function, std::vector<double *> parameter_blocks,
                          std::vector<int> marg_para_index)
                : cost_function_(std::move(cost_function))
                , loss_function_(std::move(loss_function))
                , parameter_blocks_(std::move(parameter_blocks))
                , marg_para_index_(std::move(marg_para_index)) {
        }

        void Evaluate() {
            residuals_.resize(cost_function_->num_residuals());

            std::vector<int> block_sizes = cost_function_->parameter_block_sizes();
            auto raw_jacobians           = new double *[block_sizes.size()];
            jacobians_.resize(block_sizes.size());

            for (int i = 0; i < static_cast<int>(block_sizes.size()); i++) {
                jacobians_[i].resize(cost_function_->num_residuals(), block_sizes[i]);
                raw_jacobians[i] = jacobians_[i].data();
            }
            cost_function_->Evaluate(parameter_blocks_.data(), residuals_.data(), raw_jacobians);

            delete[] raw_jacobians;

            if (loss_function_) {
                // 鲁棒核函数调整, 参考ceres/internal/ceres/corrector.cc
                double residual_scaling, alpha_sq_norm;

                double sq_norm, rho[3];

                sq_norm = residuals_.squaredNorm();
                loss_function_->Evaluate(sq_norm, rho);

                double sqrt_rho1 = sqrt(rho[1]);

                if ((sq_norm == 0.0) || (rho[2] <= 0.0)) {
                    residual_scaling = sqrt_rho1;
                    alpha_sq_norm    = 0.0;
                } else {
                    // 0.5 *  alpha^2 - alpha - rho'' / rho' *  z'z = 0
                    const double D     = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
                    const double alpha = 1.0 - sqrt(D);
                    residual_scaling   = sqrt_rho1 / (1 - alpha);
                    alpha_sq_norm      = alpha / sq_norm;
                }

                for (size_t i = 0; i < parameter_blocks_.size(); i++) {
                    // J = sqrt_rho1 * (J - alpha_sq_norm * r* (r.transpose() * J))
                    jacobians_[i] =
                            sqrt_rho1 * (jacobians_[i] - alpha_sq_norm * residuals_ * (residuals_.transpose() * jacobians_[i]));
                }
                residuals_ *= residual_scaling;
            }
        }

        const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &jacobians() {
            return jacobians_;
        }

        const std::vector<int> &parameterBlockSizes() {
            return cost_function_->parameter_block_sizes();
        }

        const std::vector<double *> &parameterBlocks() {
            return parameter_blocks_;
        }

        const Eigen::VectorXd &residuals() {
            return residuals_;
        }

        const std::vector<int> &marginalizationParametersIndex() {
            return marg_para_index_;
        }

    private:
        std::shared_ptr<ceres::CostFunction> cost_function_;
        std::shared_ptr<ceres::LossFunction> loss_function_;

        std::vector<double *> parameter_blocks_;

        std::vector<int> marg_para_index_;

        std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobians_;
        Eigen::VectorXd residuals_;
    };

    //边迹化
    class cMarginalizationInfo {

    public:
        cMarginalizationInfo() = default;

        ~cMarginalizationInfo() {
            for (auto &block : parameter_block_data_)
                delete[] block.second;
        }

        bool isValid() const {
            return isvalid_;
        }

        static int localSize(int size) {
            return size == POSE_GLOBAL_SIZE ? POSE_LOCAL_SIZE : size;
        }

        static int globalSize(int size) {
            return size == POSE_LOCAL_SIZE ? POSE_GLOBAL_SIZE : size;
        }

        void addResidualBlockInfo(const std::shared_ptr<cResidualBlockInfo> &blockinfo) {
            factors_.push_back(blockinfo);

            const auto &parameter_blocks = blockinfo->parameterBlocks();
            const auto &block_sizes      = blockinfo->parameterBlockSizes();

            for (size_t k = 0; k < parameter_blocks.size(); k++) {
                parameter_block_size_[reinterpret_cast<long>(parameter_blocks[k])] = block_sizes[k];
            }

            // 被边缘化的参数, 先加入表中以进行后续的排序
            for (int index : blockinfo->marginalizationParametersIndex()) {
                parameter_block_index_[reinterpret_cast<long>(parameter_blocks[index])] = 0;
            }
        }

        bool marginalization() {

            // 对边缘化的参数和保留的参数按照local size分配索引, 边缘化参数位于前端
            if (!updateParameterBlocksIndex()) {
                isvalid_ = false;
                // 释放内存
                releaseMemory();

                return false;
            }


            // 计算每个残差块参数, 进行参数内存拷贝
            preMarginalization();

            // 构造增量线性方程
            constructEquation();

            // Schur消元
            schurElimination();

            // 求解线性化雅克比和残差
            linearization();

            // 释放内存
            releaseMemory();

            return true;
        }

        std::vector<double *> getParamterBlocks(std::unordered_map<long, double *> &address) {
            std::vector<double *> remained_block_addr;

            remained_block_data_.clear();
            remained_block_index_.clear();
            remained_block_size_.clear();

            for (const auto &block : parameter_block_index_) {
                // 保留的参数
                if (block.second >= marginalized_size_) {
                    remained_block_data_.push_back(parameter_block_data_[block.first]);
                    remained_block_size_.push_back(parameter_block_size_[block.first]);
                    remained_block_index_.push_back(parameter_block_index_[block.first]);
                    remained_block_addr.push_back(address[block.first]);
                }
            }

            return remained_block_addr;
        }

        const Eigen::MatrixXd &linearizedJacobians() {
            return linearized_jacobians_;
        }

        const Eigen::VectorXd &linearizedResiduals() {
            return linearized_residuals_;
        }

        int marginalizedSize() const {
            return marginalized_size_;
        }

        int remainedSize() const {
            return remained_size_;
        }

        const std::vector<int> &remainedBlockSize() {
            return remained_block_size_;
        }

        const std::vector<int> &remainedBlockIndex() {
            return remained_block_index_;
        }

        const std::vector<double *> &remainedBlockData() {
            return remained_block_data_;
        }

    private:
        // 线性化
        void linearization() {
            // SVD分解求解雅克比, Hp = J^T * J = V * S^{1/2} * S^{1/2} * V^T
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(Hp_);
            Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > EPS).select(saes2.eigenvalues().array(), 0));
            Eigen::VectorXd S_inv =
                    Eigen::VectorXd((saes2.eigenvalues().array() > EPS).select(saes2.eigenvalues().array().inverse(), 0));

            Eigen::VectorXd S_sqrt     = S.cwiseSqrt();
            Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();

            // J0 = S^{1/2} * V^T
            linearized_jacobians_ = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
            // e0 = -{J0^T}^{-1} * bp = - S^{-1/2} * V^T * bp
            linearized_residuals_ = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * -bp_;
        }

        // Schur消元, 求解 Hp * dx_r = bp
        void schurElimination() {
            // H0 * dx = b0
            Eigen::MatrixXd Hmm = 0.5 * (H0_.block(0, 0, marginalized_size_, marginalized_size_) +
                                         H0_.block(0, 0, marginalized_size_, marginalized_size_).transpose());
            Eigen::MatrixXd Hmr = H0_.block(0, marginalized_size_, marginalized_size_, remained_size_);
            Eigen::MatrixXd Hrm = H0_.block(marginalized_size_, 0, remained_size_, marginalized_size_);
            Eigen::MatrixXd Hrr = H0_.block(marginalized_size_, marginalized_size_, remained_size_, remained_size_);
            Eigen::VectorXd bmm = b0_.segment(0, marginalized_size_);
            Eigen::VectorXd brr = b0_.segment(marginalized_size_, remained_size_);

            // SVD分解Amm求逆
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Hmm);
            Eigen::MatrixXd Hmm_inv =
                    saes.eigenvectors() *
                    Eigen::VectorXd((saes.eigenvalues().array() > EPS).select(saes.eigenvalues().array().inverse(), 0))
                            .asDiagonal() *
                    saes.eigenvectors().transpose();

            // Hp = Hrr - Hrm * Hmm^-1 * Hmr
            Hp_ = Hrr - Hrm * Hmm_inv * Hmr;
            // bp = br - Hrm * Hmm^-1 * bm
            bp_ = brr - Hrm * Hmm_inv * bmm;
        }

        // 构造增量方程 H * dx = b, 计算 H 和 b
        void constructEquation() {
            H0_ = Eigen::MatrixXd::Zero(local_size_, local_size_);
            b0_ = Eigen::VectorXd::Zero(local_size_);

            for (const auto &factor : factors_) {
                for (size_t i = 0; i < factor->parameterBlocks().size(); i++) {
                    int row0 = parameter_block_index_[reinterpret_cast<long>(factor->parameterBlocks()[i])];
                    int rows = parameter_block_size_[reinterpret_cast<long>(factor->parameterBlocks()[i])];
                    rows     = (rows == POSE_GLOBAL_SIZE) ? POSE_LOCAL_SIZE : rows;

                    Eigen::MatrixXd jacobian_i = factor->jacobians()[i].leftCols(rows);
                    for (size_t j = i; j < factor->parameterBlocks().size(); ++j) {
                        int col0 = parameter_block_index_[reinterpret_cast<long>(factor->parameterBlocks()[j])];
                        int cols = parameter_block_size_[reinterpret_cast<long>(factor->parameterBlocks()[j])];
                        cols     = (cols == POSE_GLOBAL_SIZE) ? POSE_LOCAL_SIZE : cols;

                        Eigen::MatrixXd jacobian_j = factor->jacobians()[j].leftCols(cols);

                        // H = J^T * J
                        if (i == j) {
                            // Hmm, Hrr
                            H0_.block(row0, col0, rows, cols) += jacobian_i.transpose() * jacobian_j;
                        } else {
                            // Hmr, Hrm = Hmr^T
                            H0_.block(row0, col0, rows, cols) += jacobian_i.transpose() * jacobian_j;
                            H0_.block(col0, row0, cols, rows) = H0_.block(row0, col0, rows, cols).transpose();
                        }
                    }
                    // b = - J^T * e
                    b0_.segment(row0, rows) -= jacobian_i.transpose() * factor->residuals();
                }
            }
        }

        bool updateParameterBlocksIndex() {
            int index = 0;
            // 只有被边缘化的参数预先加入了表
            for (auto &block : parameter_block_index_) {
                block.second = index;
                index += localSize(parameter_block_size_[block.first]);
            }
            marginalized_size_ = index;

            // 加入保留的参数, 分配索引
            for (const auto &block : parameter_block_size_) {
                if (parameter_block_index_.find(block.first) == parameter_block_index_.end()) {
                    parameter_block_index_[block.first] = index;
                    index += localSize(block.second);
                }
            }
            remained_size_ = index - marginalized_size_;

            local_size_ = index;

            return marginalized_size_ > 0;
        }

        // 边缘化预处理, 评估每个残差块, 拷贝参数
        void preMarginalization() {
            for (const auto &factor : factors_) {
                factor->Evaluate();

                std::vector<int> block_sizes = factor->parameterBlockSizes();
                for (size_t k = 0; k < block_sizes.size(); k++) {
                    long addr = reinterpret_cast<long>(factor->parameterBlocks()[k]);
                    int size  = block_sizes[k];

                    // 拷贝参数块数据
                    if (parameter_block_data_.find(addr) == parameter_block_data_.end()) {
                        auto *data = new double[size];
                        memcpy(data, factor->parameterBlocks()[k], sizeof(double) * size);
                        parameter_block_data_[addr] = data;
                    }
                }
            }
        }

        void releaseMemory() {
            // 释放因子所占有的内存, 尤其是边缘化因子及其占有的边缘化信息数据结构
                 factors_.clear();
               // factors_.resize(0);
               // vector<std::shared_ptr<cResidualBlockInfo>>().swap(factors_);
        }

    private:
        // 增量线性方程参数
        Eigen::MatrixXd H0_, Hp_;
        Eigen::VectorXd b0_, bp_;

        // 以内存地址为key的无序表

        // 存放参数块的global size
        std::unordered_map<long, int> parameter_block_size_;
        // 存放参数块索引, 待边缘化参数索引在前, 保留参数索引在后, 用于构造边缘化 H * dx = b
        std::unordered_map<long, int> parameter_block_index_;
        // 存放参数快数据指针
        std::unordered_map<long, double *> parameter_block_data_;

        // 保留的参数
        std::vector<int> remained_block_size_;  // global size
        std::vector<int> remained_block_index_; // local size
        std::vector<double *> remained_block_data_;

        // local size in total
        int marginalized_size_{0};
        int remained_size_{0};
        int local_size_{0};

        // 边缘化参数相关的残差块
        std::vector<std::shared_ptr<cResidualBlockInfo>> factors_;

        const double EPS = 1e-8;

        // 边缘化求解的残差和雅克比
        Eigen::MatrixXd linearized_jacobians_;
        Eigen::VectorXd linearized_residuals_;

        // 若无待边缘化参数, 则无效
        bool isvalid_{true};
    };

    //边迹化残差
    class cMarginalizationFactor : public ceres::CostFunction {

    public:
        cMarginalizationFactor() = delete;
        explicit cMarginalizationFactor(std::shared_ptr<cMarginalizationInfo> marg_info)
                : marg_info_(std::move(marg_info)) {

            // 给定每个参数块数据大小
            for (auto size : marg_info_->remainedBlockSize()) {
                mutable_parameter_block_sizes()->push_back(size);
            }

            // 残差大小
            set_num_residuals(marg_info_->remainedSize());
        }

        bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {
            int marginalizaed_size = marg_info_->marginalizedSize();
            int remained_size      = marg_info_->remainedSize();

            const vector<int> &remained_block_index     = marg_info_->remainedBlockIndex();
            const vector<int> &remained_block_size      = marg_info_->remainedBlockSize();
            const vector<double *> &remained_block_data = marg_info_->remainedBlockData();

            Eigen::VectorXd dx(remained_size);
            for (size_t i = 0; i < remained_block_size.size(); i++) {
                int size  = remained_block_size[i];
                int index = remained_block_index[i] - marginalizaed_size;

                Eigen::VectorXd x  = Eigen::Map<const Eigen::VectorXd>(parameters[i], size);
                Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(remained_block_data[i], size);

                // dx = x - x0
                if (size == POSE_GLOBAL_SIZE) {
                    Eigen::Quaterniond dq(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() *
                                          Eigen::Quaterniond(x(6), x(3), x(4), x(5)));

                    dx.segment(index, 3)     = x.head<3>() - x0.head<3>();
                    dx.segment(index + 3, 3) = 2.0 * dq.vec();
                    if (dq.w() < 0) {
                        dx.segment<3>(index + 3) = -2.0 * dq.vec();
                    }
                } else {
                    dx.segment(index, size) = x - x0;
                }
            }

            // e = e0 + J0 * dx
            Eigen::Map<Eigen::VectorXd>(residuals, remained_size) =
                    marg_info_->linearizedResiduals() + marg_info_->linearizedJacobians() * dx;

            if (jacobians) {

                for (size_t i = 0; i < remained_block_size.size(); i++) {
                    if (jacobians[i]) {
                        int size       = remained_block_size[i];
                        int index      = remained_block_index[i] - marginalizaed_size;
                        int local_size = marg_info_->localSize(size);

                        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(
                                jacobians[i], remained_size, size);

                        // J = J0
                        jacobian.setZero();
                        jacobian.leftCols(local_size) = marg_info_->linearizedJacobians().middleCols(index, local_size);
                    }
                }
            }

            return true;
        }

    private:
        std::shared_ptr<cMarginalizationInfo> marg_info_;
    };

    //位姿参数化
    class cPoseParameterization : public ceres::LocalParameterization {
        // 四元数定义顺序为, x, y, z, w
        // Quaternion order (x, y, z, w)

    public:
        bool Plus(const double *x, const double *delta, double *x_plus_delta) const override {
            Eigen::Map<const Eigen::Vector3d> _p(x);
            Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

            Eigen::Map<const Eigen::Vector3d> dp(delta);

            Eigen::Quaterniond dq = RotationVector2Quaternion(Eigen::Map<const Eigen::Vector3d>(delta + 3));

            Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
            Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

            p = _p + dp;
            q = (_q * dq).normalized();

            return true;
        }

        bool ComputeJacobian(const double *x, double *jacobian) const override {
            Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> j(jacobian);
            j.topRows<6>().setIdentity();
            j.bottomRows<1>().setZero();

            return true;
        }

        int GlobalSize() const override {
            return 7;
        }

        int LocalSize() const override {
            return 6;
        }
    };

	// GNSS Solution Factor
    class cGnssFactor_SOL : public ceres::SizedCostFunction<3, 7> {

    public:
        explicit cGnssFactor_SOL(tSolInfoUnit gnss, Vector3d lever)
                : gnss_(std::move(gnss))
                , lever_(std::move(lever)) {
        }

        void updateGnssState(const tSolInfoUnit &gnss) {
            gnss_ = gnss;
        }

        bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {
            Vector3d p{parameters[0][0], parameters[0][1], parameters[0][2]};
            Quaterniond q{parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]};


            Eigen::Map<Eigen::Matrix<double, 3, 1>> error(residuals);

            error = p + q.toRotationMatrix() * lever_ - gnss_.pos;

            Matrix3d weight = Matrix3d::Zero();
            weight(0, 0)    = 1.0 / gnss_.q_pos[0];
            weight(1, 1)    = 1.0 / gnss_.q_pos[1];
            weight(2, 2)    = 1.0 / gnss_.q_pos[2];

            error = weight * error;



            if (jacobians) {
                if (jacobians[0]) {
                    Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose(jacobians[0]);
                    jacobian_pose.setZero();

                    jacobian_pose.block<3, 3>(0, 0) = Matrix3d::Identity();
                    jacobian_pose.block<3, 3>(0, 3) = -q.toRotationMatrix() * VectorSkew(lever_);

                    jacobian_pose = weight * jacobian_pose;
                    cout<<"jacobian pos : "<<jacobian_pose<<endl;
                }
            }else{
              //  cout<<"JACOBIAN NONE"<<endl;
            }

            return true;
        }

    private:
        tSolInfoUnit gnss_;
        Vector3d lever_;
    };
}




#endif //FAST_LIO_IGCOUPLED_SOLVER_H
