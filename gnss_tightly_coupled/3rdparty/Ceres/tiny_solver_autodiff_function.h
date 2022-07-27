// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2017 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: mierle@gmail.com (Keir Mierle)
//
// WARNING WARNING WARNING
// WARNING WARNING WARNING  Tiny solver is experimental and will change.
// WARNING WARNING WARNING

#ifndef CERES_PUBLIC_TINY_SOLVER_AUTODIFF_FUNCTION_H_
#define CERES_PUBLIC_TINY_SOLVER_AUTODIFF_FUNCTION_H_

#include "Eigen/Core"

#include "Ceres/jet.h"
#include "Ceres/types.h"  // For kImpossibleValue.

namespace ceres {

// An adapter around autodiff-style CostFunctors to enable easier use of
// TinySolver. See the example below showing how to use it:
//
//   // Same as an autodiff cost functor, but taking only 1 parameter.
//   struct MyFunctor {
//     template<typename T>
//     bool operator()(const T* const parameters, T* residuals) const {
//       const T& x = parameters[0];
//       const T& y = parameters[1];
//       const T& z = parameters[2];
//       residuals[0] = x + 2.*y + 4.*z;
//       residuals[1] = y * z;
//       return true;
//     }
//   };
//
//   typedef TinySolverAutoDiffFunction<MyFunctor, 2, 3>
//       AutoDiffFunction;
//
//   MyFunctor my_functor;
//   AutoDiffFunction f(my_functor);
//
//   Vec3 x = ...;
//   TinySolver<AutoDiffFunction> solver;
//   solver.Solve(f, &x);
//
// WARNING: The cost function adapter is not thread safe.
template<typename CostFunctor,
         int kNumResiduals,
         int kNumParameters,
         typename T = double>
class TinySolverAutoDiffFunction {
 public:
   TinySolverAutoDiffFunction(const CostFunctor& cost_functor)
     : cost_functor_(cost_functor) {}

  typedef T Scalar;
  enum {
    NUM_PARAMETERS = kNumParameters,
    NUM_RESIDUALS = kNumResiduals,
  };

  // This is similar to AutoDiff::Differentiate(), but since there is only one
  // parameter block it is easier to inline to avoid overhead.
  bool operator()(const T* parameters,
                  T* residuals,
                  T* jacobian) const {
    if (jacobian == NULL) {
      // No jacobian requested, so just directly call the cost function with
      // doubles, skipping jets and derivatives.
      return cost_functor_(parameters, residuals);
    }
    // Initialize the input jets with passed parameters.
    for (int i = 0; i < kNumParameters; ++i) {
      jet_parameters_[i].a = parameters[i];  // Scalar part.
      jet_parameters_[i].v.setZero();        // Derivative part.
      jet_parameters_[i].v[i] = T(1.0);
    }

    // Initialize the output jets such that we can detect user errors.
    for (int i = 0; i < kNumResiduals; ++i) {
      jet_residuals_[i].a = kImpossibleValue;
      jet_residuals_[i].v.setConstant(kImpossibleValue);
    }

    // Execute the cost function, but with jets to find the derivative.
    if (!cost_functor_(jet_parameters_, jet_residuals_)) {
      return false;
    }

    // Copy the jacobian out of the derivative part of the residual jets.
    Eigen::Map<Eigen::Matrix<T,
                             kNumResiduals,
                             kNumParameters> > jacobian_matrix(jacobian);
    for (int r = 0; r < kNumResiduals; ++r) {
      residuals[r] = jet_residuals_[r].a;
      // Note that while this looks like a fast vectorized write, in practice it
      // unfortunately thrashes the cache since the writes to the column-major
      // jacobian are strided (e.g. rows are non-contiguous).
      jacobian_matrix.row(r) = jet_residuals_[r].v;
    }
    return true;
  }

 private:
  const CostFunctor& cost_functor_;

  // To evaluate the cost function with jets, temporary storage is needed. These
  // are the buffers that are used during evaluation; parameters for the input,
  // and jet_residuals_ are where the final cost and derivatives end up.
  //
  // Since this buffer is used for evaluation, the adapter is not thread safe.
  mutable Jet<T, kNumParameters> jet_parameters_[kNumParameters];
  mutable Jet<T, kNumParameters> jet_residuals_[kNumResiduals];
};

}  // namespace ceres

#endif  // CERES_PUBLIC_TINY_SOLVER_AUTODIFF_FUNCTION_H_
