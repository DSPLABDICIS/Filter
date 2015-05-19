/*
 * predictModel.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: luzdora
 */

#include "predictModel.hpp"

using namespace boost_incl;

//using namespace jblas;
using namespace filter;


/*
 * class LinearPredictModel
 */

LinearPredictModel::LinearPredictModel(mat const& F_,
                                       sym_mat const& Q_) :
  BasePredictModel(F_.size1())
{
  //JFR_INVARIANT(F.size1() == sizeState() && F.size2() == sizeState(),
  //              "LinearPredictModel::LinearPredictModel: size of F does not match");
  // JFR_INVARIANT(Q.size1() == sizeState() && Q.size2() == sizeState(),
  //              "LinearPredictModel::LinearPredictModel: size of Q does not match");
  F.assign(F_);
  Q.assign(Q_);
}

/*
 * class LinearCommandPredictModel
 */

LinearCommandPredictModel::LinearCommandPredictModel(const mat& F_,
						     const mat& G_,
						     const sym_mat& Q_) :
  BaseCommandPredictModel(F_.size1(), G_.size2())
{
  //JFR_INVARIANT(F.size1() == sizeState() && F.size2() == sizeState(),
   //             "LinearCommandPredictModel::LinearCommandPredictModel: size of F does not match");
  //JFR_INVARIANT(Q.size1() == sizeState() && Q.size2() == sizeState(),
   //             "LinearCommandPredictModel::LinearCommandPredictModel: size of Q does not match");
  //JFR_INVARIANT(G.size1() == sizeState() && G.size2() == sizeCommand(),
  //              "LinearCommandPredictModel::LinearCommandPredictModel: size of Q does not match");
  F.assign(F_);
  Q.assign(Q_);
  G.assign(G_);
}
