/*
 * observeModel.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: luzdora
 */


#include "observeModel.hpp"

//using namespace jblas;
using namespace filter;
using namespace boost::numeric::ublas;
/*
 * class BaseObserveModel
 */

BaseObserveModel::BaseObserveModel(std::size_t sizeObs_, std::size_t sizeState_) :
  p_sizeObs(sizeObs_),
  p_sizeState(sizeState_),
  p_sizeInnovation(sizeObs_),
  p_sizePrediction(sizeObs_),
  p_isCorrelated(false),
  R(sizeObs_, sizeObs_),
  z(sizeObs_)
{}

BaseObserveModel::BaseObserveModel(std::size_t sizeObs_, std::size_t sizeState_, std::size_t sizeInnovation_, std::size_t sizePrediction_) :
  p_sizeObs(sizeObs_),
  p_sizeState(sizeState_),
  p_sizeInnovation(sizeInnovation_),
  p_sizePrediction(sizePrediction_),
  p_isCorrelated(false),
  R(sizeInnovation_, sizeInnovation_),
  z(sizeInnovation_)
{
}

BaseObserveModel::~BaseObserveModel() {}

void BaseObserveModel::setCorrelatedR(sym_mat const& R_) {
  if (R_.size1() == sizeObs())
	      std::cout << "BaseObserveModel setCorrelatedR: size of R_ does not match";
  p_isCorrelated = true;
  R.assign(R_);
}

void BaseObserveModel::setUncorrelatedR(vec const& R_) {
  if(R_.size() == sizeObs())
	    std::cout << "BaseObserveModel setUncorrelatedR: size of R_ does not match";
  R.clear();
  p_isCorrelated = false;
  for (std::size_t i=0 ; i < sizeObs() ; ++i) {
    R(i,i) = R_(i);
  }
}

/*
 * class LinearObserveModel
 */


LinearObserveModel::LinearObserveModel(std::size_t sizeObs_, std::size_t sizeState_) :
  BaseObserveModel(sizeObs_, sizeState_),
  H(sizeObs_, sizeState_)
{}

LinearObserveModel::LinearObserveModel(mat const& H_) :
  BaseObserveModel(H_.size1(), H_.size2()),
  H(H_)
{}

LinearObserveModel:: ~LinearObserveModel() {}

vec const& LinearObserveModel::predictObservation(const vec& x_) const
{
  if(x_.size() == sizeState())
            std::cout <<  "LinearObserveModel::predictObservation: size of x_ does not match";
  z.assign(prod(H, x_));
  return z;
}

