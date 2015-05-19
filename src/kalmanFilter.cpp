/*
 * kalmanFilter.cpp
 *
 *  Created on: Feb 6, 2013
 *      Author: luzdora
 */


#include <cassert>

#include "boost/version.hpp"
//#include "jmath/ublasExtra.hpp"
#include "kalmanFilter.hpp"
//#include "filterException.hpp"

using namespace boost_incl;
//using namespace jblas;
using namespace filter;
using namespace boost::numeric::ublas;
/*
 * BaseKalmanFilter class
 */

BaseKalmanFilter::BaseKalmanFilter(std::size_t size_) :
  covUpdateType(STANDARD), x(size_), P(size_, size_),
  p_curTime(0,0,0,0),
  K(size_,1),
  KH_tmp(size_, size_),
  consistencyCheckLevel(CONSISTENCY_NONE),
  chi2Threshold()
{
  x.clear();
  P.clear();
}

BaseKalmanFilter::~BaseKalmanFilter() {}

void BaseKalmanFilter::setCovUpdateType(CovUpdateType covUpdateType_) {
  covUpdateType = covUpdateType_;
}

BaseKalmanFilter::CovUpdateType BaseKalmanFilter::getCovUpdateType() const
{
  return covUpdateType;
}


void BaseKalmanFilter::clear() {
  x.clear();
  P.clear();
}

void BaseKalmanFilter::init(const vec& state_, const sym_mat& cov_, boost::posix_time::time_duration const& curTime) {
//  JFR_PRECOND(state_.size() == cov_.size1() && state_.size()==cov_.size2(),
//              "BaseKalmanFilter::init: size of init state vector and covariance matrix does not match");

  x.assign(state_);
  P.assign(cov_);
  p_curTime = curTime;
}

void BaseKalmanFilter::init(const vec& state_, const sym_mat& cov_) {
  init(state_, cov_, boost::posix_time::microseconds(0));
}

void BaseKalmanFilter::init(boost::posix_time::time_duration const& curTime) {
  x.clear();
  P.clear();
  p_curTime = curTime;
}

/*
void BaseKalmanFilter::writeLogHeader(kernel::DataLogger& log) const
{
  log.writeComment("filter: Kalman filter");
  log.writeLegend("current time (ms)");
  std::ostringstream os;
  for (std::size_t i = 0 ; i < x.size() ; ++i) {
    os << "x" << i << " ";
  }
  log.writeComment("state vector: ");
  log.writeLegendTokens(os.str());
  log.writeComment("standard deviations: ");
  log.writeLegendTokens(os.str());
}

void BaseKalmanFilter::writeLogData(kernel::DataLogger& log) const
{
  log.writeData(p_curTime.total_milliseconds());

  // state vector
  log.writeDataVector(x);

  // standard deviations
  for (std::size_t i = 0 ; i < x.size() ; ++i) {
    log.writeData(sqrt(P(i,i)));
  }
}
*/

//=========================================================

/*
 * KalmanFilter class
 */

KalmanFilter::KalmanFilter(std::size_t size_) :
  BaseKalmanFilter(size_) {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::predict(const LinearPredictModel& m_) {
  //JFR_PRECOND(m_.sizeState() == x.size(),
  //            "KalmanFilter::predict: size of predict model does not match.");

  //std::cout<<std::endl<<x(0)<<"  "<<x(1)<<"  "<<x(2)<<"  "<<x(3)<<"  "<<std::endl;
  x = prod(m_.F,x);
  //std::cout<<m_.F(0,0)<<m_.F(0,1)<<m_.F(0,2)<<m_.F(0,3)<<m_.F(1,0)<<m_.F(1,1)<<m_.F(1,2)<<m_.F(1,3)<<std::endl;
  //std::cout<<x(0)<<"  "<<x(1)<<"  "<<x(2)<<"  "<<x(3)<<"  "<<std::endl;

  P = prod(m_.F, mat(prod(P,trans(m_.F)))) + m_.Q;
  p_curTime += m_.getTimeStep();
}

void KalmanFilter::predict(LinearCommandPredictModel& m, vec const& u) {
//  JFR_PRECOND(m.sizeState() == x.size(),
//              "KalmanFilter::predict: size of predict model does not match.");
//  JFR_PRECOND( u.size() == m.sizeCommand(),
//	       "KalmanFilter::predict: size of u does not match.");

  m.computeUCov(u);

  x = prod(m.F, x) + prod(m.G, u);
  P = prod(m.F, mat(prod(P, trans(m.F)))) + prod(m.G, mat(prod(m.getUCov(), trans(m.G)))) + m.Q;
  p_curTime += m.getTimeStep();
}

void KalmanFilter::update(LinearObserveModel& m_, vec& z_) {
 // JFR_PRECOND(m_.sizeObs() == z_.size() && m_.sizeState() == x.size(),
 //             "KalmanFilter::update: size of m_ does not match.");

  std::size_t sizeState = x.size();
  std::size_t sizeObs = z_.size();

  switch(m_.isCorrelated()) {
  case true:
   // JFR_RUN_TIME("KalmanFilter::update: correlated noise not implemented yet");
	  std::cout << "KalmanFilter::update: correlated noise not implemented yet";
    break;

  case false:
	  mat z (sizeObs, sizeObs);
	       for (std::size_t j=0; j<sizeObs; j++)
	       z(j,0)= z_(j);

    for (std::size_t i = 0 ; i<sizeObs ; i++) {
      mat_range hr(m_.H, range(i,i+1), range(0,sizeState));
     // vec_range zr(z_, 2);

      mvec_range zr(z, range(i,i+1), range(0,1));

      updateLinear(hr, zr, m_.getR()(i,i));
    }
    break;
  }
}

/*void KalmanFilter::updateLinear(const mat_range& H_,
				const vector_range<vec>& z_,
				double R_)
				*/
void KalmanFilter::updateLinear(const mat_range& H_,
				const mvec_range & z_,
				double R_)
{

  vec aux = prod(H_,x);
  vec y = z_ - aux;
 // vec_range y;

  //y = z_;

  sym_mat S ( prod<sym_mat>(H_, prod<mat>(P, trans(H_))) );
  S(0,0) += R_;
  sym_mat Sinv(1,1);
  Sinv(0,0) = 1/S(0,0);

  checkConsistency(y, Sinv);

  K.assign( prod(P , trans(H_)*Sinv(0,0)) );

  x.plus_assign( prod(K, y) );

  switch(covUpdateType) {
  case SIMPLE:
    // FIX_ME
    P -= prod(mat(prod(K,H_)), P);
    break;
  case STANDARD:
    // Bar-Shalom p294
    P.minus_assign(prod(K, mat(prod(S, trans(K)))));
    break;
  case JOSEPH:
    KH_tmp.assign(identity_mat(x.size()) - prod(K, H_));
    P = prod(KH_tmp, mat(prod(P, trans(KH_tmp))))
      + R_ * prod(K, trans(K));
    break;
  }

}

//int main()
//{
//	std::cout << "Pues este es main" << std::endl;
//	return 0;
//}

/*
std::ostream& filter::operator<<(std::ostream& s, const BlockExtendedKalmanFilter& f_) {
  s << (*f_.x_r) << std::endl;
  s << (*f_.P_r) << std::endl;
  return s;
}
*/
