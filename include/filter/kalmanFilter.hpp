/*
 * kalmanFilter.hpp
 *
 *  Created on: Feb 6, 2013
 *      Author: luzdora
 */

#ifndef KALMANFILTER_HPP_
#define KALMANFILTER_HPP_
#include <iostream>
#include "predictModel.hpp"
#include "observeModel.hpp"
#include "../../../KLT/include/boost_incl.hpp"
#include "boost/date_time/posix_time/posix_time_types.hpp"

 namespace filter {

 using namespace boost_incl;
 class BaseKalmanFilter {

    public:

      /// Different kind of covariance update (Bar-Shalom p294)
      enum CovUpdateType {SIMPLE,   /**< very naive update, from the books */
                          STANDARD, /**< still numerical errors in substraction */
                          JOSEPH    /**< Joseph form, less sensitive
                                       to round-off error, preserve
                                       positive definiteness,
                                       computationnaly more
                                       expensive */
      };

      /// Different level of consistency check
      enum ConsistencyCheckLevel {CONSISTENCY_NONE,     ///< no check
				  CONSISTENCY_WARNING,  ///< send a warning to the output
				  CONSISTENCY_EXCEPTION ///< throw a filter exception with id=INCONSISTENT_UPDATE
      };

    protected:

      /// the covariance update method
      CovUpdateType covUpdateType;

      /// state vector of the filter
      vec x;

      /// covariance matrix of the state vector
      sym_mat P;

      /// filter current time.
      boost::posix_time::time_duration p_curTime;

      /// Kalman gain for mono dimension update
     mat K;

      /// preallocated matrix to speed-up computation (used by Joseph update, and simple update)
     mat KH_tmp;

      /// whether the consistency is checked
      ConsistencyCheckLevel consistencyCheckLevel;
      /// consistency threshold \sa checkConsistency()
      double chi2Threshold;

      /** Check the consistency of an update. Innovation \a y is a
       * centered gaussian distribution with covariance S.
       */
      template<class Vec>
      void checkConsistency(Vec const& y, sym_mat const& Sinv)
      {
	if (consistencyCheckLevel == CONSISTENCY_NONE)
	  return;

//	JFR_PRECOND(Sinv.size1() == y.size() && Sinv.size2() == y.size(),
//		    "BaseKalmanFilter::checkConsistency: invalid size of S");

//	using namespace boost_incl;
	double ySiy = inner_prod(y, prod(Sinv, y) );
	checkConsistency(y, Sinv, ySiy);
      }

      void checkConsistency(double y, double Sinv) {
	if (consistencyCheckLevel == CONSISTENCY_NONE)
	  return;

	double ySiy = y * Sinv * y;
	checkConsistency(y, Sinv, ySiy);
      }

      template<typename YType, typename SinvType>
      void checkConsistency(YType y, SinvType Sinv, double ySiy) {
	if (ySiy > chi2Threshold) {
	  switch(consistencyCheckLevel) {
	  case CONSISTENCY_NONE:
	    // shut up gcc warning
	    break;
	  case CONSISTENCY_WARNING:
	    std::cout << "ySiy=" << ySiy << " threshold=" << chi2Threshold << " y=" << y << "LaSinv=" << Sinv  << ")"<< std::endl;
	    break;
	/*  case CONSISTENCY_EXCEPTION:
	    std::ostringstream s;
	    s << "ySiy=" << ySiy << " threshold=" << chi2Threshold
	      << " y=" << y << " Sinv=" << Sinv  << ")";
	    throw InconsistentUpdateException(ySiy, chi2Threshold, s.str(), __FILE__, __LINE__);
	    break; */
	  default: break;

	  }
	}
      }

      // logging method
    //  void writeLogHeader(jafar::kernel::DataLogger& log) const;
    //  void writeLogData(jafar::kernel::DataLogger& log) const;

    public:

       /** Create a kalman filter.
       * @param size_ The size of the state
       */
      BaseKalmanFilter(std::size_t size_);

      virtual ~BaseKalmanFilter();

      vec const& getX() const {return x;};
      vec& getX() {return x;};
      sym_mat const& getP() const {return P;};
      sym_mat& getP() {return P;};

      /// clear the state vector and the covariance matrix
      void clear();

      boost::posix_time::time_duration const& getCurrentTime() const {return p_curTime;}

      void setCovUpdateType(CovUpdateType covUpdateType_);
      CovUpdateType getCovUpdateType() const;

      void setupConsistencyCheck(ConsistencyCheckLevel level_,
				 double chi2Threshold_)
      {
	consistencyCheckLevel = level_; chi2Threshold = chi2Threshold_;
      }

      ConsistencyCheckLevel getConsistencyCheckLevel() const
      {
        return consistencyCheckLevel;
      }

      double getChi2Threshold() const
      {
        return chi2Threshold;
      }


      /** Init the filter with \a state_ and covariance \a cov_ and \c curTime.
       */
      //void init(const vec& state_, const sym_mat& cov_, boost::posix_time::time_duration const& curTime);

      void init(const vec& state_, const sym_mat& cov_, boost::posix_time::time_duration const& curTime);

      /** Init the filter with \a state_ and covariance \a cov_,
       * current time is set to 0.
       */

      void init(const vec& state_, const sym_mat& cov_);

      /** Init the filter at \c curTime, state is set to zero.
       */
      void init(boost::posix_time::time_duration const& curTime);

      /** @return size of the filter state.
       */
      virtual std::size_t sizeState() const {return x.size();};

      //      void resize(std::size_t size_);

    }; // class BaseKalmanFilter


 class KalmanFilter : public BaseKalmanFilter {

   public :

     /** Create a kalman filter.
      * @param size_ dimension of the state
      */
     KalmanFilter(std::size_t size_);

     virtual ~KalmanFilter();

     /** This function computes the prediction step.
      * @param m_ linear predict model
      */
     void predict(const LinearPredictModel& m_);

#ifndef SWIG
     /** This function computes the prediction step.
      *
      * @param m_ linear command predict model
      */
     void predict(LinearCommandPredictModel& m_, vec const& u_);
#endif // SWIG

     /** This function computes the observation step.
      * @param m_ linear observation model
      * @param z_ measure
      */
     void update(LinearObserveModel& m_, vec& z_);

   protected :

     /** This function computes the observation step (observation is given). In the case of a scalar observation.
      * @param H_ observation matrix (only one row)
      * @param z_ measure
      * @param R_ measure noise variance
      */
     void updateLinear(const mat_range& H_,
                       const mvec_range& z_, double R_);

   }; // class KalmanFilter

 }


#endif /* KALMANFILTER_HPP_ */
