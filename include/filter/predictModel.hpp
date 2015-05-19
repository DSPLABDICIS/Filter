/*
 * predictModel.hpp
 *
 *  Created on: Feb 8, 2013
 *      Author: luzdora
 */

#ifndef PREDICTMODEL_HPP_
#define PREDICTMODEL_HPP_

#include "../../../KLT/include/boost_incl.hpp"


#include "boost/date_time/posix_time/posix_time_types.hpp"

//#include "jmath/jblas.hpp"

  namespace filter {

  using namespace boost_incl;
    class KalmanFilter;
  //  class ExtendedKalmanFilter;
   // class BlockExtendedKalmanFilter;

    /** Base class for all prediction models.
     *
     * \ingroup filter
     */
    class BasePredictModel {

    private:

      std::size_t m_sizeState;

    protected:

      /// timeStep in seconds (default value is 1.0)
      boost::posix_time::time_duration m_timeStep;

    public:

      /// linear(ized) prediction function
      mat F;

      /// process noise covariance
      sym_mat Q;

      BasePredictModel(std::size_t sizeState) :
	m_sizeState(sizeState),
	m_timeStep(0,0,1,0),
	F(sizeState,sizeState),
	Q(sizeState,sizeState)
      {}

      virtual ~BasePredictModel() {}

      std::size_t sizeState() const {return m_sizeState;};

      boost::posix_time::time_duration const& getTimeStep() const {return m_timeStep;};

      /** Change value of timeStep.
       */
      virtual void setTimeStep(boost::posix_time::time_duration const& timeStep) {
	m_timeStep = timeStep;
      }

    }; // class BasePredictModel

    class BaseCommandPredictModel : public BasePredictModel {

       private:
         std::size_t m_sizeCommand;

       protected:
         /// covariance of the command
         sym_mat m_uCov;

       public:

         /// linear(ized) command function
         mat G;

         BaseCommandPredictModel(std::size_t sizeState, std::size_t sizeCommand) :
   	BasePredictModel(sizeState),
   	m_sizeCommand(sizeCommand),
   	m_uCov(sizeCommand, sizeCommand),
   	G(sizeState, sizeCommand)
         {
   	m_uCov.clear();
         }

         virtual ~BaseCommandPredictModel() {}

         std::size_t sizeCommand() const {return m_sizeCommand;}

         void setUCov(sym_mat const& uCov) {
   	m_uCov.assign(uCov);
         }

        sym_mat const& getUCov() const {return m_uCov;}

         /// compute the covariance of the command, by default do nothing.
         virtual void computeUCov(vec const& u) {}

       };



    /** Linear prediction model.
     * x' = F*x + v
     *
     * \ingroup filter
     */
    class LinearPredictModel : public BasePredictModel {

    public:

      LinearPredictModel(std::size_t sizeState) :
	BasePredictModel(sizeState) {}

      /** This constructor is usefull in case \a F and \a Q are
       * constant.
       */
      LinearPredictModel(const mat& F_,
                         const sym_mat& Q_);

      virtual ~LinearPredictModel() {}

    }; // class LinearPredictModel

    /** Linear prediction model with command.
     * x' = F*x + v + G*(u+w)
     *
     * \ingroup filter
     */

    class LinearCommandPredictModel : public BaseCommandPredictModel {

    public:

      LinearCommandPredictModel(std::size_t sizeState, std::size_t sizeCommand) :
	BaseCommandPredictModel(sizeState, sizeCommand) {}

      LinearCommandPredictModel(const mat& F_,
				const mat& G_,
				const sym_mat& Q_);

      virtual ~LinearCommandPredictModel() {}

    };

  } // namespace filter

#endif /* PREDICTMODEL_HPP_ */
