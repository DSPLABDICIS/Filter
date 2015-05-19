/*
 * observeModel.hpp
 *
 *  Created on: Feb 8, 2013
 *      Author: luzdora
 */

#ifndef OBSERVEMODEL_HPP_
#define OBSERVEMODEL_HPP_

#include "../../../KLT/include/boost_incl.hpp"


 namespace filter {

 using namespace boost_incl;
    /** Base class for observation models.
     *
     * \ingroup filter
     */
    class BaseObserveModel {

    public:

      BaseObserveModel(std::size_t sizeObs_, std::size_t sizeState_);
      BaseObserveModel(std::size_t sizeObs_, std::size_t sizeState_, std::size_t sizeInnovation_, std::size_t sizePrediction_);

      virtual ~BaseObserveModel();

      inline std::size_t sizeObs() const {return p_sizeObs;}
      inline std::size_t sizeState() const {return p_sizeState;}
      inline std::size_t sizeInnovation() const {return p_sizeInnovation;}
      inline std::size_t sizePrediction() const {return p_sizePrediction;}
      inline bool isCorrelated() const {return p_isCorrelated;}

      virtual sym_mat const& getR() const {return R;}

      /// Set a constant noise model
      void setCorrelatedR(sym_mat const& R_);

      /// Set a constant noise model
      void setUncorrelatedR(vec const& R_);

//       /** Compute observation noise covariance. A default
//        * implementation is provided, it does nothing (constant model).
//        */
//       virtual void computeR(jblas::vec const& z_);

    protected:
      std::size_t p_sizeObs;
      std::size_t p_sizeState;
      std::size_t p_sizeInnovation;
      std::size_t p_sizePrediction;

      bool p_isCorrelated;

      /// Covariance of observation noise
      sym_mat R;

      /// pre-allocated vector for predicted obervation or innovation.
      mutable vec z;

    }; // class BaseObserveModel

    //======================================================

    /** Linear observation model.
     *
     * \ingroup filter
     */
    class LinearObserveModel : public BaseObserveModel {

    public:

      /// linear observation function
      mat H;

      LinearObserveModel(std::size_t sizeObs_, std::size_t sizeState_);
      LinearObserveModel(mat const& H_);



      virtual vec const& predictObservation(const vec& x_) const;


      //       virtual jblas::vec const& computeInnovation(const jblas::vec& z_,
      //					   const jblas::vec& x_) const;
      ~LinearObserveModel();

    }; // class LinearObserveModel


  } // namespace filter

#endif /* OBSERVEMODEL_HPP_ */
