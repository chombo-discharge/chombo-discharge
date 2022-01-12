/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaJSONFunctions.H
  @brief  Declaration of various functions for usage with CD_CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaJSONFunctions_H
#define CD_CdrPlasmaJSONFunctions_H

// Std includes
#include <functional>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Parent function that encapsulates a function which returns a quantity f = f(E,N). 
      @details Users must override the operator function to return a quantity which is a function of E and N.
    */
    class CdrPlasmaJSONFunctionEN {
    public:

      /*!
	@brief Constructor. 
      */
      CdrPlasmaJSONFunctionEN() = default;

      /*!
	@brief Destructor. 
      */      
      virtual ~CdrPlasmaJSONFunctionEN() = default;

      /*!
	@brief Evaluation function
      */
      virtual Real operator()(const RealVect a_E, const Real a_N) const = 0;
    };

    /*!
      @brief Implementation of CdrPlasmaJSONFunctionEN which returns a function A * |E|^B / N^C
    */
    class CdrPlasmaJSONFunctionENA : public CdrPlasmaJSONFunctionEN {
    public:

      /*!
	@brief Constructor
	@param[in] a_A A-coefficient
	@param[in] a_B B-coefficient
	@param[in] a_C C-coefficient
      */
      CdrPlasmaJSONFunctionENA(const Real a_A, const Real a_B, const Real a_C);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaJSONFunctionENA();

      /*!
	@brief Implementation of parent function. 
      */
      Real operator()(const RealVect a_E, const Real a_N) const override final;

    protected:
      Real m_A;
      Real m_B;
      Real m_C;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
