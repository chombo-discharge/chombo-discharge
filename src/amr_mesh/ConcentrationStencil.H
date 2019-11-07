/*!
  @file   ConcentrationStencil.H
  @author Robert Marskar
  @date   2019
*/

#ifndef _CONCENTRATIONSTENCIL_H_
#define _CONCENTRATIONSTENCIL_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "BaseIVFAB.H"
#include "EBLevelRedist.H"

/*!
   The ConcentrationStencil class holds the stencil for mass injection from regular cells
   and into cut cells in order to avoid the appearance of negative densities
  */
class ConcentrationStencil {
public:

  /*!
     Default constructor.  User must subsequently call define().
  */
  ConcentrationStencil();

  /*!
    @brief Deconstructor
  */
  ~ConcentrationStencil();


  /*!
    @brief Defining constructor - calls define
  */
  ConcentrationStencil(const DisjointBoxLayout& a_dbl,
                const EBISLayout&        a_ebisl,
                const ProblemDomain&     a_domain,
                const int&               a_redistRadius);

  /*
    @brief Define fucntion
  */
  void define(const DisjointBoxLayout& a_dbl,
              const EBISLayout&        a_ebisl,
              const ProblemDomain&     a_domain,
              const int&               a_redistRadius,
              bool                     a_do2DStencil = false);

  /*!
     @brief Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  /*!
    @brief Modify the weights
  */
  void resetWeights(const LevelData<EBCellFAB>& a_modifier, const int& a_ivar);

  /*!
     @brief Returns the redistribution stencil at every irregular
     point in input  Box associated with  this DataIndex.
  */
  const BaseIVFAB<VoFStencil >& operator[](const DataIndex& datInd) const;

  /*!
     @brief Deep copy This copy operation detaches copier from copyee.
   */
  void deepCopy(const ConcentrationStencil& a_stenin);
  
  /*!
    @brief Get redistribution radius
  */
  int getRedistRadius() const;

protected:

  LayoutData<BaseIVFAB<VoFStencil > > m_stencil;
  LayoutData<BaseIVFAB<VoFStencil > > m_volsten;

  DisjointBoxLayout m_grids;
  ProblemDomain     m_domain;
  EBISLayout        m_ebisl;

  bool m_isDefined;
  bool m_hasDefaultWeights;
  bool m_alwaysRedist;

  int  m_redistRadius;

private:
  //internal use
  void computePointStencil(VoFStencil&      a_stencil,
                           const VolIndex&  a_srcVoF,
                           const DataIndex& a_datInd,
                           const bool&      a_do2DStencil);

  //disallowed for all the usual reasons
  ConcentrationStencil(const ConcentrationStencil& ebcin)
  {
    MayDay::Error("rdsi 2 invalid operator");
  }
  void operator=(const ConcentrationStencil& fabin)
  {
    MayDay::Error("rdsi 3 invalid operator");
  }
};

#endif
