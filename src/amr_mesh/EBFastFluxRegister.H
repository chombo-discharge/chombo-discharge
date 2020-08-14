/*!
  @file   EBFastFluxRegister.H
  @brief  Declaration of a new EBFluxRegister that uses the optimizations in EBFasterFR
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _EBFastFluxRegister_H_
#define _EBFastFluxRegister_H_

#include "EBFasterFR.H"

class EBCoarToFineRedist;
class EBCoarToCoarRedist;

/*!
  @brief Class that does refluxing and interaction with redistribution classes. 
*/
class EBFastFluxRegister : public EBFasterFR {
public:

  
  /*!
    @brief Default constructor.  Leaves object undefined.
  */
  EBFastFluxRegister();

  /*!
    @brief Full constructor.  Calls the define function which creates
    a levels worth of flux registers.
  */
  EBFastFluxRegister(const DisjointBoxLayout& a_dblFine,
		     const DisjointBoxLayout& a_dblCoar,
		     const EBISLayout&        a_ebislFine,
		     const EBISLayout&        a_ebislCoar,
		     const Box&               a_domainCoar,
		     const int&               a_nref,
		     const int&               a_nvar,
		     const EBIndexSpace*      ebisPtr,
		     const bool               a_forceNoEBCF = false);

  /*!
    @brief Dtor
  */
  virtual ~EBFastFluxRegister();

  /*!
    @brief Full define function.  Creates a levels worth of flux registers.
  */
  void define(const DisjointBoxLayout& a_dblFine,
	      const DisjointBoxLayout& a_dblCoar,
	      const EBISLayout&        a_ebislFine,
	      const EBISLayout&        a_ebislCoar,
	      const ProblemDomain&     a_domainCoar,
	      const int&               a_nref,
	      const int&               a_nvar,
	      const EBIndexSpace*      ebisPtr,
	      const bool               a_forceNoEBCF = false);

  /*!
    @brief Increments the register with data from coarseFlux, multiplied by scale.
    Increment coarse buffer with the area-weighted sum
    of the fluxes on the faces.
    buf += -scale*sum(areaFrac*flux)
    CoarseFlux must contain the coarse fluxes in the dir direction
    for the grid   m_coarseLayout[coarseDataIndex].
  */
  void incrementCoarseRegular(const EBFaceFAB&   a_coarseFlux,
			      const Real&        a_scale,
			      const DataIndex&   a_coarseDataIndex,
			      const Interval&    a_variables,
			      const int&         a_dir){
    for (SideIterator sit; sit.ok(); ++sit){
      incrementCoarRegul(a_coarseFlux,
			 a_scale,
			 a_coarseDataIndex,
			 a_variables,
			 a_dir, sit());
    }
  }

  /*!
    @brief Increment coarse
  */
  void incrementCoarse(const EBFaceFAB&   a_coarseFlux,
                       const Real&        a_scale,
                       const DataIndex&   a_coarseDataIndex,
                       const Interval&    a_varSrc,
                       const Interval&    a_varDst,
                       const int&         a_dir){
    CH_assert(a_varSrc==a_varDst);
    for (SideIterator sit; sit.ok(); ++sit){
      incrementCoarseBoth(a_coarseFlux, a_scale, a_coarseDataIndex, a_varSrc, a_dir, sit());
    }
  }


  /*!
    @brief increments the register with data from coarseFlux, multiplied by scale.
    Increment coarse buffer with the area-weighted sum
    of the fluxes on the faces.
    buf += -scale*sum(areaFrac*flux)
    CoarseFlux must contain the coarse fluxes in the dir direction
    for the grid   m_coarseLayout[coarseDataIndex].
  */
  void incrementCoarseIrregular(const BaseIFFAB<Real>& a_coarseFlux,
				const Real&            a_scale,
				const DataIndex&       a_coarseDataIndex,
				const Interval&        a_variables,
				const int&             a_dir);

  /*!
    @brief For when you have an EBFaceFAB lying around instead of a BaseIFFAB
  */
  void incrementCoarseIrregular(const EBFaceFAB&       a_coarFlux,
                                const Real&            a_scale,
                                const DataIndex&       a_coarDatInd,
                                const Interval&        a_variables,
                                const int&             a_dir){
    for (SideIterator sit; sit.ok(); ++sit){
      incrementCoarIrreg(a_coarFlux,
			 a_scale,
			 a_coarDatInd,
			 a_variables,
			 a_dir,
			 sit());
    }
  }

  /*!
    @brief increments the register with data from fineFlux, multiplied by scale.
    fineFlux must contain the coarse fluxes in the dir direction
    for the grid m_fineLayout[fineDataIndex].
    This operation is local.
  */
  void incrementFineRegular(const EBFaceFAB&      a_fineFlux,
			    const Real&           a_scale,
			    const DataIndex&      a_fineDataIndex,
			    const Interval&       a_variables,
			    const int&            a_dir,
			    const Side::LoHiSide& a_sd){
    incrementFineRegul(a_fineFlux,
                       a_scale,
                       a_fineDataIndex,
                       a_variables,
                       a_dir,
                       a_sd);
  }

  /*!
    @brief Increment fine
  */
  void incrementFine(const EBFaceFAB&      a_fineFlux,
		     const Real&           a_scale,
		     const DataIndex&      a_fineDataIndex,
		     const Interval&       a_varSrc,
		     const Interval&       a_varDst,
		     const int&            a_dir){
    CH_assert(a_varSrc == a_varDst);
    for (SideIterator sit; sit.ok(); ++sit){
      incrementFineBoth(a_fineFlux, a_scale, a_fineDataIndex, a_varSrc, a_dir, sit());
    }
  }

  /*!
    @brief increments the register with data from fineFlux, multiplied by scale.
    fineFlux must contain the coarse fluxes in the dir direction
    for the grid m_fineLayout[fineDataIndex].
    This operation is local.
  */
  void incrementFineIrregular(const BaseIFFAB<Real>& a_fineFlux,
			      const Real&            a_scale,
			      const DataIndex&       a_fineDataIndex,
			      const Interval&        a_variables,
			      const int&             a_dir,
			      const Side::LoHiSide&  a_sd);

  /*!
    @brief For when you have a EBFaceFAB lying around instead of a BaseIFFAB
  */
  void incrementFineIrregular(const EBFaceFAB&       a_fineFlux,
			      const Real&            a_scale,
			      const DataIndex&       a_fineDataIndex,
			      const Interval&        a_variables,
			      const int&             a_dir,
			      const Side::LoHiSide&  a_sd){
    incrementFineIrreg(a_fineFlux,
                       a_scale,
                       a_fineDataIndex,
                       a_variables,
                       a_dir,
                       a_sd);
  }

  /*!
     @brief Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToFineRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);

  /*!
     @brief Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToCoarRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);


  /*!
    @brief Check if object is defined or not
  */
  bool isDefined() const;

protected:


private:

  /*!
    @brief Returns true if a_src contains any cells
  */
  bool copyBIFFToEBFF(EBFaceFAB&             a_dst,
		      const BaseIFFAB<Real>& a_src,
		      const Box            & a_box,
		      const EBISBox&         a_ebisBox);

  //for all the usual reasons,
  //there is no copy constructor for this class
  //and there is no operator= for this class
  //(there is no rule six).
  void operator= (const EBFastFluxRegister&){
    MayDay::Error("EBFastFluxRegister::invalid operator");
  }
  
  EBFastFluxRegister(const EBFastFluxRegister&){
    MayDay::Error("EBFastFluxRegister::invalid operator");
  }
};

#endif

