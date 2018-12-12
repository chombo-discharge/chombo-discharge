/*!
  @file   DomainFluxIFFAB.H
  @brief  Declaration of a BaseIFFAB wrapper that holds domain fluxes 
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _DOMAINFLUXIFFAB_
#define _DOMAINFLUXIFFAB_

#include <BaseIFFAB.H>
#include <ProblemDomain.H>
#include <EBISLayout.H>

/*!
  @brief Class for holding fluxes on domain boundaries
  @details The internal functionality of this class is to use two SpaceDim-sized tuples that hold the fluxes in each direction
*/
class DomainFluxIFFAB {
public:

  static int preAllocatable(){
    return 1; // Symmetric messaging
  }
  
  DomainFluxIFFAB();
  DomainFluxIFFAB(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_ncomp);
  ~DomainFluxIFFAB();

  // Inline this shit
  const Box& box() const{ return m_box; }
  int nComp() const { return m_nComp; }
  const ProblemDomain& get_domain() const { return m_domain;}
  EBISBox getEBISBox() const {return m_ebisbox;}
  void define(const DomainFluxIFFAB& a_copy){define(a_copy.get_domain(), a_copy.getEBISBox(), a_copy.box(), a_copy.nComp());}

  // Define function
  void define(const ProblemDomain& a_domain, const EBISBox& a_ebisbox, const Box& a_box, const int a_ncomp);

  // Get functions
  BaseIFFAB<Real>&       operator()(const int a_dir, const Side::LoHiSide a_side);
  const BaseIFFAB<Real>& operator()(const int a_dir, const Side::LoHiSide a_side) const;

  // Linearization functions
  int  size(const Box& R, const Interval& comps) const;
  void linearOut(void* buf, const Box& R, const Interval& comps) const;
  void linearIn(void* buf, const Box& R, const Interval& comps);

  /** Copies from a subsection of one box into another.
      Assumes the boxes are both in the same index space, and that
      box R is completely contained in both the src and destination
      boxes.  This function equired by BoxLayoutData */
  void copy(const Box& Rfrom,
            const Interval& Cdest,
            const Box& Rto,
            const DomainFluxIFFAB& src,
            const Interval& Csrc);
  
  
protected:

  ProblemDomain m_domain;
  Box     m_box;
  EBISBox m_ebisbox;
  int     m_nComp;
  bool m_isDefined;

  Tuple<BaseIFFAB<Real>*, CH_SPACEDIM> m_flux_lo;
  Tuple<BaseIFFAB<Real>*, CH_SPACEDIM> m_flux_hi;

  void setDefaultValues();
  void clear();
};

#endif
