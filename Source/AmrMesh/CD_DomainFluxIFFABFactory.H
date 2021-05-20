/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DomainFluxIFFABFactory.H
  @brief  Declaration of factory class for DomainFluxIFFAB
  @author Robert Marskar
*/

#ifndef CD_DomainFluxIFFABFactory_ 
#define  CD_DomainFluxIFFABFactory_ 

// Chombo includes
#include <Vector.H>
#include <ProblemDomain.H>
#include <EBISLayout.H>
#include <DataIndex.H>

// Our includes
#include <CD_DomainFluxIFFAB.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Factory class for DomainFluxIFFAB
*/
class DomainFluxIFFABFactory : public DataFactory<DomainFluxIFFAB> {
public:

  DomainFluxIFFABFactory(const EBISLayout& a_ebisl, const ProblemDomain& a_domain);
  ~DomainFluxIFFABFactory();

  virtual DomainFluxIFFAB* create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const;

private:
  EBISLayout    m_ebisl;
  ProblemDomain m_domain;
};
  
#include <CD_NamespaceFooter.H>

#endif
