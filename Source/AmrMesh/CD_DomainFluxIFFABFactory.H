/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DomainFluxIFFABFactory.H
  @brief  Declaration of factory class for DomainFluxIFFAB
  @author Robert Marskar
*/

#ifndef CD_DomainFluxIFFABFactory_
#define CD_DomainFluxIFFABFactory_

// Chombo includes
#include <Vector.H>
#include <ProblemDomain.H>
#include <EBISLayout.H>
#include <DataIndex.H>

// Our includes
#include <CD_DomainFluxIFFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making DomainFluxIFFAB
*/
class DomainFluxIFFABFactory : public DataFactory<DomainFluxIFFAB>
{
public:
  /*!
    @brief Full constructor
    @param[in] a_ebisl  EBISLayout
    @param[in] a_domain Computational domain
  */
  DomainFluxIFFABFactory(const EBISLayout& a_ebisl, const ProblemDomain& a_domain);

  /*!
    @brief Destructor (does nothing)
  */
  ~DomainFluxIFFABFactory();

  /*!
    @brief Factory method
    @param[in] a_box   Input region
    @param[in] a_nComp Number of components
    @param[in] a_dit   Grid index
    @returns Returns a new DomainIFFluxFAB defined over the input region
  */
  virtual DomainFluxIFFAB*
  create(const Box& a_box, int a_nComp, const DataIndex& a_dit) const;

protected:
  /*!
    @brief EBIS layout
  */
  EBISLayout m_ebisl;

  /*!
    @brief Domain
  */
  ProblemDomain m_domain;
};

#include <CD_NamespaceFooter.H>

#endif
