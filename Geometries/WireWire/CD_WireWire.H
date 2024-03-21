/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_WireWire.H
  @brief  Declaration of a Wire-wire geometry
  @author Robert Marskar
*/

#ifndef CD_WireWire_H
#define CD_WireWire_H

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief This is a class that defines a wire-wire geometry
*/
class WireWire : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor
  */
  WireWire();

  /*!
    @brief Destructor
  */
  virtual ~WireWire();

protected:
  /*!
    @brief Add another cylinder electrode
    @param[in] a_pp Parameter parser object for the wire. 
  */
  virtual std::pair<std::shared_ptr<EBGeometry::ImplicitFunction<Real>>, Real>
  addElectrode(ParmParse& a_pp);

  /*!
    @brief Add another cylinder dielectric
    @param[in] a_pp Parameter parser object for the wire. 
  */
  virtual std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
  addDielectric(ParmParse& a_pp);
};

#include <CD_NamespaceFooter.H>

#endif
