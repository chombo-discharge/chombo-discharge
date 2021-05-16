/*!
  @file cdr_gdnv.H
  @brief Godunov solver for CDR equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_GDNV_
#define _CDR_GDNV_

#include "cdr_tga.H"

#include <EBAdvectLevelIntegrator.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Godunov implementation for advection. 
*/
class cdr_gdnv : public cdr_tga {
public:

  /*!
    @brief Constructor
  */
  cdr_gdnv();

  /*!
    @brief Destructor
  */
  virtual ~cdr_gdnv();

  /*!
    @brief Parse options
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions() override;

  /*!
    @brief Query number of ghost cells
  */
  virtual int queryGhost() const override;

  /*!
    @brief Allocate internal stuff
  */
  virtual void allocateInternals() override;

  /*!
    @brief Average velocities to faces
  */
  virtual void averageVelocityToFaces() override;

  /*!
    @brief Call parent function: Then fix no influx
  */
  virtual void averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity) override;

protected:

  /*!
    @brief Godunov face extrapolation method for advection
  */
  virtual void advect_to_faces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override;

  /*!
    @brief Level advect integrator
  */
  Vector<RefCountedPtr<EBAdvectLevelIntegrator> > m_level_advect;

  /*!
    @brief Turn on/off slope limiting
  */
  bool m_slopelim;

  /*!
    @brief Parses slope limiter options
  */
  virtual void parse_slopelim();     // Parses slope limiter settings
};
#include "CD_NamespaceFooter.H"

#endif
