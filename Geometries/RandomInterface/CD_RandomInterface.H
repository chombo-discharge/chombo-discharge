/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RandomInterface.H
  @brief  Random interface geometry
  @author Robert Marskar
*/

#ifndef CD_RandomInterface_H
#define CD_RandomInterface_H

#include <CD_EBGeometryIF.H>
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that represents two random interface smooshed together
*/
class RandomInterface : public ComputationalGeometry
{
public:
  /*!
    @brief Constructor.
  */
  RandomInterface() noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~RandomInterface() noexcept;

protected:
  /*!
    @brief Alias for 3D vector type
  */
  using Vec3 = EBGeometry::Vec3T<Real>;

  /*!
    @brief Clamped noisy plane signed distance function. 
  */
  template <class URNG>
  class ClampedNoisePlane : public EBGeometry::SignedDistanceFunction<Real>
  {
  public:
    /*!
      @brief Disallowed weak construction
    */
    ClampedNoisePlane() = delete;

    /*!
      @brief Full constructor
    */
    ClampedNoisePlane(const Vec3 a_point,
                      const Vec3 a_normal,
                      const Vec3 a_clampLo,
                      const Vec3 a_clampHi,
                      const Vec3 a_clampDx,
                      const Vec3 a_noiseFrequency,
                      const Real a_noiseAmplitude,
                      const Real a_noisePersistence,
                      const int  a_noiseOctaves,
                      const URNG a_rng) noexcept;

    /*!
      @brief Signed distance function
      @param[in] a_point Input point
    */
    virtual Real
    signedDistance(const Vec3& a_point) const noexcept override;

  protected:
    /*!
      @brief Point on plane
    */
    Vec3 m_point;

    /*!
      @brief Plane normal vector
    */
    Vec3 m_normal;

    /*!
      @brief Low-end clamping
    */
    Vec3 m_clampLo;

    /*!
      @brief High-end clamping
    */
    Vec3 m_clampHi;

    /*!
      @brief Clamping length
    */
    Vec3 m_clampDx;

    /*!
      @brief Random number generator
    */
    URNG m_rng;

    /*!
      @brief Base implicit function
    */
    std::shared_ptr<EBGeometry::PlaneSDF<Real>> m_plane;

    /*!
      @brief Noise function
    */
    std::shared_ptr<EBGeometry::PerlinSDF<Real>> m_perlin;
  };
};

#include <CD_NamespaceFooter.H>

#endif
