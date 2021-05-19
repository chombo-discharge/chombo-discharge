/*!
  @file   rod_dielectric.H
  @brief  Definition of a electrode rod and dielectric dielectric geometry
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ROD_DIELECTRIC_
#define _ROD_DIELECTRIC_

#include <CD_ComputationalGeometry.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Class that defines an electrode rod and a dielectric (some shapes supported)
*/
class rod_dielectric : public ComputationalGeometry {
public:

  rod_dielectric();
  ~rod_dielectric();

protected:

  void define_electrode();
  void define_insulator();

  RefCountedPtr<BaseIF> get_box();
  RefCountedPtr<BaseIF> get_plane();
  RefCountedPtr<BaseIF> get_sphere();
  RefCountedPtr<BaseIF> get_perlin_box();
};
#include "CD_NamespaceFooter.H"

#endif
