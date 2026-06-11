/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ScanShop.H
  @brief  Declare a class geometry generation (ScanShop) that uses the signed distance function for optimizations. 
  @author Robert Marskar
*/

#ifndef CD_SCANSHOP_H
#define CD_SCANSHOP_H

// Chombo includes
#include <GeometryShop.H>
#include <GeometryService.H>

// Our includes
#include <CD_Timer.H>
#include <CD_BoxSorting.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Derived version of GeometryShop that load balances based on the cut cells themselves. 
  @details The internals of this class are basically the same as GeometryShop, with the exception that this class will generate an internal map of the computational
  domain on each level. This map is used to load balance the EB generation, most often leading to orders of magnitude speedup compared to Chombo's approach. 
*/
class ScanShop : public GeometryShop
{
public:
  /**
    @brief Disallowed default constructor
  */
  ScanShop() = delete;

  /**
    @brief Constructor.
    @param[in] a_localGeom     The implicit function for the geometry
    @param[in] a_verbosity     Verbosity
    @param[in] a_dx            Resolution on the finest level
    @param[in] a_probLo        Physical coordinates of the lower-left corner of simulation domain
    @param[in] a_finestDomain  Finest grid level
    @param[in] a_scanLevel     Scan level, i.e. on which level to initiate the load balancing sequence.
    @param[in] a_ebGhost       Number of EB ghost cells
    @param[in] a_thrshdVoF     Threshold for GeometryShop
  */
  ScanShop(const BaseIF&        a_localGeom,
           int                  a_verbosity,
           Real                 a_dx,
           const RealVect&      a_probLo,
           const ProblemDomain& a_finestDomain,
           const ProblemDomain& a_scanLevel,
           int                  a_ebGhost   = 4,
           Real                 a_thrshdVoF = 1.0e-16);

  /**
    @brief Destructor (does nothing)
  */
  ~ScanShop() override;

  /**
    @brief Disallow copy construction.
  */
  ScanShop(const ScanShop&) = delete;

  /**
    @brief Disallow copy assignment.
  */
  ScanShop&
  operator=(const ScanShop&) = delete;

  /**
    @brief Disallow move construction.
  */
  ScanShop(ScanShop&&) = delete;

  /**
    @brief Disallow move assignment.
  */
  ScanShop&
  operator=(ScanShop&&) = delete;

  using GeometryShop::InsideOutside;
  using GeometryShop::isCovered;
  using GeometryShop::isRegular;

  /**
    @brief If doing profiling, set the output name.
    @param[in] a_fileName File name for profiling output
  */
  void
  setProfileFileName(std::string a_fileName);

  /**
    @brief This grid generation method is called by EBISLevel when using distributed data. 
    @param[in] a_domain Problem domain on level
    @param[out] a_grids Load-balanced grids. 
    @param[in]  a_maxGridSize Maximum grid size. 
    @param[in]  a_maxIrregGridSize Maximum irregular grid size. 
    @details This function will create an internal "map" of the cut-cell geometry and return grids that are load-balanced with respect to the cut cells.
  */
  void
  makeGrids(const ProblemDomain& a_domain,
            DisjointBoxLayout&   a_grids,
            const int&           a_maxGridSize,
            const int&           a_maxIrregGridSize) override;

  /**
    @brief This method is called by EBSILevel when using distributed data. 
    @details This allows EBGraph boxes to be designed covered/regular or irregular without actually calling the implicit function. Since this class
    uses recursion froma deeper level, we are immediately able to say that a refined box does not contain cut-cells if the coarse box didn't either. 
    @param[in] a_region Region for inside/outside testing
    @param[in] a_domain Problem domain on level
    @param[in] a_probLo Lower-left corner of physical domain
    @param[in] a_dx     Grid resolution on level
    @param[in] a_dit    Grid index
    @return Return value
  */
  GeometryService::InOut
  InsideOutside(const Box&           a_region,
                const ProblemDomain& a_domain,
                const RealVect&      a_probLo,
                const Real&          a_dx,
                const DataIndex&     a_dit) const override;

  /**
    @brief Override of GeometryShop's fillGraph.
    @details This just wraps a timer around the base version to spot load imbalance.
    @param[in,out] a_regIrregCovered Regular/Covered/Irregular cells
    @param[in,out] a_nodes           Nodes
    @param[in]    a_validRegion     Grid region
    @param[in]    a_ghostRegion     Grid region, including ghost cells
    @param[in]    a_domain          Domain
    @param[in]    a_probLo          Lower-left corner of computational domain
    @param[in]    a_dx              Grid resolution
    @param[in]    a_di              Grid index
  */
  void
  fillGraph(BaseFab<int>&        a_regIrregCovered,
            Vector<IrregNode>&   a_nodes,
            const Box&           a_validRegion,
            const Box&           a_ghostRegion,
            const ProblemDomain& a_domain,
            const RealVect&      a_probLo,
            const Real&          a_dx,
            const DataIndex&     a_di) const override;

protected:
  /**
    @brief Set output file name (if doing profiling)
  */
  std::string m_fileName;

  /**
    @brief For arranging boxes in space when we load balance
  */
  BoxSorting m_boxSorting;

  /**
    @brief Timer for when we use run-time profiling
  */
  mutable Timer m_timer;

  /**
    @brief Bool for run-time profiling of ScanShop
  */
  bool m_profile;

  /**
    @brief Scan level where we first begin to break up boxes. This is relative the EBIS level. 
  */
  int m_scanLevel;

  /**
    @brief Maximum number of ghost cells near the EB
  */
  int m_ebGhost;

  /**
    @brief Lower-left corner of simulation domain. 
  */
  RealVect m_probLo;

  /**
    @brief Implicit function used to generate geometries. 
  */
  const BaseIF* m_baseIF;

  /**
    @brief Check if scan level has been built
  */
  bool m_hasScanLevel;

  /**
    @brief Level resolutions. Note that index 0 is the finest level
  */
  Vector<Real> m_dx;

  /**
    @brief Level domains. Note that index 0 is the finest level
  */
  Vector<ProblemDomain> m_domains;

  /**
    @brief Grids on each level, used for generating the EBIS information. 
  */
  Vector<DisjointBoxLayout> m_grids;

  /**
    @brief Box maps. 
    @details BoxType is a essentially just an integer on each patch which is =0,1,2 depending on whether or not the box is covered/regular/cut. 
  */
  Vector<RefCountedPtr<LayoutData<GeometryService::InOut>>> m_boxMap;

  /**
    @brief Check if level has been scanned. 0 = false, != 0 implies true
  */
  std::vector<bool> m_hasThisLevel;

  /**
    @brief Create the problem domain and resolutions
    @param[in] a_dx           Grid resolution
    @param[in] a_probLo       Lower-left corner of computational domain
    @param[in] a_finestDomain Finest AMR domain
    @param[in] a_scanLevel    Domain where we initiate the load balanced recursion process. 
  */
  void
  makeDomains(Real                 a_dx,
              const RealVect&      a_probLo,
              const ProblemDomain& a_finestDomain,
              const ProblemDomain& a_scanLevel);

  /**
    @brief Refine levels recursively down to the finest level
    @param[in] a_coarserLevel Coarser level
    @param[in] a_maxGridSize Maximum patch size
  */
  void
  buildFinerLevels(int a_coarserLevel, int a_maxGridSize);

  /**
    @brief Build the "coarse level" where we initiate the recursion process.
    @details This level must be defined in order to recurse towards deeper levels. Which is why this is ALWAYS called exact ONCE before buildFinerLevels builds
    the finer levels.
    @param[in] a_level       Level index for the coarse level
    @param[in] a_maxGridSize Maximum patch size
  */
  void
  buildCoarseLevel(int a_level, int a_maxGridSize);

  /**
    @brief Check if every point in input box is regular
    @param[in] a_box    Cell-centered box
    @param[in] a_probLo Lower-left corner of simulation domain
    @param[in] a_dx     Grid resolution
  */
  inline bool
  isRegular(Box a_box, const RealVect& a_probLo, Real a_dx) const;

  /**
    @brief Check if every point in box is covered
    @param[in] a_box    Cell-centered box
    @param[in] a_probLo Lower-left corner of simulation domain
    @param[in] a_dx     Grid resolution
    @return True if covered, false otherwise
  */
  inline bool
  isCovered(Box a_box, const RealVect& a_probLo, Real a_dx) const;

  /**
    @brief Sort boxes lexicographically. 
    @details A strange but true thing that is necessary because DisjointBoxlayout sorts the boxes under the hood
    @return Sorted boxes and types
    @param[in,out] a_boxes Description
    @param[in,out] a_types Description
  */
  inline std::vector<std::pair<Box, int>>
  getSortedBoxesAndTypes(const Vector<Box>& a_boxes, const Vector<int>& a_types) const;

  /**
    @brief Define the "box map" on a specified level
    @param[in,out] coveredBoxes Description
    @param[in,out] regularBoxes Description
    @param[in,out] cutCellBoxes Description
    @param[in] a_level Description
  */
  void
  defineLevel(Vector<Box>& coveredBoxes, Vector<Box>& regularBoxes, Vector<Box>& cutCellBoxes, int a_level);
};

#include <CD_NamespaceFooter.H>

#include <CD_ScanShopImplem.H>

#endif
