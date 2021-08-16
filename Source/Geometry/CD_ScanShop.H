/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief  CD_ScanShop.H
  @brief  Declare a class geometry generation (ScanShop) that uses the signed distance function for optimizations. 
  @author Robert Marskar
*/

#ifndef CD_ScanShop_H
#define CD_ScanShop_H

// Chombo includes
#include <GeometryShop.H>
#include <GeometryService.H>

// Our includes
#include <CD_BoxType.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Derived version of GeometryShop that load balances based on the cut cells themselves. 
  @details The internals of this class are basically the same as GeometryShop, with the exception that this class will generate an internal map of the computational
  domain on each level. This map is used to load balance the EB generation, most often leading to orders of magnitude speedup compared to Chombo's approach. 
*/
class ScanShop : public GeometryShop {
public:

  /*!
    @brief Disallowed default constructor
  */
  ScanShop() = delete;
  
  /*!
    @brief Constructor. 
    @param[in] a_geometry      The implicit function for the geometry
    @param[in] a_verbosity     Verbosity
    @param[in] a_finestDx      Resolution on the finest level
    @param[in] a_probLo        Physical coordinates of the lower-left corner of simulation domain
    @param[in] a_finestDomain  Finest grid level
    @param[in] a_scanLevel     Scan level, i.e. on which level to initiate the load balancing sequence. 
    @param[in] a_thsdhVoF      Threshold for GeometryShop
  */
  ScanShop(const BaseIF&       a_geometry,
	   const int           a_verbosity,
	   const Real          a_finestDx,
	   const RealVect      a_probLo,
	   const ProblemDomain a_finestDomain,
	   const ProblemDomain a_scanLevel,
	   const int           a_maxGhostEB = 4,
	   const Real          a_thrshdVoF = 1.0e-16);

  /*!
    @brief Destructor (does nothing)
  */
  ~ScanShop();

  /*!
    @brief This grid generation method is called by EBISLevel when using distributed data. 
    @param[in] a_domain Problem domain on level
    @param[out] a_grids Load-balanced grids. 
    @param[in]  a_maxGridSize Maximum grid size. 
    @param[in]  a_maxIrregGridSize Maximum irregular grid size. 
    @details This function will create an internal "map" of the cut-cell geometry and return grids that are load-balanced with respect to the cut cells.
  */
  void makeGrids(const ProblemDomain&      a_domain,
		 DisjointBoxLayout&        a_grids,
		 const int&                a_maxGridSize,
		 const int&                a_maxIrregGridSize) override;
  
  /*!
    @brief This method is called by EBSILevel when using distributed data. 
    @details This allows EBGraph boxes to be designed covered/regular or irregular without actually calling the implicit function. Since this class
    uses recursion froma deeper level, we are immediately able to say that a refined box does not contain cut-cells if the coarse box didn't either. 
    @param[in] a_region Region for inside/outside testing
    @param[in] a_domain Problem domain on level
    @param[in] a_probLo Lower-left corner of physical domain
    @param[in] a_dx     Grid resolution on level
    @param[in] a_dit    Grid index
  */
  GeometryService::InOut InsideOutside(const Box&           a_region,
				       const ProblemDomain& a_domain,
				       const RealVect&      a_probLo,
				       const Real&          a_dx,
				       const DataIndex&     a_dit) const override;


protected:

  /*!
    @brief Scan level where we first begin to break up boxes. This is relative the EBIS level. 
  */
  int m_scanLevel;

  /*!
    @brief Maximum number of ghost cells near the EB
  */
  int m_maxGhostEB;

  /*!
    @brief Lower-left corner of simulation domain. 
  */
  RealVect m_probLo;

  /*!
    @brief Implicit function used to generate geometrys. 
  */
  const BaseIF* m_baseIF;

  /*!
    @brief Check if scan level has been built
  */
  bool m_hasScanLevel;

  /*!
    @brief Level resolutions. Note that index 0 is the finest level
  */
  Vector<Real> m_dx;

  /*!
    @brief Level domains. Note that index 0 is the finest level
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Grids on each level, used for generating the EBIS information. 
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Box maps. 
    @details BoxType is a essentially just an integer on each patch which is =0,1,2 depending on whether or not the box is covered/regular/cut. 
  */
  Vector<RefCountedPtr<LevelData<BoxType> > > m_boxMaps;

  /*!
    @brief Check if level has been scanned. 0 = false, != 0 implies true
  */
  Vector<int> m_hasThisLevel;

  /*!
    @brief Create the problem domain and resolutions
    @param[in] a_dx           Grid resolution
    @param[in] a_probLo       Lower-left corner of computational domain
    @param[in] a_finestDomain Finest AMR domain
    @param[in] a_scanLevel    Domain where we initiate the load balanced recusion process. 
  */
  void makeDomains(const Real a_dx, const RealVect a_probLo, const ProblemDomain a_finestDomain, const ProblemDomain a_scanLevel);

  /*!
    @brief Refine levels recursively down to the finest level
    @param[in] a_coarserLevel Coarser level
    @param[in] a_maxGridSize Maximum patch size
  */
  void buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize);

  /*!
    @brief Build the "coarse level" where we initiate the recursion process.
    @details This level must be defined in order to recurse towards deeper levels. Which is why this is ALWAYS called exact ONCE before buildFinerLevels builds
    the finer levels.
    @param[in] a_coarserLevel Coarser level
    @param[in] a_maxGridSize Maximum patch size
  */
  void buildCoarseLevel(const int a_finerLevel, const int a_maxGridSize);

  /*!
    @brief Prints number of boxes of each type on this level
    @param[in] a_level Grid level
    @details This will print (to std::cout) the number of regular/covered/cutcell boxes on each grid level
  */
  void printNumBoxesLevel(const int a_level) const;

  /*!
    @brief Gather boxes in parallel.
    @param[inout] a_boxes Boxes to gather. 
    @details On input, a_boxes contains the local boxes. On output, boxes have been gather on every rank
  */
  void gatherBoxesParallel(Vector<Box>& a_boxes) const;

  /*!
    @brief Check if every point in input box is regular
    @param[in] a_box    Cell-centered box
    @param[in] a_probLo Lower-left corner of simulation domain
    @param[in] a_dx     Grid resolution
  */
  bool isRegular(const Box a_box, const RealVect a_probLo, const Real a_dx) const;

  /*!
    @brief Check if every point in box is covered
    @param[in] a_box    Cell-centered box
    @param[in] a_probLo Lower-left corner of simulation domain
    @param[in] a_dx     Grid resolution
  */
  bool isCovered(const Box a_box, const RealVect a_probLo, const Real a_dx) const;
};

#include <CD_NamespaceFooter.H>

#endif
