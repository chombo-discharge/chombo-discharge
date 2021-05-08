/*!
  @brief  ScanShop.H
  @brief  Version of GeometryShop that also may use a bounding box for location EB volumes
  @author Robert Marskar
  @date   2019
*/

#ifndef _SCANSHOP_
#define _SCANSHOP_

#include <GeometryShop.H>
#include <GeometryService.H>

#include "BoxType.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Derived version of GeometryShop that load balances based on the cut cells themselves. The internals of this class
  is basically the same as GeometryShop, with the exception that this class will generate an internal map of the computational
  domain on each level. This map is used to load balance the EB graph generation. 
*/
class ScanShop : public GeometryShop {
public:

  static bool s_irregularBalance; // Balance with irregular grid 
  static bool s_recursive;         // Use recursive box division in order to lighten the load on Copier
  static int  s_grow;              // Used to grow boxes when testing for real/covered/other. Ideally, should equal EBIS ghosts 

  /*!
    @brief Constructor. 
    @param[in] a_localGeometry The implicit function for the geometry
    @param[in] a_verbosity     Verbosity
    @param[in] a_dx            Resolution on the finest level
    @param[in] a_origin        Domain origin
    @param[in] a_finestdomain  Finest grid level
    @param[in] a_scanLevel     Scan level
    @param[in] a_thsdhVoF      Threshold for GeometryShop
  */
  ScanShop(const BaseIF&       a_localGeom,
	   const int           a_verbosity,
	   const Real          a_dx,
	   const RealVect      a_origin,
	   const ProblemDomain a_finestDomain,
	   const ProblemDomain a_scanLevel,
	   const Real          a_thrshdVoF = 1.0e-16);
  
  ~ScanShop();

  /*!
    @brief This grid generation method is called by EBISLevel when using distributed data. This function will scan
    an internal "map" of the cut-cell geometry and return grids that are load-balanced with respect to the cut cells.
  */
  void makeGrids(const ProblemDomain&      a_domain,
		 DisjointBoxLayout&        a_grids,
		 const int&                a_maxGridSize,
		 const int&                a_maxIrregGridSize);
  
  /*!
    @brief This method is called by EBSILevel when using distributed data. It allows EBGraph boxes to be designed covered/regular
    or irregular without actually calling the implicti function. 
  */
  GeometryService::InOut InsideOutside(const Box&           a_region,
				       const ProblemDomain& a_domain,
				       const RealVect&      a_origin,
				       const Real&          a_dx,
				       const DataIndex&     a_dit) const;


protected:


  /*!
    @brief Scan level where we first begin to break up boxes
  */
  int m_scanLevel;

  /*!
    @brief Origin
  */
  RealVect m_origin;

  /*!
    @brief BaseIF function
  */
  const BaseIF* m_baseif;

  /*!
    @brief Check if scan level has been built
  */
  bool m_hasScanLevel;

  /*!
    @brief Level resolutions. Index 0 is the finest level
  */
  Vector<Real> m_dx;         // Index 0 is the finest EBIS level

  /*!
    @brief Level domains. Index 0 is the finest level
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief Box maps
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Box maps
  */
  Vector<RefCountedPtr<LevelData<BoxType> > > m_boxMaps;

  /*!
    @brief Check if level has been scanned. 0 = false, != 0 implies true
  */
  Vector<int> m_hasThisLevel;

  /*!
    @brief Create the problem domain and resolutions
  */
  void makeDomains(const Real a_dx, const RealVect a_origin, const ProblemDomain a_finestDomain, const ProblemDomain a_scanLevel);

  /*!
    @brief Refine levels recursively down to the finest level
  */
  void buildFinerLevels(const int a_coarserLevel, const int a_maxGridSize);

  /*!
    @brief Refine levels recursively down to the finest level
  */
  void buildCoarseLevel(const int a_finerLevel, const int a_maxGridSize);

  /*!
    @brief Prints number of boxes of each type on this level
  */
  void printNumBoxesLevel(const int a_level);

  /*!
    @brief Gather boxes in parallel.
    @details On input, a_boxes contains the local boxes. On output, boxes have been gather on every rank
  */
  void gatherBoxesParallel(Vector<Box>& a_boxes);

  /*!
    @brief Check if every point in box is regular
  */
  bool isRegular(const Box a_box, const RealVect m_origin, const Real a_dx);

  /*!
    @brief Check if every point in box is covered
  */
  bool isCovered(const Box a_box, const RealVect m_origin, const Real a_dx);
};
#include "CD_NamespaceFooter.H"

#endif
