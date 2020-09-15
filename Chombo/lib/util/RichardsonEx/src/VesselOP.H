#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VESSELOP_H_
#define _VESSELOP_H_

#include "LevelOP.H"

/*
  This class is necessary for using RE on grids that implement AMR. Single
  level grids can simply use the LevelOP class directly. This class contains
  LevelOP' for each level being read into the algorithm.
*/


class VesselOP
{
public:

  /// Constructor
  VesselOP();

  /// Destructor
  ~VesselOP();
  
  /// Definition of the vessel
  void define(string    a_fileName,
              const int a_numLevel,
              const int a_refRatio,
              const int a_gridOrder,
              const int a_comp,
              const int a_oofA,
              const int a_verbosity);

  /// Added a new LevelOP for each level in the vessel
  LevelOP* newLevel(const int a_lev);

  /// Calculate the error between the RE solution and the data
  void calcError(const DisjointBoxLayout&    a_solDBL,
                 const LevelData<FArrayBox>& a_solLD);

  /// Write the HDF5 file of the error
  void writeData(string a_fileName);

  /// Used within writeData for each level
  void writeLevel(const int   a_lev,
                  HDF5Handle& a_handle);

  /// Computes the norms 
  void ComputeNorms();

  /// Returns the norms
  Vector<Real> returnNorms();


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Vector<LevelOP*> m_levelOP;           ///< Vector containing pointers to the 
                                        ///< LevelOPs
  Vector<IntVect> m_amrrefRatio;        ///< Vector of the refinement ratios
  int m_refRatio;                       ///< Ratio between solution grids, 
                                        ///< not levels
  int m_gridOrder;                      ///< Which grid is being read in,
                                        ///< 2 - coarse, 1 - medium, 0 - fine
  int m_numLevel;                       ///< Number of levels in the solution
  int m_oofA;                           ///< Expected order of accuracy
  int m_comp;                           ///< Component number to do the RE
  int m_verbosity;                      ///< Verbosity
  Vector<Real> m_lnorms;                ///< Vector of the norms
                                        ///< 0 - L infinity norm
                                        ///< 1 - L 1 norm
                                        ///< 2 - L 2 norm
                                        ///< on the coarsest level
  string m_fileName;
};

/*--------------------------------------------------------------------*/
// Construct the Vessel object
/** 
 *//*-----------------------------------------------------------------*/

inline
VesselOP::VesselOP()
  :
  m_gridOrder(0),
  m_numLevel(0),
  m_oofA(0),
  m_comp(0),
  m_verbosity(0)
{
  m_fileName = string("start");
  m_levelOP.resize(0);
  m_lnorms.resize(3);
}

/*--------------------------------------------------------------------*/
// Destruct the Vessel object
/** 
 *//*-----------------------------------------------------------------*/

inline
VesselOP::~VesselOP()
{
  for(int lev = 0; lev < m_numLevel; ++lev)
    {
      if(m_levelOP[lev] != NULL)
        {
          delete m_levelOP[lev];
          m_levelOP[lev] = NULL;
        }
    }
}

/*--------------------------------------------------------------------*/
// Define data structures in the class
/** \param[in]  a_fileName     Name of the file being read in
 *  \param[in]  a_numLevel     Number of levels for this file
 *  \param[in]  a_gridOrder    Which grid is being read in.  This
 *                             determines which order the solution is
 *                             2 - coarse
 *                             1 - medium
 *                             0 - fine
 *  \param[in]  a_comp         The component to be used for RE
 *  \param[in]  a_oofA         Expected order of accuracy
 *  \param[in]  a_verbosity    
 *//*-----------------------------------------------------------------*/

inline
void
VesselOP::define(string    a_fileName,
                 const int a_numLevel,
                 const int a_refRatio,
                 const int a_gridOrder,
                 const int a_comp,
                 const int a_oofA,
                 const int a_verbosity)
{
  m_fileName = a_fileName;
  m_numLevel = a_numLevel;
  m_amrrefRatio.resize(m_numLevel+1);
  m_refRatio = a_refRatio;
  m_gridOrder = a_gridOrder;  
  m_comp = a_comp;
  m_oofA = a_oofA;

  m_verbosity = a_verbosity;

  m_levelOP.resize(m_numLevel+1);

  if(m_verbosity > 3)
    {
      pout() << "Defining Vessel " << m_gridOrder << endl;
    }
  /// Define a LevelOP for each level in the solution
  for(int lev = 0; lev <= m_numLevel; ++lev)
    {
      m_levelOP[lev] = newLevel(lev);
      m_amrrefRatio[lev] = m_levelOP[lev]->getAMRRefRatio();

      /// This is the ratio between the current level and the fine
      /// single level solution
      int refRate = (pow(m_amrrefRatio[lev][0],(m_numLevel-lev))*
                     pow(m_refRatio,m_gridOrder));
      m_levelOP[lev]->setSolRefRate(refRate);
    }
}

/*--------------------------------------------------------------------*/
// Allocate the memory for each level in the solution data
/** \param[in]  lev  Current level 
 *//*-----------------------------------------------------------------*/

inline
LevelOP* VesselOP::newLevel(const int a_lev)
{
  LevelOP* lvlPtr = new LevelOP(m_comp,
                                m_oofA,
                                m_verbosity);
  lvlPtr->define(m_fileName,
                 a_lev,
                 true);
  if(m_verbosity > 3)
    {
      pout() << "Vessel level defined: " << a_lev << endl;
    }
  return (static_cast <LevelOP*> (lvlPtr));
}

/*--------------------------------------------------------------------*/
// Solves for the error
/** \param[in]  a_pfDBL Disjoint box layout of the fine grid
 *  \param[in]  a_pfLD  Level data of the fine grid RE solution
 *//*-----------------------------------------------------------------*/

inline void
VesselOP::calcError(const DisjointBoxLayout&    a_solDBL,
                    const LevelData<FArrayBox>& a_solLD)
{
  if(m_verbosity > 3)
    {
      pout() << "Calculating error for " << m_gridOrder << endl;
    }
  for(int lev = 0; lev <= m_numLevel; ++lev)
    {
      m_levelOP[lev]->defineCoarseAverage(a_solDBL,
                                          a_solLD);
      m_levelOP[lev]->differenceSolve();
    }
}

/*--------------------------------------------------------------------*/
// Write the data to an HDF5 file
/** \param[in]  a_fileName Name of the HDF5 file
 *//*-----------------------------------------------------------------*/

inline void
VesselOP::writeData(string a_fileName)
{
  if(m_verbosity > 3)
    {
      pout() << "Writing error HDF5 file for " << m_gridOrder << endl;
    }
  const int numComp = m_levelOP[0]->diffLD().nComp();
  HDF5Handle handle(a_fileName,HDF5Handle::CREATE);

  HDF5HeaderData header1;
  header1.m_int ["ref_ratio"] = m_amrrefRatio[0][0];
  header1.m_int ["max_level"] = m_numLevel;
  header1.m_int ["num_levels"] = m_numLevel+1;
  header1.m_int ["iteration"] = 1;
  header1.m_real ["time"] = 1.;
  header1.writeToFile(handle);
  HDF5HeaderData header2;
  header2.m_int ["num_components"] = numComp;
  char compStr[30];
  for(int comp = 0; comp != numComp; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header2.m_string[compStr] = "Error" + std::to_string(comp);
    }
  header2.writeToFile(handle);
  for(int lev = 0; lev <= m_numLevel; ++lev)
    {
      writeLevel(lev,handle);
    }
  handle.close();
}

/*--------------------------------------------------------------------*/
// Write the level to the HDF5 file
/** \param[in]  a_lev       Current level being read in 
 *  \param[in]  a_handle    Handle for the file
 *//*-----------------------------------------------------------------*/

inline void
VesselOP::writeLevel(const int   a_lev,
                     HDF5Handle& a_handle)
{
  LevelData<FArrayBox> tempwritedata; 
  
  tempwritedata.define(m_levelOP[a_lev]->diffLD());

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",a_lev);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  Real dx = m_levelOP[a_lev]->gridDX(0);

  Box dbox = m_levelOP[a_lev]->boxdomLD();

  header.m_int ["ref_ratio"]   = m_amrrefRatio[a_lev][0];
  header.m_real["dx"]          = dx;
  header.m_real["dt"]          = 1;
  header.m_real["time"]        = 1;
  header.m_box ["prob_domain"] = dbox;

  // Write the header for this level
  header.writeToFile(a_handle);


  // Write the data for this level
  const int numComp = m_levelOP[a_lev]->diffLD().nComp();
  int numPlotVar = numComp;
  LevelData<FArrayBox> plotData(tempwritedata.getBoxes(), numPlotVar);

  // first copy data to plot data holder
  m_levelOP[a_lev]->diffLD().copyTo(tempwritedata.interval(), 
                                  plotData, tempwritedata.interval());

  write(a_handle,tempwritedata.boxLayout());
  write(a_handle,plotData,"data");
}

/*--------------------------------------------------------------------*/
// Solve for the norms
/** 
 *//*-----------------------------------------------------------------*/

inline void
VesselOP::ComputeNorms()
{
  Interval icomps(m_comp-1,m_comp-1);

  Vector<LevelData<FArrayBox>* > phi(m_numLevel+1);
  
  RealVect crsedx = m_levelOP[0]->gridDX();

  for(int lev = 0; lev <= m_numLevel; ++lev)
    {
      phi[lev] = &(m_levelOP[lev]->diffLD());
    }
  for(int no = 0; no != 3; ++no)
    {
      m_lnorms[no] = computeNorm(phi,
                                 m_amrrefRatio,
                                 crsedx,
                                 icomps,
                                 no);
    }
}

/*--------------------------------------------------------------------*/
// Output the norms
/** 
 *  \param[out] m_lnorms     A vector of the norms
 *//*-----------------------------------------------------------------*/

inline Vector<Real>
VesselOP::returnNorms()
{
  return m_lnorms;
}

#endif  /* ! defined _VESSELOP_H_ */
