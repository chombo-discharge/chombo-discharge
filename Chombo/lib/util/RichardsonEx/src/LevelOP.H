#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELOP_H_
#define _LEVELOP_H_

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "Vector.H"

#include "DisjointBoxLayout.H"
#include "BoxLayout.H"
#include "LevelData.H"
#include "PatchOpIF_F.H"
#include "DebugOut.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "PatchOpIF_F.H"

/*
This class is contains all the data necessary for each level,
as well as the operations necessary to complete the RE procedure.
*/

class LevelOP
{
public:

  /// Constructor
  LevelOP(const int a_compT,
          const int a_oofA,
          const int a_verbosity);

  /// Define data structures in the class
  int define(string    a_fileName,
             const int a_level,
             bool      a_wAMR);

  /// Average the refined solution data or the refined RE solution
  /// to the current level
  void defineCoarseAverage(const DisjointBoxLayout&    a_pfDBL,
                           const LevelData<FArrayBox>& a_pfLD);

  /// Write the error data for this level to an HDF5 file
  void writePlotLevel(HDF5Handle& a_handle);

  /// Write the error data to an HDF5 file
  void writeData(string a_fileName);

  /// Find the RE solution between the grids
  void solRE();

  /// Find difference between current solution and solution member data
  void differenceSolutionSolve();

  /// Find difference between the current solution and the interpolated 
  /// member data
  void differenceSolve();

  /// Return the dx in only one direction
  Real gridDX(const int dir);

  /// Return the dx vector
  RealVect& gridDX();

  /// Return the DBL for the current level
  DisjointBoxLayout& DBL();

  /// Return the solution for the current level
  LevelData<FArrayBox>& LD();

  /// Return the RE solution for this data
  LevelData<FArrayBox>& solLD();

  /// Return the error for this data
  LevelData<FArrayBox>& diffLD();

  /// Return the error for this data
  Box& boxdomLD();

  /// Return the amr refinement ratios
  IntVect& getAMRRefRatio();

  /// Set the refinement between the current solution and the 
  /// single grid fine solution
  void setSolRefRate(const int a_solRefRate);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  LevelData<FArrayBox> m_ld;            ///< Level data
  LevelData<FArrayBox> m_interpLD;      ///< Level data averaged to
                                        ///< coarser grid: for single level
                                        ///< this is the higher grid LD
                                        ///< for AMR, this is the solution on
                                        ///< fine mesh
  LevelData<FArrayBox> m_reSol;         ///< LD for higher order solution
                                        ///< that comes from RE (only used
                                        ///< for single grids
  LevelData<FArrayBox> m_diffLD;        ///< Difference between U and u
  ProblemDomain m_probDomain;           ///< Problem domain
  Box m_domainBox;                      ///< Box of problem domain
  RealVect m_dx;                        ///< Mesh spacing
  int m_oofA;                           ///< Order of accuracy of solution
  int m_level;                          ///< Level if LevelData uses AMR
  int m_compT;                          ///< Component of interest
  DisjointBoxLayout m_dbl;              ///< Disjoint box layout
  IntVect m_numCell;                    ///< Number of cells in each direction
  Vector<Box> m_boxes;                  ///< Vector of boxes
  int m_verbosity;                      ///< Verbosity
  int m_solRef;                         ///< Ratio from the fine mesh
  IntVect m_amrrefRatio;                ///< Refinement ratios if AMR is used
  IntVect m_ghostCells;
  bool m_wAMR;                          ///< Does the solution use AMR
  bool m_ldDef;                         ///< True if m_ld is defined
  bool m_interpDef;                     ///< True if m_interpLD is defined
  bool m_solDef;                        ///< True if m_reSol is defined
  bool m_diffDef;                       ///< True if m_diffLD is defined
  string m_fileName;                    ///< Input file name to define m_ld
};


/*******************************************************************************
 *
 * Class LevelOP: inline member definitions
 *
 ******************************************************************************/

inline
LevelOP::LevelOP(const int a_compT,
                 const int a_oofA,
                 const int a_verbosity)
:
  m_dx(RealVect::Zero),
  m_level(0),
  m_solRef(1),
  m_wAMR(false),
  m_ldDef(false),
  m_interpDef(false),
  m_solDef(false),
  m_diffDef(false)
{
  m_ghostCells = IntVect::Zero;
  m_verbosity = a_verbosity;
  m_compT = a_compT;
  m_oofA = a_oofA;
  if(m_verbosity >= 3)
    {
      pout() << "LevelOP construction" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Define data structures in the class
/** \param[in]  a_fileName     Name of the file being read in
 *  \param[in]  a_level        AMR level being read in
 *  \param[in]  a_wAMR         Bool: if false, single level grids
 *//*-----------------------------------------------------------------*/

inline int
LevelOP::define(string    a_fileName,
                const int a_level,
                bool      a_wAMR)
{
  int ier = 0;
  m_fileName = a_fileName;
  m_wAMR = a_wAMR;
  m_ghostCells = IntVect::Zero;
  // Read in the reference file
  HDF5Handle handle(a_fileName,HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle);
  m_level = a_level;
  if(m_verbosity >= 3)
    {
      pout() << "Defining the DBL for " << a_fileName << endl;
      pout() << "At level " << a_level << endl;
    }
  Real a_dt;
  Real a_time;
  int readin = readLevel<FArrayBox>(handle,
                                    m_level,
                                    m_ld,
                                    m_dx,
                                    a_dt,
                                    a_time,
                                    m_domainBox,
                                    m_amrrefRatio);
  if(readin != 0)
    {
      pout() << "Error reading data from file" << endl;
      ier = 1;
    }

  m_dbl = m_ld.getBoxes();
  handle.close();
  if(m_verbosity >= 3)
    {
      pout() << a_fileName << " has been read" << endl;
      pout() << "Problem domain: " << m_domainBox << endl;
    }
  m_ldDef = true;
  return ier;
}

/*--------------------------------------------------------------------*/
// Solves the interpLD. For single grids, this averages the data from
// the finer solution onto the current DBL. For AMR grids, the RE solution
// found between the postfine and the fine single grid solutions are 
// averaged down to the current DBL.
/** \param[in]  a_pfDBL Disjoint box layout of the finer grid
 *  \param[in]  a_pfLD  Level data of the finer grid
 *//*-----------------------------------------------------------------*/

inline void
LevelOP::defineCoarseAverage(const DisjointBoxLayout&    a_pfDBL,
                             const LevelData<FArrayBox>& a_pfLD)
{
  CH_assert(m_ldDef);
  int numComp = a_pfLD.nComp();
  
  m_interpLD.define(m_dbl,numComp,m_ghostCells);
  CoarseAverage cam(a_pfDBL,
                    m_dbl,
                    numComp,
                    m_solRef,
                    m_ghostCells);
  cam.averageToCoarse(m_interpLD,a_pfLD);
  m_interpDef = true;
  if(m_verbosity >= 3)
    {
      pout() << m_fileName << ": Coarse average defined" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Find the RE solution between the interpolated data and the current data
/**
 *//*-----------------------------------------------------------------*/

inline void
LevelOP::solRE()
{
  CH_assert(m_ldDef);
  CH_assert(m_interpDef);
  CH_assert(!m_wAMR);
  const int numComp = m_ld.nComp();
  m_reSol.define(m_dbl,numComp,m_ghostCells);
  for(DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& coarseFab = m_ld[dit()];
      const FArrayBox& interpFab = m_interpLD[dit()];
      FArrayBox& solFab = m_reSol[dit()];

      Box fabBox = m_ld.disjointBoxLayout()[dit];
      fabBox.grow(m_ghostCells);
      for(int comp = 0; comp != numComp; ++comp)
        {
          FORT_RESOLVE(
            CHF_BOX(fabBox),
            CHF_CONST_INT(m_oofA),
            CHF_CONST_INT(m_solRef),
            CHF_CONST_FRA1(coarseFab,comp),
            CHF_CONST_FRA1(interpFab,comp),
            CHF_FRA1(solFab,comp));
        }
    }
  m_solDef = true;

  if(m_verbosity >= 3)
    {
      pout() << m_fileName << ": RE solution defined" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Find the difference between the RE solution on the single fine grid
// and the data on the current level.
/**
 *//*-----------------------------------------------------------------*/

inline void
LevelOP::differenceSolutionSolve()
{
  CH_assert(m_solDef);
  const int numComp = m_ld.nComp();
  m_diffLD.define(m_dbl,numComp,m_ghostCells);
      
  for(DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& solFab = m_reSol[dit()];
      const FArrayBox& curFab = m_ld[dit()];
      FArrayBox& diffFab = m_diffLD[dit()];

      Box fabBox = m_dbl[dit];
      fabBox.grow(m_ghostCells);
      for(int comp = 0; comp != numComp; ++comp)
        {
          FORT_DIFFSOLVE(
            CHF_BOX(fabBox),
            CHF_CONST_FRA1(curFab,comp),
            CHF_CONST_FRA1(solFab,comp),
            CHF_FRA1(diffFab,comp));
        }
    }
  m_diffDef = true;
  if(m_verbosity >= 3)
    {
      pout() << m_fileName << ": Difference defined" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Find the difference between the interpolated data on the single fine grid
// and the data on the current level.
/**
 *//*-----------------------------------------------------------------*/
inline
void
LevelOP::differenceSolve()
{
  CH_assert(m_interpDef);
  const int numComp = m_ld.nComp();
  m_diffLD.define(m_dbl,numComp,m_ghostCells);
      
  for(DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& intFab = m_interpLD[dit()];
      const FArrayBox& curFab = m_ld[dit()];
      FArrayBox& diffFab = m_diffLD[dit()];

      Box fabBox = m_dbl[dit];
      fabBox.grow(m_ghostCells);
      for(int comp = 0; comp != numComp; ++comp)
        {
          FORT_DIFFSOLVE(
            CHF_BOX(fabBox),
            CHF_CONST_FRA1(curFab,comp),
            CHF_CONST_FRA1(intFab,comp),
            CHF_FRA1(diffFab,comp));
        }
    }
  m_diffDef = true;
  if(m_verbosity >= 3)
    {
      pout() << m_fileName << ": Difference defined" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Write the level to the HDF5 file
/** \param[in]  a_handle 
 *//*-----------------------------------------------------------------*/
inline
void
LevelOP::writePlotLevel(HDF5Handle& a_handle)
{
  CH_assert(m_diffDef);

  LevelData<FArrayBox> tempwritedata;
  
  tempwritedata.define(m_diffLD);

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_amrrefRatio[0];
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = 1.;
  header.m_real["time"]        = 0.;
  header.m_box ["prob_domain"] = m_domainBox;

  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  int numPlotVar = m_ld.nComp();
  LevelData<FArrayBox> plotData(tempwritedata.getBoxes(), numPlotVar);

  // first copy data to plot data holder
  m_diffLD.copyTo(tempwritedata.interval(), plotData, 
                  tempwritedata.interval());

  write(a_handle,tempwritedata.boxLayout());
  write(a_handle,plotData,"data");
  if(m_verbosity >= 3)
    {
      pout() << m_fileName << ": File written" << endl;
    }
}

/*--------------------------------------------------------------------*/
// Write the data to an HDF5 file
/** \param[in]  a_fileName Name of the HDF5 file
 *//*-----------------------------------------------------------------*/
inline
void
LevelOP::writeData(string a_fileName)
{
  CH_assert(!m_wAMR);
  const int numComp = m_ld.nComp();
  HDF5Handle handle(a_fileName,HDF5Handle::CREATE);

  HDF5HeaderData header1;
  header1.m_int ["ref_ratio"] = m_amrrefRatio[0];
  header1.m_int ["max_level"] = 0;
  header1.m_int ["num_levels"] = 1;
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
  writePlotLevel(handle);
  handle.close();
}


/*--------------------------------------------------------------------*/
// Output the disjoint box layout associated with this grid and level
/** 
 *//*-----------------------------------------------------------------*/
inline
DisjointBoxLayout&
LevelOP::DBL()
{
  CH_assert(m_ldDef);
  return m_dbl;
}

/*--------------------------------------------------------------------*/
// Output the level data associated with this grid and level
/** 
 *//*-----------------------------------------------------------------*/
inline
LevelData<FArrayBox>&
LevelOP::LD()
{
  CH_assert(m_ldDef);
  return m_ld;
}

/*--------------------------------------------------------------------*/
// Output the RE solution associated with this grid and level (only for 
// single grids)
/** 
 *//*-----------------------------------------------------------------*/
inline
LevelData<FArrayBox>&
LevelOP::solLD()
{
  CH_assert(m_solDef);
  return m_reSol;
}

/*--------------------------------------------------------------------*/
// Output the error associated with this grid and level
/** 
 *//*-----------------------------------------------------------------*/
inline
LevelData<FArrayBox>&
LevelOP::diffLD()
{
  CH_assert(m_diffDef);
  return m_diffLD;
}

/*--------------------------------------------------------------------*/
// Output the grid spacing in a single direction
/** 
 *//*-----------------------------------------------------------------*/
inline
Real
LevelOP::gridDX(const int dir)
{
  CH_assert(m_ldDef);
  return m_dx[dir];
}


/*--------------------------------------------------------------------*/
// Output the grid spacing vector
/** 
 *//*-----------------------------------------------------------------*/
inline
RealVect&
LevelOP::gridDX()
{
  CH_assert(m_ldDef);
  return m_dx;
}


/*--------------------------------------------------------------------*/
// Output the domain box associated with this grid
/** 
 *//*-----------------------------------------------------------------*/
inline
Box&
LevelOP::boxdomLD()
{
  CH_assert(m_ldDef);
  return m_domainBox;
}

/*--------------------------------------------------------------------*/
// Return the AMR refinement ratios
/** 
 *//*-----------------------------------------------------------------*/
inline
IntVect&
LevelOP::getAMRRefRatio()
{
  CH_assert(m_ldDef);
  return m_amrrefRatio;
}
 

/*--------------------------------------------------------------------*/
// Set the refinement ratio between the current level and the finest 
// single grid
/** 
 *//*-----------------------------------------------------------------*/
inline
void
LevelOP::setSolRefRate(const int a_solRefRate)
{
  m_solRef = a_solRefRate;
}


#endif
