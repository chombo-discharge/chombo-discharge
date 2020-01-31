#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIFFUSIONPARAMS_H_
#define _DIFFUSIONPARAMS_H_

#include <string>
using std::string;

#include "IntVect.H"
#include "RealVect.H"
#include "ProblemDomain.H"

#include "UsingNamespace.H"

/// A class to hold all the solver parameters
class DiffusionParams
{
  public:
    /// Null constructor
    DiffusionParams();

    /// Print the parameters
    void print();

    /// Type of geometry to use
    string m_geometry;

    /// Diffusion constant
    Real m_diffusionConstant;

    /// Whether to use the variable coefficient solver
    bool m_useVariableCoeff;

    /// How much the cofficient swings around the mean of m_diffusionConstant
    Real m_diffusionEps;

    /// Strengths of the source and sink
    Real m_sourceScaling;
    Real m_sinkScaling;

    /// Solution at t = 0 (constant)
    Real m_initialValue;

    /// The physical domain corner and cell size (at the coarsest level)
    RealVect m_loCorner;
    Real m_dx;

    /// The number of AMR levels
    int m_numLevels;

    /// Number of cells in each direction at the coarsest level
    IntVect m_numCells;

    /// The refinement ratio between each pair of levels
    Vector<int> m_refRatio;

    /// The coarsest problem domain
    ProblemDomain m_coarsestDomain;

    /// Time step and time to end the simulation
    Real m_dt;
    Real m_endTime;

    /// How often to output and a prefix for the output name
    int m_outputInterval;
    string m_outputPrefix;

    /// Used to determine box sizes and alignments
    int m_maxBoxSize;
    int m_blockFactor;

    /// Used to generate AMR grids
    Real   m_fillRatio;
    int    m_nestingRadius;
    string m_tagType;

    /// Multigrid parameters
    int  m_mgNumCycles;
    int  m_mgNumSmooths;
    int  m_mgRelaxType;
    bool m_mgLazyRelax;
    Real m_mgToler;
    Real m_mgHangToler;
    int  m_mgIterMax;
    int  m_mgNumPrecondIter;

    /// Ghost cell information
    int m_numGhostEBISLayout;
    IntVect m_numGhostSoln;
    IntVect m_numGhostSource;

    /// Type of refluxing to use in the AMR operator
    int m_whichReflux;
};

#endif
