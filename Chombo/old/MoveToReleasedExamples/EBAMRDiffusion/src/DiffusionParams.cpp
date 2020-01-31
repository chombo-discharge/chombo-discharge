#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <vector>
using std::vector;

#include "parstream.H"
#include "ParmParse.H"

#include "DiffusionParams.H"

// Set all the solver parameters based on the values read from an input file
DiffusionParams::DiffusionParams()
{
  ParmParse pp;

  pp.get("geometry",m_geometry);

  pp.get("diffusion_constant",m_diffusionConstant);

  pp.get("use_variable_coeff", m_useVariableCoeff);
  if (m_useVariableCoeff)
    {
      pp.get("diffusion_eps",m_diffusionEps);
    }
  else
    {
      m_diffusionEps = 0.;
    }

  pp.get("source_scaling",m_sourceScaling);
  pp.get("sink_scaling",m_sinkScaling);

  pp.get("initial_value",m_initialValue);

  vector<Real> loCorner(SpaceDim);
  pp.getarr("lo_corner",loCorner,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    m_loCorner[idir] = loCorner[idir];
  }

  pp.get("dx",m_dx);

  pp.get("num_levels",m_numLevels);

  vector<int> numCells(SpaceDim);
  pp.getarr("num_cells",numCells,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    m_numCells[idir] = numCells[idir];
  }

  // Compute the indices of lo and hi end of the coarsest domain
  IntVect loEnd;
  IntVect hiEnd;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    loEnd[idir] = 0;
    hiEnd[idir] = m_numCells[idir] - 1;
  }

  m_refRatio.resize(m_numLevels);
  if (m_numLevels > 1)
  {
    vector<int> refRatio(m_numLevels-1);
    pp.getarr("ref_ratio",refRatio,0,m_numLevels-1);

    for (int ilev = 0; ilev < m_numLevels-1; ilev++)
    {
      m_refRatio[ilev] = refRatio[ilev];
    }
    m_refRatio[m_numLevels-1] = refRatio[m_numLevels-2];
  }
  else
  {
    m_refRatio[0] = 2;
  }

  // Make a box that size
  Box domainBox(loEnd,hiEnd);

  // Define a non-periodic problem domain
  m_coarsestDomain.define(domainBox);

  pp.get("dt",m_dt);
  pp.get("end_time",m_endTime);

  pp.get("output_interval",m_outputInterval);
  pp.get("output_prefix",m_outputPrefix);

  pp.get("max_box_size",m_maxBoxSize);
  pp.get("block_factor",m_blockFactor);

  pp.get("fill_ratio",m_fillRatio);
  m_nestingRadius = 2;
  pp.get("tag_type",m_tagType);

  pp.get("mg_num_cycles",m_mgNumCycles);
  pp.get("mg_num_smooths",m_mgNumSmooths);
  pp.get("mg_relax_type",m_mgRelaxType);
  pp.get("mg_lazy_relax",m_mgLazyRelax);
  pp.get("mg_toler",m_mgToler);
  pp.get("mg_hang_toler",m_mgHangToler);
  pp.get("mg_iter_max",m_mgIterMax);
  pp.get("mg_num_precond_iter",m_mgNumPrecondIter);

  m_numGhostEBISLayout = 4;

  m_numGhostSoln   = IntVect::Unit;
  m_numGhostSource = IntVect::Zero;

  pp.get("which_reflux",m_whichReflux);
}

// Print all the solver parameters
void DiffusionParams::print()
{
  pout() << "geometry = " << m_geometry << "\n";
  pout() << "\n";
  pout() << "diffusion constant = " << m_diffusionConstant << "\n";
  pout() << "\n";
  pout() << "source scaling = " << m_sourceScaling << "\n";
  pout() << "sink   scaling = " << m_sinkScaling   << "\n";
  pout() << "\n";
  pout() << "initial value = " << m_initialValue << "\n";
  pout() << "\n";
  pout() << "lo corner = " << m_loCorner << "\n";
  pout() << "\n";
  pout() << "dx = " << m_dx       << "\n";
  pout() << "\n";
  pout() << "num levels = " << m_numLevels << "\n";
  pout() << "num cells  = " << m_numCells  << "\n";
  if (m_numLevels > 1)
  {
    pout() << "ref ratio  = " << m_refRatio  << "\n";
  }
  pout() << "\n";
  pout() << "dt       = " << m_dt      << "\n";
  pout() << "end time = " << m_endTime << "\n";
  pout() << "\n";
  pout() << "output interval = " << m_outputInterval << "\n";
  pout() << "output prefix   = " << m_outputPrefix   << "\n";
  pout() << "\n";
  pout() << "max box size = " << m_maxBoxSize  << "\n";
  pout() << "block factor = " << m_blockFactor << "\n";
  pout() << "\n";
  pout() << "fill ratio = " << m_fillRatio << "\n";
  pout() << "tag type   = " << m_tagType   << "\n";
  pout() << "\n";
  pout() << "mg num cycles       = " << m_mgNumCycles      << "\n";
  pout() << "mg num smooths      = " << m_mgNumSmooths     << "\n";
  pout() << "mg relax type       = " << m_mgRelaxType      << "\n";
  pout() << "mg lazy relax       = " << m_mgLazyRelax      << "\n";
  pout() << "mg toler            = " << m_mgToler          << "\n";
  pout() << "mg hang toler       = " << m_mgHangToler      << "\n";
  pout() << "mg iter max         = " << m_mgIterMax        << "\n";
  pout() << "mg num precond iter = " << m_mgNumPrecondIter << "\n";
  pout() << "\n";
  pout() << "which reflux = " << m_whichReflux << "\n";
  pout() << "\n";
}
