#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRTracerUtils.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"

#include "NamespaceHeader.H"

void
getProblemDomain(ProblemDomain& a_domain)
{
  ParmParse pp;
  Vector<int> ncell(SpaceDim);
  pp.getarr("n_cell", ncell, 0, SpaceDim);
  IntVect hiEnd(D_DECL(ncell[0]-1,ncell[1]-1,ncell[2]-1));
  Box level0Domain(IntVect::Zero, hiEnd);

  Vector<int> v_is_periodic(SpaceDim);
  pp.getarr("periodic_bc", v_is_periodic, 0, SpaceDim);

  bool is_periodic[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++) is_periodic[idir] = (v_is_periodic[idir]==1);

  ProblemDomain prob_domain(level0Domain.smallEnd(),
                            level0Domain.bigEnd(),
                            &(is_periodic[0]));

  a_domain = prob_domain;
}

void
getAdvectionVelocityFunction(AdvectionVelocityFunction& a_velFunc)
{
  ParmParse pp;
  int ivel;
  pp.get("advection_vel_type", ivel);
  if (ivel == 0)
    {
      pout() << "constant advection velocity" << endl;
      a_velFunc = constantAdvection;
    }
  else if (ivel == 1)
    {
      pout() << "rotating advection velocity" << endl;
      a_velFunc = rotatingAdvection;
    }
  else
    {
      MayDay::Error("bogus velocity flag");
    }
}

void
getAMRLevelTracerFactory(RefCountedPtr<AMRLevelTracerFactory>&  a_fact,
			 AdvectionVelocityFunction a_advFunc)
{
  ParmParse pp;
  Real cfl = 0.4;
  pp.get("cfl",cfl);

  Real domainLength;
  pp.get("domain_length",domainLength);

  Real maxGridSize;
  pp.get("fixed_box_size", maxGridSize);  

  a_fact = RefCountedPtr<AMRLevelTracerFactory>
    (new AMRLevelTracerFactory(a_advFunc, cfl, domainLength, maxGridSize));

}
void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelTracerFactory>&         a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat)
{
  ParmParse pp;
  int max_level = 0;
  pp.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals; // (num_read_levels,1);
  pp.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  int block_factor = 1;
  pp.get("fixed_box_size", block_factor);

  int max_grid_size = 32;
  pp.get("fixed_box_size", max_grid_size);

  Real fill_ratio = 0.75;
  pp.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  pp.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  pp.get("plot_interval",plot_interval);

  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);
  AMR amr;
  a_amr.define(max_level, a_refRat,
               a_prob_domain,&(*a_fact));

  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  pp.get("grid_buffer_size",gridBufferSize);
  a_amr.gridBufferSize(gridBufferSize);

  // set output parameters
  a_amr.checkpointInterval(checkpoint_interval);
  a_amr.plotInterval(plot_interval);
  a_amr.regridIntervals(regrid_intervals);
  a_amr.maxDtGrow(max_dt_growth);
  a_amr.dtToleranceFactor(dt_tolerance_factor);
  if (pp.contains("use_subcycling"))
    {
      bool useSubcycling;
      pp.get("use_subcycling", useSubcycling);
      if (!useSubcycling)
        {
          pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
        }
      a_amr.useSubcyclingInTime(useSubcycling);
    }
  if (pp.contains("plot_prefix"))
    {
      std::string prefix;
      pp.get("plot_prefix",prefix);
      a_amr.plotPrefix(prefix);
    }

  if (pp.contains("chk_prefix"))
    {
      std::string prefix;
      pp.get("chk_prefix",prefix);
      a_amr.checkpointPrefix(prefix);
    }

  int verbosity;
  pp.get("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  a_amr.verbosity(verbosity);
}

void
setupAMRHierarchy(Vector<Vector<Box> >& a_hierarchy)
{

  a_hierarchy.clear();

  int num_levels = 3;
  for (int l = 0; l < num_levels; ++l)
  {

    Vector<Box> levelBoxes;
    
    int numPointsLevel = 32*pow(2, l);
    int loCorner = (numPointsLevel - 32) / 2;
    int hiCorner = loCorner + 32 - 1;
    IntVect lo(D_DECL(loCorner, loCorner, loCorner));
    IntVect hi(D_DECL(hiCorner, hiCorner, hiCorner));
    Box levelDomain(lo, hi);

    domainSplit(levelDomain, levelBoxes, 8);
    a_hierarchy.push_back(levelBoxes);
  }
}

void
setupAMRForFixedHierarchyRun(AMR& a_amr, 
                             const Vector<Vector<Box> >& a_hierarchy)
{
  ParmParse pp;

  if (!pp.contains("restart_file"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      //a_amr.setupForFixedHierarchyRun(a_hierarchy);
      a_amr.setupForNewAMRRun();
    }
  else
    {
      std::string restart_file;
      pp.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      a_amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrTracer restart only defined with hdf5");
#endif
    }

  }

#include "NamespaceFooter.H"
