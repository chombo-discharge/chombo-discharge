#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "EBAMRIO.H"
#include "EBIndexSpace.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "CH_Timer.H"
#include "CH_Attach.H"
#include "memusage.H"
#include "memtrack.H"
#include "MultiSphereIF.H"
#include "GeometryShop.H"
#include <string>
#include <iostream>
#include <sstream>
class PlotParams
{
public:
  string dataFileName;
  string geomFileName;
  string outpFileName;
};                  
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
////////////
void
parseCommandLine(PlotParams& a_params, int argc, char* argv[])
{
  pout() << "usage: sparseToDense -f datafilein  -e ebfilein -o fileout" << endl;
  bool dataFileFound = false;
  bool geomFileFound = false;
  bool outpFileFound = false;
  for(int iarg = 0; iarg < argc-1; iarg++)
    {
      if(strcmp(argv[iarg],"-f") == 0)
        {
          a_params.dataFileName = string(argv[iarg+1]);
          dataFileFound = true;
        }
      else if(strcmp(argv[iarg],"-e") == 0)
        {
          a_params.geomFileName = string(argv[iarg+1]);
          geomFileFound = true;
        }
      else if(strcmp(argv[iarg], "-o") == 0)
        {
          a_params.outpFileName = string(argv[iarg+1]);
          outpFileFound = true;
        }
    }
  if(!dataFileFound)
    {
      pout() << "no data filename given in inputs"  << endl;
      exit(0);
    }
  if(!geomFileFound)
    {
      pout() << "no eb geometry filename given in inputs"  << endl;
      exit(0);
    }
  if(!geomFileFound)
    {
      pout() << "no output filename given in inputs"  << endl;
      exit(0);
    }
}
/******/
void defineGeometry(const PlotParams& a_params)
{
  std::string ebis_file = a_params.geomFileName;
  pout() << " recreating  geometry from file " << ebis_file << endl;
  //define ebis anew from file input
#ifdef CH_USE_HDF5
  HDF5Handle handleIn(ebis_file, HDF5Handle::OPEN_RDONLY);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(handleIn);
  handleIn.close();
#else
  MayDay::Error("not compiled with hdf5");
#endif
}
/******/
void convertData(const PlotParams& a_params)
{
#ifdef CH_USE_HDF5
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  HDF5Handle handleIn(a_params.dataFileName, HDF5Handle::OPEN_RDONLY);
  pout() << "reading data from file"  << a_params.dataFileName << endl;
  int             numLevels;
  Vector<int>     refRat;
  ProblemDomain   domain;
  RealVect        coarseDx = RealVect::Unit;
  IntVect         ghost;
  metaDataEBFile(handleIn, numLevels, refRat, domain,  ghost);
  Vector<LevelData<EBCellFAB> *> data(numLevels, NULL);
  for(int ilev = 0; ilev < numLevels; ilev++)
    {
      data[ilev] = new LevelData<EBCellFAB>();
      readCellCentered(handleIn, ilev, ebisPtr, ghost[0], data[ilev]);
    }
  handleIn.close();

  pout() << "outputting data to a dense file that visit can read" << endl;
  //this stuff can come in through the command line if people want that
  Real dx = 1; Real dt = 1; Real time = 1;

  int nComp = data[0]->nComp();

  Vector<Real>  coveredValues(nComp, 0.);
  Vector<DisjointBoxLayout> grids(numLevels);
  Vector<string>            names(nComp);
  for(int ivar = 0; ivar < nComp; ivar++)
    {
      names[ivar] = string("var") + convertInt(ivar);
    }
  for(int ilev = 0; ilev < numLevels; ilev++)
    {
      grids[ilev] = data[ilev]->disjointBoxLayout();
    }

  //the true is replace covered
  writeEBHDF5(a_params.outpFileName, grids, data, 
              names, domain, dx, dt, time, refRat, 
              numLevels, true, coveredValues, ghost);

  for(int ilev = 0; ilev < numLevels; ilev++)
    {
      delete data[ilev];
    }
#else
  MayDay::Error("not compiled with hdf5");
#endif
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    PlotParams params;
    parseCommandLine(params, argc, argv);

    defineGeometry(params);
    
    convertData(params);
    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();
  }
  // End scoping trick
#ifdef CH_MPI
  MPI_Finalize();
#endif
}
