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
#include <cstring>
using std::ifstream;
using std::ios;
using std::string;
#include "PostScriptPlot.H"
#include "AMRIO.H"

#include "CH_HDF5.H"


#include "UsingNamespace.H"
class PlotParams
{
public:
  PlotParams()
  {
    numContour = 42;
    useGray = false;
    drawBoxes = true;
    ivar  = 0;
    normaldir = 2;
    varname = string("var");
    outputname = string("output.ps");
  }
  string filename;
  string outputname;
  int numContour;
  bool useGray;
  bool drawBoxes;
  int ivar ;
  Box subbox;
  int normaldir;
  int slicePos;
  string varname;
};                  
////////////
void
parseCommandLine(PlotParams& a_params, int argc, char* argv[])
{
  pout() << "usage: plotFromHDF5 -i filein -o fileout -n numContours -b on/off -gray on/off -v varname -iv varnum -nd normaldir -d slicePos -xlo xlo -xhi xhi -ylo ylo -yhi yhi" << endl;
  bool filefound = false;
    int xlo, xhi, ylo, yhi;
  xlo =-1;
  xhi =-1;
  ylo =-1;
  yhi =-1;
  for(int iarg = 0; iarg < argc-1; iarg++)
    {
      if(strcmp(argv[iarg],"-i") == 0)
        {
          a_params.filename = string(argv[iarg+1]);
          filefound = true;
        }
      else if(strcmp(argv[iarg], "-o") == 0)
        {
          a_params.outputname = string(argv[iarg+1]);
        }
      else if(strcmp(argv[iarg],"-v") == 0)
        {
          a_params.varname = string(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-xlo") == 0)
        {
          xlo = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-xhi") == 0)
        {
          xhi = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-ylo") == 0)
        {
          ylo = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-yhi") == 0)
        {
          yhi = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-n") == 0)
        {
          a_params.numContour = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-d") == 0)
        {
          a_params.slicePos = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-nd") == 0)
        {
          a_params.normaldir = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-iv") == 0)
        {
          a_params.ivar = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-b") == 0)
        {
          a_params.drawBoxes = (strcmp(argv[iarg+1],"off") != 0);
        }
      else  if(strcmp(argv[iarg], "-g") == 0)
        {
          a_params.useGray = (strcmp(argv[iarg+1],"on") == 0);
        }
    }
  if((xlo >= 0) && (xhi >= 0) && (ylo >= 0) && (yhi >= 0))
    {
      IntVect lovect, hivect;
      if(SpaceDim == 2)
        {
          lovect[0] = xlo;
          lovect[1] = ylo;
          hivect[0] = xhi;
          hivect[1] = yhi;
        }
      else
        {
          int idir1, idir2;
          int inormal =  a_params.normaldir;
          if(inormal == 2)
            {
              idir1 = 0;
              idir2 = 1;
            }
          else if(inormal == 1)
            {
              idir1 = 0;
              idir2 = 2;
            }
          else if(inormal == 0)
            {
              idir1 = 1;
              idir2 = 2;
            }
          lovect[inormal] = a_params.slicePos;
          hivect[inormal] = a_params.slicePos;
          lovect[idir1] = xlo;
          lovect[idir2] = ylo;
          hivect[idir1] = xhi;
          hivect[idir2] = yhi;
        }
      a_params.subbox = Box(lovect, hivect);
      pout() << " using subbox " << a_params.subbox << endl;
    }
  if(!filefound)
    {
      pout() << "no filename input"  << endl;
      exit(0);
    }
}
////////////
int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc, &a_argv);
#endif
  {   
    PlotParams params;
    parseCommandLine(params, a_argc, a_argv);
    Vector<DisjointBoxLayout>        vectGrids;
    Vector<LevelData<FArrayBox>* >   vectData;
    Box                              domain;
    Vector<int>                      refRatio;
    int                              numLevels;

    ReadAMRHierarchyHDF5(params.filename,vectGrids, 
                         vectData, domain, refRatio, numLevels);

    string contourout = string("contour.") + params.outputname;
    string pseudoout =   string("pseudo.") + params.outputname;

    PostScriptPlot::writeContourPlot(vectData, refRatio, contourout, params.ivar,
                                     params.varname, params.drawBoxes, params.normaldir,
                                     params.slicePos, params.numContour, params.subbox);
    
    PostScriptPlot::writePseudoColor(vectData, refRatio, pseudoout, params.ivar,
                                     params.varname, params.drawBoxes, params.normaldir,
                                     params.slicePos, params.useGray, params.subbox);
    
  }
#ifdef CH_MPI
    MPI_Finalize();
#endif
}
