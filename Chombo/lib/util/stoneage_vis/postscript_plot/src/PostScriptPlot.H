#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POSTSCRIPTPLOT_H_
#define _POSTSCRIPTPLOT_H_

#include "Vector.H"
#include "Box.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "EBCellFAB.H"
#include "EBAlias.H"
#include "IntVectSet.H"
#include "PiecewiseLinearFillPatch.H"
#include <string>
#include "Tuple.H"

///sloppy class to accomodate ancient, sloppy code
class PostScriptData
{
public:
  PostScriptData(){;}

  PostScriptData(Vector< LevelData<FArrayBox>* >& a_data,
                 Vector<int>                    & a_refRat,
                 const string                   & a_filename,
                 int                            & a_ivar,
                 std::string                    & a_varname,
                 bool                           & a_drawBoxes,
                 int                            & a_normalDir,
                 int                            & a_slicePosition,
                 int                            & a_numContours,
                 bool                           & a_useGray,
                 Box                            & a_subBox,
                 bool                           & a_ebflag,
                 Vector<IntVectSet>             & a_coveredCells
                 );

  Vector<LevelData<FArrayBox>* > vect_mf;
  Vector<DisjointBoxLayout> vect_dbl;

  Vector<Box> vect_box;
  Vector<int> vect_ratio;
  std::string varname;
  std::string filename;
  std::string palettename;
  int cur_var;
  int xdraw;
  int ydraw;
  Real max;
  Real min;
  Real mag;
  Real eps;
  bool drawboxes;
  int inormal;
  int idepth;
  int numcont;
  Tuple<int, 2> axisdir;
  Box lev0subbox;
  bool doingsubregion;
  int rcol[256];
  int gcol[256];
  int bcol[256];
  //eb addenda
  bool ebflag;
  Vector<IntVectSet> coveredCells;
};

//relics of a former age 
void maxxx(Real a,Real b,Real c,Real *mx);
void minnn(Real a,Real b,Real c,Real* mn);
void HSVtoRGB( Real *r, Real *g, Real *b, Real hh, Real ss, Real vv );
void RGBtoHSV( Real rr, Real gg, Real bb, Real* h, Real *s, Real *v );

void setRedGreenCM(PostScriptData *me);
void setRainbowCM(PostScriptData *me);
void setGrayCM(PostScriptData *me);
void getMaxMinMag(PostScriptData *me);
Real spainterp(Real t1, Real t2, Real t, Real v1, Real v2);

void dumpContourPlot(PostScriptData* data);
void dumpColorMap(PostScriptData* data);
void dumpPseudoColor(PostScriptData* data);
void getAxes(PostScriptData* data);

void getDataLocs(PostScriptData* me,IntVect& iv,Box& domain,
                 FArrayBox& bigstate,Tuple<Real,5>& flocs);

void  getPSCorners(PostScriptData* me,IntVect& iv, Box& domain, Real dx, Real dy,
                   Real rimin, Real rjmin,Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs);

///Plotting utility for those who cannot stand windowing systems.
/**
   So.   You have finally come to the realization that all windowing
   systems suck.   You do not want to figure out the vaguries of
   OpenGL or XLib.   You do not want to have to port QT or the latest
   flavor of the month to every supercomputer and laptop you run on.
   You cannot get visit or ChomboVis to work because they are built on
   the unstable foundations of vtk and X and spending your life
   figuring out what went wrong sounds like the seventh circle of hell.

   Yet, somehow,  you need to visualize your data.

   This class will help you generate plots (either pseudocolor or
   contour --- someday, if Allah wills, spreadsheets---) that are
   plain, old postscript files.   You can look at them with gv or
   ghostscript or whatever.   The plots are 2d---so you have to set 
   slice directions and locations in 3d.

   Even more beautifully, you do not need HDF5 to work this.
   This can be done inline.
   Right now, this stuff is strictly serial.  
*/
class PostScriptPlot
{
public:

  ///
  static void getCoveredCells(Vector<IntVectSet>& a_ivs,
                              Vector<LevelData<EBCellFAB>* >& a_data);


  ///
  ~PostScriptPlot()
  {;}


  ///
  /**
   */
  static void writeContourPlot(Vector< LevelData<FArrayBox>* >& a_data,
                               Vector<int>                    & a_refRat,
                               const string                   & a_filename,
                               int                            & a_ivar,
                               std::string                    & a_varname,
                               bool                             a_drawBoxes = true,
                               int                              a_normalDir = 2,
                               int                              a_slicePosition = -1,
                               int                              a_numContours = 42,
                               Box                              a_subBox = Box())
  {

    bool useGray = true;
    bool ebflag = false;
    Vector<IntVectSet> coveredCells(a_data.size());
    PostScriptData steamy(a_data,
                          a_refRat,
                          a_filename,    
                          a_ivar,        
                          a_varname,     
                          a_drawBoxes,
                          a_normalDir, 
                          a_slicePosition,
                          a_numContours , 
                          useGray, 
                          a_subBox,
                          ebflag,
                          coveredCells);
    dumpContourPlot(&steamy);
  }

  ///
  /**
   */
  static void writeContourPlot(Vector< LevelData<EBCellFAB>* >& a_data,
                               Vector<int>                    & a_refRat,
                               const string                   & a_filename,
                               int                            & a_ivar,
                               std::string                    & a_varname,
                               bool                             a_drawBoxes = true,
                               int                              a_normalDir = 2,
                               int                              a_slicePosition = -1,
                               int                              a_numContours = 42,
                               Box                              a_subBox = Box())
  {
    Vector<LevelData<FArrayBox>*> mfdata(a_data.size());
    for(int ilev = 0; ilev < a_data.size(); ilev++)
      {
        aliasEB(*mfdata[ilev], *a_data[ilev]);
      }
    bool useGray = false;
    bool ebflag = true;
    Vector<IntVectSet> coveredCells(a_data.size());
    getCoveredCells(coveredCells, a_data);
    PostScriptData steamy( mfdata,
                           a_refRat,
                           a_filename,    
                           a_ivar,        
                           a_varname,     
                           a_drawBoxes,
                           a_normalDir, 
                           a_slicePosition,
                           a_numContours , 
                           useGray, 
                           a_subBox,
                           ebflag,
                           coveredCells);
    dumpContourPlot(&steamy);
  }          


  ///
  /**
   */
  static void writePseudoColor(Vector< LevelData<FArrayBox>* >& a_data,
                               Vector<int>                    & a_refRat,
                               const string&                    a_filename,
                               int                              a_ivar,
                               std::string&                     a_varname,
                               bool                             a_drawBoxes = true,
                               int                              a_normalDir = 2,
                               int                              a_slicePosition = -1,
                               bool                             a_useGrayScale = false,
                               Box                              a_subBox = Box())              
  {
    int numCont = -1;
    bool ebflag = false;
    Vector<IntVectSet> coveredCells(a_data.size());

    PostScriptData steamy(a_data,
                          a_refRat,
                          a_filename,    
                          a_ivar,        
                          a_varname,     
                          a_drawBoxes,
                          a_normalDir, 
                          a_slicePosition,
                          numCont,
                          a_useGrayScale,
                          a_subBox,
                          ebflag,
                          coveredCells);

    dumpColorMap(&steamy);
    dumpPseudoColor(&steamy);
  }


  ///
  /**
   */
  static void writePseudoColor(Vector< LevelData<EBCellFAB>* >& a_data,
                               Vector<int>                    & a_refRat,
                               const string&                    a_filename,
                               int                              a_ivar,
                               std::string&                     a_varname,
                               bool                             a_drawBoxes = true,
                               int                              a_normalDir = 2,
                               int                              a_slicePosition = -1,
                               bool                             a_useGrayScale = false,
                               Box                              a_subBox = Box())              
  {
    Vector<LevelData<FArrayBox>*> mfdata(a_data.size());
    for(int ilev = 0; ilev < a_data.size(); ilev++)
      {
        aliasEB(*mfdata[ilev], *a_data[ilev]);
      }
 
    int numCont = -1;
    bool ebflag = true;
    Vector<IntVectSet> coveredCells(a_data.size());
    getCoveredCells(coveredCells, a_data);
    PostScriptData steamy(mfdata,
                          a_refRat,
                          a_filename,    
                          a_ivar,        
                          a_varname,     
                          a_drawBoxes,
                          a_normalDir, 
                          a_slicePosition,
                          numCont,
                          a_useGrayScale,
                          a_subBox,
                          ebflag,
                          coveredCells);

    dumpColorMap(&steamy);
    dumpPseudoColor(&steamy);
  }

};




#endif
