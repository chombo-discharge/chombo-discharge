#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DATATYPE_H_
#define _DATATYPE_H_
#
/* This structure contains information relevant to your program.
 * You should fill it in with information that you need.
 *
 */
#include <X11/Intrinsic.h>
#include "Vector.H"
#include "Box.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "PiecewiseLinearFillPatch.H"
#include <string>
#include "Tuple.H"
using std::string;
/* define's */

/*
#define X_SIZE 512
#define Y_SIZE 512
#define INFOB 125
#define NUMCONT 30
*/

typedef LevelData<FArrayBox> MultiFab;

enum FileType
{
  AMRBinary, EBAMRASCII, SingleFAB, XDangFAB, AMRASCII
};

enum PlotType
{
  ContourPlot, ColorMap
};

typedef struct datatype
{
  Vector<Box> vect_box;
  Vector<int> vect_ratio;
  Vector<LevelData<FArrayBox>* > vect_mf;
  Vector<string> vect_char;
  Vector<DisjointBoxLayout> vect_dbl;
  int cur_var;
  int xdraw;
  int ydraw;
  int zdraw;
  char filename[250];
  PlotType plottype;
  int width;
  int height;
  int infob;
  Real max;
  Real min;
  Real mag;
  Real eps;
  bool drawboxes;
  int inormal;
  int idepth;
  int numcont;
  char* filein[1000];
  Tuple<int, 2> axisdir;
  FileType filetype;
  Widget window;
  Widget drawwidget;
  int num_windows;
  int nfiles;
  Box lev0subbox;
  Box mousebox;
  bool doingsubregion;
  IntVect mouse_start_pix;
  IntVect mouse_end_pix;
  unsigned char rcol[256];
  unsigned char gcol[256];
  unsigned char bcol[256];
  char filecolors[250];
  bool hascolorfile;
  //eb addenda
  bool ebflag;
  Vector<IntVectSet> coveredCells;
  Vector<IntVectSet> multiValuedCells;
} datatype;

#endif
