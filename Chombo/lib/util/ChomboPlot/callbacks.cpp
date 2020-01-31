#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctype.h>
#include <cmath>
#include <string>
#include "SPACE.H"

#include "libsx.h"
#include "libsx_private.h"
#include "REAL.H"
#include "chomboPlot.h"
#include "datatype.h"
#include "callbacks.h"
#include "LoHiSide.H"
#include "DebugDump.H"
#include "DatasetClient.H"
#include "AMRIO.H"
#include "PiecewiseLinearFillPatch.H"
using std::string;

int g_xsize;
int g_ysize;
const int g_infob = 125;
const int g_paletteb = 200;
const int g_numcont = 30;
const int g_minlen = 200;
/************************************************/
void do_contourplot(Widget w,void *data)
{
  datatype *me = (datatype *)data;
  me->plottype = ContourPlot;
  /* no need to redisplay --- no picture on small window
     redisplay(NULL,me->width,me->height,me);
  */
}
/************************************************/
void do_colormap(Widget w,void *data)
{
  datatype *me = (datatype *)data;
  me->plottype = ColorMap;
  /* no need to redisplay --- no picture on small window
     redisplay(NULL,me->width,me->height,me);
  */
}
/*****************/
void setDefaultCM(datatype *me)
{
  for (int i = 0; i< 256;i++)
    {
      me->rcol[i] = i;
      me->gcol[i] = i;
      me->bcol[i] = i;

    }
  me->rcol[0] = 0;
  me->gcol[0] = 0;
  me->bcol[0] = 0;
}
/*****************/
void setColorMap(void *data)
{
  datatype *me = (datatype *)data;
  int isize = 256;
  SetMyColorMap(isize,me->rcol,me->gcol,me->bcol);
}
/************************************************/
Real roundInt(Real input)
{
  int retvalue = (int)(input+0.5);
  return ((Real) retvalue);
}
/************************************************/
Real roundDown(Real input)
{
  int retvalue = (int)(input);
  return ((Real) retvalue);
}
/************************************************/
Real roundUp(Real input)
{
  int retvalue = (int)(input+0.99);
  return ((Real) retvalue);
}

// -------------------------------------------------------------
void
DrawDataSheet(LevelData<FArrayBox>* levelDataInPtr)
{
  const LevelData<FArrayBox>& levelDataIn = *levelDataInPtr;
  LevelData<FArrayBox> levelDataLoc(levelDataIn.getBoxes(),
                                    levelDataIn.nComp(),
                                    IntVect::Zero);
  Interval interv = levelDataIn.interval();
  levelDataIn.copyTo(interv, levelDataLoc, interv);
  MultiArrayViewFab(&levelDataLoc);
}

void keypress(Widget w, char* input, int upordown,  void *data)
{
  //if i do not check for this, it will do everrything twice
  //once on the press, once on the release
  if (upordown == 1)
    {
      char charin = *input;
      if ((charin == 'L') || (charin == 'l'))
        {
          load( w, data);
        }
      else if ((charin == 'v') || (charin == 'V'))
        {
          getVar( w, data);
        }
      else if ((charin == 'p') || (charin == 'P'))
        {
          dumpPS( w, data);
        }
      else if ((charin == 'S') || (charin == 's'))
        {
          subregion( w, data);
        }
      else if ((charin == 'n') || (charin == 'N'))
        {
          getNumCont( w, data);
        }
      else if ((charin == 't') || (charin == 'T'))
        {
          drawboxes( w, data);
        }
      else if ((charin == 'i') || (charin == 'I'))
        {
          normaldir( w, data);
        }
      else if ((charin == 'c') || (charin == 'C'))
        {
          close( w, data);
        }
      else if ((charin == 'D') || (charin == 'd'))
        {
          makeDataSheet(w,data);
        }
      else if ((charin == 'Q') || (charin == 'q'))
        {
          quit( w, data);
        }
      else
        {
          //          cout << "i am confused about the key you pressed =="
          //               << charin << endl;
        }
    }

}
void makeDataSheet(Widget w,void *data)
{
  datatype* me = (datatype *) data;
  char* filename;

  filename = GetString("\nEnter Level NUMBER for DataSheet\n","");

  if (filename)
    {
      int itest =  atoi(filename);
      int numlevels = me->vect_mf.size();
      LevelData<FArrayBox>* displevel = NULL;
      if ((itest >= 0) && (itest < numlevels))
        {
          displevel = me->vect_mf[itest];
        }
      else if (itest < 0)
        {
          GetYesNo("negative level number entered - displaying level 0\n");
          displevel = me->vect_mf[0];
        }
      else
        {
          GetYesNo("level entered > maxlevel - displaying level 0\n");
          displevel = me->vect_mf[0];
        }
      DrawDataSheet(displevel);
      redisplay(NULL,me->width, me->height,me);
    }
  else
    {
      /*psyche*/
    }

}
/************************************/
void subregion(Widget w,void *data)
{
  datatype* me = (datatype *) data;
  if (!me->mousebox.isEmpty())
    {
      datatype* menew = new datatype;
      datatypeequate(menew, me);
      menew->doingsubregion = true;
      menew->lev0subbox = me->mousebox;
      menew->num_windows++;
      menew->window = MakeWindow("ChomboPlot", NULL, NONEXCLUSIVE_WINDOW);

      if (menew->plottype == ContourPlot)
        makecontourplotter(menew);
      else if (menew->plottype == ColorMap)
        makecolormapper(menew);
      else
        {
          MayDay::Warning("resetting plottype to contourplot");
          menew->plottype = ContourPlot;
          makecontourplotter(menew);
        }
    }
}

/************************************/
void quit(Widget w,void *data)
{
  cleanUp((datatype*) data);
  exit(0);
}
/************************************/
void close(Widget w,void *arg)
{
  datatype *wi=(datatype *)arg;

  wi->num_windows--;

  SetCurrentWindow(XtParent(XtParent(w)));
  CloseWindow();
}
/************************************/
void drawboxes(Widget w,void *data)
{
  datatype *me = (datatype *) data;
  if (me->drawboxes)
    me->drawboxes = false;
  else
    me->drawboxes = true;
  redisplay(NULL,me->width, me->height,me);
}
/************************************/
void cleanUp(datatype *me)
{

  for (int ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
    if (me->vect_mf[ilev] != NULL) delete me->vect_mf[ilev];
  me->vect_box.resize(0);
  me->vect_ratio.resize(0);
  me->vect_mf.resize(0);
  me->vect_char.resize(0);
  me->vect_dbl.resize(0);
}
/************************************/
//fill ghost cells in interior
/************************************/
void makeCFStuff(datatype *data)
{
  int ncomps = data->vect_mf[0]->nComp();
  int interpRad = 1;
  data->vect_mf[0]->exchange (data->vect_mf[0]->interval());
  for (int  ilev = 1; ilev < int(data->vect_mf.size()); ilev++)
    {
      PiecewiseLinearFillPatch filler(data->vect_dbl[ilev],
                                      data->vect_dbl[ilev-1],
                                      ncomps,
                                      data->vect_box[ilev-1],
                                      data->vect_ratio[ilev-1],
                                      interpRad);

      LevelData<FArrayBox>& mfcur = *data->vect_mf[ilev];
      LevelData<FArrayBox>& uccur = *data->vect_mf[ilev-1];

      Real timeCoeff = 0.5;
      filler.fillInterp(mfcur, uccur, uccur,
                        timeCoeff, 0, 0, mfcur.nComp());
      mfcur.exchange(mfcur.interval());
    }
}
/************************************/
void file_callback(Widget w, char *str, int index, void *data)
{
  datatype *me = (datatype *)data;

  if (str)
    {
      datatype* menew = new datatype;
      datatypeequate(menew, me);
      menew->num_windows++;
      menew->window = MakeWindow("ChomboPlot", NULL, NONEXCLUSIVE_WINDOW);

      int status = fillFromFile(menew, str);
      if (status != 0)
        {
          delete menew;
        }
      else
        {
          makeCFStuff(menew);

          if (menew->plottype == ContourPlot)
            makecontourplotter(menew);
          else if (menew->plottype == ColorMap)
            makecolormapper(menew);
          else
            {
              MayDay::Warning("resetting plottype to ContourPlot");
              menew->plottype = ContourPlot;
              makecontourplotter(menew);
            }
        }
    }
  else
    {
    }
}
void makecontourplotter(void* data)
{

  datatype* me = (datatype *)data;

  Widget w0 = MakeButton("Load",   load,     me);
  Widget w2 = MakeButton("Variable", getVar,  me);
  Widget w21= MakeButton("PSDump", dumpPS,    me);
  Widget w22= MakeButton("DataSheet", makeDataSheet, me);
  Widget w3 = MakeButton("Close",    close,     me);
  Widget w10 = MakeButton("Subregion",  subregion,     me);

  Widget w4 = MakeDrawArea(me->width,me->height, redisplay, me);
  me->drawwidget = w4;
  Widget w5 = MakeScrollList(me->filein, 175, me->height, file_callback, me);
  Widget w6 = MakeLabel("ChomboPlot");
  Widget w7 = MakeButton("Number of Contours", getNumCont,  me);
  Widget w8 = MakeButton("Toggle Drawboxes",    drawboxes,     me);

  SetWidgetPos(w0, PLACE_UNDER, w6, NO_CARE, NULL);
  SetWidgetPos(w2, PLACE_UNDER, w0, NO_CARE, NULL);
  SetWidgetPos(w21, PLACE_UNDER, w2, NO_CARE, NULL);
  SetWidgetPos(w22, PLACE_UNDER, w21, NO_CARE, NULL);
  SetWidgetPos(w7, PLACE_UNDER, w22, NO_CARE, NULL);
  SetWidgetPos(w8, PLACE_UNDER, w7, NO_CARE, NULL);
  SetWidgetPos(w3, PLACE_UNDER, w8, NO_CARE, NULL);
  SetWidgetPos(w10, PLACE_UNDER, w3, NO_CARE, NULL);
#if (CH_SPACEDIM==3)
  Widget w11 = MakeButton("INormal DIR",  normaldir,     me);
  SetWidgetPos(w11, PLACE_UNDER, w10, NO_CARE, NULL);
  Widget w12 = MakeButton("sLice Position",  slicepos,     me);
  SetWidgetPos(w12, PLACE_UNDER, w11, NO_CARE, NULL);
#endif
  SetWidgetPos(w5, PLACE_RIGHT, w7, NO_CARE, NULL);
  SetWidgetPos(w4, PLACE_RIGHT, w5, NO_CARE, NULL);

  /* button up for this draw area */
  SetButtonDownCB(w4, button_down);
  SetButtonUpCB(w4,   button_up);
  SetKeypressCB(w4,   keypress);

  GetStandardColors();

#if (CH_SPACEDIM==3)
  SetFgColor(w11,  BLACK);
  SetBgColor(w11,  WHITE);
  SetFgColor(w12,  BLACK);
  SetBgColor(w12,  WHITE);
#endif
  SetFgColor(w10,  BLACK);
  SetFgColor(w0,  BLACK);
  SetFgColor(w2,  BLACK);
  SetFgColor(w21, BLACK);
  SetFgColor(w3,  BLACK);
  SetFgColor(w4,  BLACK);
  SetFgColor(w5,  BLACK);
  SetFgColor(w6,  BLACK);
  SetFgColor(w7,  BLACK);
  SetFgColor(w8,  BLACK);
  SetFgColor(w22,  BLACK);

  SetBgColor(w22,  WHITE);
  SetBgColor(w10,  WHITE);
  SetBgColor(w0,  WHITE);
  SetBgColor(w2,  WHITE);
  SetBgColor(w21, WHITE);
  SetBgColor(w3,  WHITE);
  SetBgColor(w5,  WHITE);
  SetBgColor(w6,  WHITE);
  SetBgColor(w7,  WHITE);
  SetBgColor(w8,  WHITE);
  ShowDisplay();
}

void makecolormapper(void* data)
{

  datatype* me = (datatype *)data;

  Widget w0 = MakeButton("Load",   load,     me);
  Widget w2 = MakeButton("Variable", getVar,  me);
  Widget w21= MakeButton("PSDump", dumpPSCM,    me);
  Widget w22= MakeButton("DataSheet", makeDataSheet, me);
  Widget w3 = MakeButton("Close",    close,     me);
  Widget w10 = MakeButton("Subregion",  subregion,     me);

  Widget w4 = MakeDrawArea(me->width + g_paletteb, me->height, redisplay, me);
  me->drawwidget = w4;
  XFont draw_font;
  draw_font = GetFont("fixed");
  draw_font = GetFont("9x15");
  Widget w5 = MakeScrollList(me->filein, 175, me->height, file_callback, me);
  Widget w6 = MakeLabel("ChomboPlot");
  //Widget w7 = MakeButton("Change Color Map", loadColorMap,  me);
  Widget w7 = MakeMenu("Change Color Map");
  Widget w71 = MakeMenuItem(w7, "Greyscale", loadGreyscale, me);
  Widget w72 = MakeMenuItem(w7, "Color", loadColorScale, me);
  Widget w73 = MakeMenuItem(w7, "Load Custom Map", loadColorMap,  me);
  Widget w8 = MakeButton("Toggle DrawBoxes",    drawboxes,     me);

  SetWidgetPos(w0, PLACE_UNDER, w6, NO_CARE, NULL);
  SetWidgetPos(w2, PLACE_UNDER, w0, NO_CARE, NULL);
  SetWidgetPos(w21, PLACE_UNDER, w2, NO_CARE, NULL);
  SetWidgetPos(w22, PLACE_UNDER, w21, NO_CARE, NULL);
  SetWidgetPos(w7, PLACE_UNDER, w22, NO_CARE, NULL);
  SetWidgetPos(w8, PLACE_UNDER, w7, NO_CARE, NULL);
  SetWidgetPos(w3, PLACE_UNDER, w8, NO_CARE, NULL);
  SetWidgetPos(w10, PLACE_UNDER, w3, NO_CARE, NULL);
#if (CH_SPACEDIM==3)
  Widget w11 = MakeButton("INormal DIR",  normaldir,     me);
  SetWidgetPos(w11, PLACE_UNDER, w10, NO_CARE, NULL);
  Widget w12 = MakeButton("sLice Position",  slicepos,     me);
  SetWidgetPos(w12, PLACE_UNDER, w11, NO_CARE, NULL);
#endif
  SetWidgetPos(w5, PLACE_RIGHT, w7, NO_CARE, NULL);
  SetWidgetPos(w4, PLACE_RIGHT, w5, NO_CARE, NULL);

  /* button up for this draw area */
  SetButtonDownCB(w4, button_down);
  SetButtonUpCB(w4,   button_up);
  SetKeypressCB(w4,   keypress);

  GetStandardColors();

#if (CH_SPACEDIM==3)
  SetFgColor(w11,  BLACK);
  SetBgColor(w11,  WHITE);
  SetFgColor(w12,  BLACK);
  SetBgColor(w12,  WHITE);
#endif
  SetFgColor(w10,  BLACK);
  SetFgColor(w0,  BLACK);
  SetFgColor(w2,  BLACK);
  SetFgColor(w21, BLACK);
  SetFgColor(w3,  BLACK);
  SetFgColor(w4,  BLACK);
  SetFgColor(w5,  BLACK);
  SetFgColor(w6,  BLACK);
  SetFgColor(w7,  BLACK);
  SetFgColor(w71,  BLACK);
  SetFgColor(w72,  BLACK);
  SetFgColor(w73,  BLACK);
  SetFgColor(w8,  BLACK);
  SetFgColor(w22,  BLACK);

  SetBgColor(w22,  WHITE);
  SetBgColor(w10,  WHITE);
  SetBgColor(w0,  WHITE);
  SetBgColor(w2,  WHITE);
  SetBgColor(w21, WHITE);
  SetBgColor(w3,  WHITE);
  SetBgColor(w5,  WHITE);
  SetBgColor(w6,  WHITE);
  SetBgColor(w7,  WHITE);
  SetBgColor(w71,  WHITE);
  SetBgColor(w72,  WHITE);
  SetBgColor(w73,  WHITE);
  SetBgColor(w8,  WHITE);
  ShowDisplay();
}
/*****/
void datatypeinit(datatype* data, int argc, char **argv)
{
  data->vect_box.resize(0);
  data->vect_char.resize(0);
  data->vect_ratio.resize(0);
  data->vect_mf.resize(0);
  data->xdraw = g_xsize;
  data->ydraw = g_ysize;
  data->width = g_xsize;
  data->numcont = g_numcont;
  data->height = g_ysize + g_infob;
  data->drawboxes = true;
  data->inormal = 2;
  data->idepth = 0;
  data->filetype = AMRBinary;
  data->plottype = ContourPlot;
#if 1
  int argpos = 1;
  if ( (argc > 1) && !strcmp(argv[1],"-ebamrascii"))
    {
      data->filetype = EBAMRASCII;
      ++argpos;
    }
  if ( (argc > 1) && !strcmp(argv[1],"-sf"))
    {
      data->filetype = AMRBinary;
      ++argpos;
    }
  if ( (argc > 1) && !strcmp(argv[1],"-xdang"))
    {
      data->filetype = XDangFAB;
      ++argpos;
    }
  if ( (argc > 1) && !strcmp(argv[1],"-amrascii"))
    {
      data->filetype = AMRASCII;
      ++argpos;
    }
  if ( (argc > 1) && !strcmp(argv[1],"-fab"))
    {
      data->filetype = SingleFAB;
      ++argpos;
    }
#endif
  data->drawwidget = NULL;
  data->hascolorfile = false;
  //eb addenda
  data->ebflag = false;
  data->coveredCells.resize(0);
  data->multiValuedCells.resize(0);
  setDefaultCM(data);
}

void resetWidths(datatype* me)
{
  ///hmmmbegin
  Box bbase = me->vect_box[0];
  if (me->doingsubregion)
    bbase = me->lev0subbox;
  int xlen = bbase.size(me->axisdir[0]);
  int ylen = bbase.size(me->axisdir[1]);
  int magiclen = g_xsize/2;
  while ((xlen <= magiclen)&&(ylen <= magiclen))
    {
      xlen *= 2;
      ylen *= 2;
    }
  xlen = Max(xlen, g_minlen);
  ylen = Max(ylen, g_minlen);
  xlen = Min(xlen, g_xsize);
  ylen = Min(ylen, g_ysize);
  me->xdraw = xlen;
  me->ydraw = ylen;
  me->height= ylen + g_infob;
  me->width = xlen;
  me->infob = g_infob;

  ////hmmmend
}
void datatypeequate(datatype* dtto, datatype* dtfrom)
{
  dtto->vect_box = dtfrom->vect_box;
  dtto->vect_ratio = dtfrom->vect_ratio;
  dtto->vect_mf = dtfrom->vect_mf;
  dtto->vect_char = dtfrom->vect_char;
  dtto->cur_var=dtfrom->cur_var;
  dtto->xdraw=dtfrom->xdraw;
  dtto->ydraw=dtfrom->ydraw;
  dtto->zdraw=dtfrom->zdraw;
  dtto->plottype=dtfrom->plottype;
  dtto->width=dtfrom->width;
  dtto->height=dtfrom->height;
  dtto->infob=dtfrom->infob;
  dtto->max=dtfrom->max;
  dtto->min=dtfrom->min;
  dtto->mag=dtfrom->mag;
  dtto->eps=dtfrom->eps;
  dtto->drawboxes=dtfrom->drawboxes;
  dtto->inormal=dtfrom->inormal;
  dtto->idepth=dtfrom->idepth;
  dtto->numcont=dtfrom->numcont;
  dtto->axisdir= dtfrom->axisdir;
  dtto->filetype= dtfrom->filetype;

  int nfiles = dtfrom->nfiles;
  for (int i = 0; i < nfiles; i++)
    {
      dtto->filein[i] = new char[250];
      strcpy(dtto->filein[i],dtfrom->filein[i]);
    }
  strcpy(dtto->filename,dtfrom->filename);
  dtto->filein[nfiles] = NULL;
  dtto->nfiles = nfiles;
  dtto->lev0subbox = dtfrom->lev0subbox;
  dtto->mousebox = dtfrom->mousebox;
  dtto->doingsubregion = dtfrom->doingsubregion;
  dtto->mouse_start_pix = dtfrom->mouse_start_pix;
  dtto->mouse_end_pix = dtfrom->mouse_end_pix;
  dtto->hascolorfile = dtfrom->hascolorfile;
  if (dtfrom->hascolorfile)
    {
      strcpy(dtto->filecolors,dtfrom->filecolors);
      readColorsFromFile(dtto);
    }
  else
    {
      setDefaultCM(dtto);
    }
  //eb addenda
  dtto->ebflag = dtfrom->ebflag;
  dtto->coveredCells = dtfrom->coveredCells;
  dtto->multiValuedCells = dtfrom->multiValuedCells;
}
void readColorsFromFile(datatype* me)
{
  fstream is;
  is.open(me->filecolors, ios::in);
  if ( !is.is_open() )
    {
      MayDay::Warning ( "fillFromFile: file open failed" );
      return;
    }
  for (int i =0; i < 256; i++)
    {
      int ir, ig, ib;
      is  >> ir >> ig >> ib;
      me->rcol[i] = ir;
      me->gcol[i] = ig;
      me->bcol[i] = ib;
    }
  me->rcol[0] = 0;
  me->gcol[0] = 0;
  me->bcol[0] = 0;

  me->hascolorfile = true;
  is.close();
}

//
// return code:
// 0 = success
// -1 = file open failed
// -2 readHierarchy failed
int
fillFromFile(datatype* me, char* filename)
{
  strcpy(me->filename,filename);
  me->cur_var  = 0;

  fstream is;
  is.open(filename, ios::in);
  if ( !is.is_open() )
    {
      MayDay::Warning ( "fillFromFile: file open failed" );
      return (-1);
    }
  cleanUp(me);

  string filestring(filename);
  if ((me->filetype ==  AMRBinary)||
      (me->filetype ==  AMRASCII) ||
      (me->filetype ==   EBAMRASCII))
    {
      int status = readHierarchy(me, filestring);
      if (status != 0)
        {
          MayDay::Warning ( "fillFromFile: readHierarchy failed" );
          return (-2);
        }
    }
  else if (me->filetype ==  SingleFAB || me->filetype == XDangFAB)
    {
      readSingleFAB(me, is);
    }
  else
    {
      cerr << "illegal file type " << me->filetype << endl;
      abort();
    }
  int nratios = me->vect_ratio.size();
  int ilastratio = 1;
  if (nratios > 0)
    {
      ilastratio = me->vect_ratio[nratios-1];
    }
  me->vect_ratio.push_back(ilastratio);
  //reset stuff that should not be continuous from file to file
  me->lev0subbox = me->vect_box[0];
  me->mousebox = me->vect_box[0];
  me->mouse_start_pix = IntVect::Zero;
  me->mouse_end_pix = IntVect::Zero;
  getAxes(me);
  resetWidths(me);
  me->doingsubregion = false;

  return (0);
}
/***********************************/
//
// return codes:
// 0 = success
// -1 = ReadAMRHierarchy failed
int
readHierarchy(datatype* me, const string& filename)
{
  int numlevels;
  Box basedom;
  Real dx, dt, time;
  if (me->filetype ==  AMRBinary)
    {
#ifdef CH_USE_HDF5
      int status =
        ReadAMRHierarchyHDF5(filename,
                             me->vect_dbl,
                             me->vect_mf,
                             me->vect_char,
                             basedom,
                             dx, dt, time,
                             me->vect_ratio,
                             numlevels,
                             IntVect::Unit);
      if ( status != 0 )
        {
          MayDay::Warning ( "readHierarchy: ReadAMRHierarchyHDF5 failed" );
          return (-1);
        }
#else
      MayDay::Error("Attempted to use HDF filetype without compiling with HDF");
#endif
    }
  else if (me->filetype ==  AMRASCII)
    {
      int status =
        ReadAMRHierarchyASCII(filename,
                              me->vect_dbl,
                              me->vect_mf,
                              me->vect_char,
                              basedom,
                              dx, dt, time,
                              me->vect_ratio,
                              numlevels,
                              IntVect::Unit);
      if ( status != 0 )
        {
          MayDay::Warning ( "readHierarchy: ReadAMRHierarchyASCII failed" );
          return (-1);
        }
      me->ebflag = false;
    }
  else if (me->filetype ==  EBAMRASCII)
    {
      int status =
        ReadEBAMRASCII(filename,
                       me->vect_dbl,
                       me->vect_mf,
                       me->vect_char,
                       basedom,
                       me->coveredCells,
                       me->multiValuedCells,
                       me->vect_ratio,
                       numlevels,
                       IntVect::Unit);
      me->ebflag = true;
      if ( status != 0 )
        {
          MayDay::Warning ( "readHierarchy: ReadEBAMRASCII failed" );
          return (-1);
        }
    }
  else
    {
      MayDay::Warning("readHierarchy: unknown filetype");
      return (-1);
    }
  me->vect_box.resize(numlevels, basedom);
  for (int ilev = 1; ilev < numlevels; ilev++)
    me->vect_box[ilev] = refine(me->vect_box[ilev-1], me->vect_ratio[ilev-1]);

  return (0);
}
/***************************/
void
readSingleFAB(datatype* me, istream& is)
{
  me->ebflag = false;
  FArrayBox temp;
  if (me->filetype ==  SingleFAB )
    {
      readAFabASCII(is, temp);
    }
  else if ( me->filetype == XDangFAB)
    {
      readXDangFAB( me, is, temp);

    }
  else
    {
      cerr << "got to readirregfab in error. " << endl;
      cerr << " needs to be -fab or -xdang option" << endl;
      abort();
    }
  int nvar = temp.nComp();
  CH_assert(nvar > 0);
  Box fabbox = temp.box();
  CH_assert(!fabbox.isEmpty());
  DisjointBoxLayout ba;
  ba.addBox(fabbox, 0);
  ba.close();
  me->vect_mf.resize(1, NULL);
  me->vect_mf[0] =
    new LevelData<FArrayBox>(ba, nvar, IntVect::Unit);

  LevelData<FArrayBox>& mf = *me->vect_mf[0];
  DataIterator dit = mf.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    mf[dit()].copy(temp);
  me->vect_box.resize(1,fabbox);
  me->vect_ratio.resize(1,2);
  for (int ivar = 0; ivar < nvar; ivar++)
    {
      char name[250];
      sprintf(name, "fab%d", ivar);
      string astname(name);
      me->vect_char.push_back(astname);
    }
}
void
readXDangFAB(datatype* me,istream& is, FArrayBox& fabin)
{
  int nx,ny,nvar;
#if (CH_SPACEDIM==2)
  is >> nx >> ny >> nvar;
  IntVect ivlo(0,0);
  IntVect ivhi(nx-1,ny-1);
#elif (CH_SPACEDIM==3)
  int nz;
  is >> nx >> ny >> nz >> nvar;
  IntVect ivlo(0,0,0);
  IntVect ivhi(nx-1,ny-1,nz-1);
#else
  invalid spacedim;
#endif
  Box domain(ivlo,ivhi);
  fabin.resize(domain,nvar);

#if (CH_SPACEDIM==3)
  for (int k =0; k < nz; k++)
    {
#endif
      for (int j =0; j < ny; j++)
        {
          for (int i =0; i < nx; i++)
            {
              for (int ivar =0; ivar < nvar; ivar++)
                {
#if (CH_SPACEDIM==3)
                  int iloc = i + j*nx + k*nx*ny + ivar*nx*ny*nz;
#else
                  int iloc = i + j*nx + ivar*nx*ny;
#endif
                  Real* dataloc = fabin.dataPtr() + iloc;
                  Real rtemp;
                  is >> rtemp;
                  *dataloc = rtemp;
                }
            }
        }
#if (CH_SPACEDIM==3)
    }
#endif

}
void load(Widget w, void *data)
{
  datatype *meold = (datatype *)data;

  char* filename;

  filename = GetString("\nEnter filename to load\n","");

  if (filename)
    {
      datatype* menew = new datatype;
      datatypeequate(menew, meold);
      menew->num_windows++;
      menew->window = MakeWindow("ChomboPlot", NULL, NONEXCLUSIVE_WINDOW);

      int status = fillFromFile(menew, filename);
      if (status != 0)
        {
          delete menew;
        }
      else
        {
          // need to redisplay since the picture was changed
          makeCFStuff(menew);
          if (menew->plottype == ContourPlot)
            makecontourplotter(menew);
          else if (menew->plottype == ColorMap)
            makecolormapper(menew);
          else
            {
              MayDay::Warning("resetting plottype to ColorMap");
              menew->plottype = ContourPlot;
              makecontourplotter(menew);
            }
        }
    }
  else
    {
      /*psyche*/
    }
}
void loadColorMap(Widget w, void *data)
{
  datatype *meold = (datatype *)data;

  char* filename;
  filename = GetString("\nEnter Color Map Filename\n","");

  if (filename)
    {
      strcpy(meold->filecolors, filename);
      readColorsFromFile(meold);
      redisplay(NULL,meold->width, meold->height,meold);
    }
  else
    {
      /*psyche*/
    }
}

void loadColorScale(Widget w, void *data)
{
  datatype *meold = (datatype *)data;

  // this is a hacked way to do this
#include "defaultColors.H"

  meold->rcol[0] = 0;
  meold->gcol[0] = 0;
  meold->bcol[0] = 0;

  redisplay(NULL,meold->width, meold->height,meold);

}

void loadGreyscale(Widget w, void *data)
{
  datatype *meold = (datatype *)data;
  for (int i = 0; i< 256;i++)
    {
      meold->rcol[i] = 255-i;
      meold->gcol[i] = 255-i;
      meold->bcol[i] = 255-i;

    }
  meold->rcol[0] = 0;
  meold->gcol[0] = 0;
  meold->bcol[0] = 0;

  redisplay(NULL,meold->width, meold->height,meold);

}

/***************************/
void getVar(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char* filename;

  filename = GetString("\nEnter Variable NUMBER \n","");

  if (filename)
    {
      int itest =  atoi(filename);
      if ((itest >= 0) && (itest < me->vect_char.size()))
        {
          me->cur_var = itest;
          redisplay(NULL,me->width, me->height,me);
        }
      else
        {
          GetYesNo("bogus # entered --setting current to 0\n");
          me->cur_var = 0;
        }
    }
  else
    {
      /*psyche*/
    }
}
void normaldir(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char* filename;

  filename = GetString("\nEnter normal direction NUMBER (0,1,2) \n","");

  if (filename)
    {
      int itest =  atoi(filename);
      if ((itest >= 0) && (itest < CH_SPACEDIM))
        {
          me->inormal = itest;
        }
      else
        {
          GetYesNo("bogus # entered --setting normal to 2\n");
          me->inormal = 2;
        }
      getAxes(me);
      redisplay(NULL,me->width, me->height,me);
    }
  else
    {
      /*psyche*/
    }
}
void slicepos(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char* filename;

  filename = GetString("\nEnter slice position NUMBER on coarsest grid\n","");

  if (filename)
    {
      int ibiggest = me->vect_box[0].size(me->inormal);
      int itest =  atoi(filename);
      if ((itest >= 0) && (itest < ibiggest))
        {
          me->idepth = itest;
        }
      else
        {
          GetYesNo("bogus # entered --setting depth to 0\n");
          me->idepth = 0;
        }
      redisplay(NULL,me->width, me->height,me);
    }
  else
    {
      cout << "I did not grok the fullness of that filename\n";
    }
}
/***************************/
void getNumCont(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char* contour_string;

  contour_string = GetString("\nEnter Number of Contours \n","");

  int maxcont = 1000;
  if (contour_string)
    {
      int itest =  atoi(contour_string);
      if ((itest >= 0) && (itest <= maxcont))
        {
          me->numcont = itest;
          redisplay(NULL,me->width, me->height,me);
        }
      else if (itest >= 0)
        {
          GetYesNo("num_cont > maxcont --setting contours to maxcont\n");
          me->cur_var = maxcont;
        }
      else
        {
          GetYesNo("num_cont < 0 --setting contours to 30\n");
          me->cur_var = 30;
        }
    }
  else
    {
      /*psyche*/
    }
}
/************************************************/
void do_con_plot_cart(Widget w,void *data)
{
  datatype *me = (datatype *)data;
  /* need to redisplay since the picture was changed */
  redisplay(NULL,me->width,me->height,me);
}

/**************************************************/
/*
 * The following is the redisplay code for the
 * drawing area widget.
 *
 * Each time it needs to be redisplayed
 * (either because it the window was resized or because
 * it was obscured), this function gets called.
 */
/**************************************************/
void redisplay(Widget w, int width, int height, void *data)
{
  datatype *me = (datatype *)data;

  SetDrawArea(me->drawwidget);
  SetColor(WHITE);
  ClearDrawArea();
  draw_stuff(width, height,me);
}

/*************************************************/
void draw_stuff(int width, int height, void *data)
{
  datatype *me = (datatype *)data;
  if (me->vect_mf.size() == 0)
    {
      draw_text(width, height,me);
    }
  else if (me->plottype == ContourPlot)
    {
      FreeAllColors();
      GetStandardColors();
      con_plot_cart(width, height,me);
    }
  else if (me->plottype == ColorMap)
    {
      FreeAllColors();
      GetAllColors();
      colorMap(width, height,me);
    }
  else
    {
      FreeAllColors();
      GetStandardColors();
      draw_text(width, height,me);
    }
}

/*************************************************/
void getAxes(datatype *me)
{
  //do spacedim-dependent stuff
  Box bbase = me->vect_box[0];
  if (me->doingsubregion)
    bbase = me->lev0subbox;

#if (CH_SPACEDIM == 2)
  me->axisdir[0] = 0;
  me->axisdir[1] = 1;
#elif (CH_SPACEDIM == 3)
  int idir1, idir2;
  CH_assert(me->inormal >= 0 && me->inormal < SpaceDim);
  if (me->inormal == 2)
    {
      idir1 = 0;
      idir2 = 1;
    }
  else if (me->inormal == 1)
    {
      idir1 = 0;
      idir2 = 2;
    }
  else if (me->inormal == 0)
    {
      idir1 = 1;
      idir2 = 2;
    }
  if (bbase.size(idir1) >= bbase.size(idir2))
    {
      me->axisdir[0] = idir1;
      me->axisdir[1] = idir2;
    }
  else
    {
      me->axisdir[0] = idir2;
      me->axisdir[1] = idir1;
    }

#else
#error unimplemented CH_SPACEDIM
#endif
}
/********************************/
void getMaxMinMag(datatype *me)
{
  CH_assert((me->cur_var >= 0) &&
         (me->cur_var < me->vect_char.size()));
  Real fmax = -7.77e10;
  Real fmin =  7.77e10;
  Real fmag = -7.77e10;
  //find max, min, mag

  for (int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
    {
      const LevelData<FArrayBox>& mfcur = *(me->vect_mf[ilev]);
      DataIterator dit = mfcur.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& fabcur = mfcur[dit()];
          Box subbox = grow(fabcur.box(), -1);
          fmax = Max(fmax, fabcur.max(subbox, me->cur_var));
          fmin = Min(fmin, fabcur.min(subbox, me->cur_var));
        }
    }
  fmag = fmax - fmin;
  me->eps = 1.0e-10;
  if (fmag < me->eps)
    {
      fmag = 1.0;
    }
  me->max =fmax;
  me->min =fmin;
  me->mag =fmag;
}
/* useful function which linearly interpolates between v1 and v2
   (which exist at T1, T2) to valmid (which lives at time)
   --- returns valmid*/
Real spainterp(Real t1, Real t2, Real t, Real v1, Real v2)
{
  Real TCOld;
  Real TCNew;
  Real Time;
  Real valold;
  Real valnew;
  if (t2 > t1)
    {
      TCOld = t1;
      TCNew = t2;
      valold = v1;
      valnew = v2;
      Time = t;
    }
  else
    {
      TCOld = t2;
      TCNew = t1;
      valold = v2;
      valnew = v1;
      Time = t;
    }
  Real tcdiff = 1.0;
  Real tfdiff = 0.0;
  if ((TCNew - TCOld) > 1.0e-12)
    {
      tcdiff = TCNew -TCOld;
      tfdiff = Time -TCOld;
    }
  Real valmid = valold + tfdiff*(valnew-valold)/tcdiff;
  return valmid;

}
/***********************************************/
void dumpPS(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char fileout[300];

  sprintf(fileout,"%s.%s.cp.eps",me->filename,me->vect_char[me->cur_var].c_str());

  ofstream ps_strm(fileout);

  /* draw stuff */
  /* obligatory postscript begining */
  Real papersize = 8*72;
  Real rxdraw = me->xdraw;
  Real rydraw = me->ydraw;
  Real rimin = 0;
  Real rjmin = 0;
  Real rimax = papersize-36;
  Real rjmax = papersize-36;
  Real maxside = Max(rxdraw,rydraw);
  Real xfac = rxdraw/maxside;
  Real yfac = rydraw/maxside;
  rimax = rimax*xfac;
  rjmax = rjmax*yfac;

  rimax += 36;
  rjmax += 36;
  rimin = 36;
  rjmin = 36;

  int imax = int(rimax);
  int jmax = int(rjmax);
  int imin = int(rimin);
  int jmin = int(rjmin);

  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin << " " << imax << " " << jmax << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/st {stroke} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;
  ps_strm << "/cp {closepath} def" << endl;

  /* draw box around plot */
  if (!me->doingsubregion)
    {
      ps_strm << "3 setlinewidth" << endl;
      ps_strm << imin << " " << jmin << " mt" << endl;
      ps_strm << imin << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmin << " lt" << endl;
      ps_strm << "cp" << endl << "st" << endl;
    }

  /* plot */
  /* upper left corner is origin---ergo weirdness */

  Real rwidth = rimax - rimin;
  Real rheight = rjmax - rjmin;
  getMaxMinMag(me);
  getAxes(me);
  Real mindiff = 2.*me->eps;
  int numcont = me->numcont;
  Real diff = me->max - me->min;
  Real df = diff/numcont;
#if (CH_SPACEDIM == 3)
  int idepth = me->idepth;
#endif
  if (diff > mindiff)
    {
      for (int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
        {
#if (CH_SPACEDIM == 3)
          if (ilev > 0)
            idepth *= me->vect_ratio[ilev-1];
#endif
          //gets set in getridxpix
          Box domain;
          //gets set in getridxpix
          Real ridxpix;

          getridxpix(domain, ridxpix, rwidth, rheight, data, ilev);
          CH_assert(!domain.isEmpty());

          LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
          //ghost cells already filled
          DataIterator dit = mfcur.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& bigfab = mfcur[dit()];
              Box bigbox = bigfab.box();
              Box curbox = grow(bigbox, -1);

#if (CH_SPACEDIM == 3)
              Box threedbox = curbox;
              IntVect loiv = threedbox.smallEnd();
              IntVect hiiv = threedbox.bigEnd();
              loiv[me->inormal]= idepth;
              hiiv[me->inormal]= idepth;
              curbox= Box(loiv,hiiv);
              curbox &= threedbox;
#endif
              curbox &= domain;
              if (!curbox.isEmpty())
                {
                  //draw area in which to do contours
                  IntVect ivcurlo = curbox.smallEnd();
                  IntVect ivcurhi = curbox.bigEnd();
                  IntVect ivdomlo = domain.smallEnd();
                  ivcurlo -= ivdomlo;
                  ivcurhi -= ivdomlo;
                  ivcurhi += IntVect::Unit;

                  Tuple<int, 2> axdir = me->axisdir;
                  Real rilo = rimin+ivcurlo[axdir[0]]*ridxpix;
                  Real rjlo = rjmin+ivcurlo[axdir[1]]*ridxpix;
                  Real rihi = rimin+ivcurhi[axdir[0]]*ridxpix;
                  Real rjhi = rjmin+ivcurhi[axdir[1]]*ridxpix;
                  Real bwidth = rihi - rilo;
                  Real bhight = rjhi - rjlo;
                  //sets color to white
                  ps_strm << "1.0 setgray" << endl;
                  ps_strm << rilo   << " " << rjlo << " "
                          << bwidth << " " << bhight << " rf" << endl;
                  ps_strm << "0.5 setlinewidth" << endl;
                  //sets color to black
                  ps_strm << "0.0 setgray" << endl;

                  //stuff is not normalized to one
                  Box interiorbox = curbox;
                  for (int nc = 0; nc < numcont; nc++)
                    {
                      Real fcont = me->min + df*(nc + 0.5);
                      BoxIterator bit(interiorbox);
                      for (bit.begin(); bit.ok(); bit.next())
                        {
                          IntVect ivint = bit();
                          bool coveredCell = false;
                          //i am carefully nesting these so as to not
                          //kill performance in the case of no eb
                          if (me->ebflag)
                            {
                              if (me->coveredCells[ilev].contains(ivint))
                                coveredCell = true;
                            }
                          Vector<Tuple<Real, 2> > vecpoints;
                          Tuple<Real,5> xlocs;
                          Tuple<Real,5> ylocs;
                          getPSCorners(me,ivint,domain,
                                       ridxpix,rimin,rjmin,
                                       xlocs,ylocs);
                          if (!coveredCell)
                            {
                              Tuple<Real,5> flocs;
                              getDataLocs(me,ivint,domain,bigfab,
                                          flocs);
                              for (int iloc = 0; iloc <= 3; iloc++)
                                {
                                  Real fint = flocs[iloc];
                                  Real xint = xlocs[iloc];
                                  Real yint = ylocs[iloc];
                                  Real fout = flocs[iloc+1];
                                  Real xout = xlocs[iloc+1];
                                  Real yout = ylocs[iloc+1];
                                  if (
                                     ((fint > fcont)&& (fout < fcont))
                                     ||((fint < fcont)&& (fout > fcont))
                                     )
                                    {
                                      Real xpt =  roundInt(
                                                           spainterp(fint, fout, fcont, xint, xout)
                                                           );
                                      Real ypt =  roundInt(
                                                           spainterp(fint, fout, fcont, yint, yout)
                                                           );
                                      Tuple<Real, 2> tupp;
                                      tupp[0] = xpt;
                                      tupp[1] = ypt;
                                      vecpoints.push_back(tupp);
                                    }
                                }
                              if (vecpoints.size() > 0)
                                {
                                  Tuple<Real, 2> pt1 = vecpoints[0];
                                  ps_strm << pt1[0] << " " << pt1[1] << " mt"
                                          << endl;
                                  for (int ivec=0;ivec < int(vecpoints.size()); ivec++)
                                    {
                                      Tuple<Real, 2> pt2 = vecpoints[ivec];
                                      ps_strm << pt2[0] << " " << pt2[1] << " lt"
                                              << endl;
                                    }
                                  ps_strm << "st" << endl;
                                }
                            } //if (!covered)
                          else
                            {
                              //draw black box if cell is covered
                              Real rcolor = 0;
                              ps_strm << setprecision(3)
                                      << rcolor << "   "
                                      << rcolor << "   "
                                      << rcolor << "   "
                                      << "rgb" << endl;

                              Real rilo = roundInt(xlocs[0]);
                              Real rjlo = roundInt(ylocs[0]);

                              ps_strm << rilo    << "   "  <<  rjlo << "   "
                                      << roundUp(ridxpix) << "   "  << roundUp(ridxpix) << endl;
                              ps_strm << "rf"  << endl;
                              //ps_strm << "st" << endl;
                            }
                        }
                    }
                  //sets color to black
                  ps_strm << "0.0 setgray" << endl;
                  if (me->drawboxes)
                    {

                      ps_strm << "1 setlinewidth" << endl;
                      ps_strm << rilo << " " << rjlo << " mt" << endl;
                      ps_strm << rilo << " " << rjhi << " lt" << endl;
                      ps_strm << rihi << " " << rjhi << " lt" << endl;
                      ps_strm << rihi << " " << rjlo << " lt" << endl;
                      ps_strm << "cp" << endl
                              << "st" << endl;
                    }

                }
            }
        }
    }

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
}
/***********************************************/
void dumpPSCM(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char fileout[300];

  sprintf(fileout,"%s.%s.cm.eps",me->filename,me->vect_char[me->cur_var].c_str());

  ofstream ps_strm(fileout);

  /* draw stuff */
  /* obligatory postscript begining */
  Real papersize = 8*72;
  Real rxdraw = me->xdraw;
  Real rydraw = me->ydraw;
  Real rimin = 0;
  Real rjmin = 0;
  Real rimax = papersize-36;
  Real rjmax = papersize-36;
  Real maxside = Max(rxdraw,rydraw);
  Real xfac = rxdraw/maxside;
  Real yfac = rydraw/maxside;
  rimax = rimax*xfac;
  rjmax = rjmax*yfac;

  rimax += 36;
  rjmax += 36;
  rimin = 36;
  rjmin = 36;

  int imax = int(rimax);
  int jmax = int(rjmax);
  int imin = int(rimin);
  int jmin = int(rjmin);

  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin << " " << imax << " " << jmax << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;
  ps_strm << "/st {stroke} def" << endl;
  ps_strm << "/cp {closepath} def" << endl;

  /* draw box around plot */
  if (!me->doingsubregion)
    {
      ps_strm << "3 setlinewidth" << endl;
      ps_strm << imin << " " << jmin << " mt" << endl;
      ps_strm << imin << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmin << " lt" << endl;
      ps_strm << "cp" << endl << "st" << endl;
    }

  /* plot */
  /* upper left corner is origin---ergo weirdness */

  Real rwidth = rimax - rimin;
  Real rheight = rjmax - rjmin;
  getMaxMinMag(me);
  getAxes(me);
  Real mindiff = 2.*me->eps;
  int numcont = me->numcont;
  Real diff = me->max - me->min;
#if (CH_SPACEDIM == 3)
  int idepth = me->idepth;
#endif
  if (diff > mindiff)
    {
      for (int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
        {
#if (CH_SPACEDIM == 3)
          if (ilev > 0)
            idepth *= me->vect_ratio[ilev-1];
#endif
          //gets set in getridxpix
          Box domain;
          //gets set in getridxpix
          Real ridxpix;

          getridxpix(domain, ridxpix, rwidth, rheight, data, ilev);
          CH_assert(!domain.isEmpty());

          LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
          //ghost cells already filled
          DataIterator dit = mfcur.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& bigfab = mfcur[dit()];
              Box bigbox = bigfab.box();
              Box curbox = grow(bigbox, -1);

#if (CH_SPACEDIM == 3)
              Box threedbox = curbox;
              IntVect loiv = threedbox.smallEnd();
              IntVect hiiv = threedbox.bigEnd();
              loiv[me->inormal]= idepth;
              hiiv[me->inormal]= idepth;
              curbox= Box(loiv,hiiv);
              curbox &= threedbox;
#endif
              curbox &= domain;
              if (!curbox.isEmpty())
                {
                  //draw area in which to do contours
                  IntVect ivcurlo = curbox.smallEnd();
                  IntVect ivcurhi = curbox.bigEnd();
                  IntVect ivdomlo = domain.smallEnd();
                  ivcurlo -= ivdomlo;
                  ivcurhi -= ivdomlo;
                  ivcurhi += IntVect::Unit;

                  Tuple<int, 2> axdir = me->axisdir;
                  Real rilo = rimin+ivcurlo[axdir[0]]*ridxpix;
                  Real rjlo = rjmin+ivcurlo[axdir[1]]*ridxpix;
                  Real rihi = rimin+ivcurhi[axdir[0]]*ridxpix;
                  Real rjhi = rjmin+ivcurhi[axdir[1]]*ridxpix;
                  Real bwidth = rihi - rilo;
                  Real bhight = rjhi - rjlo;
                  //sets color to white
                  ps_strm << "1.0 setgray" << endl;
                  ps_strm << rilo   << " " << rjlo << " "
                          << bwidth << " " << bhight << " rf" << endl;
                  ps_strm << "0.5 setlinewidth" << endl;
                  //sets color to black
                  ps_strm << "0.0 setgray" << endl;

                  //stuff is not normalized to one
                  Box interiorbox = curbox;
                  for (int nc = 0; nc < numcont; nc++)
                    {
                      BoxIterator bit(interiorbox);
                      for (bit.begin(); bit.ok(); bit.next())
                        {
                          IntVect ivint = bit();
                          bool coveredCell = false;
                          //i am carefully nesting these so as to not
                          //kill performance in the case of no eb
                          if (me->ebflag)
                            {
                              if (me->coveredCells[ilev].contains(ivint))
                                coveredCell = true;
                            }
                          Vector<Tuple<Real, 2> > vecpoints;
                          Tuple<Real,5> xlocs;
                          Tuple<Real,5> ylocs;
                          getPSCorners(me,ivint,domain,
                                       ridxpix,rimin,rjmin,
                                       xlocs,ylocs);

                          Real fval = bigfab(ivint, me->cur_var);
                          Real fmag = me->mag;
                          Real fmin = me->min;
                          Real fscale = (fval-fmin)/fmag;

                          int ir, ig, ib;
                          if (coveredCell)
                            {
                              ir = 0;
                              ig = 0;
                              ib = 0;
                            }
                          else
                            {
                              int icolor = (int)(253.*fscale);
                              icolor =abs(icolor)+2;
                              if ((icolor>256)||(icolor<0)) icolor = 0;
                              ir = me->rcol[icolor];
                              ig = me->gcol[icolor];
                              ib = me->bcol[icolor];
                            }
                          Real rcolor = ((Real) ir)/256.;
                          Real gcolor = ((Real) ig)/256.;
                          Real bcolor = ((Real) ib)/256.;
                          ps_strm << setprecision(3)
                                  << rcolor << "   "
                                  << gcolor << "   "
                                  << bcolor << "   "
                                  << "rgb" << endl;

                          Real rilo = roundInt(xlocs[0]);
                          //strangeness with coordinate system means
                          // high is low,  war is peace
                          //freedom is slavery and ignorance is strength
                          Real rjlo = roundInt(ylocs[0]);

                          ps_strm << rilo    << "   "  <<  rjlo << "   "
                                  << roundUp(ridxpix) << "   "  << roundUp(ridxpix) << endl;
                          ps_strm << "rf"  << endl;

                        }
                    }
                  //sets color to black
                  ps_strm << "0.0 setgray" << endl;

                  if (me->drawboxes)
                    {
                      /* draw box around plot */
                      ps_strm << "1 setlinewidth" << endl;
                      ps_strm << rilo << " " << rjlo << " mt" << endl;
                      ps_strm << rilo << " " << rjhi << " lt" << endl;
                      ps_strm << rihi << " " << rjhi << " lt" << endl;
                      ps_strm << rihi << " " << rjlo << " lt" << endl;
                      ps_strm << "cp" << endl
                              << "st" << endl;
                    }
                }
            }
        }
    }

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
  dumpPSPalCM(w, data);
}
/***********************************************/
void dumpPSPalCM(Widget w, void *data)
{
  datatype *me = (datatype *)data;
  char fileout[300];

  sprintf(fileout,"%s.%s.pal.eps",me->filename,me->vect_char[me->cur_var].c_str());

  ofstream ps_strm(fileout);

  /* draw stuff */
  /* obligatory postscript begining */
  Real papersize = 8*72;
  Real rxdraw = g_paletteb;
  Real rydraw = me->ydraw;
  Real rimin = 0;
  Real rjmin = 0;
  Real rimax = papersize-36;
  Real rjmax = papersize-36;
  Real maxside = Max(rxdraw,rydraw);
  //Real xfac = .5 * rxdraw/maxside;
  Real xfac = 0.25;
  Real yfac = .5 * rydraw/maxside;
  rimax *= xfac;
  rjmax *= yfac;
  /*
    rimax += 36;
    rjmax += 36;
    rimin += 36;
    rjmin += 36;
  */

  int imax = int(rimax);
  int jmax = int(rjmax);
  int imin = int(rimin);
  int jmin = int(rjmin);

  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  // below, 12 is the font size
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin-12 << " " << imax << " " << jmax+12 << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/Helvetica 12 selectfont" << endl;
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;

  getMaxMinMag(me);
  Real fmin = me->min;
  Real fmag = me->mag;
  int ioffset = 20;

  Real  rbarheight = rjmax - rjmin;

  Real rdycol = rbarheight/256.0;
  Real rdxcol = ioffset;
  Real rilo = rimin;
  Real rjlo = rjmin;
  // the following loop starts at 1, not 0, to avoid displaying color 0,
  // which is always black regardless of the colormap.
  for (int ivec = 1; ivec < 256 ; ivec++)
    {
      Real rcolor = ((Real) me->rcol[ivec])/256.;
      Real gcolor = ((Real) me->gcol[ivec])/256.;
      Real bcolor = ((Real) me->bcol[ivec])/256.;
      ps_strm << rcolor << "   "
              << gcolor << "   "
              << bcolor << "   "
              << "rgb" << endl;
      ps_strm << rilo << " " << rjlo << " mt" << endl;
      ps_strm << rilo    << "   "  <<  rjlo << "   "
              << rdxcol << "   "  << rdycol << endl;
      ps_strm << "rf"  << endl;
      rjlo += rdycol;
    }

  int num_lines = 8;
  Real line_spacing = rbarheight/num_lines;
  Real dval = fmag/num_lines;
  Real value = fmin;
  rilo = rimin + 2*rdxcol;
  // 5 is about half the height of the characters
  rjlo = rjmin - 5;

  ps_strm << "0 setgray" << endl;

  for (int ivec = 0; ivec <= num_lines; ivec++)
    {
      ps_strm << rilo << " " << rjlo << " mt" << endl;
      ps_strm << "( " << value << " ) show" << endl;
      value +=  dval;
      rjlo  += line_spacing;
    }

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
}

void con_plot_cart(int width, int height,void *data)
{
  datatype *me = static_cast<datatype *>(data);
  SetColor(WHITE);
  DrawFilledBox(0,0,me->width,me->height);
  getMaxMinMag(me);
  getAxes(me);

  /* upper left corner is origin---ergo weirdness */
  Real mindiff = 2.*me->eps;
  int numcont = me->numcont;
  Real diff = me->max - me->min;
  Real df = diff/numcont;
#if (CH_SPACEDIM == 3)
  int idepth = me->idepth;
#endif
  if (diff > mindiff)
    {
      for (int  ilev = 0; ilev < me->vect_mf.size(); ilev++)
        {
#if (CH_SPACEDIM == 3)
          if (ilev > 0)
            idepth *= me->vect_ratio[ilev-1];
#endif
          Box domain= me->vect_box[ilev];
          Real ridxpix = -1.0;
          Real rxdraw = me->xdraw;
          Real rydraw = me->ydraw;
          getridxpix(domain, ridxpix, rxdraw, rydraw, data, ilev);
          CH_assert(!domain.isEmpty());
          CH_assert(ridxpix >= 0);

          LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
          DataIterator dit = mfcur.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& bigfab = mfcur[dit()];
              Box bigbox = bigfab.box();
              Box curbox = grow(bigbox, -1);

#if (CH_SPACEDIM == 3)
              Box threedbox = curbox;
              IntVect loiv = threedbox.smallEnd();
              IntVect hiiv = threedbox.bigEnd();
              loiv[me->inormal]= idepth;
              hiiv[me->inormal]= idepth;
              curbox= Box(loiv,hiiv);
              curbox &= threedbox;
#endif
              //this culls out stuff
              curbox &= domain;
              if (!curbox.isEmpty())
                {
                  //draw area to do contours in
                  IntVect ivcurlo = curbox.smallEnd();
                  IntVect ivcurhi = curbox.bigEnd();
                  IntVect ivdomlo = domain.smallEnd();
                  ivcurlo -= ivdomlo;
                  ivcurhi -= ivdomlo;
                  ivcurhi += IntVect::Unit;

                  Tuple<int, 2> axdir = me->axisdir;
                  Real rilo = ivcurlo[axdir[0]]*ridxpix;
                  Real rjlo = ivcurlo[axdir[1]]*ridxpix;
                  Real rihi = ivcurhi[axdir[0]]*ridxpix;
                  Real rjhi = ivcurhi[axdir[1]]*ridxpix;

                  int ilo = (int)(roundInt(rilo));
                  int jlo = (int)(roundInt(rjlo));
                  int ihi = (int)(roundInt(rihi));
                  int jhi = (int)(roundInt(rjhi));
                  //origin weirdness
                  jlo = me->height - jlo;
                  jhi = me->height - jhi;
                  int bwidth = ihi - ilo;
                  int bhight = jhi - jlo;

                  //                int ivar = me->cur_var;
                  SetColor(WHITE);
                  DrawFilledBox(ilo,jlo,bwidth,bhight);
                  Box interiorbox = curbox;
                  for (int nc = 0; nc < numcont; nc++)
                    {
                      Real fcont = me->min + df*(nc + 0.5);
                      //stuff is not normalized to one
                      BoxIterator bit(interiorbox);
                      for (bit.begin(); bit.ok(); bit.next())
                        {
                          IntVect ivint = bit();
                          bool coveredCell = false;
                          //i am carefully nesting these so as to not
                          //kill performance in the case of no eb
                          if (me->ebflag)
                            {
                              if (me->coveredCells[ilev].contains(ivint))
                                coveredCell = true;
                            }
                          Vector<Tuple<int, 2> > vecpoints;

                          Tuple<Real,5> xlocs;
                          Tuple<Real,5> ylocs;
                          getCornerLocs(me,ivint,domain,ridxpix,
                                        xlocs,ylocs);
                          if (!coveredCell)
                            {
                              Tuple<Real,5> flocs;
                              getDataLocs(me,ivint,domain,bigfab, flocs);
                              for (int iloc = 0; iloc <= 3; iloc++)
                                {
                                  Real fint = flocs[iloc];
                                  Real xint = xlocs[iloc];
                                  Real yint = ylocs[iloc];
                                  Real fout = flocs[iloc+1];
                                  Real xout = xlocs[iloc+1];
                                  Real yout = ylocs[iloc+1];
                                  if (
                                     ((fint > fcont)&& (fout < fcont))
                                     ||((fint < fcont)&& (fout > fcont))
                                     )
                                    {
                                      Real xpt =  roundInt(
                                                           spainterp(fint, fout, fcont, xint, xout)
                                                           );
                                      Real ypt =  roundInt(
                                                           spainterp(fint, fout, fcont, yint, yout)
                                                           );
                                      Tuple<int, 2> tupp;
                                      tupp[0] = int(xpt);
                                      tupp[1] = int(ypt);
                                      vecpoints.push_back(tupp);
                                    }
                                }
                              //now draw lines at this location
                              SetColor(BLACK);
                              SetLineWidth(1);
                              int vecsize = vecpoints.size();
                              for (int ivec=0; ivec < vecsize-1; ivec++)
                                {
                                  Tuple<int, 2> pt1 = vecpoints[ivec];
                                  Tuple<int, 2> pt2 = vecpoints[ivec+1];
                                  DrawLine(pt1[0],pt1[1],pt2[0],pt2[1]);
                                }
                            }// end if (!coveredCell)
                          else
                            {
                              //draw black box for covered cell
                              Real rilo = roundDown(xlocs[0]);
                              Real rjlo = roundDown(ylocs[1]);
                              int iloloc = int(rilo);
                              int jloloc = int(rjlo);
                              Real ridpixround = roundUp(ridxpix);
                              int idxpix = (int) ridpixround;
                              SetColor(BLACK);
                              DrawFilledBox(iloloc,jloloc,idxpix,idxpix);
                            }
                        }//end loop over cells
                      if (me->drawboxes || ilev == 0)
                        {
                          if (ilev == 0)
                            {
                              SetColor(BLACK);
                            }
                          else if ( (ilev-1) % 3 == 0)
                            {
                              SetColor(BLUE);
                            }
                          else if ( (ilev-1) % 3 == 1)
                            {
                              SetColor(RED);
                            }
                          else if ( (ilev-1) % 3 == 2)
                            {
                              SetColor(GREEN);
                            }
                          SetLineWidth(1);
                          //draw box around area
                          DrawBox(ilo,jlo,bwidth,bhight);
                        }
                    } //end loop over contours
                } //end if (!curbox.isEmpty())
            }//end loop over grids
        } //end loop over levels
    }//end if (diff > mindiff)
  else
    {
      SetColor(WHITE);
      DrawFilledBox(0,0,width,height);
      SetColor(BLACK);
      char string[250];
      sprintf(string,"Constant = %e",me->max);
      DrawText(string, (width/2)-50, height/2);
    }
  infobfull( width, height, data);
}

void colorMap(int width, int height,void *data)
{
  setColorMap(data);
  datatype *me = static_cast<datatype *>(data);
  SetColor(WHITE);
  DrawFilledBox(0,0,me->width,me->height);
  getMaxMinMag(me);
  getAxes(me);

  /* upper left corner is origin---ergo weirdness */
  Real mindiff = 2.*me->eps;
  Real diff = me->max - me->min;
#if (CH_SPACEDIM == 3)
  int idepth = me->idepth;
#endif
  if (diff > mindiff)
    {
      for (int  ilev = 0; ilev < me->vect_mf.size(); ilev++)
        {
#if (CH_SPACEDIM == 3)
          if (ilev > 0)
            idepth *= me->vect_ratio[ilev-1];
#endif
          Box domain= me->vect_box[ilev];
          Real ridxpix = -1.0;
          Real rxdraw = me->xdraw;
          Real rydraw = me->ydraw;
          getridxpix(domain, ridxpix, rxdraw, rydraw, data, ilev);
          CH_assert(!domain.isEmpty());
          CH_assert(ridxpix >= 0);

          LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
          DataIterator dit = mfcur.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& bigfab = mfcur[dit()];
              Box bigbox = bigfab.box();
              Box curbox = grow(bigbox, -1);

#if (CH_SPACEDIM == 3)
              Box threedbox = curbox;
              IntVect loiv = threedbox.smallEnd();
              IntVect hiiv = threedbox.bigEnd();
              loiv[me->inormal]= idepth;
              hiiv[me->inormal]= idepth;
              curbox= Box(loiv,hiiv);
              curbox &= threedbox;
#endif
              //this culls out stuff
              curbox &= domain;
              if (!curbox.isEmpty())
                {
                  //draw area to do contours in
                  IntVect ivcurlo = curbox.smallEnd();
                  IntVect ivcurhi = curbox.bigEnd();
                  IntVect ivdomlo = domain.smallEnd();
                  ivcurlo -= ivdomlo;
                  ivcurhi -= ivdomlo;
                  ivcurhi += IntVect::Unit;

                  Tuple<int, 2> axdir = me->axisdir;
                  Real rilo = ivcurlo[axdir[0]]*ridxpix;
                  Real rjlo = ivcurlo[axdir[1]]*ridxpix;
                  Real rihi = ivcurhi[axdir[0]]*ridxpix;
                  Real rjhi = ivcurhi[axdir[1]]*ridxpix;

                  int ilo = (int)(roundInt(rilo));
                  int jlo = (int)(roundInt(rjlo));
                  int ihi = (int)(roundInt(rihi));
                  int jhi = (int)(roundInt(rjhi));
                  //origin weirdness
                  jlo = me->height - jlo;
                  jhi = me->height - jhi;
                  int bwidth = ihi - ilo;
                  int bhight = jhi - jlo;

                  SetColor(0);
                  DrawFilledBox(ilo,jlo,bwidth,bhight);
                  Box interiorbox = curbox;
                  //stuff is not normalized to one
                  BoxIterator bit(interiorbox);
                  for (bit.begin(); bit.ok(); bit.next())
                    {

                      IntVect ivint = bit();
                      bool coveredCell = false;
                      //i am carefully nesting these so as to not
                      //kill performance in the case of no eb
                      if (me->ebflag)
                        {
                          if (me->coveredCells[ilev].contains(ivint))
                            coveredCell = true;
                        }

                      //in case of embedded boundary, do nothing in case
                      //where volume fraction is zero
                      //vol frac is always last.
                      Tuple<Real,5> xlocs;
                      Tuple<Real,5> ylocs;
                      getCornerLocs(me,ivint,domain,ridxpix,
                                    xlocs,ylocs);
                      Real rilo = roundDown(xlocs[0]);
                      //strangeness with coordinate system means
                      // high is low,  war is peace
                      //freedom is slavery and ignorance is strength
                      //mad?  you think me mad?
                      //if anything, my senses have sharpened
                      Real rjlo = roundDown(ylocs[1]);
                      int iloloc = int(rilo);
                      int jloloc = int(rjlo);

                      if (coveredCell)
                        {
                          SetColor(BLACK);
                        }
                      else
                        {
                          Real fval = bigfab(ivint, me->cur_var);
                          Real fmag = me->mag;
                          Real fmin = me->min;
                          Real fscale = (fval-fmin)/fmag;

                          int icolor = (int)(253.*fscale);
                          icolor =abs(icolor)+2;
                          if ((icolor>256)||(icolor<0)) icolor =0;
                          SetColor(icolor);
                        }
                      Real ridpixround = roundUp(ridxpix);
                      int idxpix = (int) ridpixround;
                      DrawFilledBox(iloloc,jloloc,idxpix,idxpix);
                    }
                  if (me->drawboxes || ilev == 0)
                    {
                      if (ilev == 0)
                        {
                          SetColor(BLACK);
                        }
                      else if ( (ilev-1) % 3 == 0)
                        {
                          SetColor(BLUE);
                        }
                      else if ( (ilev-1) % 3 == 1)
                        {
                          SetColor(RED);
                        }
                      else if ( (ilev-1) % 3 == 2)
                        {
                          SetColor(GREEN);
                        }
                      SetLineWidth(1);

                      //draw box around area
                      DrawBox(ilo,jlo,bwidth,bhight);
                      //label box with level number if just drawing levels
                    }
                }
            }
        }
    }
  else
    {
      SetColor(WHITE);
      DrawFilledBox(0,0,width,height);
      SetColor(BLACK);
      char string[250];
      sprintf(string,"Constant = %e",me->max);
      DrawText(string, (width/2)-50, height/2);
    }
  infobfull( width, height, data);
  drawPalette( width, height, data);
}

void drawPalette(int width, int height, void *data)
{
  setColorMap(data);
  datatype *me = static_cast<datatype *>(data);
  SetColor(WHITE);
  DrawFilledBox(me->width,0,g_paletteb,me->height);
  getMaxMinMag(me);

  Real fmin = me->min;
  Real fmag = me->mag;
  int ioffset = 20;

  Real  rbarheight = me->height - 2*ioffset;

  Real rdycol = rbarheight/256.;
  Real rdxcol = ioffset;
  int idxcol = (int) rdxcol;
  int idycol = (int) roundUp(rdycol);
  int ilow = ioffset + me->width;
  Real rjlow = ioffset;
  // the following loop starts at 1, not 0, to avoid displaying color 0,
  // which is always black regardless of the colormap.
  for (int ivec = 1; ivec < 256 ; ivec++)
    {
      SetColor(ivec);
      int jlow = int(roundUp(rjlow));
      DrawFilledBox(ilow,jlow,idxcol,idycol);
      rjlow += rdycol;
    }

  int num_lines = 8;
  Real line_spacing = me->height/(num_lines+1);
  Real dval = fmag/num_lines;
  Real value = fmin;
  ilow = 3*ioffset+me->width;
  rjlow = ioffset;

  char str[200];
  for (int ivec = 0; ivec <= num_lines; ivec++)
    {
      sprintf(str, "%7.5e",value);
      int jlow = int(rjlow);
      DrawText(str, ilow, jlow);
      value +=  dval;
      rjlow += line_spacing;
    }
}

void  getCornerLocs(datatype* me,IntVect& iv, Box& domain, Real ridxpix,
                    Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs)
{
  Real dx = ridxpix;
  //start at lower left, go counter-clockwise
  IntVect ivdomlo = domain.smallEnd();
  int ivd = ivdomlo[me->axisdir[0]];
  int jvd = ivdomlo[me->axisdir[1]];
  int im = iv[me->axisdir[0]];
  int jm = iv[me->axisdir[1]];
  Real height = me->height;
  Real yl = height - (jm-jvd)*dx;
  Real yh = height - (jm+1-jvd)*dx;
  Real xl = (im-ivd)*dx;
  Real xh = (im-ivd+1)*dx;

  xlocs[0] = xl;
  ylocs[0] = yl;
  xlocs[1] = xl;
  ylocs[1] = yh;
  xlocs[2] = xh;
  ylocs[2] = yh;
  xlocs[3] = xh;
  ylocs[3] = yl;
  xlocs[4] = xlocs[0];
  ylocs[4] = ylocs[0];
}
void  getPSCorners(datatype* me,IntVect& iv, Box& domain, Real ridxpix,
                   Real rimin, Real rjmin,Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs)
{
  Real dx = ridxpix;
  //start at lower left, go counter-clockwise
  int im = iv[me->axisdir[0]];
  int jm = iv[me->axisdir[1]];
  IntVect ivdomlo = domain.smallEnd();
  int ivd = ivdomlo[me->axisdir[0]];
  int jvd = ivdomlo[me->axisdir[1]];
  Real yl = rjmin+(jm-jvd)*dx;
  Real yh = rjmin+(jm-jvd+1)*dx;
  Real xl = rimin+(im-ivd)*dx;
  Real xh = rimin+(im-ivd+1)*dx;

  xlocs[0] = xl;
  ylocs[0] = yl;
  xlocs[1] = xl;
  ylocs[1] = yh;
  xlocs[2] = xh;
  ylocs[2] = yh;
  xlocs[3] = xh;
  ylocs[3] = yl;
  xlocs[4] = xlocs[0];
  ylocs[4] = ylocs[0];
}
void getDataLocs(datatype* me,IntVect& iv,Box& domain,
                 FArrayBox& bigstate,Tuple<Real,5>& flocs)
{
  int ivar = me->cur_var;
  SideIterator sit;
  int idir0 = me->axisdir[0];
  int idir1 = me->axisdir[1];
  IntVect iv0 = BASISV(idir0);
  IntVect iv1 = BASISV(idir1);
  Real fmid = bigstate(iv, ivar);
  int nf = 1;
  //lower left corner of cell
  Tuple<IntVect, 3> outiv;

  Real fave = fmid;
  outiv[0] = iv - iv0 - iv1;
  outiv[1] = iv - iv0;
  outiv[2] = iv - iv1;
  for (int i = 0; i < 3; i++)
    {
      if (domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[0] = fave/nf;

  //upper left corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv - iv0 + iv1;
  outiv[1] = iv - iv0;
  outiv[2] = iv + iv1;
  for (int i = 0; i < 3; i++)
    {
      if (domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[1] = fave/nf;

  //upper right corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv + iv0 + iv1;
  outiv[1] = iv + iv0;
  outiv[2] = iv + iv1;
  for (int i = 0; i < 3; i++)
    {
      if (domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[2] = fave/nf;

  //lower right corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv + iv0 - iv1;
  outiv[1] = iv + iv0;
  outiv[2] = iv - iv1;
  for (int i = 0; i < 3; i++)
    {
      if (domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[3] = fave/nf;
  flocs[4] = flocs[0];
}
/***********************************************/
//
// write text in the info bar
void infobfull(int width, int height,void *data)
{
  datatype *me = (datatype *)data;
  int iwx,iwy,infob;
  char string[800];

  iwx = me->width;
  iwy = me->height;
  infob = me->infob;

  XFont draw_font;
  draw_font = GetFont("fixed");
  draw_font = GetFont("9x15");
  SetWidgetFont(me->drawwidget, draw_font);

  SetColor(WHITE);
  DrawFilledBox(0,0,iwx,infob);
  SetColor(BLACK);
  SetLineWidth(2);
  DrawLine(0,infob,iwx,infob);
  DrawLine(0,iwy+infob,iwx,iwy+infob);

  const int line_spacing = 15;

  DrawText(me->filename, 0, line_spacing);
  sprintf(string, "var = %s",me->vect_char[me->cur_var].c_str());
  DrawText(string, 0, 2*line_spacing);

  sprintf(string, "max = %e",me->max);
  DrawText(string, 0, 3*line_spacing);

  sprintf(string, "min = %e",me->min);
  DrawText(string, 0, 4*line_spacing);

#if CH_SPACEDIM == 2
  sprintf(string, "base grid = %d x %d",
          me->vect_box[0].size(0),
          me->vect_box[0].size(1) );
#endif
#if CH_SPACEDIM == 3
  sprintf(string, "base grid = %d x %d x %d",
          me->vect_box[0].size(0),
          me->vect_box[0].size(1),
          me->vect_box[0].size(2) );
#endif
  DrawText(string, 0, 5*line_spacing);

  const int num_levels = me->vect_box.size ();
  if ( num_levels > 1)
    {
      sprintf ( string,
                "numLevels = %d", num_levels );
      DrawText(string, 0, 6*line_spacing);

      sprintf ( string,
                "refRatios =");
      std::string stdstring(string);
      for (int level = 0; level < num_levels - 1; ++level)
        {
          sprintf ( string, " %d ", me->vect_ratio[level] );
          stdstring += string;
        }
      strcpy ( string, stdstring.c_str () );
    }
  else
    {
      sprintf ( string,
                "number of levels = 1" );
    }
  DrawText(string, 0, 7*line_spacing);
}
/***********************************/
void button_down(Widget w, int which_button, int x, int y, void *data)
{
  datatype *me=(datatype *) data;
  me->mouse_start_pix = IntVect::Zero;
  //set pixel mouse start to x and y
  me->mouse_start_pix[me->axisdir[0]] = x;
  me->mouse_start_pix[me->axisdir[1]] = y;

}
/***********************************/
/* this routine changes me->mousebox */
/*which becomes me->lev0subbox for new window */
/*if subregion is clicked. */
/*****/
void button_up(Widget w, int which_button, int x, int y, void *data)
{
  datatype *me=(datatype *) data;
  getAxes(me);

  int istart = me->mouse_start_pix[me->axisdir[0]];
  int jstart = me->mouse_start_pix[me->axisdir[1]];

  int ilo = Min(x, istart);
  int jlo = Min(y, jstart);
  int bwidth = Abs(x-istart);
  int bhight = Abs(y-jstart);
  SetColor(RED);
  SetLineWidth(1);
  DrawBox(ilo,jlo,bwidth,bhight);

  Real rxdraw = me->xdraw;
  Real rydraw = me->ydraw;
  //gets set in getridxpix
  Box domain;
  //gets set in getridxpix
  Real ridxpix;
  getridxpix(domain, ridxpix, rxdraw, rydraw, data, 0);
  CH_assert(!domain.isEmpty());

  int  idxpix = int(roundInt(ridxpix));
  //now transform to cell space
  if (idxpix > 0)
    {
      //origin weirdness
      int yloc0 = (me->height - jstart)/idxpix;
      int xloc0 = istart/idxpix;
      int yloc1 = (me->height - y)/idxpix;
      int xloc1 = x/idxpix;
      int ilo = Min(xloc0, xloc1);
      int jlo = Min(yloc0, yloc1);
      int ihi = Max(xloc0-1, xloc1-1);
      int jhi = Max(yloc0-1, yloc1-1);
      IntVect loend = IntVect::Zero;
      IntVect hiend = IntVect::Zero;
      loend[me->axisdir[0]] = ilo;
      loend[me->axisdir[1]] = jlo;
      hiend[me->axisdir[0]] = ihi;
      hiend[me->axisdir[1]] = jhi;
      me->mousebox = Box(loend, hiend);
      me->mousebox &= me->vect_box[0];
    }
}
//returns domain and ridxpix
void getridxpix(Box& domain, Real& ridxpix,
                Real xwidth, Real ywidth, void *data, int ilev)
{
  datatype *me=(datatype *) data;
  domain=me->vect_box[ilev];
  if ((me->doingsubregion)&&(!me->lev0subbox.isEmpty()))
    {
      if (domain.contains(me->lev0subbox))
        {
          domain = me->lev0subbox;
          for (int  iref = 0; iref < ilev; iref++)
            {
              domain.refine(me->vect_ratio[iref]);
            }
        }
    }

  Real dxt = xwidth/domain.size(me->axisdir[0]);
  Real dyt = ywidth/domain.size(me->axisdir[1]);
  ridxpix = Min(dxt, dyt);

}
/*****************/
void draw_text(int width, int height,void *data)
{
  datatype *me = (datatype *)data;
  SetColor(WHITE);
  DrawFilledBox(0,0,width,height);
  SetColor(BLACK);
  XFont draw_font;
  draw_font = GetFont("fixed");
  draw_font = GetFont("10x20");
  SetWidgetFont(me->drawwidget, draw_font);
  DrawText("need to load file", (width/2)-50, height/2);
}

/* This function sets up the display.  For any kind of a real program,
 * you'll probably want to save the values returned by the MakeXXX calls
 * so that you have a way to refer to the display objects you have
 * created (like if you have more than one drawing area, and want to
 * draw into both of them).
 */

void init_display(int argc, char **argv, datatype *me)
{
  /* parse input string */
  if (OpenDisplay(argc, argv) == FALSE)
    return;

  int xsize = DisplayWidth  ( lsx_curwin->display, lsx_curwin->screen ) * 7/8;
  int ysize = DisplayHeight ( lsx_curwin->display, lsx_curwin->screen ) * 7/8;
  g_xsize = 2;
  g_ysize = 2;
  while ( (2*g_xsize <= xsize) && (2*g_ysize <= ysize) )
    {
      g_xsize *= 2;
      g_ysize *= 2;
    }
  g_xsize = Max(g_xsize, g_minlen);
  g_ysize = Max(g_ysize, g_minlen);

  datatypeinit(me, argc, argv);
  int nfiles = 0;
  for (int i = 1; i < argc; i++)
    {
      if ( argv[i][0] != '-')
        {
          me->filein[nfiles] = new char[250];
          strcpy(me->filein[nfiles],argv[i]);
          ++nfiles;
        }
    }
  me->filein[nfiles] = NULL;
  me->nfiles = nfiles;

  Widget w0 = MakeButton("Load",   load,     me);
  Widget w3 = MakeButton("Quit",    quit,     me);
  Widget w4 = MakeMenu("Plot Type");
  Widget w41 = MakeMenuItem(w4,"Contour Plot",do_contourplot, me);
  Widget w42 = MakeMenuItem(w4,"Color Map",do_colormap, me);

  Widget w5 = MakeScrollList(me->filein, 175, 175, file_callback, me);
  Widget w6 = MakeLabel("ChomboPlot");

  SetWidgetPos(w0, PLACE_RIGHT, w6, NO_CARE, NULL);
  SetWidgetPos(w4, PLACE_RIGHT, w0, NO_CARE, NULL);
  SetWidgetPos(w3, PLACE_RIGHT, w4, NO_CARE, NULL);
  SetWidgetPos(w5, PLACE_RIGHT, w3, NO_CARE, NULL);
  GetStandardColors();
  SetFgColor(w5,  BLACK);
  SetBgColor(w5,  WHITE);
  SetFgColor(w4,  BLACK);
  SetBgColor(w4,  WHITE);
  SetFgColor(w6,  BLACK);
  SetBgColor(w6,  WHITE);
  SetFgColor(w3,  BLACK);
  SetBgColor(w3,  WHITE);
  SetFgColor(w0,  BLACK);
  SetBgColor(w0,  WHITE);
  SetFgColor(w41,  BLACK);
  SetBgColor(w41,  WHITE);
  SetFgColor(w42,  BLACK);
  SetBgColor(w42,  WHITE);

  /**
**/

  /* This call actually causes the whole thing to be displayed on the
   * screen.  You have to call this function before doing any drawing
   * into the window.
   */
  ShowDisplay();
}
