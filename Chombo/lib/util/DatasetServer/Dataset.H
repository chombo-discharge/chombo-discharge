#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// -------------------------------------------------------------
// Dataset.H
// -------------------------------------------------------------
#ifndef _DATASET_H_
#define _DATASET_H_
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <strstream>

extern "C"
{
#include <Xm/Xm.h>
#include <Xm/MainW.h>
#include <Xm/RowColumn.h>
#include <Xm/DrawingA.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/SelectioB.h>
#include <Xm/ScrollBar.h>
#include <Xm/ScrolledW.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/Scale.h>
#include <Xm/Label.h>

#include <unistd.h>
}

#include "Box.H"
// ParmParse is included for the PP_Array and PP_String classes
#include "ParmParse.H"
#include "FArrayBox.H"
#include "LayoutData.H"

const int MAXSTRINGWIDTH = 64;

enum DSPlane
{
  DS_XYPLANE = 0, DS_XZPLANE, DS_YZPLANE
};
enum DSDir
{
  DS_XDIR = 0, DS_YDIR, DS_ZDIR
};

class StringLoc;

class Dataset
{
public:
  //Dataset(BaseFab<Real> *fab);
  //Dataset(unsigned long *fab);  // to enter the address of a fab
  //Dataset(Real *realptr, int xlo, int ylo, int xhi, int yhi);
  Dataset(Widget parent, BaseFab<Real> *fab, const PP_String &formatstring,
          const PP_String &datalabel, bool is_irregular);
  Dataset(Widget parent,
          LayoutData<BaseFab<Real> >* layoutdata,
          const PP_String &formatstring,
          const PP_String &datalabel);
  virtual ~Dataset();

private:
  Widget         wParent, wDatasetTopLevel, wDatasetForm;
  Widget         wDatasetTools, wScrollArea, wPixArea;
  Widget         wFormat, wColorButton, wQuitButton;
  Widget         wHScrollBar, wVScrollBar;
  Widget         wPlaneRadioBox, wXYPlane, wXZPlane, wYZPlane;
  Widget         wSliceScale, wComponentScale, wElementScale, wIrregScale;
  Widget         wSliceLabel, wComponentLabel, wElementLabel, wIrregLabel;
  Widget         wGlobalMin, wGlobalMax, wSliceMin,  wSliceMax;
  Widget         wGlobalBadFPLabel, wSliceBadFPLabel;
  XtAppContext dataapp;
  Arg            args[20];
  Box            dataBox, sliceBox;
  BaseFab<Real> *sliceFab, *originalFab;
  LayoutData<BaseFab<Real> >*   originalMulti;
  unsigned int pixSizeX, pixSizeY;
  unsigned int dataAreaWidth, dataAreaHeight;
  unsigned int dataItemWidth, dataItemHeight;
  unsigned int characterWidth;
  unsigned int indexWidth, indexHeight;
  int            hStringOffset, vStringOffset;
  bool   stringOk, showColor;
  int            numStrings;
  char   buffer[BUFSIZ], dataString[MAXSTRINGWIDTH];
  PP_String      formatString, dataLabel;
  StringLoc   *dataStringArray, *hIndexArray, *vIndexArray;
  int    dragging, drags;
  PP_Array<Real>         dataMin, dataMax;
  PP_Array<bool>  minMaxCalculated;
  PP_Array<bool>         containsBadFloats;  // one for each fab component
  Real   sliceMin, sliceMax;
  int          maxComp;
  int          maxIrreg;
  Real  *dataPoint, *numIrreg;
  PP_String      hAxisString, vAxisString;
  int            hScrollBarPos, vScrollBarPos;
  int            hIndexAreaHeight, hIndexAreaStart, hIndexAreaEnd;
  int            vIndexAreaWidth,  vIndexAreaStart, vIndexAreaEnd;
  int            xh, yv;
  int            currentElement;
  PP_Array<int>  currentFabComponent;  // one per LayoutData element
  PP_Array<PP_Array<int> >       currentSlice; // one Array per LayoutData element
  PP_Array<DSPlane>      currentPlane; // one per LayoutData element
  int          currentIrreg;
  Dimension      width, height;
  DSDir          hDir, vDir, sDir;
  int            fpClass;
  bool   isMulti;
  bool         isIrregular;

  static PP_List<Dataset *> datasetList;
  static bool   topLevelCreated;

  // X baggage
  Display     *display;
  Visual      *visual;
  Screen      *screen;
  XVisualInfo  visualInfo;
  GC           gc;
  Window       root, dataWindow;
  int          screenNumber;
  int          blackIndex, whiteIndex;

  Dataset();
  void Init(Widget parent, const PP_String &formatstring,
            const PP_String &datalabel);
  bool Render(BaseFab<Real> *fab);
  bool Render(LayoutData<BaseFab<Real> >* layoutdata, int element);
  bool Render(Real *data, Real* num_irreg,
              const Box &databox, Real datamin, Real datamax);
  bool CalculateMinMax(Real *data, Real* num_irreg, int numpoints,
                       Real &fmin, Real &fmax, bool &badfloats);
  //PP_String FPClassString(Real fpNumber);
  //bool BadFPClass(Real fpNumber);

  void DoRaise();
  void DoExpose(Widget w, bool fromExpose);
  static void CBColorButton(Widget, XtPointer, XtPointer);
  static void CBQuitButton(Widget, XtPointer, XtPointer);
  static void CBReadString(Widget, XtPointer, XtPointer);
  static void CBPixInput(Widget, XtPointer, XtPointer);
  static void CBEndScrolling(Widget, XtPointer, XtPointer);
  static void CBScrolling(Widget, XtPointer, XtPointer);
  static void CBDoExposeDataset(Widget, XtPointer, XEvent *, Boolean *);
  static void CBResize(Widget, XtPointer, XEvent *, Boolean *);
  static void CBChangePlane(Widget, XtPointer, XtPointer);
  static void CBChangeSlice(Widget, XtPointer, XtPointer);
  static void CBChangeComponent(Widget, XtPointer, XtPointer);
  static void CBChangeElement(Widget, XtPointer, XtPointer);
  static void CBChangeIrreg(Widget, XtPointer, XtPointer);

  void DoColorButton(Widget w, bool onoff);
  void DoQuitButton();
  void DoReadString(Widget, XmSelectionBoxCallbackStruct *);
  void DoPixInput(XmDrawingAreaCallbackStruct *);
  void DoChangePlane(Widget w);
  void DoChangeSlice(int newslice);
  void DoChangeComponent(int newcomponent);
  void DoChangeElement(int newelement);
  void DoChangeIrreg(int newirreg);
  void DoResize(Widget w);
  void SetScrollIncrements();

  void DrawIndicies(const Box &tempBox);
  void DrawGrid(int startX, int startY, int finishX, int finishY,
                int gridspacingX, int gridspacingY,
                int foregroundIndex, int backgroundIndex);
};

class StringLoc
{
public:
  StringLoc()
  {
  }

  ~StringLoc()
  {
  }

  char ds[MAXSTRINGWIDTH];
  int xloc, yloc, color, olflag, dslen;
};
#endif
