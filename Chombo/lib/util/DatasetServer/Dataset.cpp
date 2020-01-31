#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// -------------------------------------------------------------------
//  Dataset.cpp
// -------------------------------------------------------------------
#include "ArrayViewData.H"
#include "Dataset.H"
#include "Misc.H"

#include <iostream>
using std::ends;

#include <strstream>
using std::ostrstream;

#include <csignal>
#include <cstdlib>
#include <cfloat>

extern "C"
{
#include <sys/socket.h>
#include <sys/errno.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>
}

// i dunno if this is still necessary...
// #include <cstring>
// #ifdef CH_Solaris
// #  include <strings.h>
// #endif

#define SHOWVAL(val)                      \
  {                                       \
    cout << #val << " = " << val << endl; \
  }

#ifdef CH_USE_FLOAT
#  ifdef CH_Solaris
#    define DS_FP_CLASS fpclass
#  elif (defined(CH_IRIX) || defined(CH_IRIX64))
#    define DS_FP_CLASS fp_class_f
#  else
#    define DS_FP_CLASS fp_classf
#  endif
#  define DS_REAL_MIN FLT_MIN
#  define DS_REAL_MAX FLT_MAX
#else
#  ifdef CH_Solaris
#    define DS_FP_CLASS fpclass
#  elif (defined(CH_IRIX) || defined(CH_IRIX64))
#    define DS_FP_CLASS fp_class_d
#  else
#    define DS_FP_CLASS fp_class
#  endif
#  define DS_REAL_MIN DBL_MIN
#  define DS_REAL_MAX DBL_MAX
#endif

//const int CHARACTERWIDTH  = 14;  //NOTE: not used
const int DATAITEMHEIGHT  = 22;
const int MAXINDEXCHARS   = 4;
const int WOFFSET         = 4;

static char *sliceLabels[] =
{
  "Z Slice:", "Y Slice:", "X Slice:"
};
static char componentLabel[] = "nVar: ";
static char elementLabel[] = "Fab: ";
static char irregLabel[] = "Node: ";
#if (CH_SPACEDIM == 3)
static char *planeLabels[] =
{
  "XY", "XZ", "YZ"
};
#endif

bool Dataset::topLevelCreated = false;

// -------------------------------------------------------------------
// -------------------------------------------------------------------
Dataset::Dataset()
{
  cerr << "_in Dataset()" << endl;
  exit(-1);
}

// -------------------------------------------------------------------
Dataset::Dataset(Widget parent,
                 BaseFab<Real>* fab,
                 const PP_String& formatstring,
                 const PP_String& datalabel,
                 bool is_irregular)
{
  isMulti = false;
  currentElement = 0;
  originalMulti = NULL;
  originalFab = fab;
  currentFabComponent.resize(1);
  currentSlice.resize(1);
  currentSlice[0].resize(CH_SPACEDIM);
  currentPlane.resize(1);
  containsBadFloats.resize(1);
  isIrregular = is_irregular;
  currentIrreg = 0;
  maxIrreg = 1;
  maxComp = fab->nComp();
  Init(parent, formatstring, datalabel);
  Render(originalFab);
}

// -------------------------------------------------------------------

Dataset::Dataset(Widget parent,
                 LayoutData<BaseFab<Real> >* layoutdata,
                 const PP_String& formatstring,
                 const PP_String& datalabel)
{
  int i;
  isMulti = true;
  currentElement = 0;
  originalMulti = layoutdata;
  const int num_elem = layoutdata->boxLayout().size();
  ArrayViewData data(layoutdata);
  originalFab = &data[currentElement];  // start with the first element
  currentFabComponent.resize(num_elem);
  currentSlice.resize(num_elem);
  for (i = 0; i < currentSlice.length(); i++)
  {
    currentSlice[i].resize(CH_SPACEDIM);
#   if (CH_SPACEDIM == 2)
    sDir = DS_XDIR;  // arbitrarily
    currentSlice[i][DS_XDIR] = 0;
    currentSlice[i][DS_YDIR] = 0;
#   else
    currentSlice[i][DS_XDIR] = data[i].box().smallEnd(DS_XDIR);
    currentSlice[i][DS_YDIR] = data[i].box().smallEnd(DS_YDIR);
    currentSlice[i][DS_ZDIR] = data[i].box().smallEnd(DS_ZDIR);
#   endif

  }
  currentPlane.resize(num_elem);
  containsBadFloats.resize(num_elem);
  isIrregular = false;
  currentIrreg = 0;
  maxIrreg = 1;
  maxComp = originalFab->nComp();

  Init(parent, formatstring, datalabel);
  Render(originalFab);
}

// -------------------------------------------------------------------
void
Dataset::Init(Widget parent,
              const PP_String& formatstring,
              const PP_String& datalabel)
{

  // X window initialization
  display = XtDisplay(parent);
  root = RootWindow(display, DefaultScreen(display));
  screen = XtScreen(parent);
  screenNumber = DefaultScreen(display);
  visual = visualInfo.visual;
  XFontStruct* font = XLoadQueryFont(display, "10x20");
  if (font != NULL)
  {
    dataItemHeight = font->ascent + font->descent + 7;
    characterWidth = font->max_bounds.width + 4;
    XGCValues values;
    values.foreground = WhitePixel(display,screenNumber);
    values.background = BlackPixel(display,screenNumber);
    values.font = font->fid;
    gc = XCreateGC(display, RootWindow(display,screenNumber),
                   (GCForeground|GCBackground|GCFont), &values);
  }
  else
  {
    dataItemHeight = 22;
    characterWidth = 13;
    gc = screen->default_gc;
  }
  blackIndex = BlackPixel(display, screenNumber);
  whiteIndex = WhitePixel(display, screenNumber);

  indexWidth  = MAXINDEXCHARS * characterWidth;
  indexHeight = DATAITEMHEIGHT + 7;

  int i;
  stringOk = true;
  showColor = false;
  wParent = parent;
  dataLabel = datalabel;
  formatString = formatstring;
  hStringOffset = 12;
  vStringOffset = -4;
  char *format, tempFormatString[32];
  strcpy(tempFormatString, formatString.c_str());
  XmString sFormatString = XmStringCreateSimple(tempFormatString);
  XmStringGetLtoR(sFormatString, XmSTRING_DEFAULT_CHARSET, &format);
  XmStringFree(sFormatString);

  hDir = DS_XDIR;
  vDir = DS_YDIR;
  sDir = DS_ZDIR;
  hAxisString = "x";
  vAxisString = "y";
  currentFabComponent[currentElement] = 0;

  dataMin.resize(originalFab->nComp());
  dataMax.resize(originalFab->nComp());
  minMaxCalculated.resize(originalFab->nComp(), false);
  containsBadFloats.resize(originalFab->nComp(), true);

  CalculateMinMax(originalFab->dataPtr(currentFabComponent[currentElement]),
                  originalFab->dataPtr(originalFab->nComp()-1),
                  originalFab->box().numPts(),
                  dataMin[currentFabComponent[currentElement]],
                  dataMax[currentFabComponent[currentElement]],
                  containsBadFloats[currentFabComponent[currentElement]]);

  minMaxCalculated[currentFabComponent[currentElement]] = true;

#if (CH_SPACEDIM == 2)
  sDir = DS_XDIR;  // arbitrarily
  currentSlice[currentElement][DS_XDIR] = 0;
  currentSlice[currentElement][DS_YDIR] = 0;
#else
  currentSlice[currentElement][DS_XDIR] = originalFab->box().smallEnd(DS_XDIR);
  currentSlice[currentElement][DS_YDIR] = originalFab->box().smallEnd(DS_YDIR);
  currentSlice[currentElement][DS_ZDIR] = originalFab->box().smallEnd(DS_ZDIR);
#endif

  pixSizeX = 1;
  pixSizeY = 1;

  // ************************************************ Dataset Window
  wDatasetTopLevel = XtVaCreatePopupShell(dataLabel.c_str(),
                                          topLevelShellWidgetClass,
                                          wParent,
                                          XmNwidth,     850,
                                          XmNheight,    500,
                                          NULL);

  wDatasetForm = XtVaCreateManagedWidget("datasetform",
                                         xmFormWidgetClass, wDatasetTopLevel,
                                         NULL);

  wDatasetTools = XtVaCreateManagedWidget("datasettools",
                                          xmFormWidgetClass, wDatasetForm,
                                          XmNtopAttachment,       XmATTACH_FORM,
                                          XmNleftAttachment,      XmATTACH_FORM,
                                          XmNrightAttachment,     XmATTACH_FORM,
#if (CH_SPACEDIM == 2)
                                          XmNheight,              40,
#else
                                          XmNheight,              80,
#endif
                                          NULL);

  // ************************************************ Format Text Field
  i=0;
  XtSetArg(args[i], XmNtopAttachment, XmATTACH_FORM);      i++;
  XtSetArg(args[i], XmNtopOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNleftAttachment, XmATTACH_FORM);      i++;
  XtSetArg(args[i], XmNleftOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNvalue, formatString.c_str());      i++;
  XtSetArg(args[i], XmNcolumns, 7);      i++;
  wFormat = XtCreateManagedWidget("format", xmTextFieldWidgetClass,
                                  wDatasetTools, args, i);
  XtAddCallback(wFormat, XmNactivateCallback, &Dataset::CBReadString,
                (XtPointer) this);
  Dimension bHeight;
  XtVaGetValues(wFormat, XmNheight, &bHeight, NULL);

  // ************************************************ Component Label
  XmString sComponentLabelString = XmStringCreateSimple(componentLabel);
  wComponentLabel = XtVaCreateManagedWidget("Component", xmLabelWidgetClass,
                                            wDatasetTools,
                                            XmNlabelString, sComponentLabelString,
                                            XmNtopAttachment, XmATTACH_WIDGET,
                                            XmNtopWidget, wFormat,
                                            XmNtopOffset, WOFFSET + 18,
                                            XmNleftAttachment, XmATTACH_FORM,
                                            XmNleftOffset, WOFFSET,
                                            NULL);
  XmStringFree(sComponentLabelString);

  // ************************************************ Component Scale
  wComponentScale = XtVaCreateManagedWidget("Component", xmScaleWidgetClass,
                                            wDatasetTools,
                                            //XtVaTypedArg, XmNtitleString, XmRString, "Component", 10,
#if 0
                                            XmNminimum, 0,
                                            XmNmaximum, 1,
                                            XmNvalue, 0,
                                            XmNscaleMultiple, 1,
#endif
                                            XmNshowValue, true,
                                            XmNwidth, 750,
                                            XmNorientation, XmHORIZONTAL,
                                            XmNtopAttachment, XmATTACH_WIDGET,
                                            XmNtopWidget, wFormat,
                                            XmNtopOffset, WOFFSET,
                                            XmNleftAttachment, XmATTACH_WIDGET,
                                            XmNleftWidget, wComponentLabel,
                                            NULL);

  XtAddCallback(wComponentScale, XmNvalueChangedCallback,
                &Dataset::CBChangeComponent, (XtPointer) this);

  // ************************************************ Element Label
  XmString sElementLabelString = XmStringCreateSimple(elementLabel);
  wElementLabel = XtVaCreateManagedWidget("Element", xmLabelWidgetClass,
                                          wDatasetTools,
                                          XmNlabelString, sElementLabelString,
                                          XmNtopAttachment, XmATTACH_WIDGET,
                                          XmNtopWidget, wComponentLabel,
                                          XmNtopOffset, WOFFSET + 14,
                                          XmNleftAttachment, XmATTACH_FORM,
                                          XmNleftOffset, WOFFSET,
                                          NULL);
  XmStringFree(sElementLabelString);

  // ************************************************ Element Scale
  wElementScale = XtVaCreateManagedWidget("Element", xmScaleWidgetClass,
                                          wDatasetTools,
#if 0
                                          XmNminimum, 0,
                                          XmNmaximum, 1,
                                          XmNvalue, 0,
                                          XmNscaleMultiple, 1,
#endif
                                          XmNshowValue, true,
                                          XmNwidth, 750,
                                          XmNorientation, XmHORIZONTAL,
                                          XmNtopAttachment, XmATTACH_WIDGET,
                                          XmNtopWidget, wComponentLabel,
                                          XmNleftAttachment, XmATTACH_WIDGET,
                                          XmNleftWidget, wElementLabel,
                                          XmNleftOffset, WOFFSET + 2,
                                          NULL);

  XtAddCallback(wElementScale, XmNvalueChangedCallback,
                &Dataset::CBChangeElement, (XtPointer) this);

  // ************************************************ Irregular Node Label
  XmString sIrregLabelString = XmStringCreateSimple(irregLabel);
  wIrregLabel = XtVaCreateManagedWidget("Irreg", xmLabelWidgetClass,
                                        wDatasetTools,
                                        XmNlabelString, sIrregLabelString,
                                        XmNtopAttachment, XmATTACH_FORM,
                                        XmNtopOffset, WOFFSET + 14,
                                        XmNleftAttachment, XmATTACH_WIDGET,
                                        XmNleftWidget, wFormat,
                                        XmNleftOffset, WOFFSET,
                                        NULL);
  XmStringFree(sIrregLabelString);

  // ************************************************ Irregular Node Scale
  wIrregScale = XtVaCreateManagedWidget("Irreg", xmScaleWidgetClass,
                                        wDatasetTools,
#if 0
                                        XmNminimum, 0,
                                        XmNmaximum, 1,
                                        XmNvalue, 0,
                                        XmNscaleMultiple, 1,
#endif
                                        XmNshowValue, true,
                                        XmNwidth, 400,
                                        XmNorientation, XmHORIZONTAL,
                                        XmNtopAttachment, XmATTACH_FORM,
                                        XmNleftAttachment, XmATTACH_WIDGET,
                                        XmNleftWidget, wIrregLabel,
                                        NULL);

  XtAddCallback(wIrregScale, XmNvalueChangedCallback,
                &Dataset::CBChangeIrreg, (XtPointer) this);

  // ************************************************ Color Button
  i=0;
  XtSetArg(args[i], XmNtopAttachment, XmATTACH_FORM);      i++;
  XtSetArg(args[i], XmNtopOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNleftAttachment, XmATTACH_WIDGET);      i++;
  XtSetArg(args[i], XmNleftWidget, wFormat);      i++;
  XtSetArg(args[i], XmNleftOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNheight, bHeight);      i++;
  wColorButton = XmCreateToggleButton(wDatasetTools, "Color", args, i);
  ////XtAddCallback(wColorButton, XmNactivateCallback, &Dataset::CBColorButton,
  ////(XtPointer) this);

#if (CH_SPACEDIM == 3)
  // ************************************************ Slice Label
  XmString sSliceLabelString = XmStringCreateSimple(sliceLabels[DS_XYPLANE]);
  wSliceLabel = XtVaCreateManagedWidget("Slice", xmLabelWidgetClass,
                                        wDatasetTools,
                                        XmNlabelString, sSliceLabelString,
                                        XmNtopAttachment, XmATTACH_WIDGET,
                                        XmNtopWidget, wElementLabel,
                                        XmNtopOffset, WOFFSET + 14,
                                        XmNleftAttachment, XmATTACH_FORM,
                                        XmNleftOffset, WOFFSET,
                                        NULL);
  XmStringFree(sSliceLabelString);

  // ************************************************ Slice Scale
  wSliceScale = XtVaCreateManagedWidget("Slice", xmScaleWidgetClass,
                                        wDatasetTools,
                                        //XtVaTypedArg, XmNtitleString, XmRString, "X Slice", 8,
#if 0
                                        XmNscaleMultiple, 1,
#endif
                                        XmNshowValue, true,
                                        XmNwidth, 600,
                                        XmNorientation, XmHORIZONTAL,
                                        XmNtopAttachment, XmATTACH_WIDGET,
                                        XmNtopWidget, wElementLabel,
                                        XmNleftAttachment, XmATTACH_WIDGET,
                                        XmNleftWidget, wSliceLabel,
                                        XmNleftOffset, WOFFSET,
                                        NULL);
  XtAddCallback(wSliceScale, XmNvalueChangedCallback, &Dataset::CBChangeSlice,
                (XtPointer) this);

  // ************************************************ Plane Radio Buttons
  i=0;
  //XtSetArg(args[i], XmNbottomAttachment, XmATTACH_FORM);  i++;
  XtSetArg(args[i], XmNtopAttachment, XmATTACH_WIDGET);   i++;
  XtSetArg(args[i], XmNtopWidget, wElementLabel);         i++;
  XtSetArg(args[i], XmNtopOffset, 3 * WOFFSET);           i++;
  XtSetArg(args[i], XmNleftAttachment, XmATTACH_WIDGET);  i++;
  XtSetArg(args[i], XmNleftWidget, wSliceScale);          i++;
  XtSetArg(args[i], XmNleftOffset, WOFFSET + 14);         i++;
  XtSetArg(args[i], XmNspacing, 0);                       i++;
  XtSetArg(args[i], XmNnumColumns, 3);                    i++;
  wPlaneRadioBox = XmCreateRadioBox(wDatasetTools, "Plane", args, i);
  wXYPlane = XtVaCreateManagedWidget(planeLabels[DS_XYPLANE],
                                     xmToggleButtonGadgetClass,
                                     wPlaneRadioBox,
                                     XmNset, true,
                                     NULL);
  XtAddCallback(wXYPlane, XmNvalueChangedCallback, &Dataset::CBChangePlane,
                (XtPointer) this);
  wXZPlane = XtVaCreateManagedWidget(planeLabels[DS_XZPLANE],
                                     xmToggleButtonGadgetClass,
                                     wPlaneRadioBox, NULL);
  XtAddCallback(wXZPlane, XmNvalueChangedCallback, &Dataset::CBChangePlane,
                (XtPointer) this);
  wYZPlane = XtVaCreateManagedWidget(planeLabels[DS_YZPLANE],
                                     xmToggleButtonGadgetClass,
                                     wPlaneRadioBox, NULL);
  XtAddCallback(wYZPlane, XmNvalueChangedCallback, &Dataset::CBChangePlane,
                (XtPointer) this);
#endif

  // ************************************************ Close Button
  i=0;
  XtSetArg(args[i], XmNtopAttachment, XmATTACH_FORM);      i++;
  XtSetArg(args[i], XmNtopOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNrightAttachment, XmATTACH_FORM);      i++;
  XtSetArg(args[i], XmNrightOffset, WOFFSET);      i++;
  XtSetArg(args[i], XmNheight, bHeight);      i++;
  wQuitButton = XmCreatePushButton(wDatasetTools, "Close", args, i);
  XtAddCallback(wQuitButton, XmNactivateCallback, &Dataset::CBQuitButton,
                (XtPointer) this);

  // ************************************************ Global Max Label
  char globalMaxString[MAXSTRINGWIDTH];
  char globalMaxValueString[MAXSTRINGWIDTH];
  sprintf(globalMaxValueString, formatString.c_str(),
          dataMax[currentFabComponent[currentElement]]);
  sprintf(globalMaxString, "Max:  %s", globalMaxValueString);
  XmString sGlobalMaxString = XmStringCreateSimple(globalMaxString);
  wGlobalMax = XtVaCreateManagedWidget("gmax", xmLabelWidgetClass,
                                       wDatasetTools,
                                       XmNlabelString, sGlobalMaxString,
                                       XmNtopAttachment, XmATTACH_FORM,
                                       XmNtopOffset, WOFFSET,
                                       XmNrightAttachment, XmATTACH_WIDGET,
                                       XmNrightWidget, wQuitButton,
                                       XmNrightOffset, 4 * WOFFSET,
                                       NULL);
  XmStringFree(sGlobalMaxString);

  // ************************************************ Global Min Label
  char globalMinString[MAXSTRINGWIDTH];
  char globalMinValueString[MAXSTRINGWIDTH];
  sprintf(globalMinValueString, formatString.c_str(),
          dataMin[currentFabComponent[currentElement]]);
  sprintf(globalMinString, "Min:  %s", globalMinValueString);
  XmString sGlobalMinString = XmStringCreateSimple(globalMinString);
  wGlobalMin = XtVaCreateManagedWidget("gmin", xmLabelWidgetClass,
                                       wDatasetTools,
                                       XmNlabelString, sGlobalMinString,
                                       XmNtopAttachment, XmATTACH_FORM,
                                       XmNtopOffset, WOFFSET,
                                       XmNrightAttachment, XmATTACH_WIDGET,
                                       XmNrightWidget, wGlobalMax,
                                       XmNrightOffset, 2 * WOFFSET,
                                       NULL);
  XmStringFree(sGlobalMinString);

  // ************************************************ Bad Global Float Label
  XmString sGlobalBadFPString = XmStringCreateSimple("*** Bad Float Value ***");
  wGlobalBadFPLabel = XtVaCreateManagedWidget("globalbadfp", xmLabelWidgetClass,
                                              wDatasetTools,
                                              XmNlabelString, sGlobalBadFPString,
                                              XmNtopAttachment, XmATTACH_FORM,
                                              XmNtopOffset, WOFFSET,
                                              XmNrightAttachment, XmATTACH_WIDGET,
                                              XmNrightWidget, wGlobalMin,
                                              XmNrightOffset, 2 * WOFFSET,
                                              NULL);
  XmStringFree(sGlobalBadFPString);

#if (CH_SPACEDIM == 3)
  // ************************************************ Slice Max Label
  char sliceMaxString[MAXSTRINGWIDTH];
  char sliceMaxValueString[MAXSTRINGWIDTH];
  sprintf(sliceMaxValueString, formatString.c_str(),
          dataMax[currentFabComponent[currentElement]]);  // set to global value for now
  sprintf(sliceMaxString, "Max:  %s", sliceMaxValueString);
  XmString sSliceMaxString = XmStringCreateSimple(sliceMaxString);
  wSliceMax = XtVaCreateManagedWidget("smax", xmLabelWidgetClass,
                                      wDatasetTools,
                                      XmNlabelString, sSliceMaxString,
                                      XmNtopAttachment, XmATTACH_WIDGET,
                                      XmNtopWidget, wGlobalMax,
                                      XmNtopOffset, WOFFSET,
                                      XmNrightAttachment, XmATTACH_WIDGET,
                                      XmNrightWidget, wQuitButton,
                                      XmNrightOffset, 4 * WOFFSET,
                                      NULL);
  XmStringFree(sSliceMaxString);

  // ************************************************ Slice Min Label
  char sliceMinString[MAXSTRINGWIDTH];
  char sliceMinValueString[MAXSTRINGWIDTH];
  sprintf(sliceMinValueString, formatString.c_str(),
          dataMin[currentFabComponent[currentElement]]);  // set to global value for now
  sprintf(sliceMinString, "Slice Min:  %s", sliceMinValueString);
  XmString sSliceMinString = XmStringCreateSimple(sliceMinString);
  wSliceMin = XtVaCreateManagedWidget("smin", xmLabelWidgetClass,
                                      wDatasetTools,
                                      XmNlabelString, sSliceMinString,
                                      XmNtopAttachment, XmATTACH_WIDGET,
                                      XmNtopWidget, wGlobalMin,
                                      XmNtopOffset, WOFFSET,
                                      XmNrightAttachment, XmATTACH_WIDGET,
                                      XmNrightWidget, wSliceMax,
                                      XmNrightOffset, 2 * WOFFSET,
                                      NULL);
  XmStringFree(sSliceMinString);

  // ************************************************ Bad Slice Float Label
  XmString sSliceBadFPString = XmStringCreateSimple("*** Bad Float Value ***");
  wSliceBadFPLabel = XtVaCreateManagedWidget("globalbadfp", xmLabelWidgetClass,
                                             wDatasetTools,
                                             XmNlabelString, sSliceBadFPString,
                                             XmNtopAttachment, XmATTACH_WIDGET,
                                             XmNtopWidget, wGlobalMin,
                                             XmNtopOffset, WOFFSET,
                                             XmNrightAttachment, XmATTACH_WIDGET,
                                             XmNrightWidget, wSliceMin,
                                             XmNrightOffset, 2 * WOFFSET,
                                             NULL);
  XmStringFree(sSliceBadFPString);

#endif

  // ************************************************ data area
  wScrollArea = XtVaCreateManagedWidget("scrollArea",
                                        xmScrolledWindowWidgetClass,
                                        wDatasetForm,
                                        XmNtopAttachment,       XmATTACH_WIDGET,
                                        XmNtopWidget,           wDatasetTools,
                                        XmNleftAttachment,      XmATTACH_FORM,
                                        XmNrightAttachment,     XmATTACH_FORM,
                                        XmNbottomAttachment,    XmATTACH_FORM,
                                        XmNscrollingPolicy,     XmAUTOMATIC,
                                        NULL);

  String trans =
    "<Btn1Motion>: DrawingAreaInput() ManagerGadgetButtonMotion() \n\
         <Btn1Down>: DrawingAreaInput() ManagerGadgetButtonMotion() \n\
         <Btn1Up>: DrawingAreaInput() ManagerGadgetButtonMotion()";

  // ************************************************ drawing area
  wPixArea = XtVaCreateManagedWidget("pixArea", xmDrawingAreaWidgetClass,
                                     wScrollArea,
                                     XmNtranslations,   XtParseTranslationTable(trans),
                                     XmNwidth,          pixSizeX,
                                     XmNheight,         pixSizeY,
                                     XmNbackground,             blackIndex,
                                     NULL);
  XtAddCallback(wPixArea, XmNinputCallback, &Dataset::CBPixInput, (XtPointer) this);
  XtVaSetValues(wScrollArea, XmNworkWindow, wPixArea, NULL);

  XtManageChild(wScrollArea);
  XtManageChild(wPixArea);
  //XtManageChild(wColorButton);
#if (CH_SPACEDIM == 3)
  XtManageChild(wPlaneRadioBox);
  XtUnmanageChild(wSliceBadFPLabel);
#endif
  XtUnmanageChild(wGlobalBadFPLabel);
  XtManageChild(wQuitButton);

  XtPopup(wDatasetTopLevel, XtGrabNone);

  dataWindow = XtWindow(wPixArea);

  XtAddEventHandler(wPixArea, ExposureMask, false,
                    &Dataset::CBDoExposeDataset, (XtPointer) this);

  XtAddEventHandler(wScrollArea, StructureNotifyMask, false,
                    &Dataset::CBResize, (XtPointer) this);

  XtVaGetValues(wScrollArea,
                XmNhorizontalScrollBar, &wHScrollBar,
                XmNverticalScrollBar,   &wVScrollBar,
                NULL);

  XtAddCallback(wHScrollBar, XmNdragCallback,
                &Dataset::CBScrolling, (XtPointer) this);
  XtAddCallback(wHScrollBar, XmNvalueChangedCallback,
                &Dataset::CBEndScrolling, (XtPointer) this);
  XtAddCallback(wVScrollBar, XmNdragCallback,
                &Dataset::CBScrolling, (XtPointer) this);
  XtAddCallback(wVScrollBar, XmNvalueChangedCallback,
                &Dataset::CBEndScrolling, (XtPointer) this);

  dragging = false;
  drags = 0;
}  // end Dataset::Init

// -------------------------------------------------------------------
Dataset::~Dataset()
{
#if (CH_SPACEDIM == 3)
  delete sliceFab;
#endif
  delete [] dataStringArray;
  delete [] hIndexArray;
  delete [] vIndexArray;
  XtDestroyWidget(wDatasetTopLevel);
  exit(0);
}

// -------------------------------------------------------------------
bool
Dataset::Render(BaseFab<Real>* fab)
{

  if (isMulti)
  {
    // ---------------------- set the element scale
    const int num_elem = originalMulti->boxLayout().size();
    if (num_elem == 1)
    {
      XtUnmanageChild(wElementScale);
      char tempLabel[32];
      sprintf(tempLabel, "%s  0", elementLabel);
      XmString sElementLabelString = XmStringCreateSimple(tempLabel);
      XtVaSetValues(wElementLabel,
                    XmNlabelString, sElementLabelString,
                    NULL);
      XmStringFree(sElementLabelString);
    }
    else
    {
      XtVaSetValues(wElementScale,
                    XmNminimum, 0,
                    XmNmaximum, num_elem - 1,
                    XmNvalue, currentElement,
                    XmNscaleMultiple, 1,
                    NULL);
    }
  }
  else
  {
    XtUnmanageChild(wElementLabel);
    XtUnmanageChild(wElementScale);
  }

#if (CH_SPACEDIM == 2)
  sliceFab = fab;
#endif

#if (CH_SPACEDIM == 3)
  // ---------------------- set the slice scale
  if (fab->size()[sDir] == 1)
  {
    // a 3d fab one cell thick in the slicedir
    XtUnmanageChild(wSliceScale);
  }
  else
  {
    XtManageChild(wSliceScale);
    XtVaSetValues(wSliceScale,
                  XmNminimum, fab->box().smallEnd(sDir),
                  XmNmaximum, fab->box().bigEnd(sDir),
                  XmNvalue, currentSlice[currentElement][sDir],
                  XmNscaleMultiple, 1,
                  NULL);
  }

  // ---------------------- create a fab for the slice and copy data into it
  sliceBox = fab->box();
  sliceBox.setSmall(sDir, currentSlice[currentElement][sDir]);
  sliceBox.setBig(sDir,   currentSlice[currentElement][sDir]);
  sliceFab = new BaseFab<Real>(sliceBox, fab->nComp());
  sliceFab->copy(*fab, sliceBox);

#endif

// figure out the max number of components and max number of nodes
// for irregular data
  int slice_max_irreg = 1;
  if (isIrregular)
  {
    Real* num_irreg = fab->dataPtr(fab->nComp()-1);
    for (int n = 0; n < fab->box().numPts(); ++n)
    {
      maxIrreg = Max(maxIrreg, int(num_irreg[n]));
    }
    maxComp = (fab->nComp() - 1) / maxIrreg;
    num_irreg = sliceFab->dataPtr(sliceFab->nComp()-1);
    for (int n = 0; n < sliceFab->box().numPts(); ++n)
    {
      slice_max_irreg = Max(slice_max_irreg, int(num_irreg[n]));
    }
  }
  else
  {
    XtUnmanageChild(wIrregLabel);
    XtUnmanageChild(wIrregScale);
  }

  // ---------------------- set the component scale
  if (maxComp == 1)
  {
    XtUnmanageChild(wComponentScale);
    char tempLabel[32];
    sprintf(tempLabel, "%s  0", componentLabel);
    XmString sComponentLabelString = XmStringCreateSimple(tempLabel);
    XtVaSetValues(wComponentLabel,
                  XmNlabelString, sComponentLabelString,
                  NULL);
    XmStringFree(sComponentLabelString);
  }
  else
  {
    XtVaSetValues(wComponentScale,
                  XmNminimum, 0,
                  XmNmaximum, (maxComp - 1),
                  XmNvalue, currentFabComponent[currentElement],
                  XmNscaleMultiple, 1,
                  NULL);
  }

  // ---------------------- set the irregular nodes scale
  if (!isIrregular || (slice_max_irreg == 1))
  {
    XtUnmanageChild(wIrregScale);
    char tempLabel[32];
    sprintf(tempLabel, "%s  0", irregLabel);
    XmString sIrregLabelString = XmStringCreateSimple(tempLabel);
    XtVaSetValues(wIrregLabel,
                  XmNlabelString, sIrregLabelString,
                  NULL);
    XmStringFree(sIrregLabelString);
  }
  else
  {
    XtVaSetValues(wIrregScale,
                  XmNminimum, 0,
                  XmNmaximum, (slice_max_irreg-1),
                  XmNvalue, currentIrreg,
                  XmNscaleMultiple, 1,
                  NULL);
  }

  int comp = maxIrreg*currentFabComponent[currentElement] + currentIrreg;
  Real* slice_data = sliceFab->dataPtr(comp);
  return (Render(slice_data,
                 sliceFab->dataPtr(sliceFab->nComp()-1),
                 sliceFab->box(),
                 dataMin[currentFabComponent[currentElement]],
                 dataMax[currentFabComponent[currentElement]]));
}  // end Render

// -------------------------------------------------------------------
bool
Dataset::Render(Real* data,
                Real* num_irreg,
                const Box& databox,
                Real datamin,
                Real datamax)
{
  int c, d, ddl, ns;
  Box tempBox;
  char title[256];

  dataPoint = data;
  dataBox   = databox;
  dataMin[currentFabComponent[currentElement]] = datamin;
  dataMax[currentFabComponent[currentElement]] = datamax;
  numIrreg = num_irreg;

  ostrstream titlestream(title, sizeof(title));
  titlestream << dataLabel << "   dataptr = " << originalFab->dataPtr()
              << "   box = " << originalFab->box() << ends;
  XtVaSetValues(wDatasetTopLevel, XmNtitle, title, NULL);

  bool badSliceFloats = containsBadFloats[currentFabComponent[currentElement]];
#if (CH_SPACEDIM == 3)
  CalculateMinMax(data,
                  num_irreg,
                  databox.numPts(),
                  sliceMin,
                  sliceMax,
                  badSliceFloats);
  if (badSliceFloats)
  {
    XtManageChild(wSliceBadFPLabel);
  }
  else
  {
    XtUnmanageChild(wSliceBadFPLabel);
  }
#endif
  if (containsBadFloats[currentFabComponent[currentElement]])
  {
    XtManageChild(wGlobalBadFPLabel);
  }
  else
  {
    XtUnmanageChild(wGlobalBadFPLabel);
  }

  // ------------ find largest data width and count # of data strings
  int largestWidth = 0;
  if (badSliceFloats)
  {
    for (d = 0; d < dataBox.size(vDir); d++)
    {
      // dont use numPts here
      ddl = d * dataBox.size(hDir);
      for (c = 0; c < dataBox.size(hDir); c++)
      {
          sprintf(dataString, formatString.c_str(), dataPoint[c+ddl]);
        largestWidth = Max((int)strlen(dataString), largestWidth);
      }
    }
  }
  else
  {
    for (d = 0; d < dataBox.size(vDir); d++)
    {
      // dont use numPts here
      ddl = d * dataBox.size(hDir);
      for (c = 0; c < dataBox.size(hDir); c++)
      {
        sprintf(dataString, formatString.c_str(), dataPoint[c+ddl]);
        largestWidth = Max((int)strlen(dataString), largestWidth);
      }
    }
  }

  // ------------ change the min and max labels
  char globalMaxString[MAXSTRINGWIDTH];
  char globalMaxValueString[MAXSTRINGWIDTH];
  sprintf(globalMaxValueString, formatString.c_str(),
          dataMax[currentFabComponent[currentElement]]);
  sprintf(globalMaxString, "Max:  %s", globalMaxValueString);
  XmString sGlobalMaxString = XmStringCreateSimple(globalMaxString);
  XtVaSetValues(wGlobalMax,
                XmNlabelString, sGlobalMaxString,
                NULL);
  XmStringFree(sGlobalMaxString);

  char globalMinString[MAXSTRINGWIDTH];
  char globalMinValueString[MAXSTRINGWIDTH];
  sprintf(globalMinValueString, formatString.c_str(),
          dataMin[currentFabComponent[currentElement]]);
  sprintf(globalMinString, "Min:  %s", globalMinValueString);
  XmString sGlobalMinString = XmStringCreateSimple(globalMinString);
  XtVaSetValues(wGlobalMin,
                XmNlabelString, sGlobalMinString,
                NULL);
  XmStringFree(sGlobalMinString);

#if (CH_SPACEDIM == 3)
  // ------------ change the slice min and max labels
  char sliceMaxString[MAXSTRINGWIDTH];
  char sliceMaxValueString[MAXSTRINGWIDTH];
  sprintf(sliceMaxValueString, formatString.c_str(), sliceMax);
  sprintf(sliceMaxString, "Max:  %s", sliceMaxValueString);
  XmString sSliceMaxString = XmStringCreateSimple(sliceMaxString);
  XtVaSetValues(wSliceMax,
                XmNlabelString, sSliceMaxString,
                NULL);
  XmStringFree(sSliceMaxString);

  char sliceMinString[MAXSTRINGWIDTH];
  char sliceMinValueString[MAXSTRINGWIDTH];
  sprintf(sliceMinValueString, formatString.c_str(), sliceMin);
  sprintf(sliceMinString, "Slice Min:  %s", sliceMinValueString);
  XmString sSliceMinString = XmStringCreateSimple(sliceMinString);
  XtVaSetValues(wSliceMin,
                XmNlabelString, sSliceMinString,
                NULL);
  XmStringFree(sSliceMinString);
#endif

  dataItemWidth  = largestWidth * characterWidth;
//  dataItemHeight = DATAITEMHEIGHT;

  int stringCount = dataBox.size(vDir) * dataBox.size(hDir);  // dont use numpts
  numStrings = stringCount;
  dataStringArray = new StringLoc[numStrings];
  if (dataStringArray == NULL)
  {
    cout << "Error in Dataset::Render:  out of memory" << endl;
    return false;
  }
  for (ns = 0; ns < numStrings; ns++)
  {
    dataStringArray[ns].olflag = false;
  }

  // ------------ determine size of data area
  dataAreaWidth  = dataBox.size(hDir) * dataItemWidth;
  dataAreaHeight = dataBox.size(vDir) * dataItemHeight;
  pixSizeX = dataAreaWidth  + indexWidth;
  pixSizeY = dataAreaHeight + indexHeight;

  // create StringLoc array and define color scheme
  XtVaSetValues(wPixArea,
                XmNwidth,       pixSizeX,
                XmNheight,      pixSizeY,
                NULL);

  stringCount = 0;
  int colorSlots = 256;
  int csm1 = colorSlots - 1;
  Real globalDiff = dataMax[currentFabComponent[currentElement]] -
    dataMin[currentFabComponent[currentElement]];
  if (globalDiff == 0.0)
  {
    globalDiff = 1.0;  // so we dont divide by zero
  }
  int paletteStart = 0, paletteEnd = 255;

  tempBox = dataBox;
  tempBox.shift(hDir, -dataBox.smallEnd(hDir));
  tempBox.shift(vDir, -dataBox.smallEnd(vDir));

  // ------------ fill the data string array
  if (isIrregular)
  {
    for (d = 0; d < dataBox.size(vDir); d++)
    {
      ddl = d * dataBox.size(hDir);
      for (c = 0; c < dataBox.size(hDir); c++)
      {
//// cells that do not have irregular data are left blank
          if (numIrreg[c+ddl] > 0.0)
          {
            PP_String local_format_string = "";
            for (int i = 1; i < numIrreg[c+ddl]; ++i)
            {
              local_format_string += "(";
            }
            local_format_string += " ";
            local_format_string += formatString;
            local_format_string += " ";
            for (int i = 1; i < numIrreg[c+ddl]; ++i)
            {
              local_format_string += ")";
            }

            sprintf(dataString, local_format_string.c_str(), dataPoint[c+ddl]);
          }
          else
          {
            sprintf(dataString, "%s", "");
          }
          if (dataPoint[c+ddl] > dataMax[currentFabComponent[currentElement]])
          {
            dataStringArray[stringCount].color = paletteEnd;    // clip
          }
          else if (dataPoint[c+ddl] < dataMin[currentFabComponent[currentElement]])
          {
            dataStringArray[stringCount].color = paletteStart;  // clip
          }
          else
          {
            dataStringArray[stringCount].color = (int) (((dataPoint[c+ddl] -
                                                          dataMin[currentFabComponent[currentElement]]) / globalDiff) * csm1 );
            dataStringArray[stringCount].color += paletteStart;
          }

        // generalize this for flipping axes
        dataStringArray[stringCount].xloc =
          ((tempBox.smallEnd(hDir) + c ) * dataItemWidth) + hStringOffset;

        dataStringArray[stringCount].yloc = (dataAreaHeight - 1) -
          ((tempBox.smallEnd(vDir) + d) * dataItemHeight) + vStringOffset;

        strcpy(dataStringArray[stringCount].ds, dataString);
        dataStringArray[stringCount].dslen = strlen(dataString);

        stringCount++;

      }  // end for (c...)
    }  // end for (d...)
  }
  else
  {
    for (d = 0; d < dataBox.size(vDir); d++)
    {
      ddl = d * dataBox.size(hDir);
      for (c = 0; c < dataBox.size(hDir); c++)
      {

          sprintf(dataString, formatString.c_str(), dataPoint[c+ddl]);
          if (dataPoint[c+ddl] > dataMax[currentFabComponent[currentElement]])
          {
            dataStringArray[stringCount].color = paletteEnd;    // clip
          }
          else if (dataPoint[c+ddl] < dataMin[currentFabComponent[currentElement]])
          {
            dataStringArray[stringCount].color = paletteStart;  // clip
          }
          else
          {
            dataStringArray[stringCount].color = (int) (((dataPoint[c+ddl] -
                                                          dataMin[currentFabComponent[currentElement]]) / globalDiff) * csm1 );
            dataStringArray[stringCount].color += paletteStart;
          }

        // generalize this for flipping axes
        dataStringArray[stringCount].xloc =
          ((tempBox.smallEnd(hDir) + c ) * dataItemWidth) + hStringOffset;

        dataStringArray[stringCount].yloc = (dataAreaHeight - 1) -
          ((tempBox.smallEnd(vDir) + d) * dataItemHeight) + vStringOffset;

        strcpy(dataStringArray[stringCount].ds, dataString);
        dataStringArray[stringCount].dslen = strlen(dataString);

        stringCount++;

      }  // end for (c...)
    }  // end for (d...)
  }

  // fill the box index arrays
  hIndexArray = new StringLoc[dataBox.size(hDir)];
  vIndexArray = new StringLoc[dataBox.size(vDir)];
  // horizontal
  for (d = 0; d < dataBox.size(hDir); d++)
  {
    sprintf(dataString, "%d", d + dataBox.smallEnd(hDir));
    hIndexArray[d].color = blackIndex;
    hIndexArray[d].xloc = ((tempBox.smallEnd(hDir) + d) * dataItemWidth) +
      hStringOffset;
    hIndexArray[d].yloc = 0;  // find this dynamically when drawing
    strcpy(hIndexArray[d].ds, dataString);
    hIndexArray[d].dslen = strlen(dataString);
    hIndexArray[d].olflag = false;
  }  // end for (d...)

  // vertical
  for (d = 0; d < dataBox.size(vDir); d++)
  {
    sprintf(dataString, "%d", d + dataBox.smallEnd(vDir));
    vIndexArray[d].color = blackIndex;
    vIndexArray[d].xloc = 0;  // find this dynamically when drawing
    vIndexArray[d].yloc = (dataAreaHeight - 1) -
      ((tempBox.smallEnd(vDir) + d) * dataItemHeight) + vStringOffset;
    strcpy(vIndexArray[d].ds, dataString);
    vIndexArray[d].dslen = strlen(dataString);
    vIndexArray[d].olflag = false;
  }  // end for (d...)

  SetScrollIncrements();

  return true;
}  // end Dataset::Render

// -------------------------------------------------------------------
void
Dataset::DrawGrid(int startX, int startY, int finishX, int finishY,
                  int gridspacingX, int gridspacingY,
                  int foregroundIndex, int backgroundIndex)
{
  int i;

  XSetForeground(display, gc, foregroundIndex);
  XSetBackground(display, gc, backgroundIndex);

  XDrawLine(display, dataWindow, gc, startX+1, startY, startX+1, finishY);
  for (i = startX; i <= finishX; i += gridspacingX)
  {
    XDrawLine(display, dataWindow, gc, i, startY, i, finishY);
  }
  XDrawLine(display, dataWindow, gc, finishX-1, startY, finishX-1, finishY);

  XDrawLine(display, dataWindow, gc, startX, startY+1, finishX, startY+1);
  for (i = startY; i <= finishY; i += gridspacingY)
  {
    XDrawLine(display, dataWindow, gc, startX, i, finishX, i);
  }
  XDrawLine(display, dataWindow, gc, startX, finishY-1, finishX, finishY-1);
}  // end DrawGrid(...)

// -------------------------------------------------------------------
void
Dataset::CBQuitButton(Widget, XtPointer client_data, XtPointer)
{
  Dataset *obj = (Dataset*) client_data;
  obj->DoQuitButton();
}

// -------------------------------------------------------------------
void
Dataset::CBColorButton(Widget w, XtPointer client_data, XtPointer)
{
  Dataset* obj = (Dataset *) client_data;
  obj->DoColorButton(w, XmToggleButtonGetState(w));
}

// -------------------------------------------------------------------
void
Dataset::DoColorButton(Widget w, bool onoff)
{
  showColor = onoff;
  DoExpose(w, false);
}

// -------------------------------------------------------------------
void
Dataset::CBPixInput(Widget, XtPointer client_data, XtPointer call_data)
{
  Dataset* obj = (Dataset *) client_data;
  obj->DoPixInput((XmDrawingAreaCallbackStruct *) call_data);
}

// -------------------------------------------------------------------
void
Dataset::DoPixInput(XmDrawingAreaCallbackStruct* cbs)
{
}

// -------------------------------------------------------------------
void
Dataset::CBReadString(Widget w, XtPointer client_data, XtPointer call_data)
{
  Dataset* obj = (Dataset *) client_data;
  obj->DoReadString(w, (XmSelectionBoxCallbackStruct *) call_data);
}

// -------------------------------------------------------------------
void
Dataset::DoReadString(Widget w, XmSelectionBoxCallbackStruct *)
{
  char tempString[32];
  strcpy(tempString, XmTextFieldGetString(w));
  if (tempString[0] != '%')
  {
    formatString  = "%";
    formatString += tempString;
  }
  else
  {
    formatString = tempString;
  }
  // unexhaustive string check to prevent errors
  stringOk = true;
  for (int i=0; i < formatString.length(); i++)
  {
    if (formatString[i] == 's' || formatString[i] == 'u' ||
       formatString[i] == 'p' || formatString[i] == 'w')
    {
      stringOk = false;
    }
  }
  if (stringOk)
  {
    Render(dataPoint, numIrreg, dataBox,
           dataMin[currentFabComponent[currentElement]],
           dataMax[currentFabComponent[currentElement]]);
  }
  DoExpose(w, false);
}

// -------------------------------------------------------------------
void
Dataset::DoQuitButton()
{
  // call registered quit function
  //pltAppPtr->QuitDataset();
  delete this;
}

// -------------------------------------------------------------------
void
Dataset::DoRaise()
{
  XtPopup(wDatasetTopLevel, XtGrabNone);
  XMapRaised(display, XtWindow(wDatasetTopLevel));
}

// -------------------------------------------------------------------
void
Dataset::DoExpose(Widget w, bool fromExpose)
{
  if (fromExpose && drags)
  {
    drags--;
    return;
  }

  unsigned int stringCount;
  Box tempBox;
  int xloc, yloc;

  int hSliderSize, hIncrement, hPageIncrement;
  int vSliderSize, vIncrement, vPageIncrement;
  XmScrollBarGetValues(wHScrollBar, &hScrollBarPos, &hSliderSize,
                       &hIncrement, &hPageIncrement);
  XmScrollBarGetValues(wVScrollBar, &vScrollBarPos, &vSliderSize,
                       &vIncrement, &vPageIncrement);

  Dimension scrollAreaSpacing;
  XtVaGetValues(wScrollArea,
                XmNwidth, &width,
                XmNheight, &height,
                XmNspacing, &scrollAreaSpacing,
                NULL);

  xh = (int) hScrollBarPos - dataItemWidth;
  yv = (int) vScrollBarPos - dataItemWidth;

  int hScrollBarBuffer = 32;
  int vScrollBarBuffer = 32;
  if ((int)pixSizeY == vSliderSize)
  {
    // the vertical scroll bar is not visible
    vScrollBarBuffer = scrollAreaSpacing;
  }
  if ((int)pixSizeX == hSliderSize)
  {
    // the horizontal scroll bar is not visible
    hScrollBarBuffer = scrollAreaSpacing;
  }

  hIndexAreaHeight = indexHeight;
  hIndexAreaEnd    = Min((int) pixSizeY, vScrollBarPos + height - hScrollBarBuffer);
  hIndexAreaStart  = hIndexAreaEnd + 1 - hIndexAreaHeight;
  vIndexAreaWidth  = indexWidth;
  vIndexAreaEnd    = Min((int) pixSizeX, hScrollBarPos + width - vScrollBarBuffer);
  vIndexAreaStart  = vIndexAreaEnd + 1 - vIndexAreaWidth;

  tempBox = dataBox;
  tempBox.shift(hDir, -dataBox.smallEnd(hDir));
  tempBox.shift(vDir, -dataBox.smallEnd(vDir));

  XClearArea(display, dataWindow, 0, 0, dataAreaWidth, dataAreaHeight, false);

  // the starting location for DrawGrid can be optimized.
  DrawGrid(0, 0, vIndexAreaStart, hIndexAreaStart, dataItemWidth, dataItemHeight,
           whiteIndex, blackIndex);

  DrawIndicies(tempBox);

  if (dragging)
  {
    return;
  }

  // draw the data strings
  if (showColor)
  {
    cerr << "showColor not implemented." << endl;
    exit(-2);
  }
  else
  {
    // ! showColor
    XSetForeground(display, gc, whiteIndex);
    for (stringCount=0; (int)stringCount<numStrings; stringCount++)
    {
      xloc = dataStringArray[stringCount].xloc;
      yloc = dataStringArray[stringCount].yloc;
      if (dataStringArray[stringCount].olflag == false  &&
         (xloc > xh) && (xloc < vIndexAreaStart) &&
         (yloc > yv) && (yloc < hIndexAreaEnd))
      {
        XDrawString(display, dataWindow, gc, xloc, yloc,
                    dataStringArray[stringCount].ds, dataStringArray[stringCount].dslen);
      }
    }  // end for (...)
  }  // end if (showColor)

  DrawIndicies(tempBox);
  SetScrollIncrements();

}  // end DoExpose

// -------------------------------------------------------------------
void
Dataset::DrawIndicies(const Box& tempBox)
{
  int xloc, yloc;
  unsigned int stringCount;

  XSetForeground(display, gc, whiteIndex);
  // horizontal
  XFillRectangle(display, dataWindow, gc, hScrollBarPos,
                 hIndexAreaStart, Min((unsigned int) width,  pixSizeX),
                 hIndexAreaHeight);
  // vertical
  XFillRectangle(display, dataWindow, gc, vIndexAreaStart,
                 vScrollBarPos, vIndexAreaWidth,
                 Min((unsigned int) height, pixSizeY));
  // draw the horizontal box index grid
  DrawGrid(tempBox.smallEnd(hDir) * dataItemWidth, hIndexAreaStart,
           vIndexAreaStart, hIndexAreaEnd,
           dataItemWidth, hIndexAreaHeight, blackIndex, whiteIndex);
  // draw the vertical box index grid
  DrawGrid(vIndexAreaStart,
           dataAreaHeight - (tempBox.bigEnd(vDir)+1) * dataItemHeight,
           vIndexAreaEnd, hIndexAreaStart, vIndexAreaWidth, dataItemHeight,
           blackIndex, whiteIndex);

  XSetForeground(display, gc, blackIndex);

  // draw the corner axis labels
  XDrawLine(display, dataWindow, gc, vIndexAreaStart, hIndexAreaStart,
            vIndexAreaEnd,   hIndexAreaEnd - 1);
  // frame the corner box
  XDrawLine(display, dataWindow, gc, vIndexAreaStart,   hIndexAreaEnd - 1,
            vIndexAreaEnd,     hIndexAreaEnd - 1);
  XDrawLine(display, dataWindow, gc, vIndexAreaEnd - 1, hIndexAreaStart,
            vIndexAreaEnd - 1, hIndexAreaEnd);
  XDrawString(display, dataWindow, gc,
              vIndexAreaStart + hStringOffset,
              hIndexAreaEnd   + vStringOffset,
              hAxisString.c_str(), hAxisString.length());
  XDrawString(display, dataWindow, gc,
              vIndexAreaStart + (indexWidth/2)  + hStringOffset,
              hIndexAreaEnd   - (indexHeight/2) + vStringOffset,
              vAxisString.c_str(), vAxisString.length());

  // draw the box indicies
  // horizontal
  yloc = hIndexAreaEnd + vStringOffset;
  for (stringCount = 0; (int)stringCount < dataBox.size(hDir); stringCount++)
  {
    xloc = hIndexArray[stringCount].xloc;
    if ((xloc > xh) && (xloc < (vIndexAreaStart - (int)(indexWidth / 3))))
    {
      XDrawString(display, dataWindow, gc, xloc, yloc,
                  hIndexArray[stringCount].ds, hIndexArray[stringCount].dslen);
    }
  }  // end for (...)

  // vertical
  xloc = vIndexAreaStart + hStringOffset;
  for (stringCount = 0; (int)stringCount < dataBox.size(vDir); stringCount++)
  {
    yloc = vIndexArray[stringCount].yloc;
    if ((yloc > yv) && (yloc < (hIndexAreaStart)))
    {
      XDrawString(display, dataWindow, gc, xloc, yloc,
                  vIndexArray[stringCount].ds, vIndexArray[stringCount].dslen);
    }
  }  // end for (...)
}  // end DrawIndicies

// -------------------------------------------------------------------
void
Dataset::CBDoExposeDataset(Widget w, XtPointer client_data,
                           XEvent *, Boolean *)
{
  XEvent nextEvent;
  Dataset *dset = (Dataset *) client_data;
  while (XCheckTypedWindowEvent(dset->display,
                               XtWindow(dset->wPixArea), Expose, &nextEvent))
  {
    if (dset->drags)
    {
      dset->drags--;
    }
  }
  dset->DoExpose(w, true);
}

// -------------------------------------------------------------------
void
Dataset::CBResize(Widget w, XtPointer client_data, XEvent *, Boolean *)
{
  Dataset* dset = (Dataset *) client_data;
  dset->DoResize(w);
}

// -------------------------------------------------------------------
void
Dataset::DoResize(Widget w)
{
  DoExpose(w, false);
}

// -------------------------------------------------------------------
void
Dataset::SetScrollIncrements()
{
  Dimension scrollWidth, scrollHeight;
  Dimension hPageIncrement, vPageIncrement;

  XtVaGetValues(wScrollArea,
                XmNwidth, &scrollWidth,
                XmNheight, &scrollHeight,
                NULL);

  hPageIncrement = dataItemWidth *
    ((scrollWidth - dataItemWidth - indexWidth - 32) / dataItemWidth);

  vPageIncrement = dataItemHeight *
    ((scrollHeight - dataItemHeight - indexHeight - 32) / dataItemHeight);

  if (hPageIncrement < 1 || hPageIncrement > 64000)
  {
    hPageIncrement = dataItemWidth;
  }
  if (vPageIncrement < 1 || vPageIncrement > 64000)
  {
    vPageIncrement = dataItemHeight;
  }

  XtVaSetValues(wHScrollBar,
                XmNincrement, dataItemWidth,
                XmNpageIncrement, hPageIncrement,
                NULL);
  XtVaSetValues(wVScrollBar,
                XmNincrement, dataItemHeight,
                XmNpageIncrement, vPageIncrement,
                NULL);
}

// -------------------------------------------------------------------
void
Dataset::CBScrolling(Widget w, XtPointer client_data, XtPointer)
{
  Dataset* dset = (Dataset *) client_data;
  dset->dragging = true;
  dset->drags++;
  dset->DoExpose(w, false);
}

// -------------------------------------------------------------------
void
Dataset::CBEndScrolling(Widget w, XtPointer client_data, XtPointer)
{
  Dataset* dset = (Dataset *) client_data;
  dset->dragging = false;
  dset->DoExpose(w, false);
}

// -------------------------------------------------------------------
void
Dataset::CBChangePlane(Widget w, XtPointer client_data, XtPointer state)
{
  Dataset* dset = (Dataset *) client_data;
  if (((XmToggleButtonCallbackStruct *) state)->set)
  {
    dset->DoChangePlane(w);
  }
}

// -------------------------------------------------------------------
void
Dataset::DoChangePlane(Widget w)
{
  if (w == wXYPlane)
  {
    currentPlane[currentElement] = DS_XYPLANE;
    hAxisString = "x";
    vAxisString = "y";
    hDir = DS_XDIR;
    vDir = DS_YDIR;
    sDir = DS_ZDIR;
  }
  else if (w == wXZPlane)
  {
    currentPlane[currentElement] = DS_XZPLANE;
    hAxisString = "x";
    vAxisString = "z";
    hDir = DS_XDIR;
    vDir = DS_ZDIR;
    sDir = DS_YDIR;
  }
  else if (w == wYZPlane)
  {
    currentPlane[currentElement] = DS_YZPLANE;
    hAxisString = "y";
    vAxisString = "z";
    hDir = DS_YDIR;
    vDir = DS_ZDIR;
    sDir = DS_XDIR;
  }
  else
  {
    cerr << "Error in DoChangePlane:  bad widget = " << w << endl;
    exit(-4);
  }

  char tempLabel[32];
  if (originalFab->size()[sDir] == 1)
  {
    // a 3d fab one cell thick in the slicedir
    sprintf(tempLabel, "%s  %d", sliceLabels[currentPlane[currentElement]],
            originalFab->box().smallEnd(sDir));
  }
  else
  {
    sprintf(tempLabel, "%s", sliceLabels[currentPlane[currentElement]]);
  }
  XmString sSliceLabelString = XmStringCreateSimple(tempLabel);
  XtVaSetValues(wSliceLabel,
                XmNlabelString, sSliceLabelString,
                NULL);
  XmStringFree(sSliceLabelString);

  delete sliceFab;
  Render(originalFab);
  DoExpose(w, false);
}  // end DoChangePlane

// -------------------------------------------------------------------
void
Dataset::CBChangeSlice(Widget w, XtPointer client_data, XtPointer cbs)
{
  Dataset* dset = (Dataset *) client_data;
  dset->DoChangeSlice(((XmScaleCallbackStruct *) cbs)->value);
}

// -------------------------------------------------------------------
void
Dataset::DoChangeSlice(int newslice)
{
  currentSlice[currentElement][sDir] = newslice;
  delete sliceFab;
  Render(originalFab);
  DoExpose(wDatasetTopLevel, false);
}

// -------------------------------------------------------------------
void
Dataset::CBChangeComponent(Widget w, XtPointer client_data, XtPointer cbs)
{
  Dataset* dset = (Dataset *) client_data;
  dset->DoChangeComponent(((XmScaleCallbackStruct *) cbs)->value);
}

// -------------------------------------------------------------------
void
Dataset::DoChangeComponent(int newcomponent)
{
  currentFabComponent[currentElement] = newcomponent;

  if ( ! minMaxCalculated[currentFabComponent[currentElement]])
  {
    CalculateMinMax(originalFab->dataPtr(currentFabComponent[currentElement]),
                    originalFab->dataPtr(originalFab->nComp()-1),
                    originalFab->box().numPts(),
                    dataMin[currentFabComponent[currentElement]],
                    dataMax[currentFabComponent[currentElement]],
                    containsBadFloats[currentFabComponent[currentElement]]);
    minMaxCalculated[currentFabComponent[currentElement]] = true;
  }
  int comp = maxIrreg*currentFabComponent[currentElement] + currentIrreg;
  Real* slice_data = sliceFab->dataPtr(comp);
  Render(slice_data, numIrreg, sliceFab->box(),
         dataMin[currentFabComponent[currentElement]],
         dataMax[currentFabComponent[currentElement]]);
  DoExpose(wDatasetTopLevel, false);
}

// -------------------------------------------------------------------
void
Dataset::CBChangeElement(Widget w, XtPointer client_data, XtPointer cbs)
{
  Dataset* dset = (Dataset *) client_data;
  dset->DoChangeElement(((XmScaleCallbackStruct *) cbs)->value);
}

// -------------------------------------------------------------------
void
Dataset::DoChangeElement(int newelement)
{
  if (newelement == currentElement)
  {
    return;
  }
  ArrayViewData data(originalMulti);
  currentFabComponent[newelement] = currentFabComponent[currentElement];
  currentElement = newelement;
  originalFab = &(data[currentElement]);

  for (int i = 0; i < minMaxCalculated.length(); i++)
  {
    minMaxCalculated[i] = false;
  }

  CalculateMinMax(originalFab->dataPtr(currentFabComponent[currentElement]),
                  originalFab->dataPtr(originalFab->nComp()-1),
                  originalFab->box().numPts(),
                  dataMin[currentFabComponent[currentElement]],
                  dataMax[currentFabComponent[currentElement]],
                  containsBadFloats[currentFabComponent[currentElement]]);
  minMaxCalculated[currentFabComponent[currentElement]] = true;

#if (CH_SPACEDIM == 3)
  delete sliceFab;
#endif
  Render(originalFab);
  DoExpose(wDatasetTopLevel, false);
}

// -------------------------------------------------------------------
void
Dataset::CBChangeIrreg(Widget w, XtPointer client_data, XtPointer cbs)
{
  Dataset *dset = (Dataset *) client_data;
  dset->DoChangeIrreg(((XmScaleCallbackStruct *) cbs)->value);
}

// -------------------------------------------------------------------
void
Dataset::DoChangeIrreg(int newirreg)
{
  if (newirreg == currentIrreg)
  {
    return;
  }
  currentIrreg = newirreg;

  int comp = maxIrreg*currentFabComponent[currentElement] + currentIrreg;
  Real* slice_data = sliceFab->dataPtr(comp);
  Render(slice_data, numIrreg, sliceFab->box(),
         dataMin[currentFabComponent[currentElement]],
         dataMax[currentFabComponent[currentElement]]);
  DoExpose(wDatasetTopLevel, false);
}

// -------------------------------------------------------------------
bool
Dataset::CalculateMinMax(Real* data, Real* num_irreg,
                         int numpoints, Real& fmin, Real& fmax,
                         bool& badfloats)
{
  bool minMaxSet = false;  // was a valid value set for fmin and fmax
  fmin =  DS_REAL_MAX;
  fmax = -DS_REAL_MAX;
  badfloats = false;

  for (int n = 0; n < numpoints; n++)
  {
      if (!isIrregular || (num_irreg[n] > 0.0) )
      {
        fmin = Min(fmin, data[n]);
        fmax = Max(fmax, data[n]);
        minMaxSet = true;
      }
  }  // end for
  return minMaxSet;
}
