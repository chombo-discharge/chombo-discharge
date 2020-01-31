/* DrawingA.c: The DrawingArea Widget Methods */

/* Copyright 1990, David Nedde
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without fee
 * is granted provided that the above copyright notice appears in all copies.
 * It is provided "as is" without express or implied warranty.
 */

#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <X11/CoreP.h>
#include <X11/Xaw/SimpleP.h>
#include "DrawingAP.h"

static void     Initialize(DrawingAreaWidget request,
                           DrawingAreaWidget newYuk);
static void     Destroy(DrawingAreaWidget  w);
static void     Redisplay(DrawingAreaWidget w,XEvent* event,Region region);
static void     input_draw(DrawingAreaWidget w,
                           XEvent *event,
                           char *args[],
                           int  n_args);
static void     motion_draw(
                            DrawingAreaWidget w,
                            XEvent* event,
                            char* args[],
                            int  n_args);
static void     resize_draw(
                            DrawingAreaWidget w,
                            XEvent* event,
                            char* args[],
                            int   n_args);

static char defaultTranslations[] = "<BtnDown>: input() \n <BtnUp>: input() \n <KeyDown>: input() \n <KeyUp>: input() \n <Motion>: motion() \n <Configure>: resize()";
static XtActionsRec actionsList[] =
{
  {
    "input",  (XtActionProc)input_draw
  },
  {
    "motion", (XtActionProc)motion_draw
  },
  {
    "resize", (XtActionProc)resize_draw
  },
};

/* Default instance record values */
static XtResource resources[] =
{
  {
    XtNexposeCallback, XtCCallback, XtRCallback, sizeof(caddr_t),
     XtOffset(DrawingAreaWidget, drawing_area.expose_callback),
     XtRCallback, NULL
  },
  {
    XtNinputCallback, XtCCallback, XtRCallback, sizeof(caddr_t),
     XtOffset(DrawingAreaWidget, drawing_area.input_callback),
     XtRCallback, NULL
  },
  {
    XtNmotionCallback, XtCCallback, XtRCallback, sizeof(caddr_t),
     XtOffset(DrawingAreaWidget, drawing_area.motion_callback),
     XtRCallback, NULL
  },
  {
    XtNresizeCallback, XtCCallback, XtRCallback, sizeof(caddr_t),
     XtOffset(DrawingAreaWidget, drawing_area.resize_callback),
     XtRCallback, NULL
  },
};

DrawingAreaClassRec drawingAreaClassRec =
{
  /* CoreClassPart */
{
   ( WidgetClass    )    &simpleClassRec,       /* superclass             */
   ( String         )    "DrawingArea",                 /* class_name             */
   ( Cardinal       )    sizeof(DrawingAreaRec),                /* size                   */
   ( XtProc         )    NULL,                          /* class_initialize       */
   ( XtWidgetClassProc)  NULL,                          /* class_part_initialize  */
   ( XtEnum         )    FALSE,                         /* class_inited           */
   ( XtInitProc     )    Initialize,                            /* initialize             */
   ( XtArgsProc     )    NULL,                          /* initialize_hook        */
   ( XtRealizeProc  )    XtInheritRealize,                      /* realize                */
   ( XtActionList   )    actionsList,                   /* actions                */
   ( Cardinal       )    XtNumber(actionsList),         /* num_actions            */
   ( XtResourceList )    resources,                             /* resources              */
   ( Cardinal       )    XtNumber(resources),           /* resource_count         */
   ( XrmClass       )    NULLQUARK,                             /* xrm_class              */
   ( Boolean        )    FALSE,                         /* compress_motion        */
   ( XtEnum         )    FALSE,                         /* compress_exposure      */
   ( Boolean        )    TRUE,                          /* compress_enterleave    */
   ( Boolean        )    FALSE,                         /* visible_interest       */
   ( XtWidgetProc   )    Destroy,                               /* destroy                */
   ( XtWidgetProc   )    NULL,                          /* resize                 */
   ( XtExposeProc   )    Redisplay,                             /* expose                 */
   ( XtSetValuesFunc)    NULL,                          /* set_values             */
   ( XtArgsFunc     )    NULL,                          /* set_values_hook        */
   ( XtAlmostProc   )    XtInheritSetValuesAlmost,              /* set_values_almost      */
   ( XtArgsProc     )    NULL,                          /* get_values_hook        */
   ( XtAcceptFocusProc)  NULL,                          /* accept_focus           */
   ( XtVersionType  )    XtVersion,                             /* version                */
   ( XtPointer      )    NULL,                          /* callback_private       */
   ( String         )    defaultTranslations,           /* tm_table               */
   ( XtGeometryHandler)  XtInheritQueryGeometry,                /* query_geometry         */
   ( XtStringProc   )    XtInheritDisplayAccelerator,   /* display_accelerator    */
   ( XtPointer      )    NULL                           /* extension              */
  },  /* CoreClass fields initialization */
  {
    /* change_sensitive         */      XtInheritChangeSensitive
  },  /* SimpleClass fields initialization */
  {
    0,                                     /* field not used    */
  },  /* DrawingAreaClass fields initialization */
};

WidgetClass drawingAreaWidgetClass = (WidgetClass)&drawingAreaClassRec;

static void Initialize(DrawingAreaWidget  request, DrawingAreaWidget newYuk)
{
  if (request->core.width == 0)
    newYuk->core.width = 100;
  if (request->core.height == 0)
    newYuk->core.height = 100;
}

static void Destroy(DrawingAreaWidget  w)
{
  XtRemoveAllCallbacks((Widget)w, XtNexposeCallback);
  XtRemoveAllCallbacks((Widget)w, XtNinputCallback);
  XtRemoveAllCallbacks((Widget)w, XtNmotionCallback);
  XtRemoveAllCallbacks((Widget)w, XtNresizeCallback);
}

/* Invoke expose callbacks */
static void Redisplay(DrawingAreaWidget w,XEvent* event,Region region)
{
  XawDrawingAreaCallbackStruct cb;

  cb.reason = XawCR_EXPOSE;
  cb.event  = event;
  cb.window = XtWindow(w);
  XtCallCallbacks((Widget)w, XtNexposeCallback, (char *)&cb);
}

/* Invoke resize callbacks */
static void resize_draw(
                        DrawingAreaWidget w,
                        XEvent           *event,
                        char             *args[],
                        int               n_args)
{
  XawDrawingAreaCallbackStruct cb;

  cb.reason = XawCR_RESIZE;
  cb.event  = event;
  cb.window = XtWindow(w);
  XtCallCallbacks((Widget)w, XtNresizeCallback, (char *)&cb);
}

/* Invoke input callbacks */
static void input_draw(
                       DrawingAreaWidget w,
                       XEvent            *event,
                       char              *args[],
                       int                n_args)
{
  XawDrawingAreaCallbackStruct cb;

  cb.reason = XawCR_INPUT;
  cb.event  = event;
  cb.window = XtWindow(w);
  XtCallCallbacks((Widget)w, XtNinputCallback, (char *)&cb);
}

/* Invoke motion callbacks */
static void motion_draw(
                        DrawingAreaWidget w,
                        XEvent           *event,
                        char             *args[],
                        int               n_args)
{
  XawDrawingAreaCallbackStruct cb;

  cb.reason = XawCR_MOTION;
  cb.event  = event;
  cb.window = XtWindow(w);
  XtCallCallbacks((Widget)w, XtNmotionCallback, (char *)&cb);
}
