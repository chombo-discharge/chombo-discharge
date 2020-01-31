/* DrawingArea Private header file */

/* Copyright 1990, David Nedde
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without fee
 * is granted provided that the above copyright notice appears in all copies.
 * It is provided "as is" without express or implied warranty.
 */

#ifndef _DRAWINGAP_H_
#define _DRAWINGAP_H_

#include "DrawingA.h"
#ifdef X11_R3
#include <X11/SimpleP.h>
#else
#include <X11/Xaw/SimpleP.h>
#endif

/* The drawing area's contribution to the class record */
typedef struct _DrawingAreaClassPart
{
  int ignore;
} DrawingAreaClassPart;

/* Drawing area's full class record */
typedef struct _DrawingAreaClassRec
{
    CoreClassPart       core_class;
    SimpleClassPart     simple_class;
    DrawingAreaClassPart drawing_area;
} DrawingAreaClassRec;

extern DrawingAreaClassRec drawingAreaClassRec;

/* Resources added and status of drawing area widget */
typedef struct _XsDrawingAreaPart
{
  /* Resources */
  XtCallbackList        expose_callback;
  XtCallbackList        input_callback;
  XtCallbackList        motion_callback;
  XtCallbackList        resize_callback;
} DrawingAreaPart;

/* Drawing area's instance record */
typedef struct _DrawingAreaRec
{
    CorePart         core;
    SimplePart       simple;
    DrawingAreaPart  drawing_area;
} DrawingAreaRec;

#endif /* _DRAWINGAP_H_ */
