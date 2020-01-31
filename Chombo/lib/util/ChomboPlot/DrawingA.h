/* DrawingA.h - Public Header file */

/* Copyright 1990, David Nedde
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without fee
 * is granted provided that the above copyright notice appears in all copies.
 * It is provided "as is" without express or implied warranty.
 */

/* Define widget's class pointer and strings used to specify resources */

#ifndef _DRAWINGA_H_
#define _DRAWINGA_H_

#define XADCS XawDrawingAreaCallbackStruct

/* Resources ADDED to label widget:

 Name                Class              RepType         Default Value
 ----                -----              -------         -------------
 exposeCallback      Callback           Pointer         NULL
 inputCallback       Callback           Pointer         NULL
 motionCallback      Callback           Pointer         NULL
 resizeCallback      Callback           Pointer         NULL
*/

extern WidgetClass drawingAreaWidgetClass;

typedef struct _DrawingAreaClassRec *DrawingAreaWidgetClass;
typedef struct _DrawingAreaRec      *DrawingAreaWidget;

/* Resource strings */
#define XtNexposeCallback       "exposeCallback"
#define XtNinputCallback        "inputCallback"
#define XtNmotionCallback       "motionCallback"
#define XtNresizeCallback       "resizeCallback"

typedef struct _XawDrawingAreaCallbackStruct
{
  int     reason;
  XEvent *event;
  Window  window;
} XawDrawingAreaCallbackStruct;

/* Reasons */
#define XawCR_EXPOSE 1
#define XawCR_INPUT  2
#define XawCR_MOTION 3
#define XawCR_RESIZE 4

#endif /* _DRAWINGAREA_H_ */
