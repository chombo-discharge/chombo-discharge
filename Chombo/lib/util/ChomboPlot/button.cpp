/* This file contains routines that manage creating and manipulating
 * button and label widgets.  Button widgets can be clicked on, label
 * widgets just display some text.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include <X11/Xaw/Command.h>
#include "libsx.h"
#include "libsx_private.h"

extern WindowState *lsx_curwin;   /* global handle for the current window */

/*
 * Command Button Creation routine.
 */

Widget MakeButton(char *txt, ButtonCB func, void *data)
{
  int    n = 0;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget button;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  n = 0;
  if (txt)
   {
     XtSetArg(wargs[n], XtNlabel, txt);                   n++;
   }

  button = XtCreateManagedWidget("button", commandWidgetClass,
                                 lsx_curwin->form_widget,wargs,n);
  if (button == NULL)
    return NULL;

  if (func)
    XtAddCallback(button, XtNcallback, (XtCallbackProc)func, data);

  return button;
}    /* end of MakeButton() */

/*
 * Text Label Creation.
 */

Widget MakeLabel(char *txt)
{
  int    n = 0;
  int    bg_color = -1;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget label;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  n = 0;
  if (txt)
   {
     XtSetArg(wargs[n], XtNlabel, txt);               n++;
   }

  label = XtCreateManagedWidget("label", labelWidgetClass,
                                lsx_curwin->form_widget,wargs,n);
  if (label == NULL)
    return NULL;

  /* this little contortion here is to make sure there is no
   * border around the label (else it looks exactly like a command
   * button, and that's confusing)
   */

  n = 0;
  XtSetArg(wargs[n], XtNbackground, &bg_color);       n++;
  XtGetValues(label, wargs, n);

  n = 0;
  XtSetArg(wargs[n], XtNborder, bg_color);            n++;
  XtSetValues(label, wargs, n);

  return label;
}                    /* end of MakeLabel() */

void SetLabel(Widget w, char *txt)
{
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  n = 0;
  XtSetArg(wargs[n], XtNlabel, txt);                 n++;

  XtSetValues(w, wargs, n);
}
