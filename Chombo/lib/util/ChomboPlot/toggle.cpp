/* This file contains routines that manage creating toggle widgets (that
 * is, widgets that have state, either on or off) and radio group widgets
 * (radio groups are a set of mutually exclusive toggles).
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include <X11/Xaw/Toggle.h>
#include "libsx.h"
#include "libsx_private.h"

extern WindowState *lsx_curwin;  /* global handle for the current window. */

/*
 * Toggle Button Creation/Manipulation routines.
 */

Widget MakeToggle(char *txt, int state, Widget w, ButtonCB func, void *d)
{
  int    n = 0;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget toggle;
  char   *name = "toggle";

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  n = 0;
  if (txt)
   {
     XtSetArg(wargs[n], XtNlabel, txt);                   n++;
   }
  XtSetArg(wargs[n], XtNstate, (Boolean)state);           n++;
  if (w)
   {
     char *s = XtName(w);

     /*
      * If the widget we were given is not part of a radio group, the
      * user obviously is screwed up, so just dump out.
      */
     if (s==NULL || (strcmp(s,"toggle") != 0 && strcmp(s,"radio_group") != 0))
       return NULL;

     XtSetArg(wargs[n], XtNradioGroup, w);                n++;
     name = "radio_group";
   }

  toggle = XtCreateManagedWidget(name, toggleWidgetClass,
                                 lsx_curwin->form_widget,wargs,n);
  if (toggle == NULL)
    return NULL;

  if (w)
   {
     void *data;

     n = 0;
     XtSetArg(wargs[n], XtNradioData, toggle);            n++;
     XtSetValues(toggle, wargs, n);

     /*
      * Check if the radioData resource is set for the widget w.  If
      * it isn't, set it here, so it is set for all widgets that
      * are part of the radioGroup.  This only gets executed to
      * patch up the first toggle of a radio group that never had
      * its radioData resource set.
      */
     n=0;
     XtSetArg(wargs[n], XtNradioData,  &data);            n++;
     XtGetValues(w, wargs, n);
     if (data != w) /* if we already set the radioData field, then data == w */
      {
        n=0;
        XtSetArg(wargs[n], XtNradioData, w);              n++;
        XtSetValues(w, wargs, n);
      }
   }

  if (func)
    XtAddCallback(toggle, XtNcallback, (XtCallbackProc)func, d);

  return toggle;
}    /* end of MakeToggle() */

void SetToggleState(Widget w, int state)
{
  int n=0;
  Arg wargs[5];
  void *data;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  n=0;
  XtSetArg(wargs[n], XtNradioData, &data);                n++;
  XtGetValues(w, wargs, n);

  if (data != w)     /* only happens if we're not in a radio group */
   {
     Boolean old_state;

     n=0;
     XtSetArg(wargs[n], XtNstate, &old_state);            n++;
     XtGetValues(w, wargs, n);

     if (old_state == state)   /* no state change */
       return;

     n=0;
     XtSetArg(wargs[n], XtNstate, (Boolean)state);        n++;
     XtSetValues(w, wargs, n);

     /* make sure callbacks get called since we changed state */
     XtCallCallbacks(w, XtNcallback, NULL);
   }
  else    /* widget is part of a radio group, so do things right */
   {
     Widget cur_widg;

     cur_widg = (Widget)XawToggleGetCurrent(w);
     if ((cur_widg != w && state == FALSE) || (cur_widg == w && state == TRUE))
       return;  /* early out because there is no state change */

     if (state == TRUE)
       XawToggleSetCurrent(w, w);
     else if (state == FALSE)      /* current widget is on, turn it off */
       XawToggleUnsetCurrent(w);
   }
}

int GetToggleState(Widget w)
{
  Boolean state=0;
  int n=0;
  Arg wargs[2];

  if (w == NULL)
    return FALSE;

  n=0;
  XtSetArg(wargs[n], XtNstate, &state);      n++;
  XtGetValues(w, wargs, n);

  return (int)state;
}
