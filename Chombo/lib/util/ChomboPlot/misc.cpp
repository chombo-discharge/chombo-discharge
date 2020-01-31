/* This file contains miscellaneous routines that apply to all kinds
 * of widgets.  Things like setting and getting the color of a widget,
 * its font, position, etc.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include "libsx.h"
#include "libsx_private.h"

extern WindowState *lsx_curwin;  /* global handle to the current window */

/*
 * Miscellaneous functions that allow customization of widgets and
 * other sundry things.
 */

void SetFgColor(Widget w, int color)
{
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  if ((di=libsx_find_draw_info(w)) != NULL)
   {
     Display *d=XtDisplay(w);

     di->foreground = color;

     if (di->mask != 0xffffffff)
       XSetPlaneMask(d, di->drawgc,  di->foreground ^ di->background);
     else
       XSetForeground(d, di->drawgc, color);

     return;
   }

  n = 0;
  XtSetArg(wargs[n], XtNforeground, color);          n++;

  XtSetValues(w, wargs, n);
}

void SetBgColor(Widget w, int color)
{
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */
  Widget tmp;
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL || w == NULL)
     return;

  if ((di=libsx_find_draw_info(w)) != NULL)
   {
     Display *d=XtDisplay(w);

     XSetBackground(d, di->drawgc, color);
     XSetWindowBackground(d, XtWindow(w), color);
     di->background = color;

     if (di->mask != 0xffffffff)
       XSetPlaneMask(d, di->drawgc, di->foreground ^ di->background);

     return;
   }

  tmp = XtParent(w);
  if (tmp != lsx_curwin->form_widget)
   {
     if (XtNameToWidget(tmp, "menu_item"))
       w = tmp;
   }

  n = 0;
  XtSetArg(wargs[n], XtNbackground, color);          n++;

  XtSetValues(w, wargs, n);
}

int GetFgColor(Widget w)
{
  int    color;
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return -1;

  if ((di=libsx_find_draw_info(w)) != NULL)
   {
     return di->foreground;
   }

  n = 0;
  XtSetArg(wargs[n], XtNforeground, &color);         n++;

  XtGetValues(w, wargs, n);

  return color;
}

int GetBgColor(Widget w)
{
  int    color;
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */
  Widget tmp;
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return -1;

  if ((di=libsx_find_draw_info(w)) != NULL)
   {
     return di->background;
   }

  tmp = XtParent(w);
  if (tmp != lsx_curwin->form_widget)
   {
     if (XtNameToWidget(tmp, "menu_item"))
       w = tmp;
   }

  n = 0;
  XtSetArg(wargs[n], XtNbackground, &color);         n++;

  XtGetValues(w, wargs, n);

  return color;
}

void SetBorderColor(Widget w, int color)
{
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  n = 0;
  XtSetArg(wargs[n], XtNborder, color);              n++;

  XtSetValues(w, wargs, n);
}

void SetWidgetState(Widget w, int state)
{
  int    n = 0;
  Arg    wargs[1];              /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  n = 0;
  XtSetArg(wargs[n], XtNsensitive, state);           n++;

  XtSetValues(w, wargs, n);
}

int GetWidgetState(Widget w)
{
  int    n = 0, state;
  Arg    wargs[1];              /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return 0;

  n = 0;
  XtSetArg(wargs[n], XtNsensitive, &state);          n++;

  XtGetValues(w, wargs, n);

  return state;
}

void SetWidgetBitmap(Widget w, char *data, int width, int height)
{
  Pixmap pm;
  Display *d;
  Arg wargs[3];
  int  n=0;

  if (lsx_curwin->display == NULL || w == NULL)
    return;

  d = XtDisplay(w);

  pm = XCreateBitmapFromData(d, DefaultRootWindow(d), data, width, height);
  if (pm == None)
    return;

  n=0;
  XtSetArg(wargs[n], XtNbitmap, pm);    n++;
  XtSetValues(w, wargs, n);
}

void AttachEdge(Widget w, int edge, int attach_to)
{
  char *edge_name;
  static char *edges[]    =
  {
    XtNleft,
    XtNright,
    XtNtop,
    XtNbottom
  };
  static int   attached[] =
  {
    XtChainLeft,
    XtChainRight,
    XtChainTop,
    XtChainBottom
  };
  int   a;
  Arg wargs[5];
  int n=0;

  if (w == NULL || edge < 0 || edge > BOTTOM_EDGE
      || attach_to < 0 || attach_to > ATTACH_BOTTOM)
    return;

  edge_name = edges[edge];
  a         = attached[attach_to];

  n=0;
  XtSetArg(wargs[n], edge_name, a);     n++;

  XtSetValues(w, wargs, n);
}

void SetWidgetPos(Widget w, int where1, Widget from1, int where2, Widget from2)
{
  int  n = 0;
  Arg  wargs[5];                /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  /*
   * Don't want to do this for menu item widgets
   */
  if (XtName(w) && strcmp(XtName(w), "menu_item") == 0)
    return;

  /*
   * This if statement handles the case that the widget we were passed
   * was a List widget.  The reason we use the parent of the widget
   * instead of the widget we were given is because we really want to
   * set the position of the viewport widget that is the parent of the
   * List widget (because when we create a List widget, its parent is
   * a Viewport widget, not the lsx_curwin->form_widget like everyone else.
   *
   * The extra check for the name of the widget not being "form" is to
   * allow proper setting of multiple form widgets.  In the case that
   * we are setting multiple form widgets, the parent of the widget
   * we are setting (and those it is relative to) will not be
   * lsx_curwin->form_widget.  When this is the case, we just want to
   * set the widget itself, not the parent.  Basically this is an
   * exception to the previous paragraph (and it should be the only
   * one, I think).
   *
   * Kind of hackish.  Oh well...
   */
  if (XtParent(w) != lsx_curwin->form_widget && strcmp(XtName(w), "form") != 0)
    w = XtParent(w);

  /*
   * We also change the from1 and from2 fields just as above just so
   * that positioning relative to a list widget works correctly.
   *
   * If we just used the list widget, we'd use its full size for positioning
   * even though the size of the viewport widget is all that's really visible.
   */
  if (from1 && XtParent(from1) != lsx_curwin->form_widget
      && strcmp(XtName(from1), "form") != 0)
    from1 = XtParent(from1);

  if (from2 && XtParent(from2) != lsx_curwin->form_widget
      && strcmp(XtName(from2), "form") != 0)
    from2 = XtParent(from2);

  /*
   * Here we do the first half of the positioning.
   */
  if (where1 == PLACE_RIGHT && from1)
   {
     XtSetArg(wargs[n], XtNfromHoriz, from1);              n++;
   }
  else if (where1 == PLACE_UNDER && from1)
   {
     XtSetArg(wargs[n], XtNfromVert,  from1);              n++;
   }

  /*
   * Now do the second half of the positioning
   */
  if (where2 == PLACE_RIGHT && from2)
   {
     XtSetArg(wargs[n], XtNfromHoriz, from2);              n++;
   }
  else if (where2 == PLACE_UNDER && from2)
   {
     XtSetArg(wargs[n], XtNfromVert,  from2);              n++;
   }

  if (n)                      /* if any values were actually set */
    XtSetValues(w, wargs, n);
}
