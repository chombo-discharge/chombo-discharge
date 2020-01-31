/* This file contains routines that create menus and menu items within
 * those menus.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include <X11/Xaw/MenuButton.h>
#include <X11/Xaw/SimpleMenu.h>
#include <X11/Xaw/SmeBSB.h>
#include "libsx.h"
#include "libsx_private.h"
#include "check_mark.h"

extern WindowState *lsx_curwin;  /* global handle to the current window */

/*
 * Menu functions.
 */
Widget MakeMenu(char *name)
{
  int    n = 0;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget button, menu=NULL;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0,NULL) == 0)
    return NULL;
  if (name == NULL)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNlabel, name);                     n++;

  /*
   * We create two widgets, the simpleMenu widget is a child of the menu button
   * widget.  We return the reference to the button widget however so that
   * way the SetXXColor() and SetWidgetFont() functions work as expect
   * (that is they change the menu button).
   *
   * The MakeMenuItem() function is aware of this, and gets the child of
   * the menu button for creation of the actual menu items.
   *
   */
  button = XtCreateManagedWidget("menuButton", menuButtonWidgetClass,
                                 lsx_curwin->form_widget,  wargs, n);

  if (button)
    menu = XtCreatePopupShell("menu", simpleMenuWidgetClass, button,
                              NULL, 0);

  if (menu == NULL)
   {
     XtDestroyWidget(button);
     button = NULL;
   }

  return button;
}

Widget MakeMenuItem(Widget parent, char *name, ButtonCB func, void *arg)
{
  int    n = 0;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget item, menu;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0,NULL) == 0)
    return NULL;
  if (parent == NULL)
    return NULL;

  /*
   * We get the "menu" widget child of the parent widget because the
   * parent is really a reference to the menu button widget, not the
   * popup-shell widget we really want.  See the above comment in
   * MakeMenu().
   */
  if ((menu = XtNameToWidget(parent, "menu")) == NULL)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNlabel,      name);                n++;
  XtSetArg(wargs[n], XtNleftMargin, check_mark_width);    n++;

  item = XtCreateManagedWidget("menu_item", smeBSBObjectClass, menu, wargs, n);
  if (item == NULL)
    return NULL;

  if (func)
    XtAddCallback(item, XtNcallback, (XtCallbackProc)func, arg);

  return item;
}

void SetMenuItemChecked(Widget w, int state)
{
  int n=0;
  Arg wargs[5];
  Pixmap mark;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  if (lsx_curwin->check_mark == None)  /* create the check mark bitmap */
   {
     Display *d=XtDisplay(XtParent(w));

     mark = XCreateBitmapFromData(d, DefaultRootWindow(d),
                                  (char *)check_mark_bits,
                                  check_mark_width, check_mark_height);
     if (mark == None)
       return;

     lsx_curwin->check_mark = mark;
   }
  else
    mark = lsx_curwin->check_mark;

  /*
   * Now set the check mark.
   */
  n=0;
  XtSetArg(wargs[n], XtNleftMargin, 16);          n++;
  if (state)
    {
      /* checkmark on */
      XtSetArg(wargs[n], XtNleftBitmap, mark);    n++;
    }
  else
    {
      /* checkmark off */
      XtSetArg(wargs[n], XtNleftBitmap, None);    n++;
    }

  XtSetValues(w, wargs, n);
}

int GetMenuItemChecked(Widget w)
{
  int n=0;
  Arg wargs[5];
  Pixmap mark;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return FALSE;

  n=0;
  XtSetArg(wargs[n], XtNleftBitmap, &mark);    n++;
  XtGetValues(w, wargs, n);

  if (mark == None)
    return FALSE;
  else
    return TRUE;
}
