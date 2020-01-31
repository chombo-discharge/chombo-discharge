/* This file contains routines to manipulate a scrolling list widget.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include <X11/Xaw/Viewport.h>
#include <X11/Xaw/List.h>
#include "libsx.h"
#include "libsx_private.h"

extern WindowState *lsx_curwin;   /* global handle to current window */

/*
 * this structure maintains some internal state information about each
 * scrolled list widget.
 */
typedef struct ListInfo
{
  Widget w;
  void (*func)(Widget w, char *str, int index, void *data);
  void *data;
  struct ListInfo *next;
} ListInfo;

static ListInfo *scroll_lists = NULL;

/*
 * List Widget Creation Routines and stuff.
 */

static void list_callback(Widget w, XtPointer data, XtPointer call_data)
{
  ListInfo *li = (ListInfo *)data;
  XawListReturnStruct *list = (XawListReturnStruct *)call_data;

  if (li->func)
    li->func(w, list->string, list->list_index, li->data);
}

Widget MakeScrollList(char **item_list, int width, int height,
                      ListCB func, void *data)
{
  int    n = 0;
  Arg    wargs[10];             /* Used to set widget resources */
  Widget list, vport;
  ListInfo *li;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNwidth,  width);             n++;
  XtSetArg(wargs[n], XtNheight, height);            n++;
  XtSetArg(wargs[n], XtNallowVert, True);           n++;
  XtSetArg(wargs[n], XtNallowHoriz, True);          n++;
  XtSetArg(wargs[n], XtNuseBottom, True);           n++;

  vport = XtCreateManagedWidget("vport", viewportWidgetClass,
                               lsx_curwin->form_widget,wargs,n);
  if (vport == NULL)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNlist,   item_list);         n++;
  XtSetArg(wargs[n], XtNverticalList, True);        n++;
  XtSetArg(wargs[n], XtNforceColumns, True);        n++;
  XtSetArg(wargs[n], XtNdefaultColumns, 1);         n++;
  XtSetArg(wargs[n], XtNborderWidth, 1);            n++;

  /*
   * Here we create the list widget and make it the child of the
   * viewport widget so that the viewport will properly handle scrolling
   * it and all that jazz.
   */
  list = XtCreateManagedWidget("list", listWidgetClass,
                               vport,wargs,n);
  if (list == NULL)
   {
     XtDestroyWidget(vport);
     return NULL;
   }

  li = (ListInfo *)malloc(sizeof(ListInfo));
  if (li == NULL)
   {
     XtDestroyWidget(list);
     XtDestroyWidget(vport);
     return NULL;
   }

  li->func = func;
  li->data = data;
  li->w    = list;

  li->next = scroll_lists;
  scroll_lists = li;

  if (func)
    XtAddCallback(list, XtNcallback, list_callback, li);

  return list;
}    /* end of MakeScrollList() */

void SetCurrentListItem(Widget w, int list_index)
{
  if (w && list_index >= 0)
    XawListHighlight(w, list_index);
}

int GetCurrentListItem(Widget w)
{
  XawListReturnStruct *item;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return -1;

  item = XawListShowCurrent(w);
  if (item == NULL)
    return -1;

  return item->list_index;
}

void ChangeScrollList(Widget w, char **new_list)
{
  if (lsx_curwin->toplevel && w && new_list)
    XawListChange(w, new_list, -1, -1, TRUE);
}
