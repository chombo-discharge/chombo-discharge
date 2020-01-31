/*
 *  [ This file is blatantly borrowed from the pixmap distribution. ]
 *  [ It was written by the fellows below, and they disclaim all    ]
 *  [ warranties, expressed or implied, in this software.           ]
 *  [ As if anyone cares about that...                              ]
 *
 * Copyright 1991 Lionel Mallet
 *
 * Author:  Davor Matic, MIT X Consortium
 */

#include <stdio.h>
#include <stdlib.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xos.h>
#include <X11/Shell.h>
#include <X11/Xatom.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/AsciiText.h>

#include "Dialog.h"

#define min(x, y)                     (((x) < (y)) ? (x) : (y))
#define max(x, y)                     (((x) > (y)) ? (x) : (y))

//static int selected;
static union
{
  int isel;
  XtPointer xtpsel;
} selected;

static DialogButton dialog_buttons[] =
{
    {
      "Yes", Yes
    },
    {
      "No", No
    },
    {
      "Okay", Okay
    },
    {
      "Abort", Abort
    },
    {
      "Cancel", Cancel
    },
    {
      "Retry", Retry
    },
};

static void SetSelected(
                        Widget w,
                        XtPointer client_data,
                        XtPointer call_data)
{
  selected.xtpsel = client_data;
}

/*
 * Can't make this static because we need to be able to access the
 * name over in display.c to properly set the resources we want.
 */
void SetOkay(Widget w)
{
  SetSelected(w, (XtPointer)( Okay | Yes ), (XtPointer)(NULL));
}

Dialog CreateDialog(Widget top_widget, char *name, int options)
{
  int i;
  Dialog popup;

  if ((popup = (Dialog) XtMalloc(sizeof(_Dialog))) == NULL)
    return NULL;

  popup->top_widget = top_widget;
  popup->shell_widget = XtCreatePopupShell(name,
                                           transientShellWidgetClass,
                                           top_widget, NULL, 0);
  popup->dialog_widget = XtCreateManagedWidget("dialog",
                                               dialogWidgetClass,
                                               popup->shell_widget,
                                               NULL, 0);

  for (i = 0; i < XtNumber(dialog_buttons); i++)
    if (options & dialog_buttons[i].flag)
      XawDialogAddButton(popup->dialog_widget,
                         dialog_buttons[i].name,
                         SetSelected, (XtPointer)dialog_buttons[i].flag);

  popup->options = options;
  return popup;
}

void FreeDialog(Dialog dialog)
{
  if (dialog == NULL)
    return;

  XtDestroyWidget(dialog->shell_widget);
  XtFree((char *)dialog);
}

void PositionPopup(Widget shell_widget)
{
  int n;
  Arg wargs[10];
  Position popup_x, popup_y, top_x, top_y;
  Dimension popup_width, popup_height, top_width, top_height, border_width;

  n = 0;
  XtSetArg(wargs[n], XtNx, &top_x); n++;
  XtSetArg(wargs[n], XtNy, &top_y); n++;
  XtSetArg(wargs[n], XtNwidth, &top_width); n++;
  XtSetArg(wargs[n], XtNheight, &top_height); n++;
  XtGetValues(XtParent(shell_widget), wargs, n);

  n = 0;
  XtSetArg(wargs[n], XtNwidth, &popup_width); n++;
  XtSetArg(wargs[n], XtNheight, &popup_height); n++;
  XtSetArg(wargs[n], XtNborderWidth, &border_width); n++;
  XtGetValues(shell_widget, wargs, n);

  popup_x = max(0,
                min(top_x + ((Position)top_width - (Position)popup_width) / 2,
                    (Position)DisplayWidth(XtDisplay(shell_widget),
                               DefaultScreen(XtDisplay(shell_widget))) -
                    (Position)popup_width - 2 * (Position)border_width));
  popup_y = max(0,
                min(top_y+((Position)top_height - (Position)popup_height) / 2,
                    (Position)DisplayHeight(XtDisplay(shell_widget),
                               DefaultScreen(XtDisplay(shell_widget))) -
                    (Position)popup_height - 2 * (Position)border_width));
  n = 0;
  XtSetArg(wargs[n], XtNx, popup_x); n++;
  XtSetArg(wargs[n], XtNy, popup_y); n++;
  XtSetValues(shell_widget, wargs, n);
}

int PopupDialog(XtAppContext app_con, Dialog popup, char *message,
                char *suggestion, char **answer, XtGrabKind grab)
{
  int n, height = 35;
  Arg wargs[8];
  Widget value;

  n = 0;
  XtSetArg(wargs[n], XtNlabel, message); n++;
  if (suggestion)
   {
     XtSetArg(wargs[n], XtNvalue, suggestion); n++;
   }
  XtSetValues(popup->dialog_widget, wargs, n);

  /*
   * Here we get ahold of the ascii text widget for the dialog box (if it
   * exists) and we set some useful resources in it.
   */
  value = XtNameToWidget(popup->dialog_widget, "value");

  n = 0;
  XtSetArg(wargs[n], XtNresizable,        True);                n++;
  XtSetArg(wargs[n], XtNheight,           height);              n++;
  XtSetArg(wargs[n], XtNwidth,            250);                 n++;
  XtSetArg(wargs[n], XtNresize,           XawtextResizeHeight); n++;
  XtSetArg(wargs[n], XtNscrollHorizontal, XawtextScrollWhenNeeded); n++;
  if (value)
    XtSetValues(value, wargs, n);

  XtRealizeWidget(popup->shell_widget);

  PositionPopup(popup->shell_widget);

  selected.isel = Empty;

  XtPopup(popup->shell_widget, grab);

  while ((selected.isel & popup->options) == Empty)
   {
     XEvent event;

     XtAppNextEvent(app_con, &event);
     XtDispatchEvent(&event);
   }

  PopdownDialog(popup, answer);

  return (selected.isel & popup->options);
}

void PopdownDialog(Dialog popup, char **answer)
{
    if (answer)
      *answer = XawDialogGetValueString(popup->dialog_widget);

    XtPopdown(popup->shell_widget);
}
