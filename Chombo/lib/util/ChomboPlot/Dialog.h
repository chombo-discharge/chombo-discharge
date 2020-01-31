/*
 *
 *  [ This file is blatantly borrowed from the pixmap distribution. ]
 *  [ It was written by the fellows below, and they disclaim all    ]
 *  [ warranties, expressed or implied, in this software.           ]
 *  [ As if anyone cares about that...                              ]
 *
 * Copyright 1991 Lionel Mallet
 *
 * Author:  Davor Matic, MIT X Consortium
 */

#ifndef _DIALOG_H_
#define _DIALOG_H_

#define Yes    16
#define No     32
#define Empty  0
#define Okay   1
#define Abort  2
#define Cancel 4
#define Retry  8

typedef struct
{
  Widget top_widget, shell_widget, dialog_widget;
  int options;
} _Dialog, *Dialog;

typedef struct
{
    String name;
    int flag;
} DialogButton;

Dialog CreateDialog(Widget top_widget, char *name, int options);
void   FreeDialog(Dialog dialog);
void   PositionPopup(Widget shell_widget);
int    PopupDialog(XtAppContext app_con, Dialog popup, char *message,
                   char *suggestion, char **answer, XtGrabKind grab);
void   PopdownDialog(Dialog popup, char **answer);

#endif
