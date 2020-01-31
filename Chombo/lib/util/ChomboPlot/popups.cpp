/* This file contains routines that handle popping up dialog boxes.
 * They use the routines in Dialog.c to do most of the work.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include "Dialog.h"
#include "libsx.h"
#include "libsx_private.h"

extern WindowState  *lsx_curwin;    /* global handle to the current window */
extern XtAppContext  lsx_app_con;

/*
 * User input routines that take place through Dialog boxes (getting
 * a string and a simple yes/no answer).
 */

char *GetString(char *blurb, char *default_string)
{
  char *string = default_string;
  Dialog dialog = NULL;

  if (blurb == NULL || (lsx_curwin->toplevel==NULL && OpenDisplay(0,NULL)==0))
    return NULL;

  dialog = CreateDialog(lsx_curwin->toplevel, "Input Window", Okay | Cancel);

  if (dialog == NULL)   /* then there's an error */
    return NULL;

  if (string == NULL)
    string = "";

  switch(PopupDialog(lsx_app_con, dialog, blurb, string, &string,
                     XtGrabExclusive))
   {
     case Okay: /* don't have to do anything, string is already set */
                break;

     case Cancel:  string = NULL;
                   break;

     default: string = NULL;    /* shouldn't happen, but just in case */
              break;
   }  /* end of switch */

  FreeDialog(dialog);
  return string;
}

int GetYesNo(char *blurb)
{
  Dialog dialog = NULL;
  int ret;

  if (blurb == NULL || (lsx_curwin->toplevel==NULL && OpenDisplay(0,NULL)==0))
    return FALSE;

  dialog = CreateDialog(lsx_curwin->toplevel, "Input Window2", Okay|Cancel);

  if (dialog == NULL)   /* then there's an error */
    return FALSE;

  switch(PopupDialog(lsx_app_con, dialog, blurb, NULL, NULL, XtGrabExclusive))
   {
     case Okay: ret = TRUE;
                break;

     case Cancel:  ret = FALSE;
                   break;

     default: ret = FALSE;   /* unknown return from dialog, return an err */
              break;
   }  /* end of switch */

  FreeDialog(dialog);

  return ret;
}
