/* This file contains routines that manipulate single and multi-line
 * text entry widgets.
 *
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <stdlib.h>
#include "xstuff.h"
#include <X11/Xaw/AsciiText.h>
#include "libsx.h"
#include "libsx_private.h"

extern WindowState *lsx_curwin;   /* global handle to the current window */

/*
 * this structure maintains some internal state information about each
 * string entry widget.
 */
typedef struct StringInfo
{
  Widget str_widget;
  void (*func)(Widget w, char *string, void *data);
  void *user_data;

  struct StringInfo *next;
} StringInfo;

static StringInfo *string_widgets = NULL;

/*
 * String Entry Widget Creation stuff.
 */
Widget MakeStringEntry(char *txt, int size, StringCB func, void *data)
{
  static int already_done = FALSE;
  static XtTranslations trans;
  int    n = 0;
  Arg    wargs[10];             /* Used to set widget resources */
  Widget str_entry;
  StringInfo *stri;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  if (already_done == FALSE)
   {
     already_done = TRUE;
     trans = XtParseTranslationTable("#override\n\
                                      <Key>Return: done_with_text()\n\
                                      <Key>Linefeed: done_with_text()\n\
                                      Ctrl<Key>M: done_with_text()\n\
                                      Ctrl<Key>J: done_with_text()\n");
   }

  stri = (StringInfo *)malloc(sizeof(*stri));
  if (stri == NULL)
    return NULL;

  stri->func      =  func;
  stri->user_data = data;

  n = 0;
  XtSetArg(wargs[n], XtNeditType,     XawtextEdit);           n++;
  XtSetArg(wargs[n], XtNwrap,         XawtextWrapNever);      n++;
  XtSetArg(wargs[n], XtNresize,       XawtextResizeWidth);    n++;
  XtSetArg(wargs[n], XtNtranslations, trans);                 n++;
  XtSetArg(wargs[n], XtNwidth,        size);                  n++;
  if (txt)
   {
     XtSetArg(wargs[n], XtNstring,    txt);                   n++;
     XtSetArg(wargs[n], XtNinsertPosition,  strlen(txt));     n++;
   }

  str_entry = XtCreateManagedWidget("string", asciiTextWidgetClass,
                                    lsx_curwin->form_widget,wargs,n);

  if (str_entry)  /* only if we got a real widget do we bother */
   {
     stri->str_widget = str_entry;
     stri->next = string_widgets;
     string_widgets = stri;
   }
  else
    free(stri);

  return str_entry;
}                    /* end of MakeStringEntry() */

/*
 * Private internal callback for string entry widgets.
 */
void libsx_done_with_text(Widget w, XEvent *xev, String *parms,
                           Cardinal *num_parms)
{
  int    n = 0;
  Arg    wargs[10];             /* Used to get widget resources */
  char  *txt;
  StringInfo *stri;

  n = 0;
  XtSetArg(wargs[n], XtNstring,    &txt);                  n++;
  XtGetValues(w, wargs, n);

  /*
   * Find the correct ScrollInfo structure.
   */
  for (stri=string_widgets; stri; stri=stri->next)
   {
     for (; stri && stri->str_widget != w; stri=stri->next)
       ;

     if (stri)      /* didn't find it. */
       break;

     if (XtDisplay(stri->str_widget) == XtDisplay(w))  /* did find it */
       break;
   }

  if (stri == NULL)
    return;

  if (stri->func)
    stri->func(w, txt, stri->user_data);    /* call the user's function */
}

void SetStringEntry(Widget w, char *new_text)
{
  int  n = 0;
  Arg  wargs[2];                /* Used to set widget resources */

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  n = 0;
  XtSetArg(wargs[n], XtNstring, new_text);                   n++;
  XtSetValues(w, wargs, n);

  /*
   * Have to set this resource separately or else it doesn't get
   * updated in the display.  Isn't X a pain in the ass?
   *
   * (remember that with X windows, form follows malfunction)
   */
  n = 0;
  XtSetArg(wargs[n], XtNinsertPosition,  strlen(new_text));  n++;
  XtSetValues(w, wargs, n);
}

char *GetStringEntry(Widget w)
{
  int   n = 0;
  Arg   wargs[2];               /* Used to set widget resources */
  char *text;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNstring, &text);                   n++;
  XtGetValues(w, wargs, n);

  return text;
}

/*
 * Full Text Widget creation and support routines.
 */

/* forward prototype */
char *slurp_file(char *fname);

Widget MakeTextWidget(char *txt, int is_file, int editable, int w, int h)
{
  int n;
  Arg wargs[10];
  Widget text;
  char *real_txt;
  static int already_done = FALSE;
  static XtTranslations trans;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  if (already_done == FALSE)
   {
     already_done = TRUE;
     trans = XtParseTranslationTable("#override\n\
                                      <Key>Prior: previous-page()\n\
                                      <Key>Next:  next-page()\n\
                                      <Key>Home:  beginning-of-file()\n\
                                      <Key>End:   end-of-file()\n\
                                      Ctrl<Key>Up:    beginning-of-file()\n\
                                      Ctrl<Key>Down:  end-of-file()\n\
                                      Shift<Key>Up:   previous-page()\n\
                                      Shift<Key>Down: next-page()\n");
   }

  n=0;
  XtSetArg(wargs[n], XtNwidth,            w);                        n++;
  XtSetArg(wargs[n], XtNheight,           h);                        n++;
  XtSetArg(wargs[n], XtNscrollHorizontal, XawtextScrollWhenNeeded);  n++;
  XtSetArg(wargs[n], XtNscrollVertical,   XawtextScrollWhenNeeded);  n++;
  XtSetArg(wargs[n], XtNautoFill,         TRUE);                     n++;
  XtSetArg(wargs[n], XtNtranslations, trans);                        n++;

  if (is_file && txt)
   {
     real_txt = slurp_file(txt);
   }
  else
    real_txt = txt;

  if (real_txt)
    { XtSetArg(wargs[n], XtNstring,       real_txt);                 n++; }
  if (editable)
    { XtSetArg(wargs[n], XtNeditType,     XawtextEdit);              n++; }

  text = XtCreateManagedWidget("text", asciiTextWidgetClass,
                               lsx_curwin->form_widget,wargs,n);

  if (real_txt != txt && real_txt != NULL)
    free(real_txt);                         /* we're done with the buffer */

  return text;
}

void SetTextWidgetText(Widget w, char *txt, int is_file)
{
  int n;
  Arg wargs[2];
  char *real_txt;
  Widget source;

  if (lsx_curwin->toplevel == NULL || w == NULL || txt == NULL)
    return;

  source = XawTextGetSource(w);
  if (source == NULL)
    return;

  if (is_file)
   {
     real_txt = slurp_file(txt);
   }
  else
    real_txt = txt;

  n=0;
  XtSetArg(wargs[n], XtNstring, real_txt);  n++;
  XtSetValues(source, wargs, n);

  if (real_txt != txt && real_txt != NULL)
    free(real_txt);                         /* we're done with the buffer */
}

char *GetTextWidgetText(Widget w)
{
  int n;
  Arg wargs[4];
  char *text=NULL;
  Widget source;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return NULL;

  source = XawTextGetSource(w);
  if (source == NULL)
    return NULL;

  n=0;
  XtSetArg(wargs[n], XtNstring, &text);           n++;

  XtGetValues(source, wargs, n);

  return text;
}

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

char *slurp_file(char *fname)
{
  struct stat st;
  char *buff;
  int   fd, count;

  if (stat(fname, &st) < 0)
    return NULL;

  if (S_ISDIR(st.st_mode) == TRUE)   /* don't want to edit directories */
    return NULL;

  buff = (char *)malloc(sizeof(char)*(st.st_size+1));
  if (buff == NULL)
    return NULL;

  fd = open(fname, O_RDONLY);
  if (fd < 0)
   {
     free(buff);
     return NULL;
   }

  count = read(fd, buff, st.st_size);
  buff[count] = '\0';        /* null terminate */
  close(fd);

  return (buff);
}
