/* This file contains the top level routines that manipulate whole
 * windows and what not.  It's the main initialization stuff, etc.
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

/*
 * External prototypes.
 */
extern void SetOkay(Widget w);
extern void PositionPopup(Widget w);                  /* from Dialog.c */
extern void libsx_done_with_text(Widget w, XEvent *xev, String *parms,
                                 Cardinal *num_parms);/* from string_entry.c */

/* internal protos */
static Bool  is_expose_event(Display *d, XEvent *xev, char *blorg);

/*
 * Private static variables.
 */
static WindowState  *lsx_windows  = NULL,
                     empty_window =
                     {
                       0, 0, 0, NULL, NULL, NULL, NULL, 0,
                     },
                    *orig_window  = NULL,
                    *new_curwin   = NULL;

static char *app_name="";
static Display *base_display=NULL;

/*
 * Global Variables for all of libsx (but not for user applications).
 *
 * We initialize lsx_curwin to point to the static empty_window so that
 * if other functions get called before OpenDisplay(), they won't core
 * dump, and will fail gracefully instead.
 */
WindowState  *lsx_curwin = &empty_window;
XtAppContext  lsx_app_con;

/*
 * These are some internal actions and resources so that dialog boxes
 * and string entry widgets work they way god intended them to work
 * (i.e. pressing the return key in one gives the text you typed to
 * the application).
 */
static XtActionsRec actions_table[] =
{
  {
    "set-okay",        (XtActionProc) SetOkay
  },
  {
    "done_with_text", libsx_done_with_text
  }
};

static String fallback_resources[] =
{
  "*Dialog*label.resizable: True",
  "*Dialog*Text.resizable: True",
  "*Dialog.Text.translations: #override <Key>Return: set-okay()\\n\
                              <Key>Linefeed: set-okay()",
  NULL
};

/*
 *
 * This is the function that gets everything started.
 *
 */
int _OpenDisplay(int argc, char **argv, Display *dpy, Arg *wargs, int arg_cnt)
{
  static char *dummy_argv[] =
  {
    "Untitled", NULL
  };
  static int   dummy_argc  = 1;
  static int already_called = FALSE;

  if (already_called)   /* must have been called a second time */
    return FALSE;
  already_called = TRUE;

  if (argv == NULL)
    {
      argv = dummy_argv;
      argc = dummy_argc;
    }

  /*
   * First create the window state.
   */
  lsx_curwin = (WindowState *)calloc(sizeof(WindowState), 1);
  if (lsx_curwin == NULL)
    return FALSE;

  /* open the display stuff */
  if (dpy == NULL)
    lsx_curwin->toplevel = XtAppInitialize(&lsx_app_con, argv[0], NULL, 0,
                                           &argc, argv, fallback_resources,
                                           wargs, arg_cnt);
  else
    lsx_curwin->toplevel = XtAppCreateShell (argv[0], argv[0],
                                             applicationShellWidgetClass,
                                             dpy, wargs, arg_cnt);

  if (lsx_curwin->toplevel == NULL)
   {
     free(lsx_curwin);
     return FALSE;
   }

  app_name  = argv[0];   /* save this for later */

  XtAppAddActions(lsx_app_con, actions_table, XtNumber(actions_table));

  lsx_curwin->form_widget = XtCreateManagedWidget("form", formWidgetClass,
                                                  lsx_curwin->toplevel,NULL,0);

  if (lsx_curwin->form_widget == NULL)
   {
     XtDestroyWidget(lsx_curwin->toplevel);
     free(lsx_curwin);
     lsx_curwin = &empty_window;
     return FALSE;
   }
  lsx_curwin->toplevel_form = lsx_curwin->form_widget;

  lsx_curwin->next = lsx_windows;    /* link it in to the list */
  lsx_windows = lsx_curwin;

  /* save these values for later */
  lsx_curwin->display = (Display *)XtDisplay(lsx_curwin->toplevel);
  lsx_curwin->screen  = DefaultScreen(lsx_curwin->display);
  orig_window  = lsx_curwin;
  base_display = lsx_curwin->display;

  return argc;
} /* end of OpenDisplay() */

int OpenDisplay(int argc, char **argv)
{
  return _OpenDisplay(argc, argv, NULL, NULL, 0);
}

#ifdef OPENGL_SUPPORT

int OpenGLDisplay(int argc, char **argv, int *attributes)
{
  int xargc, cnt, i, retval, count;
  char **xargv;
  Arg args[5];
  Display *dpy;
  XVisualInfo *xvi;
  GLXContext cx;
  XVisualInfo   vinfo;
  XVisualInfo   *vinfo_list;    /* returned list of visuals */
  Colormap      colormap;       /* created colormap */
  Widget        top;

  /*
   * First we copy the command line arguments
   */
  xargc = argc;
  xargv = (char **)malloc(argc * sizeof (char *));
  for (i=0; i < xargc; i++)
    xargv[i] = strdup(argv[i]);

  /*
   * The following creates a _dummy_ toplevel widget so we can
   * retrieve the appropriate visual resource.
   */
  cnt = 0;
  top = XtAppInitialize(&lsx_app_con, xargv[0], NULL, 0, &xargc, xargv,
                        (String *)NULL, args, cnt);
  if (top == NULL)
    return 0;

  dpy = XtDisplay(top);

  /*
   * Check if the server supports GLX.  If not, crap out.
   */
  if (glXQueryExtension(dpy, NULL, NULL) == GL_FALSE)
   {
     XtDestroyWidget(top);
     return FALSE;
   }

  xvi = glXChooseVisual(dpy, DefaultScreen(dpy), attributes);
  if (xvi == NULL)
   {
     XtDestroyWidget(top);
     return 0;
   }

  cx  = glXCreateContext(dpy, xvi, 0, GL_TRUE);
  if (cx == NULL)
   {
     XtDestroyWidget(top);
     return 0;
   }

  cnt = 0;
  XtSetArg(args[cnt], XtNvisual, xvi->visual); cnt++;

  /*
   * Now we create an appropriate colormap.  We could
   * use a default colormap based on the class of the
   * visual; we could examine some property on the
   * rootwindow to find the right colormap; we could
   * do all sorts of things...
   */
  colormap = XCreateColormap (dpy,
                              RootWindowOfScreen(XtScreen (top)),
                              xvi->visual,
                              AllocNone);
  XtSetArg(args[cnt], XtNcolormap, colormap); cnt++;

  /*
   * Now find some information about the visual.
   */
  vinfo.visualid = XVisualIDFromVisual(xvi->visual);
  vinfo_list = XGetVisualInfo(dpy, VisualIDMask, &vinfo, &count);
  if (vinfo_list && count > 0)
   {
     XtSetArg(args[cnt], XtNdepth, vinfo_list[0].depth);
     cnt++;
     XFree((XPointer) vinfo_list);
   }

  XtDestroyWidget(top);

  /*
   * Free up the copied version of the command line arguments
   */
  for (i=0; i < xargc; i++)
    if (xargv[i])
      free(xargv[i]);
  free(xargv);

  retval = _OpenDisplay(argc, argv, dpy, args, cnt);

  if (retval > 0)
   {
     lsx_curwin->xvi        = xvi;
     lsx_curwin->gl_context = cx;
     lsx_curwin->cmap       = colormap;
   }

  return retval;
}

#endif   /* OPENGL_SUPPORT */

void ShowDisplay(void)
{
  XEvent xev;

  if (lsx_curwin->toplevel == NULL || lsx_curwin->window_shown == TRUE)
    return;

  XtRealizeWidget(lsx_curwin->toplevel);

  if (XtIsTransientShell(lsx_curwin->toplevel))  /* do popups differently */
   {
     PositionPopup(lsx_curwin->toplevel);
     XtPopup(lsx_curwin->toplevel, XtGrabExclusive);

     lsx_curwin->window  = (Window   )XtWindow(lsx_curwin->toplevel);
     lsx_curwin->window_shown = TRUE;

     return;
   }

  /*
   * wait until the window is _really_ on screen
   */
  while (!XtIsRealized(lsx_curwin->toplevel))
    ;

  /*
   * Now make sure it is really on the screen.
   */
  XPeekIfEvent(XtDisplay(lsx_curwin->toplevel), &xev, is_expose_event, NULL);

  SetDrawArea(lsx_curwin->last_draw_widget);

  lsx_curwin->window  = (Window   )XtWindow(lsx_curwin->toplevel);
  lsx_curwin->window_shown = TRUE;
}   /* end of ShowDisplay() */

static Bool is_expose_event(Display *d, XEvent *xev, char *blorg)
{
  if (xev->type == Expose)
    return TRUE;
  else
    return FALSE;
}

void MainLoop(void)
{
  if (lsx_curwin->toplevel == NULL)
    return;

  /* in case the user forgot to map the display, do it for them */
  if (lsx_curwin->window_shown == FALSE)
   {
     ShowDisplay();
     GetStandardColors();
   }

  if (XtIsTransientShell(lsx_curwin->toplevel)) /* handle popups differently */
   {
     WindowState *curwin = lsx_curwin;

     while (curwin->window_shown == TRUE)  /* while window is still open */
      {
        XEvent event;

        XtAppNextEvent(lsx_app_con, &event);
        XtDispatchEvent(&event);
      }

     /*
      * Ok, at this point the popup was just closed, so now CloseWindow()
      * stored some info for us in the global variable new_curwin (which
      * we use to change lsx_curwin (after free'ing what lsx_curwin used
      * to point to)..
      */
     free(lsx_curwin);
     lsx_curwin = new_curwin;
     new_curwin = NULL;

     return;
   }
  else
    XtAppMainLoop(lsx_app_con);
}

void SyncDisplay(void)
{
  if (lsx_curwin->display)
    XSync(lsx_curwin->display, FALSE);
}

void AddTimeOut(unsigned long interval, void (*func)(), void *data)
{
  if (lsx_curwin->toplevel && func)
    XtAppAddTimeOut(lsx_app_con, interval, (XtTimerCallbackProc)func, data);
}

void  AddReadCallback(int fd,  IOCallback func, void *data)
{
  XtInputMask mask = XtInputReadMask;

  XtAppAddInput(lsx_app_con, fd, (XtPointer)mask,
                (XtInputCallbackProc)func, data);
}

void  AddWriteCallback(int fd, IOCallback func, void *data)
{
  XtInputMask mask = XtInputWriteMask;

  XtAppAddInput(lsx_app_con, fd, (XtPointer)mask,
                (XtInputCallbackProc)func, data);
}

Widget MakeWindow(char *window_name, char *display_name, int exclusive)
{
  WindowState *win=NULL;
  Display *d=NULL;
  Arg wargs[20];
  int n=0;
  Visual *vis;
  Colormap cmap;
  char *argv[5];
  int   argc;

  if (lsx_curwin->display == NULL)   /* no other windows open yet... */
    return NULL;

  win = (WindowState *)calloc(sizeof(WindowState), 1);
  if (win == NULL)
    return NULL;

  /*
   * Setup a phony argv/argc to appease XtOpenDisplay().
   */
  if (window_name)
    argv[0] = window_name;
  else
    argv[0] = app_name;
  argv[1] = NULL;
  argc = 1;

  if (display_name != NULL)
    d = XtOpenDisplay(lsx_app_con, display_name, app_name, app_name, NULL, 0,
                      &argc, argv);
  else
    d = base_display;
  if (d == NULL)
   {
     free(win);
     return NULL;
   }

  win->display  = d;

  cmap = DefaultColormap(d, DefaultScreen(d));
  vis  = DefaultVisual(d, DefaultScreen(d));

  n=0;
  XtSetArg(wargs[n], XtNtitle,    app_name);      n++;
  XtSetArg(wargs[n], XtNiconName, app_name);      n++;
  XtSetArg(wargs[n], XtNcolormap, cmap);          n++;
  XtSetArg(wargs[n], XtNvisual,   vis);           n++;

  if (exclusive == FALSE)
   {
     win->toplevel = XtAppCreateShell(argv[0], app_name,
                                      topLevelShellWidgetClass, d, wargs, n);
   }
  else
   {
     win->toplevel = XtCreatePopupShell(argv[0], transientShellWidgetClass,
                                        lsx_curwin->toplevel, NULL, 0);
   }

  if (win->toplevel == NULL)
   {
     if (d != base_display)
       XtCloseDisplay(d);
     free(win);
     return NULL;
   }

  win->form_widget = XtCreateManagedWidget("form", formWidgetClass,
                                           win->toplevel, NULL, 0);
  if (win->form_widget == NULL)
   {
     XtDestroyWidget(win->toplevel);
     if (d != base_display)
       XtCloseDisplay(d);
     free(win);
     return NULL;
   }
  win->toplevel_form = win->form_widget;

  win->screen = DefaultScreen(win->display);

  /*
   * Now link in the new window into the window list and make it
   * the current window.
   */
  win->next   = lsx_windows;
  lsx_windows = win;
  lsx_curwin  = win;

  return win->toplevel;    /* return a handle to the user */
}

Widget MakeForm(Widget parent, int where1, Widget from1,
                int where2, Widget from2)
{
  Widget form;

  if (lsx_curwin->toplevel == NULL)
    return NULL;

  if (parent == TOP_LEVEL_FORM)
    parent = lsx_curwin->toplevel_form;
  else if (strcmp(XtName(parent), "form") != 0)  /* parent not a form widget */
    return NULL;

  form = XtCreateManagedWidget("form", formWidgetClass,
                               parent, NULL, 0);
  if (form == NULL)
    return NULL;

  SetWidgetPos(form, where1, from1, where2, from2);

  lsx_curwin->form_widget = form;

  return form;
}

void SetForm(Widget form)
{
  if (lsx_curwin->toplevel == NULL)
    return;

  if (form == TOP_LEVEL_FORM)
    lsx_curwin->form_widget = lsx_curwin->toplevel_form;
  else
    lsx_curwin->form_widget = form;
}

void SetCurrentWindow(Widget w)
{
  WindowState *win;

  if (w == ORIGINAL_WINDOW)
   {
     if (orig_window)
       lsx_curwin = orig_window;
     else if (lsx_windows)            /* hmm, don't have orig_window */
       lsx_curwin = lsx_windows;
     else
       lsx_curwin = &empty_window;    /* hmm, don't have anything left */

     SetDrawArea(lsx_curwin->last_draw_widget);
     return;
   }

  for (win=lsx_windows; win; win=win->next)
    if (win->toplevel == w && win->display == XtDisplay(w))
      break;

  if (win == NULL)
    return;

  lsx_curwin = win;
  SetDrawArea(lsx_curwin->last_draw_widget);
}

void CloseWindow(void)
{
  WindowState *tmp;
  int is_transient;

  if (lsx_curwin->toplevel == NULL)
    return;

  is_transient = XtIsTransientShell(lsx_curwin->toplevel);

  XtDestroyWidget(lsx_curwin->toplevel);

  if (lsx_curwin->display != base_display)
   {
     FreeFont(lsx_curwin->font);
     XtCloseDisplay(lsx_curwin->display);
   }

  /*
   * if the window to delete is not at the head of the list, find
   * it in the list of available windows.
   * else, just assign tmp to the head of the list.
   */
  if (lsx_curwin != lsx_windows)
    for (tmp=lsx_windows; tmp && tmp->next != lsx_curwin; tmp=tmp->next)
      ;
  else
    tmp = lsx_windows;

  if (tmp == NULL)  /* bogus window deletion... */
   {
     return;
   }
  else if (lsx_curwin == lsx_windows)   /* delete head of list */
   {
     lsx_windows = lsx_curwin->next;
     tmp = lsx_windows;
   }
  else
   {
     tmp->next = lsx_curwin->next;
   }

  if (lsx_curwin == orig_window)
    orig_window = NULL;

  if (is_transient)
   {
     lsx_curwin->window_shown = FALSE;

     /*
      * Store tmp in a global so MainLoop() can change lsx_curwin
      * to what it should be after it frees everything.
      */
     new_curwin = tmp;
   }
  else    /* it's a normal window, so get rid of it completely */
   {
     free(lsx_curwin);
     lsx_curwin = tmp;
   }

}

void Beep(void)
{
  XBell(lsx_curwin->display, 100);
}
