#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*
 * Structures and variables private to libsx.
 */

#ifndef _LIBSX_PRIVATE_H_
#define _LIBSX_PRIVATE_H_

typedef struct WindowState
{
  int      screen;
  int      window_shown;
  Window   window;
  Display *display;
  Widget   toplevel, toplevel_form, form_widget, last_draw_widget;
  int      has_standard_colors;
  int      named_colors[256];
  int      color_index;
  Colormap     cmap;
  Pixmap       check_mark;
  XFontStruct *font;

#ifdef    OPENGL_SUPPORT

  XVisualInfo *xvi;
  GLXContext  gl_context;

#endif /* OPENGL_SUPPORT */

  struct WindowState *next;
} WindowState;

extern WindowState *lsx_curwin;  /* defined in libsx.c */

typedef struct DrawInfo
{
  RedisplayCB   redisplay;
  MouseButtonCB button_down;
  MouseButtonCB button_up;
  KeyCB         keypress;
  MotionCB      motion;

  GC            drawgc;       /* Graphic Context for drawing in this widget */

  int           foreground;   /* need to save these so we know what they are */
  int           background;
  unsigned long mask;
  XFontStruct  *font;

  void        *user_data;

  Widget widget;
  struct DrawInfo *next;
} DrawInfo;

DrawInfo *libsx_find_draw_info(Widget w);  /* private internal function */

#endif
