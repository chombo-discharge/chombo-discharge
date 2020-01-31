/* This file contains routines to manage colormaps.  I've written all
 * of it except for the routines that generate the colormap data.  That
 * was taken from some code by Jeff LeBlanc who took it from some other
 * nameless soul.
 *
 *                     This code is under the GNU Copyleft.
 *
 *  Dominic Giampaolo
 *  dbg@sgi.com
 */

#include <stdio.h>
#include <math.h>
#include "xstuff.h"
#include "libsx.h"
#include "libsx_private.h"

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

/* internal prototypes */
static void grey_scale_plus(int ncells);
static void grey_scale(int ncells);
static void g_opt_2(int ncells);
static void b_opt_2(int ncells);

int WHITE  = 0,                  /* Global indicies into the color map */
    BLACK  = 0,
    RED    = 0,
    GREEN  = 0,
    BLUE   = 0,
    YELLOW = 0;

static XColor    col[256];
static int       ncells;          /* number of color cells available */

extern WindowState *lsx_curwin;

#define X_MAX_COLOR 65535

static void get_color(Colormap cmap, char *name, int *var)
{
  XColor exact, pixel_color;

  if (XAllocNamedColor(lsx_curwin->display, cmap, name, &exact, &pixel_color))
   {
     *var = pixel_color.pixel;
     lsx_curwin->named_colors[lsx_curwin->color_index++] = pixel_color.pixel;
   }
}

void GetStandardColors(void)
{
  Colormap mycmap;

  if (lsx_curwin->display == NULL || lsx_curwin->has_standard_colors)
    return;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  get_color(mycmap, "black", &BLACK);
  get_color(mycmap, "white", &WHITE);
  get_color(mycmap, "red",   &RED);
  get_color(mycmap, "green", &GREEN);
  get_color(mycmap, "blue",  &BLUE);
  get_color(mycmap, "yellow",&YELLOW);

  lsx_curwin->has_standard_colors = TRUE;
}

int GetNamedColor(char *name)
{
  Colormap mycmap;
  XColor exact, pixel_color;

  if (lsx_curwin->display == NULL)
    return -1;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  if (XAllocNamedColor(lsx_curwin->display, mycmap, name,&exact,&pixel_color))
   {
     lsx_curwin->named_colors[lsx_curwin->color_index++] = pixel_color.pixel;
     return pixel_color.pixel;
   }
  else
    return -1;
}

int GetRGBColor(int r, int g, int b)
{
  Colormap mycmap;
  XColor color;

  if (lsx_curwin->display == NULL)
    return -1;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  color.flags = DoRed | DoGreen | DoBlue;
  color.red   = (unsigned short) ((r * X_MAX_COLOR) / 256);
  color.green = (unsigned short) ((g * X_MAX_COLOR) / 256);
  color.blue  = (unsigned short) ((b * X_MAX_COLOR) / 256);

  if (XAllocColor(lsx_curwin->display, mycmap, &color))
   {
     lsx_curwin->named_colors[lsx_curwin->color_index++] = color.pixel;
     return color.pixel;
   }
  else
    return -1;
}

void FreeStandardColors(void)
{
  int i;
  Colormap mycmap;

  if (lsx_curwin->display == NULL)
    return;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  for (i=0; i < lsx_curwin->color_index; i++)
    XFreeColors(lsx_curwin->display, mycmap,
                (unsigned long *)&lsx_curwin->named_colors[i], 1, 0);

  lsx_curwin->color_index = 0;
}

int GetAllColors(void)
{
  int i;
  int mydepth;
  Visual *xv;
  Window win;

  if (lsx_curwin->display == NULL)
    return FALSE;

  if (lsx_curwin->cmap)   /* because we already have a custom cmap */
    return TRUE;

  ncells = 256;

  mydepth  = XDefaultDepth(lsx_curwin->display, lsx_curwin->screen);
  if (mydepth < 2)
    return FALSE;

  /* setup for colormap stuff */
  xv = DefaultVisual(lsx_curwin->display, lsx_curwin->screen);
  lsx_curwin->cmap = XCreateColormap(lsx_curwin->display, lsx_curwin->window,
                                     xv,AllocAll);
  if (lsx_curwin->cmap == None)
    return FALSE;

  for (i=0; i < ncells; i++)
   {
     col[i].pixel = i;
     col[i].red   = col[i].green = col[i].blue = 0; /* all black by default */
     col[i].flags = DoRed | DoGreen | DoBlue;
   }

  /*
   * Supposedly this is the correct thing to do when setting a colormap.
   *
   * That is, it sets the window colormap for the top-level shell widget
   * then it calls XSetWMColormapWindows() for the top-level shell widget
   * and the drawing area widget.
   *
   * For kicks we also call XChangeWindowAttributes() because it seems
   * like it might be a good thing to do (isn't it fun when there's 5 or
   * 6 ways to do the same thing?).  Remember that in X, form follows
   * malfunction.
   *
   * Unfortunately, the correct thing to do isn't always what works.
   * The code in the #ifdef section breaks twm rather horribly.  What
   * seems to happen is that when I call XSetWMColormapWindows() with
   * more than one id (in the same call or split across multiple calls)
   * and one of those window id's is a child of another, then twm dies.
   * On an SGI & DECstation twm gets a seg fault and dies (every single
   * time).
   *
   * So in summary, it seems that just calling XSetWindowColormap()
   * for the top-level shell widget and the drawing widget is what
   * works the most reliably across the widest number of systems and
   * window managers.
   *
   * I have personally seen this code work properly on a Sparc 10
   * w/OpenWindows version 3, a Solbourne running OpenWindows v2.0,
   * a whole schwack of SGI workstations running various versions of
   * IRIX and several different graphics systems (Elan, XS24, Reality
   * Engine, VGX, GT, and Entry Level graphics on Indigo).  It also
   * worked on a DECstation 3100 w/twm and dxwm.  I've also tested
   * this with mwm and 4Dwm on the SGI's.  Another window manager,
   * ctwm also seems to work.  That's a pretty wide range of platforms
   * (though HP is notably absent from the list and it probably doesn't
   * work there since I'm sure they're different than everyone else :^).
   */
  win = lsx_curwin->window;
  XSetWindowColormap(lsx_curwin->display, win, lsx_curwin->cmap);
  if (lsx_curwin->last_draw_widget)
    {
      win = XtWindow(lsx_curwin->last_draw_widget);
      XSetWindowColormap(lsx_curwin->display, win, lsx_curwin->cmap);
    }

#ifdef NO_BROKEN_WINDOW_MANAGERS
  winattr.colormap = lsx_curwin->cmap;
  XChangeWindowAttributes(lsx_curwin->display, win, CWColormap, &winattr);

  i = 0;
  win_ids[i] = lsx_curwin->window;
  if (lsx_curwin->last_draw_widget)
    {
      i++
      win_ids[i] = XtWindow(lsx_curwin->last_draw_widget);
    }
  XSetWMColormapWindows(lsx_curwin->display, lsx_curwin->window, win_ids, i);
#endif

  return TRUE;
}

void FreeAllColors(void)
{
  if (lsx_curwin->display == NULL)
    return;

  if (lsx_curwin->cmap == None)   /* woops, don't have a cmap to free */
    return;

  XFreeColormap(lsx_curwin->display, lsx_curwin->cmap);
  lsx_curwin->cmap = None;
}

void SetMyColorMap(int n, unsigned char *r, unsigned char *g,
                   unsigned char *b)
{
  int i;

  if (lsx_curwin->display == NULL || n < 0 || n > 256)
    return;

  if (lsx_curwin->cmap == None)
    if (GetAllColors() == FALSE)
      return;

  for (i=0; i < n; i++)
   {
     col[i].flags = DoRed | DoGreen | DoBlue;
     col[i].red   = (unsigned short)(r[i] * X_MAX_COLOR / 256);
     col[i].green = (unsigned short)(g[i] * X_MAX_COLOR / 256);
     col[i].blue  = (unsigned short)(b[i] * X_MAX_COLOR / 256);
   }

  XStoreColors(lsx_curwin->display, lsx_curwin->cmap, col, n);
}

void SetColorMap(int num)
{
  if (lsx_curwin->display == NULL)
    return;

  if (lsx_curwin->cmap == None)
    if (GetAllColors() == FALSE)
      return;

  switch (num)
   {
     case 0: grey_scale_plus(ncells);
             break;

     case 1: grey_scale(ncells);
             break;

     case 2: g_opt_2(ncells);
             break;

     case 3: b_opt_2(ncells);
             break;

     default: grey_scale_plus(ncells);
              break;
   }

  XStoreColors(lsx_curwin->display, lsx_curwin->cmap, col, ncells);
}

int GetPrivateColor(void)
{
  Colormap mycmap;
  unsigned long pixels[1], plane_masks[1];
  int result;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  result = XAllocColorCells(lsx_curwin->display, mycmap, FALSE,
                            &plane_masks[0], 0, &pixels[0], 1);

  if (result != 0)
    return pixels[0];
  else
    return -1;
}

void SetPrivateColor(int which, int r, int g, int b)
{
  Colormap mycmap;
  XColor pcol;

  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  pcol.pixel = which;
  pcol.flags = DoRed|DoGreen|DoBlue;
  pcol.red   = (short)(r * X_MAX_COLOR / 256);
  pcol.green = (short)(g * X_MAX_COLOR / 256);
  pcol.blue  = (short)(b * X_MAX_COLOR / 256);
  pcol.pad   = 0;

  XStoreColor(lsx_curwin->display, mycmap, &pcol);
}

void FreePrivateColor(int which)
{
  Colormap mycmap;
  unsigned long pixels[1];

  pixels[0] = which;
  if (lsx_curwin->cmap == None)
    mycmap = DefaultColormap(lsx_curwin->display,
                             DefaultScreen(lsx_curwin->display));
  else
    mycmap = lsx_curwin->cmap;

  XFreeColors(lsx_curwin->display, mycmap, &pixels[0], 1, 0);
}

/*****************************************************************************/

static void grey_scale_plus(int ncells)
{
  int i;
  int num;

  num = ncells - 4;
  /* gray scale */
  for (i = 0; i < num; i++)
   {
     col[i].flags = DoRed | DoGreen | DoBlue;
     col[i].red = (unsigned short)(i * X_MAX_COLOR / num);
     col[i].green = (unsigned short)(i * X_MAX_COLOR / num);
     col[i].blue = (unsigned short)(i * X_MAX_COLOR / num);
   }

  /* plus */

  /* yellow */
  col[num].flags = DoRed | DoGreen | DoBlue;
  col[num].red = (unsigned short)(X_MAX_COLOR);
  col[num].green = (unsigned short)(X_MAX_COLOR);
  col[num].blue = (unsigned short)(0);

  /* red */
  col[num+1].flags = DoRed | DoGreen | DoBlue;
  col[num+1].red = (unsigned short)(X_MAX_COLOR);
  col[num+1].green = (unsigned short)(0);
  col[num+1].blue = (unsigned short)(0);

  /* blue */
  col[num+2].flags = DoRed | DoGreen | DoBlue;
  col[num+2].red = (unsigned short)(0);
  col[num+2].green = (unsigned short)(X_MAX_COLOR);
  col[num+2].blue = (unsigned short)(0);

  /* green */
  col[num+3].flags = DoRed | DoGreen | DoBlue;
  col[num+3].red = (unsigned short)(0);
  col[num+3].green = (unsigned short)(0);
  col[num+3].blue = (unsigned short)(X_MAX_COLOR);

  BLACK  = col[0].pixel;
  WHITE  = col[num-1].pixel;
  YELLOW = col[num].pixel;
  RED    = col[num+1].pixel;
  BLUE   = col[num+2].pixel;
  GREEN  = col[num+3].pixel;

 }

/*****************************************************************************/

static void grey_scale(int ncells)
{
  int i;

  /* gray scale */
  for (i = 0; i < ncells ; i++)
   {
     col[i].flags = DoRed | DoGreen | DoBlue;
     col[i].red = (unsigned short)(i * X_MAX_COLOR / ncells);
     col[i].green = (unsigned short)(i * X_MAX_COLOR / ncells);
     col[i].blue = (unsigned short)(i * X_MAX_COLOR / ncells);
   }

  WHITE  = col[ncells-1].pixel;
  BLACK  = col[0].pixel;
  RED    = BLACK;
  GREEN  = BLACK;
  BLUE   = BLACK;
  YELLOW = BLACK;

}

/*****************************************************************************/

static void g_opt_2(int ncells)
{
  int i;
  float con;
  float x,tmp;

  con = 1.0 / ( (1-0.3455)*(1-0.3455) * (1-0.90453)*(1-0.90453) );
  for (i = 0; i < ncells ; i++)
   {
     col[i].flags = DoRed | DoGreen | DoBlue;
     /* ramp on red */
     col[i].red = (unsigned short)(i * X_MAX_COLOR / ncells);
     x = (float)(i)/(float)(ncells);

     /* double hump on green */
     tmp = con*x*(x - 0.3455)*(x - 0.3455)*(x - 0.90453)*(x - 0.90453);
     if (tmp > 1.0) tmp = 1.0;
     if (tmp < 0.0) tmp = 0.0;
     col[i].green = (unsigned short)((int)((float)(X_MAX_COLOR) * tmp));

     /* single hump on blue */
     tmp = x*(4*x - 3)*(4*x - 3);
     if (tmp > 1.0) tmp = 1.0;
     if (tmp < 0.0) tmp = 0.0;
     col[i].blue = (unsigned short)((int)((float)(X_MAX_COLOR) * tmp));
   }

  WHITE = 0;
  BLACK = 256 - ncells;
  RED = BLACK;
  GREEN = BLACK;
  BLUE = BLACK;
  YELLOW = WHITE;

}

/*****************************************************************************/

static void b_opt_2(int ncells)
{
  int i;
  float con;
  float x,tmp;

  con = 1.0 / ( (1-0.3455)*(1-0.3455) * (1-0.90453)*(1-0.90453) );
  for (i = 0; i < ncells ; i++)
   {
     col[i].flags = DoRed | DoGreen | DoBlue;
     /* ramp on red */
     col[i].red = (unsigned short)(i * X_MAX_COLOR / ncells);
     x = (float)(i)/(float)(ncells);

     /* single hump on green */
     tmp = x*(4*x - 3)*(4*x - 3);
     if (tmp > 1.0) tmp = 1.0;
     if (tmp < 0.0) tmp = 0.0;
     col[i].green = (unsigned short)((int)((float)(X_MAX_COLOR) * tmp));

     /* double hump on blue */
     tmp = con*x*(x - 0.3455)*(x - 0.3455)*(x - 0.90453)*(x - 0.90453);
     if (tmp > 1.0) tmp = 1.0;
     if (tmp < 0.0) tmp = 0.0;
     col[i].blue = (unsigned short)((int)((float)(X_MAX_COLOR) * tmp));
   }

  WHITE = 0;
  BLACK = 256 - ncells;
  RED = BLACK;
  GREEN = BLACK;
  BLUE = BLACK;
  YELLOW = WHITE;

}

/*****************************************************************************/
