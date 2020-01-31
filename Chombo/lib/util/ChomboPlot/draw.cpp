/* This file contains routines that encapsulate a drawing area widget.
 * It is slightly different in format than the other files because it
 * also has some callback wrappers that extract relevant information and
 * then pass that to the user callback functions for redisplaying,
 * keyboard input and mouse clicks/motion.
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
#include "DrawingA.h"

extern WindowState *lsx_curwin;    /* global handle for the current window */

/*
 * Internal prototypes that massage the data into something more
 * digestable by regular humans.
 */
static void  _redisplay(Widget w, void *me, XADCS *call_data);
static void  _resize(Widget w,    void *me, XADCS *call_data);
static void  _do_input(Widget w,  void *me, XADCS *call_data);
static void  _do_motion(Widget w, void *me, XADCS *call_data);
static char *TranslateKeyCode(XEvent *ev);
static GC    setup_gc(Widget w);

static DrawInfo *draw_info_head = NULL;
static DrawInfo *cur_di=NULL;   /* current drawing area info structure */
static Window    window;        /* only used below by the drawing functions. */
static GC        drawgc;
static Display  *display=NULL;

/*
 * Drawing Area Creation Routine.
 */

Widget MakeDrawArea(int width, int height, RedisplayCB redisplay, void *data)
{
  int    n = 0;
  Arg    wargs[5];              /* Used to set widget resources */
  Widget draw_widget;
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL && OpenDisplay(0, NULL) == 0)
    return NULL;

  di = (DrawInfo *)calloc(sizeof(DrawInfo), 1);
  if (di == NULL)
    return NULL;

  n = 0;
  XtSetArg(wargs[n], XtNwidth, width);          n++;
  XtSetArg(wargs[n], XtNheight,height);         n++;

  draw_widget = XtCreateManagedWidget("drawing_area", drawingAreaWidgetClass,
                                      lsx_curwin->form_widget,wargs,n);

  if (draw_widget == NULL)
   {
     free(di);
     return NULL;
   }

  di->drawgc     = setup_gc(draw_widget);
  di->foreground = BlackPixel(lsx_curwin->display, lsx_curwin->screen);
  di->background = WhitePixel(lsx_curwin->display, lsx_curwin->screen);
  di->mask       = 0xffffffff;

  di->user_data   = data;
  di->redisplay   = redisplay;

  XtAddCallback(draw_widget, XtNexposeCallback, (XtCallbackProc)_redisplay,di);
  XtAddCallback(draw_widget, XtNresizeCallback, (XtCallbackProc)_resize,   di);
  XtAddCallback(draw_widget, XtNinputCallback,  (XtCallbackProc)_do_input, di);
  XtAddCallback(draw_widget, XtNmotionCallback, (XtCallbackProc)_do_motion,di);

  lsx_curwin->last_draw_widget = draw_widget;

  di->widget = draw_widget;
  di->next = draw_info_head;
  draw_info_head = di;
  cur_di = di;

  /*
   * Make sure the font is set to something sane.
   */
  if (lsx_curwin->font == NULL)
    lsx_curwin->font = GetFont("fixed");
  SetWidgetFont(draw_widget, lsx_curwin->font);

  return draw_widget;
}

/*
 * Internal function for getting a graphics context so we can draw.
 */
static GC setup_gc(Widget w)
{
  int fore_g,back_g;      /* Fore and back ground pixels */
  GC  drawgc;

  back_g = WhitePixel(XtDisplay(w),DefaultScreen(XtDisplay(w)));
  fore_g = BlackPixel(XtDisplay(w),DefaultScreen(XtDisplay(w)));

  /* Create drawing GC */
  drawgc = XCreateGC(XtDisplay(w), DefaultRootWindow(XtDisplay(w)), 0, 0);

  XSetBackground(XtDisplay(w), drawgc, back_g);
  XSetForeground(XtDisplay(w), drawgc, fore_g);
  XSetFunction(XtDisplay(w),   drawgc, GXcopy);

  return drawgc;
} /* end of setup_gc() */

/*
 * This function searches through our list of drawing area info structures
 * and returns the one that matches the specified widget.  We have to be
 * careful and make sure that the DrawInfo structure we match is for the
 * right display.
 */
DrawInfo *libsx_find_draw_info(Widget w)
{
  DrawInfo *di;

  if (w == NULL)
    return NULL;

  for (di=draw_info_head;  di; di=di->next)
   {
     for (; di && di->widget != w; di=di->next)
       ;

     if (di == NULL)       /* we didn't find it */
       break;
     if (XtDisplay(di->widget) == XtDisplay(w)) /* then we've found it */
       break;
   }

  return di;
}

void   SetButtonDownCB(Widget w, MouseButtonCB button_down)
{
  DrawInfo *di;

  if ((di = libsx_find_draw_info(w)) == NULL)
    return;

  di->button_down = button_down;
}

void   SetButtonUpCB(Widget w, MouseButtonCB button_up)
{
  DrawInfo *di;

  if ((di = libsx_find_draw_info(w)) == NULL)
    return;

  di->button_up = button_up;
}

void   SetKeypressCB(Widget w, KeyCB keypress)
{
  DrawInfo *di;

  if ((di = libsx_find_draw_info(w)) == NULL)
    return;

  di->keypress = keypress;
}

void   SetMouseMotionCB(Widget w, MotionCB motion)
{
  DrawInfo *di;

  if ((di = libsx_find_draw_info(w)) == NULL)
    return;

  di->motion = motion;
}

void SetDrawMode(int mode)
{
  if (display == NULL)
    return;

  if (mode == SANE_XOR)
   {
     cur_di->mask = cur_di->foreground ^ cur_di->background;
     XSetForeground(display, drawgc, 0xffffffff);
     XSetBackground(display, drawgc, cur_di->background);
     XSetFunction(display,   drawgc, GXxor);
     XSetPlaneMask(display,  drawgc, cur_di->mask);
   }
  else
   {
     XSetForeground(display, drawgc, cur_di->foreground);
     XSetBackground(display, drawgc, cur_di->background);
     XSetFunction(display,   drawgc, mode);
     XSetPlaneMask(display,  drawgc, 0xffffffff);
     cur_di->mask = 0xffffffff;
   }
}

void SetLineWidth(int width)
{
  if (display && width > 0)
    XSetLineAttributes(display, drawgc, width, LineSolid, CapButt, JoinMiter);
}

void SetDrawArea(Widget w)
{
  DrawInfo *di;

  if (lsx_curwin->toplevel == NULL || w == NULL)
    return;

  if ((di=libsx_find_draw_info(w)) == NULL)  /* w isn't really a draw area */
    return;

  window  = (Window)XtWindow(w);
  drawgc  = di->drawgc;
  display = XtDisplay(w);
  cur_di  = di;

  lsx_curwin->last_draw_widget = w;

#ifdef    OPENGL_SUPPORT

  if (lsx_curwin->gl_context)
    glXMakeCurrent(display, XtWindow(w), lsx_curwin->gl_context);

#endif /* OPENGL_SUPPORT */
}

#ifdef    OPENGL_SUPPORT

void SwapBuffers(void)
{
  if (lsx_curwin == NULL || lsx_curwin->last_draw_widget == NULL)
    return;

  glXSwapBuffers(display, window);
}

#endif /* OPENGL_SUPPORT */

void GetDrawAreaSize(int *w, int *h)
{
  int n;
  Arg wargs[2];
  Dimension nwidth, nheight;

  if (lsx_curwin->toplevel == NULL || lsx_curwin->last_draw_widget == NULL
      || w == NULL || h == NULL)
    return;

  *w = *h = 0;

  n = 0;
  XtSetArg(wargs[n], XtNwidth,  &nwidth);      n++;
  XtSetArg(wargs[n], XtNheight, &nheight);     n++;

  XtGetValues(lsx_curwin->last_draw_widget, wargs, n);

  *w = nwidth;
  *h = nheight;
}

/*
 * These are the drawing area "draw" functions.  You use these functions
 * to draw into a DrawingArea widget.
 */

void ClearDrawArea(void)
{
  XClearWindow(display, window);
}

void SetColor(int color)
{
  if (cur_di == NULL || display == NULL)
    return;

  cur_di->foreground = color;

  if (cur_di->mask != 0xffffffff)
    XSetPlaneMask(display, drawgc, cur_di->foreground ^ cur_di->background);
  else
    XSetForeground(display, drawgc, color);
}

void DrawPixel(int x1, int y1)
{
  XDrawPoint(display, window, drawgc, x1, y1);
}

int GetPixel(int x1, int y1)
{
  char ch;

  GetImage(&ch, x1, y1, 1, 1);  /* gag! no other easy way to do it */

  return (int)ch;
}

void DrawLine(int x1, int y1, int x2, int y2)
{
  XDrawLine(display, window, drawgc, x1, y1, x2, y2);
}

void DrawPolyline(XPoint *points, int n)
{
  XDrawLines(display, window, drawgc, points, n, CoordModeOrigin);
}

void DrawFilledPolygon (XPoint *points, int n)
{
  XFillPolygon(display, window, drawgc, points, n, Complex, CoordModeOrigin);
}

void DrawFilledBox(int x, int y, int fwidth, int fheight)
{
  if (fwidth < 0)
   { fwidth  *= -1; x -= fwidth; }
  if (fheight < 0)
   { fheight *= -1; y -= fheight; }

  XFillRectangle(display, window, drawgc, x, y, fwidth, fheight);
}

void DrawBox(int x, int y, int bwidth, int bheight)
{
  if (bwidth < 0)
   { bwidth  *= -1; x -= bwidth; }
  if (bheight < 0)
   { bheight *= -1; y -= bheight; }

  XDrawRectangle(display, window, drawgc, x, y, bwidth, bheight);
}

void DrawText(char *string, int x, int y)
{
  XDrawImageString(display, window, drawgc, x, y, string, strlen(string));
}

void DrawArc(int x, int y, int awidth, int aheight, int angle1, int angle2)
{
  angle1 = angle1 * 64;  /* multiply by 64 because X works in 64'ths of a */
  angle2 = angle2 * 64;  /* a degree and we just want to work in degrees  */

  if (awidth < 0)
   { awidth *= -1;  x -= awidth; }
  if (aheight < 0)
   { aheight *= -1; y -= aheight; }

  XDrawArc (display, window, drawgc, x, y,
            awidth, aheight, angle1, angle2);
}

void DrawFilledArc(int x, int y, int awidth, int aheight,
                   int angle1, int angle2)
{
  angle1 = angle1 * 64;  /* multiply by 64 because X works in 64'ths of a */
  angle2 = angle2 * 64;  /* a degree and we just want to work in degrees  */

  if (awidth < 0)
   { awidth *= -1;  x -= awidth; }
  if (aheight < 0)
   { aheight *= -1; y -= aheight; }

  XFillArc (display, window, drawgc, x, y, awidth, aheight, angle1, angle2);
}

void DrawImage(char *data, int x, int y, int width, int height)
{
  XImage *xi;

  if (lsx_curwin->toplevel == NULL || data == NULL)
    return;

  xi = XCreateImage(display, DefaultVisual(display, DefaultScreen(display)),
                    8, ZPixmap, 0, data, width, height,
                    XBitmapPad(display),  width);
  if (xi == NULL)
    return;

  XPutImage(display, window, drawgc, xi, 0,0, x,y,  xi->width,xi->height);

  XFree((char *)xi);
}

void DrawBitmap(char *data, int x, int y, int width, int height)
{
  Pixmap pix;

  if (lsx_curwin->toplevel == NULL || data == NULL)
    return;

  pix = XCreatePixmapFromBitmapData(display, window, data, width, height,
                                    cur_di->foreground, cur_di->background,
                                    DefaultDepth(display, lsx_curwin->screen));
  if (pix == 0)
    return;

  XCopyArea(display, pix, window, drawgc, 0,0, width, height, x,y);

  XFreePixmap(display, pix);
}

/*
 * This function is kind of gaggy in some respects because it
 * winds up requiring twice the amount of memory really needed.
 * It would be possible to return the XImage structure directly,
 * but that kind of defeats the whole purpose of libsx in addition
 * to the fact that X packs the data in ways that might not be
 * what the user wants.  So we unpack the data and just put the
 * raw bytes of the image in the user's buffer.
 */
void GetImage(char *data, int x, int y, int width, int height)
{
  XImage *xi;
  int i,j;
  char *xi_data;

  if (lsx_curwin->toplevel == NULL || data == NULL)
    return;

  xi = XGetImage(display, window, x,y, width,height, ~0, ZPixmap);

  xi_data = xi->data;
  for (i=0; i < height; i++)
   {
     char *line_start = xi_data;

     for (j=0; j < width; j++, xi_data++, data++)
       *data = *xi_data;

     while ((xi_data - line_start) < xi->bytes_per_line)
       xi_data++;
   }

  XFree((char *)xi);
}

/*
 * Below are internal callbacks that the drawing area calls.  They in
 * turn call the user callback functions.
 */

/*
 * Internal callback routines for the drawing areas that massage
 * all the X garbage into a more digestable form.
 */

static void _redisplay(Widget w, void *data, XADCS *call_data)
{
  int new_width, new_height;
  DrawInfo *di = (DrawInfo *) data;

  if (call_data->event->xexpose.count != 0) /* Wait until last expose event */
    return;

  SetDrawArea(w);
  GetDrawAreaSize(&new_width, &new_height);   /* get the draw area size */

  if (di->redisplay)
    di->redisplay(w, new_width, new_height, di->user_data);
}

/* Called when a DrawingArea is resized.
 */
static void _resize(Widget w, void *data, XADCS *call_data)
{
  int new_width, new_height;
  DrawInfo *di = (DrawInfo *) data;

  if (call_data->event->xexpose.count != 0) /* Wait until last expose event */
    return;

  SetDrawArea(w);
  GetDrawAreaSize(&new_width, &new_height);   /* get the new draw area size */

  if (di->redisplay)
    di->redisplay(w, new_width, new_height, di->user_data);
}

/* Called when a DrawingArea has input (either mouse or keyboard).
 */
static void _do_input(Widget w, void *data, XADCS *call_data)
{
  char *input;
  DrawInfo *di = (DrawInfo *) data;

  SetDrawArea(w);
  if (call_data->event->type == ButtonPress)
   {
     if (di->button_down)
       di->button_down(w, call_data->event->xbutton.button,
                       call_data->event->xbutton.x,call_data->event->xbutton.y,
                       di->user_data);
   }
  else if (call_data->event->type == ButtonRelease)
   {
     if (di->button_up)
       di->button_up(w, call_data->event->xbutton.button,
                     call_data->event->xbutton.x, call_data->event->xbutton.y,
                     di->user_data);
   }
  else if (call_data->event->type == KeyPress)
   {
     input = TranslateKeyCode(call_data->event);

     if (input && *input != '\0' && di->keypress)
       di->keypress(w, input, 0, di->user_data);
   }
  else if (call_data->event->type == KeyRelease)
   {
     input = TranslateKeyCode(call_data->event);

     if (input && *input != '\0' && di->keypress)
       di->keypress(w, input, 1, di->user_data);
   }
}

static void _do_motion(Widget w, void *data,  XADCS *call_data)
{
  DrawInfo *di = (DrawInfo *) data;

  SetDrawArea(w);
  if (di->motion)
    di->motion(w, call_data->event->xmotion.x,  call_data->event->xmotion.y,
               di->user_data);
}

#define KEY_BUFF_SIZE 256
static char key_buff[KEY_BUFF_SIZE];

static char *TranslateKeyCode(XEvent *ev)
{
  int count;
  char *tmp;
  KeySym ks;

  if (ev)
   {
     count = XLookupString((XKeyEvent *)ev, key_buff, KEY_BUFF_SIZE, &ks,NULL);
     key_buff[count] = '\0';
     if (count == 0)
      {
        tmp = XKeysymToString(ks);
        if (tmp)
          strcpy(key_buff, tmp);
        else
          strcpy(key_buff, "");
      }

     return key_buff;
   }
  else
    return NULL;
}

#define ABS(x)     ((x < 0) ? -x : x)

#define SWAP(a,b)           \
  {                         \
    a ^= b; b ^= a; a ^= b; \
  }

void ScrollDrawArea(int dx, int dy, int x1, int y1, int x2, int y2)
{
  int w, h, x3, y3, x4, y4, _dx_, _dy_;
  Window win=window;   /* window is a static global */

  if (dx == 0 && dy == 0)
    return;

  if (display == NULL)
    return;

  if (x2 < x1) SWAP (x1,x2);
  if (y2 < y1) SWAP (y1,y2);

  _dx_ = ABS(dx);
  _dy_ = ABS(dy);

  x3 = x1 + _dx_;
  y3 = y1 + _dy_;

  x4 = x2 - _dx_ +1;
  y4 = y2 - _dy_ +1;

  w = x2 - x3 +1;
  h = y2 - y3 +1;

  if (dx <= 0)
   {
     if (dy <= 0)
      {
        XCopyArea (display,win,win,drawgc, x1, y1, w, h, x3, y3);

        if (_dy_)
          XClearArea (display,win, x1, y1, w+_dx_, _dy_, FALSE);

        if (_dx_)
          XClearArea (display,win, x1, y1, _dx_, h, FALSE);

     }
     else              /* dy > 0 */
      {
        XCopyArea (display,win,win,drawgc, x1, y3, w, h, x3, y1);

        XClearArea (display,win, x1, y4, w+_dx_, _dy_, FALSE);

        if (_dx_)
          XClearArea (display,win, x1, y1, _dx_, h, FALSE);
      }
   }
  else                 /* dx > 0 */
   {
     if (dy <= 0)
      {
        XCopyArea (display,win,win,drawgc, x3, y1, w, h, x1, y3);

        if (_dy_)
          XClearArea (display,win, x1, y1, w+_dx_, _dy_, FALSE);

        XClearArea (display,win, x4, y3, _dx_, h, FALSE);
      }
     else              /* dy > 0 */
      {
        XCopyArea (display,win,win, drawgc, x3, y3, w, h, x1, y1);

        XClearArea (display,win, x1, y4, w+_dx_, _dy_, FALSE);

        XClearArea (display,win, x4, y1, _dx_, h, FALSE);
      }
   }
}
