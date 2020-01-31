#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/* callback protos */

#ifndef _CALLBACKS_H_
#define _CALLBACKS_H_

extern void resetWidths(datatype* me);
extern void getridxpix(Box& domain, Real& idxpix,
                       Real xwidth, Real ywidth, void *data, int ilev);
extern void datatypeinit(datatype* dt,int argc, char **argv);
extern void datatypeequate(datatype* dtcopyer, datatype* dtcopyee);
extern void close(Widget w,    void *data);
extern void normaldir(Widget w,    void *data);
extern void slicepos(Widget w,    void *data);
extern void quit(Widget w,    void *data);
extern void subregion(Widget w,    void *data);
extern void keypress(Widget w, char* input, int upordown,  void *data);
extern void drawboxes(Widget w,    void *data);
extern void file_callback(Widget w, char *str, int index, void *arg);
extern void variable_callback(Widget w, char *str, int index, void *arg);
extern void cleanUp(datatype *data);
extern void load(Widget w, void *data);
extern void dumpPS(Widget w, void *data);
extern void makeDataSheet(Widget w, void *data);
extern void do_renorm(Widget w, void *data);
extern void no_renorm(Widget w, void *data);
extern void do_con_plot_cart(Widget w, void *data);
extern void button_down(Widget w, int which_button, int x, int y,void *data);
extern void button_up(Widget w, int which_button, int x, int y,void *data);
extern void getAxes(datatype *data);
extern void makeCFStuff(datatype* data);
extern void readheader(FILE*,int*,int*,int*);
extern void readstuff(FILE * p_file1,int nl0,int nl1
                     ,int nvar,float *p_x);
extern int fillFromFile(datatype* me, char* filename);
extern void redisplay(Widget w, int new_width, int new_height, void *data);
extern void  getCornerLocs(datatype* me,IntVect& iv, Box& domain,
                    Real ridxpix,
                    Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs);
extern void  getPSCorners(datatype* me,IntVect& iv, Box& domain,
                   Real ridxpix, Real imin, Real jmin,
                   Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs);
extern void getDataLocs(datatype* me,IntVect& iv,Box& domain,
                 FArrayBox& bigstate,Tuple<Real,5>& flocs);

/* drawing functions called from redisplay() */
extern void draw_stuff(int width, int height,void *data);
extern void infobfull(int width, int height,void *data);
extern void getVar(Widget w,void *data);
extern void getNumCont(Widget w,void *data);
extern void draw_text(int width, int height,void *data);
extern void renormsca(int width, int height,void *data,float* x);
extern void getMaxMinMag(datatype *data);
extern void renormvec(int width, int height,void *data,float* u,float* v);
extern void con_plot_cart(int width, int height,void *data);
extern void set_col_map(int width, int height,void *data);
extern void fill2d(void *data);
extern void readXDangFAB(datatype* me,istream& is, FArrayBox& fabin);
extern void readSingleFAB(datatype* me,istream& is);
extern int readHierarchy(datatype* me, const string& is);
extern void DrawDataSheet(LevelData<FArrayBox>* newLayoutData);
extern void do_colormap(Widget w,void *data);
extern void setColorMap(void *data);
extern void loadColorMap(Widget w, void *data);
extern void loadGreyscale(Widget w, void *data);
extern void loadColorScale(Widget w, void *data);
extern void readColorsFromFile(datatype* me);
extern void colorMap(int width, int height,void *data);
extern void dumpPSCM(Widget w, void *data);
extern void makecolormapper(void* data);
extern void makecontourplotter(void* data);
extern void drawPalette(int width, int height,void *data);
extern void displayPalette(Widget w, int width, int height, void *data);
extern void dumpPSPalCM(Widget w, void *data);

#endif
