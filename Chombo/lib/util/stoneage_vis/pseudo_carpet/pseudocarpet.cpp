#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include "REAL.H"
#include "Box.H"
#include "CH_assert.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Tuple.H"
#include "BoxIterator.H"
#include "PiecewiseLinearFillPatch.H"
#include "DebugDump.H"
#include "AMRIO.H"
using std::string;
#define WINSIZE   1024
#define NUMCOLORS 16384

class PlotParams
{
public:
  PlotParams()
  {
    useGray = false;
    drawBoxes = true;
    cur_var  = 0;
    normaldir = 2;
  }
  string filename;
  string outputname;
  bool useGray;
  bool drawBoxes;
  int cur_var ;
  Box subbox;
  int normaldir;
  int slicePos;
  string varname;
};                  

class GlobalDat
{
public:
  PlotParams params;
  Vector<DisjointBoxLayout>        vect_grids;
  Vector<LevelData<FArrayBox>* >   vect_mf;
  Box                              lev0_domain;
  Vector<int>                      vect_ratio;
  int                              numLevels;
  GLint windW;
  GLint windH;
  Real scale;
  Real shiftX, shiftY, shiftZ;
  Real angleX, angleY, angleZ;
  Tuple<int,2> axisdir;
  Real max, min, mag, eps;
  int rcol[NUMCOLORS];
  int gcol[NUMCOLORS];
  int bcol[NUMCOLORS];
  bool carpet;
  GlobalDat()
  {
    carpet = false;
    windW = WINSIZE;
    windH = WINSIZE;
    scale = 1.0;
    shiftX = 0.0; 
    shiftY = 0.0; 
    shiftZ = 0.0;
    angleX = 0.0; 
    angleY = 0.0; 
    angleZ = 0.0;
  }
};

GlobalDat fnerg;
void getMaxMinMag();
void getAxes();
void fillGhost()
{
  int nlev = fnerg.vect_mf.size();
  Box domlev = fnerg.lev0_domain;
  for(int ilev = 0; ilev < nlev ; ilev++)
    {
      bool neednewdata = false;
      IntVect ivghost = fnerg.vect_mf[ilev]->ghostVect();
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(ivghost[idir] < 1) neednewdata = true;
        }
      DisjointBoxLayout dbl = fnerg.vect_mf[ilev]->disjointBoxLayout();
      int ncomp =  (fnerg.vect_mf[ilev])->nComp();
      Interval interv(0, ncomp-1);
      if(neednewdata)
        {
          LevelData<FArrayBox>* oldData = fnerg.vect_mf[ilev];
          LevelData<FArrayBox>* newData = new LevelData<FArrayBox>(dbl, ncomp, IntVect::Unit);
          oldData->copyTo(interv, *newData, interv);
          fnerg.vect_mf[ilev] =  newData;
          //         delete oldData;
        }
      if(ilev > 0)
        {
          PiecewiseLinearFillPatch patcher(fnerg.vect_mf[ ilev  ]->disjointBoxLayout(),
                                           fnerg.vect_mf[ ilev-1]->disjointBoxLayout(),  ncomp,
                                           ProblemDomain(domlev),
                                           fnerg.vect_ratio[ilev-1], 1);
          patcher.fillInterp(*fnerg.vect_mf[ilev],*fnerg.vect_mf[ilev-1],*fnerg.vect_mf[ilev-1],1,0,0,ncomp);
        }
      domlev.refine(fnerg.vect_ratio[ilev]);
    }
}
/*************************************************/
void ButtonDown(int x,int y)
{
  //find highest level which contains location
  /* draw stuff */
  int papersize = WINSIZE; //ortho
  int imin = 0;
  int jmin = 0;
  int imax = papersize;
  int jmax = papersize;
  Real rimin = imin;
  Real rjmin = jmin;
  Real rimax = imax;
  Real rwidth  = rimax - rimin;
  Box domain = fnerg.lev0_domain;
  Real dx = rwidth/domain.size(0);
  Real dy = rwidth/domain.size(1);
  if(dx > dy) dx = dy;
  if(dy > dx) dy = dx;
  int slicePos = fnerg.params.slicePos;


  getMaxMinMag();
  getAxes();
  Real dxlev = dx;
  Box domlev = domain;
  int slicelev = slicePos;
  for(int ilev = 0; ilev < fnerg.numLevels; ilev++)
    {
      IntVect iv = slicelev*IntVect::Unit;
      iv[fnerg.axisdir[0]] = x/dxlev;
      iv[fnerg.axisdir[1]] = domlev.bigEnd()[1]- int(y/dxlev);
      for(DataIterator dit = fnerg.vect_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          Box curBox = fnerg.vect_grids[ilev][dit()];
          if(curBox.contains(iv))
            {
              pout() <<  "on level " << ilev << ", you clicked on iv = " << iv << endl;
              LevelData<FArrayBox>& dat = *fnerg.vect_mf[ilev];
              Box localBox(iv, iv);
              localBox.grow(1);
              if(SpaceDim == 3)
                localBox.grow(fnerg.params.normaldir, -1);
              localBox &= curBox;
              FArrayBox localFAB(localBox, dat.nComp());
              localFAB.copy(dat[dit()]);
              dumpFAB(&localFAB);
            }
        }
      slicelev *= fnerg.vect_ratio[ilev];
      dxlev /= fnerg.vect_ratio[ilev];
      domlev.refine(fnerg.vect_ratio[ilev]);
    }
}

/*************************************************/
void Mouse(int button,int state,int x,int y)
{
  //can also swith on GLUT_LEFT_BUTTON and so on
  if(state==GLUT_DOWN)
    ButtonDown(x,y);
  
}
/*************************************************/
void getAxes()
{
#if(CH_SPACEDIM == 2)
  fnerg.axisdir[0] = 0;
  fnerg.axisdir[1] = 1;
#elif(CH_SPACEDIM == 3)
  int idir1, idir2;
  CH_assert(fnerg.params.normaldir >= 0 && fnerg.params.normaldir < SpaceDim);
  if(fnerg.params.normaldir == 2)
    {
      idir1 = 0;
      idir2 = 1;
    }
  else if(fnerg.params.normaldir == 1)
    {
      idir1 = 0;
      idir2 = 2;
    }
  else if(fnerg.params.normaldir == 0)
    {
      idir1 = 1;
      idir2 = 2;
    }
  fnerg.axisdir[0] = idir1;
  fnerg.axisdir[1] = idir2;

#else
#error unimplemented CH_SPACEDIM
#endif
}
/********************************/
void getMaxMinMag()
{
  Real fmax = -7.77e10;
  Real fmin =  7.77e10;
  Real fmag = -7.77e10;
  //find max, min, mag

  for(int  ilev = 0; ilev < int(fnerg.vect_mf.size()); ilev++)
    {
      const LevelData<FArrayBox>& mfcur = *(fnerg.vect_mf[ilev]);
      DataIterator dit = mfcur.dataIterator();
      DisjointBoxLayout dbl = mfcur.disjointBoxLayout();
      for(dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& fabcur = mfcur[dit()];
          const Box& gridbox = dbl.get(dit());
          IntVect loiv = gridbox.smallEnd();
          IntVect hiiv = gridbox.bigEnd();
          if(SpaceDim==3)
            {
              loiv[fnerg.params.normaldir] = fnerg.params.slicePos;
              hiiv[fnerg.params.normaldir] = fnerg.params.slicePos;
            }
          Box subbox(loiv, hiiv);
          subbox &= gridbox;
          if(!subbox.isEmpty())
            {
              fmax = Max(fmax, fabcur.max(subbox, fnerg.params.cur_var));
              fmin = Min(fmin, fabcur.min(subbox, fnerg.params.cur_var));
            }
        }
    }
  fmag = fmax - fmin;
  fnerg.eps = 1.0e-10;
  if(fmag < fnerg.eps)
    {
      fmag = 1.0;
    }
  fnerg.max =fmax;
  fnerg.min =fmin;
  fnerg.mag =fmag;
  pout() << "max = " << fmax << ", min = " << fmin << ", mag = " << fmag << endl;
}
////////////
void
parseCommandLine(PlotParams& a_params, int argc, char* argv[])
{
  pout() << "usage: plotFromHDF5 -i filein -b on/off -gray on/off -iv varnum -nd normaldir -d slicePos " << endl;
  bool filefound = false;
  for(int iarg = 0; iarg < argc-1; iarg++)
    {
      if(strcmp(argv[iarg],"-i") == 0)
        {
          a_params.filename = string(argv[iarg+1]);
          filefound = true;
        }
      else if(strcmp(argv[iarg], "-o") == 0)
        {
          a_params.outputname = string(argv[iarg+1]);
        }
      else if(strcmp(argv[iarg],"-v") == 0)
        {
          a_params.varname = string(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-d") == 0)
        {
          a_params.slicePos = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-nd") == 0)
        {
          a_params.normaldir = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-iv") == 0)
        {
          a_params.cur_var = atoi(argv[iarg+1]);
        }
      else  if(strcmp(argv[iarg],"-b") == 0)
        {
          a_params.drawBoxes = (strcmp(argv[iarg+1],"off") != 0);
        }
      else  if(strcmp(argv[iarg], "-g") == 0)
        {
          a_params.useGray = (strcmp(argv[iarg+1],"on") == 0);
        }
    }
  if(!filefound)
    {                                                  // 
      pout() << "no filename input"  << endl;
      exit(0);
    }
}

/**these three never have to change*/
static void Init(void);
static void Draw(void);
static void Key(unsigned char key, int x, int y);
void
drawRect(float squad[4][3], float rr, float gg, float bb)
{
  int rcolnum = NUMCOLORS-1;
  float xpt[3];
  glColor3f(rr/rcolnum,gg/rcolnum,bb/rcolnum);
  
  glBegin(GL_POLYGON);
  xpt[2] = squad[0][2];    xpt[0] = squad[0][0];    xpt[1] = squad[0][1]; 
  glVertex3fv(xpt);
  xpt[2] = squad[1][2];    xpt[0] = squad[1][0];    xpt[1] = squad[1][1];; 
  glVertex3fv(xpt);
  xpt[2] = squad[2][2];    xpt[0] = squad[2][0];    xpt[1] = squad[2][1];
  glVertex3fv(xpt);
  xpt[2] = squad[3][2];    xpt[0] = squad[3][0];    xpt[1] = squad[3][1];
  glVertex3fv(xpt);
  glEnd();

}
void increment(int& count, Real& value, FArrayBox& data, IntVect& iv)
{
  if(data.box().contains(iv))
    {
      count++;
      Real datVal = data(iv, fnerg.params.cur_var);
      if(fnerg.mag > 1.0e-8)
        {
          value += 0.25*WINSIZE*(datVal-fnerg.min)/fnerg.mag;
        }
    }
}

Real getZVal(FArrayBox& data, IntVect iv0,IntVect  iv1,IntVect  iv2,IntVect  iv3)
{
  int  numVal = 0;
  Real datVal = 0;
  increment(numVal, datVal, data, iv0);
  increment(numVal, datVal, data, iv1);
  increment(numVal, datVal, data, iv2);
  increment(numVal, datVal, data, iv3);

  if(numVal > 1) datVal /= numVal;
  return datVal;
}
void
drawBox(const IntVect& ivdomlo,const IntVect& ivdat,const Real& dx,
        const Real& rcolor, const Real& gcolor, const Real& bcolor, int ilev,
        FArrayBox& data)
{
  //  Box domain(IntVect::Zero, 63*IntVect::Unit);
  float squad[4][3];
//  Real bsizex = Real(domain.size(0));
//  Real bsizey = Real(domain.size(1));
//
//  Real dx = Real(WINSIZE)/bsizex;
//  Real dy = Real(WINSIZE)/bsizey;
  IntVect ivdiff = ivdat - ivdomlo;
  Real x1,x2, y1, y2;
  Real zval = dx*ilev;
  x1 = dx*ivdiff[fnerg.axisdir[0]];
  x2 = x1 + dx;
  y1 = dx*ivdiff[fnerg.axisdir[1]];
  y2 = y1 + dx;
  squad[0][0] = x1;
  squad[0][1] = y1;
  squad[1][0] = x2;
  squad[1][1] = y1;
  squad[2][0] = x2; 
  squad[2][1] = y2;
  squad[3][0] = x1;
  squad[3][1] = y2;
  drawRect(squad, rcolor, gcolor, bcolor);
  if(!fnerg.carpet)
    {    
      squad[0][2] = zval;
      squad[1][2] = zval;
      squad[2][2] = zval;
      squad[3][2] = zval;
    }
  else
    {
      IntVect iv0 = BASISV(fnerg.axisdir[0]);
      IntVect iv1 = BASISV(fnerg.axisdir[1]);
      IntVect iv = ivdat;
      squad[0][2] = zval + getZVal(data, iv, iv-iv0, iv-iv1, iv-iv0-iv1);
      squad[1][2] = zval + getZVal(data, iv, iv+iv0, iv-iv1, iv+iv0-iv1);
      squad[2][2] = zval + getZVal(data, iv, iv+iv0, iv+iv1, iv+iv0+iv1);
      squad[3][2] = zval + getZVal(data, iv, iv-iv0, iv+iv1, iv-iv0+iv1);
    }
  drawRect(squad, rcolor, gcolor, bcolor);
}



#include "AMRIO.H"
#include <stdlib.h>
int main(int argc, char *argv[])
{
  parseCommandLine(fnerg.params, argc, argv);
  
  ReadAMRHierarchyHDF5(fnerg.params.filename, fnerg.vect_grids, 
                       fnerg.vect_mf, fnerg.lev0_domain, fnerg.vect_ratio, 
                       fnerg.numLevels);
  fillGhost();
  // argv++;
  // sscanf(*argv, "%d", &it);
  glutInit(&argc, argv);

  printf("key commands \n");
  printf("esc to exit  \n");
  printf("x or X moves in x  \n");
  printf("y or Y moves in y  \n");
  printf("z or Z moves in z  \n");
  printf("v or V decrements or increments varialbe  \n");
  printf("g toggles grey scale  \n");
  printf("c toggles carpet plot  \n");
  printf("0 makes the normal x\n");
  printf("1 makes the normal y\n");
  printf("2 makes the normal z\n");
  printf("r or R rotates in x  \n");
  printf("t or T rotates in y  \n");
  printf("q or Q rotates in z  \n");
  printf("P or p changes the slice position \n");

  printf("s or S scales higher or lower (zooms out or in)  \n");

  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowPosition(5, 30);
  glutInitWindowSize(fnerg.windW, fnerg.windH);
  glutCreateWindow("opendang");

  glutKeyboardFunc(Key);
  glutMouseFunc(Mouse);
  glutDisplayFunc(Draw);
  Init(); 
  glutMainLoop();
  return 0;
}
void HSVtoRGB( Real *r, Real *g, Real *b, Real hh, Real ss, Real vv )
{
  int i;
  Real f,p,q,t,fi;
  if (hh < -.5){
    *r = vv;
    *g = vv;
    *b = vv;
  }
  else{
    i = int(floor(hh));
    fi = i;
    f = hh - fi;
    p = vv*(1. - ss);
    q = vv*(1. - (ss*f));
    t = vv*(1. - ss*(1.-f));
    if (i == 0){
      *r = vv;
      *g = t;
      *b = p;
    }
    if (i == 1){
      *r = q;
      *g = vv;
      *b = p;
    }
    if (i == 2){
      *r = p;
      *g = vv;
      *b = t;
    }
    if (i == 3){
      *r = p;
      *g = q;
      *b = vv;
    }
    if (i == 4){
      *r = t;
      *g = p;
      *b = vv;
    }
    if (i == 5){
      *r = vv;
      *g = p;
      *b = q;
    }
  }
}
/*****************/
void setRainbowCM()
{
  for(int i = 0; i< NUMCOLORS;i++)
    {
      Real frac = Real(i)/Real(NUMCOLORS);;
      Real sat = 1.;
      Real val = 1;
      Real hue = 6*frac;
      Real rnumcol = NUMCOLORS-1;
      Real r, g, b;
      HSVtoRGB(&r, &g, &b, hue, sat, val);
      fnerg.rcol[i] = rnumcol*r;
      fnerg.gcol[i] = rnumcol*g;
      fnerg.bcol[i] = rnumcol*b;

    }
}
/*****************/
void setGrayCM()
{
  for(int i = 0; i< NUMCOLORS;i++)
    {
      fnerg.rcol[i] = i;
      fnerg.gcol[i] = i;
      fnerg.bcol[i] = i;
    }
}
/******/
void pseudoColorPlot()
{

  /* draw stuff */
  int papersize = WINSIZE; //ortho
  int imin = 0;
  int jmin = 0;
  int imax = papersize;
  int jmax = papersize;
  Real rimin = imin;
  Real rjmin = jmin;
  Real rimax = imax;
  Real rwidth  = rimax - rimin;

  /* plot */
  /* upper left corner is origin---ergo weirdness */

  getMaxMinMag();
  getAxes();
  if(fnerg.params.useGray)
    setGrayCM();
  else
    setRainbowCM();

  int slicePos = fnerg.params.slicePos;

  Box domain = fnerg.lev0_domain;
  Real dx = rwidth/domain.size(0);
  Real dy = rwidth/domain.size(1);
  if(dx > dy) dx = dy;
  if(dy > dx) dy = dx;
  for(int  ilev = 0; ilev < fnerg.numLevels; ilev++)
    {
      IntVect ivdomlo = domain.smallEnd();
      LevelData<FArrayBox>& mfcur = *fnerg.vect_mf[ilev];
      DisjointBoxLayout dbl = fnerg.vect_grids[ilev];
      DataIterator dit = dbl.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& bigfab = mfcur[dit()];
          Box bigbox = bigfab.box();
          Box curbox = dbl.get(dit());

#if(CH_SPACEDIM == 3)
          Box threedbox = curbox;
          IntVect loiv = threedbox.smallEnd();
          IntVect hiiv = threedbox.bigEnd();
          loiv[fnerg.params.normaldir]= slicePos;
          hiiv[fnerg.params.normaldir]= slicePos;
          curbox= Box(loiv,hiiv);
          curbox &= threedbox;
#endif
          curbox &= domain;
          if(!curbox.isEmpty())
            {
              FArrayBox smallfab(curbox, 1);
              smallfab.copy(bigfab, fnerg.params.cur_var, 0, 1);
              //stuff is not normalized to one
              for(BoxIterator bit(curbox); bit.ok(); bit.next())
                {
                  IntVect ivint = bit();

                  Real fval = smallfab(ivint, 0);
                  Real fmag = fnerg.mag;
                  Real fmin = fnerg.min;
                  Real fscale;
                  if(Abs(fmag) > 1.0e-10) 
                    fscale = (fval-fmin)/fmag;
                  else 
                    fscale = 0;
                  
                  int ir, ig, ib;
                  int icolor = (int)(NUMCOLORS*fscale);
                  if(icolor>=NUMCOLORS) icolor = NUMCOLORS-1;
                  if(icolor <0 ) icolor = 0;

                  ir = fnerg.rcol[icolor];
                  ig = fnerg.gcol[icolor];
                  ib = fnerg.bcol[icolor];

                  Real rcolor = ((Real) ir);
                  Real gcolor = ((Real) ig);
                  Real bcolor = ((Real) ib);
                  //                  pout() << "icolor r g b" << icolor  << " "<< ir << " " << ig << " " << ib << endl;
                  drawBox(ivdomlo, ivint, dx, rcolor, gcolor, bcolor, ilev,  bigfab);

                }
            } //if curbox is not empty
        } //loop over boxes

      slicePos *= fnerg.vect_ratio[ilev];
      dx /= fnerg.vect_ratio[ilev];
      dy /= fnerg.vect_ratio[ilev];
      domain.refine(fnerg.vect_ratio[ilev]);
    }//loop over levels

}
/******/
void Key(unsigned char key, int x, int y)
{
  int maxslice = 0;
  if(SpaceDim==3)
    {
      maxslice = fnerg.lev0_domain.size(fnerg.params.normaldir) - 1;
    }

  switch (key) 
    {
    case 27:
      exit(0);

    case 'x':
      fnerg.shiftX -= 20.0;
      break;
    case 'X':
      fnerg.shiftX += 20.0;
      break;
    case 'r':
      fnerg.angleX -= 2.0;
      if (fnerg.angleX < 0.0) 
        {
          fnerg.angleX = 360.0 + fnerg.angleX;
        }
      break;
    case 'R':
      fnerg.angleX += 2.0;
      if (fnerg.angleX > 360.0) 
        {
          fnerg.angleX = fnerg.angleX - 360.0;
        }
      break;
    case 't':
      fnerg.angleY -= 5.0;
      if (fnerg.angleY < 0.0) 
        {
          fnerg.angleY = 360.0 + fnerg.angleY;
        }
      break;
    case 'T':
      fnerg.angleY += 5.0;
      if (fnerg.angleY > 360.0) 
        {
          fnerg.angleY = fnerg.angleY - 360.0;
        }
      break;
    case 'q':
      fnerg.angleZ -= 5.0;
      if (fnerg.angleZ < 0.0) 
        {
          fnerg.angleZ = 360.0 + fnerg.angleZ;
        }
      break;
    case 'Q':
      fnerg.angleZ += 5.0;
      if (fnerg.angleZ > 360.0) 
        {
          fnerg.angleZ = fnerg.angleZ - 360.0;
        }
      break;
    case 'c':
      fnerg.carpet  = !fnerg.carpet;
      break;
    case 'y':
      fnerg.shiftY += 20.0;
      break;
    case 'Y':
      fnerg.shiftY -= 20.0;
      break;
    case 'v':
      fnerg.params.cur_var--;
      if(fnerg.params.cur_var < 0) fnerg.params.cur_var = 0;
      pout() << "current variable is now " << fnerg.params.cur_var<< endl;
      break;
    case 'V':
      fnerg.params.cur_var++;
      if(fnerg.params.cur_var >= fnerg.vect_mf[0]->nComp()) fnerg.params.cur_var = 0;
      pout() << "current variable is now " << fnerg.params.cur_var<< endl;
      break;
    case 'p':
      fnerg.params.slicePos--;
      if(fnerg.params.slicePos < 0) fnerg.params.slicePos = 0;
      pout() << "current slice is now " << fnerg.params.slicePos<< endl;
      break;
    case 'P':
      fnerg.params.slicePos++;
      if(fnerg.params.slicePos > maxslice) fnerg.params.slicePos = maxslice;
      pout() << "current slice is now " << fnerg.params.slicePos<< endl;
      break;
    case 'z':
      fnerg.shiftZ += 20.0;
      break;
    case 'Z':
      fnerg.shiftZ -= 20.0;
      break;
    case 'g':
      fnerg.params.useGray = !fnerg.params.useGray;
      break;
    case '0':
      fnerg.params.normaldir = 0;
      break;
    case '1':
      fnerg.params.normaldir = 1;
      break;
    case '2':
      fnerg.params.normaldir = 2;
      break;
    case 's':
      fnerg.scale -= 0.1;
      if (fnerg.scale < 0.1) 
        {
          fnerg.scale = 0.1;
        }
      break;
    case 'S':
      fnerg.scale += 0.1;
      break;
    }
  glutPostRedisplay();
}

void Draw(void)
{

  glClear( GL_COLOR_BUFFER_BIT | 
           GL_DEPTH_BUFFER_BIT | 
           GL_STENCIL_BUFFER_BIT);

  glViewport(0, 0, fnerg.windW, fnerg.windH);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glScissor(0, 0, fnerg.windW, fnerg.windH);

  glPushMatrix();

  //now in the rotated coordinate frame
  glTranslatef(fnerg.shiftX, fnerg.shiftY, fnerg.shiftZ);
  glRotatef(fnerg.angleX, 1.0, 0.0, 0.0);
  glRotatef(fnerg.angleY, 0.0, 1.0, 0.0);
  glRotatef(fnerg.angleZ, 0.0, 0.0, 1.0);
  glScalef(fnerg.scale, fnerg.scale, fnerg.scale);


  pseudoColorPlot();

  glPopMatrix();
  glPopMatrix();


  glFlush();
  //because of double buffering
  glutSwapBuffers();
}
void Init(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  Real winsize = WINSIZE;
  glOrtho(0.,winsize, 0.,winsize, -winsize, winsize);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glDepthMask(GL_FALSE);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
}

