#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctype.h>
#include <cmath>
#include <string>
#include "SPACE.H"

#include "REAL.H"
#include "PostScriptPlot.H"
#include "LoHiSide.H"
#include "DebugDump.H"
#include "AMRIO.H"
#include "PiecewiseLinearFillPatch.H"
using std::string;

int g_xsize = 1024;
int g_ysize = 1024;
const int g_infob = 125;
const int g_paletteb = 200;
const int g_numcont = 30;
const int g_minlen = 200;

/*
  So the following function goes
  from HSV values, to RGB.  "floor" is a function from the math library,
  that just removes any decimal part of a positive number (i.e. floor(5.3)
  = 5, floor(3.9) = 3).  Note that 0.<=hh<=6., 0.<=ss<=1., 0.<=vv<=1., and
  r, g, and b range from 0. to 255. (so you'll have to do a little extra
  adjustments to these numbers, if your HSV and RGB ranges are something
  else.)
*/
void
maxxx(Real a,Real b,Real c,Real *mx)
{
  if ((a >= b) && (a >= c))
    *mx = a;
  else if ((b >= a) && (b >= c))
    *mx = b;
  else
    *mx = c;
}
void
minnn(Real a,Real b,Real c,Real* mn)
{
  if ((a <= b) && (a <= c))
    *mn = a;
  else if ((b <= a) && (b <= c))
    *mn = b;
  else
    *mn = c;
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

/*Now my "colget" function as you see, goes the other way.  It calls 2
  other functions I wrote, "maxxx" and "minnn", which obtain the maximum
  and minimum respectively from among 3 numbers.  "fabs" is another math
  library function, that returns the absolute value of a Real.  Love.
*/
void RGBtoHSV( Real rr, Real gg, Real bb, Real* h, Real *s, Real *v )
{
  Real  r = rr/255.;
  Real  g = gg/255.;
  Real  b = bb/255.;
  Real maxx, minn;
  maxxx(r,g,b,&maxx);
  minnn(r,g,b,&minn);
  *v = maxx;
  if (fabs(maxx) >= 10e-10)
    *s = (maxx - minn)/maxx;
  else
    *s = 0.;
  Real delta = maxx - minn;
  if (fabs(delta) >= 10e-10){
    if (r == maxx)
      *h = (g - b)/delta;
    else if (g == maxx)
      *h = 2. + (b - r)/delta;
    else
      *h = 4. + (r - g)/delta;
    if (*h < 0.)
      *h = *h + 6.;
    if(*h > 6.)
      *h = *h - 6.;
  }
  else {
    *h = -1.;
    printf("D'OH!\n");
  }
}
/*****************/
void setRainbowCM(PostScriptData *me)
{
  for(int i = 0; i< 256;i++)
    {
      Real frac = Real(i)/Real(256);;
      Real sat = 1.;
      Real val = 1;
      Real hue = 6*frac;

      Real r, g, b;
      HSVtoRGB(&r, &g, &b, hue, sat, val);
      me->rcol[i] = 255*r;
      me->gcol[i] = 255*g;
      me->bcol[i] = 255*b;

    }
}
/*****************/
void setRedGreenCM(PostScriptData *me)
{
  for(int i = 0; i< 256;i++)
    {
      me->rcol[i] = 255-i;
      me->gcol[i] = i;
      me->bcol[i] = 0;

    }
}

/////
void
PostScriptPlot::
getCoveredCells(Vector<IntVectSet>& a_ivs,
                Vector<LevelData<EBCellFAB>* >& a_data)
{
  for(int ilev = 0; ilev < a_data.size(); ilev++)
    {
      a_ivs[ilev] = IntVectSet();
      for(DataIterator dit = a_data[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = (*a_data[ilev])[dit()].getEBISBox();
          const Box& box = a_data[ilev]->disjointBoxLayout()[dit()];
          for(BoxIterator bit(box); bit.ok(); ++bit)
            {
              if(ebisBox.isCovered(bit()))
                a_ivs[ilev] |= bit();
            }
        }
    }
}
/////
PostScriptData::
PostScriptData(Vector< LevelData<FArrayBox>* >& a_data,
               Vector<int>                    & a_refRat,
               const string                   & a_filename,
               int                            & a_ivar,
               std::string                    & a_varname,
               bool                           & a_drawBoxes,
               int                            & a_normalDir,
               int                            & a_slicePosition,
               int                            & a_numContours,
               bool                           & a_useGray,
               Box                            & a_subBox,
               bool                           & a_ebflag,
               Vector<IntVectSet>             & a_coveredCells
               )
{
  vect_mf = a_data;
  vect_ratio = a_refRat;

  vect_dbl.resize(a_data.size());
  vect_box.resize(a_data.size());
  for(int ilev =  0; ilev < a_data.size(); ilev++)
    {
      vect_dbl[ilev] = vect_mf[ilev]->disjointBoxLayout();
      vect_box[ilev] = vect_mf[ilev]->disjointBoxLayout().physDomain().domainBox();
    }

  varname  = a_varname;
  filename = a_filename;
  palettename = string("palette.");
  palettename += filename;

  cur_var = a_ivar;
  xdraw = g_xsize;
  ydraw = g_ysize;
  drawboxes = a_drawBoxes;
  inormal = a_normalDir;
  idepth  = a_slicePosition;
  numcont = a_numContours;
  lev0subbox = a_subBox;
  doingsubregion = !(lev0subbox.isEmpty());
  if(!doingsubregion)
    {
      lev0subbox = vect_box[0];
    }
  if(a_useGray)
    {
      setGrayCM(this);
    }
  else
    {
      setRainbowCM(this);
    }
  ebflag = a_ebflag;
  coveredCells = a_coveredCells;

}
/*****************/
void setGrayCM(PostScriptData *me)
{
  for(int i = 0; i< 256;i++)
    {
      me->rcol[i] = i;
      me->gcol[i] = i;
      me->bcol[i] = i;

    }
}

void doCFStuff(PostScriptData *me)
{
  Vector<LevelData<FArrayBox> *>& data = me->vect_mf;
  int ncomps = data[0]->nComp();
  int interpRad = 1;
  data[0]->exchange ();
  for(int  ilev = 1; ilev < int(data.size()); ilev++)
    {
      PiecewiseLinearFillPatch filler(data[ilev]->disjointBoxLayout(),
                                      data[ilev-1]->disjointBoxLayout(),
                                      ncomps,
                                      me->vect_box[ilev-1],
                                      me->vect_ratio[ilev-1],
                                      interpRad);

      LevelData<FArrayBox>& mfcur = *data[ilev];
      LevelData<FArrayBox>& uccur = *data[ilev-1];

      Real timeCoeff = 0.5;
      filler.fillInterp(mfcur, uccur, uccur,
                        timeCoeff, 0, 0, mfcur.nComp());
      mfcur.exchange(mfcur.interval());
    }
}

/*************************************************/
void getAxes(PostScriptData *me)
{
  //do spacedim-dependent stuff
  Box bbase = me->lev0subbox;

#if(CH_SPACEDIM == 2)
  me->axisdir[0] = 0;
  me->axisdir[1] = 1;
#elif(CH_SPACEDIM == 3)
  int idir1, idir2;
  CH_assert(me->inormal >= 0 && me->inormal < SpaceDim);
  if(me->inormal == 2)
    {
      idir1 = 0;
      idir2 = 1;
    }
  else if(me->inormal == 1)
    {
      idir1 = 0;
      idir2 = 2;
    }
  else if(me->inormal == 0)
    {
      idir1 = 1;
      idir2 = 2;
    }
  me->axisdir[0] = idir1;
  me->axisdir[1] = idir2;
#else
#error unimplemented CH_SPACEDIM
#endif
}
/********************************/
void getMaxMinMag(PostScriptData *me)
{
  Real fmax = -7.77e10;
  Real fmin =  7.77e10;
  Real fmag = -7.77e10;
  //find max, min, mag

  for(int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
    {
      const LevelData<FArrayBox>& mfcur = *(me->vect_mf[ilev]);
      DataIterator dit = mfcur.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
        {
          const FArrayBox& fabcur = mfcur[dit()];
          Box subbox = grow(fabcur.box(), -1);
          fmax = Max(fmax, fabcur.max(subbox, me->cur_var));
          fmin = Min(fmin, fabcur.min(subbox, me->cur_var));
        }
    }
  fmag = fmax - fmin;
  me->eps = 1.0e-10;
  if(fmag < me->eps)
    {
      fmag = 1.0;
    }
  me->max =fmax;
  me->min =fmin;
  me->mag =fmag;
  pout() << "max = " << fmax << ", min = " << fmin << ", mag = " << fmag << endl;
}
/* useful function which linearly interpolates between v1 and v2
   (which exist at T1, T2) to valmid (which lives at time)
   --- returns valmid*/
Real spainterp(Real t1, Real t2, Real t, Real v1, Real v2)
{
  Real TCOld;
  Real TCNew;
  Real Time;
  Real valold;
  Real valnew;
  if(t2 > t1)
    {
      TCOld = t1;
      TCNew = t2;
      valold = v1;
      valnew = v2;
      Time = t;
    }
  else
    {
      TCOld = t2;
      TCNew = t1;
      valold = v2;
      valnew = v1;
      Time = t;
    }
  Real tcdiff = 1.0;
  Real tfdiff = 0.0;
  if((TCNew - TCOld) > 1.0e-12)
    {
      tcdiff = TCNew -TCOld;
      tfdiff = Time -TCOld;
    }
  Real valmid = valold + tfdiff*(valnew-valold)/tcdiff;
  return valmid;

}
void fillGhost(PostScriptData *me)
{
  int nlev = me->vect_mf.size();
  for(int ilev = 0; ilev < nlev ; ilev++)
    {
      bool neednewdata = false;
      IntVect ivghost = me->vect_mf[ilev]->ghostVect();
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(ivghost[idir] < 1) neednewdata = true;
        }
      DisjointBoxLayout dbl = me->vect_mf[ilev]->disjointBoxLayout();
      int ncomp =  (me->vect_mf[ilev])->nComp();
      Interval interv(0, ncomp-1);
      if(neednewdata)
        {
          LevelData<FArrayBox>* oldData = me->vect_mf[ilev];
          LevelData<FArrayBox>* newData = new LevelData<FArrayBox>(dbl, ncomp, IntVect::Unit);
          oldData->copyTo(interv, *newData, interv);
          me->vect_mf[ilev] =  newData;
          //         delete oldData;
        }
      if(ilev > 0)
        {
          PiecewiseLinearFillPatch patcher(me->vect_mf[ ilev  ]->disjointBoxLayout(),
                                           me->vect_mf[ ilev-1]->disjointBoxLayout(),  ncomp,
                                           ProblemDomain(me->vect_box[ilev-1]),
                                           me->vect_ratio[ilev-1], 1);
          patcher.fillInterp(*me->vect_mf[ilev],*me->vect_mf[ilev-1],*me->vect_mf[ilev-1],1,0,0,ncomp);
        }
    }
}
/***********************************************/
void dumpContourPlot(PostScriptData *me)
{
  ofstream ps_strm(me->filename.c_str());
  fillGhost(me);

  /* draw stuff */
  /* obligatory postscript begining */
  int papersize = 8*72;
  int imin = 0;
  int jmin = 0;
  int imax = papersize;
  int jmax = papersize;
  Real rimin = imin;
  Real rjmin = jmin;
  Real rimax = imax;
  Real rwidth  = rimax - rimin;


  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin << " " << imax << " " << jmax << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/st {stroke} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rs {rectstroke} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;
  ps_strm << "/cp {closepath} def" << endl;

  /* draw box around plot */
  if(!me->doingsubregion)
    {
      ps_strm << "3 setlinewidth" << endl;
      ps_strm << imin << " " << jmin << " mt" << endl;
      ps_strm << imin << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmin << " lt" << endl;
      ps_strm << "cp" << endl << "st" << endl;
    }

  /* plot */
  /* upper left corner is origin---ergo weirdness */

  getMaxMinMag(me);
  getAxes(me);
  int idepth = me->idepth;
  Box domain = me->vect_box[0];
  if(me->doingsubregion) domain = me->lev0subbox;
  Real dx = rwidth/domain.size(0);
  Real dy = rwidth/domain.size(1);
  if(dx > dy) dx = dy;
  if(dy > dx) dy = dx;
  for(int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
    {


      LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
      DisjointBoxLayout dbl = mfcur.disjointBoxLayout();
      DataIterator dit = mfcur.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& bigfab = mfcur[dit()];
          Box bigbox = bigfab.box();
          Box curbox = dbl.get(dit());
#if(CH_SPACEDIM == 3)
          Box threedbox = curbox;
          IntVect loiv = threedbox.smallEnd();
          IntVect hiiv = threedbox.bigEnd();
          loiv[me->inormal]= idepth;
          hiiv[me->inormal]= idepth;
          curbox= Box(loiv,hiiv);
          curbox &= threedbox;
#endif
          curbox &= domain;

          if(!curbox.isEmpty())
            {

              //clear area in which to do contours
              IntVect ivcurlo = curbox.smallEnd();
              IntVect ivcurhi = curbox.bigEnd();
              IntVect ivdomlo = domain.smallEnd();
              ivcurlo -= ivdomlo;
              ivcurhi -= ivdomlo;
              ivcurhi += IntVect::Unit;

              Tuple<int, 2> axdir = me->axisdir;
              Real rilo = rimin+ivcurlo[axdir[0]]*dx;
              Real rjlo = rjmin+ivcurlo[axdir[1]]*dy;
              Real rihi = rimin+ivcurhi[axdir[0]]*dx;
              Real rjhi = rjmin+ivcurhi[axdir[1]]*dy;
              //set color to white and draw filled rectangle
              ps_strm << "1.0 setgray" << endl;
              Real length = rihi-rilo;
              Real height = rjhi-rjlo;
              //jhi because of  wacky coordinates
              ps_strm << rilo    << "   "  <<  rjhi << "   "
                      << length << "   "  << height << endl;
              ps_strm << "rf"  << endl;
              //draw a box around region if needed
              ps_strm << "0.0 setgray" << endl;
              ps_strm << "0.5 setlinewidth" << endl;
              if(me->drawboxes)
                {
                  ps_strm << rilo    << "   "  <<  rjhi << "   "
                          << length << "   "  << height << endl;
                  ps_strm << "rs"  << endl;
                }

              Real df = (me->max - me->min)/Real(me->numcont);
              //stuff is not normalized to one
              for (int nc = 0; nc < me->numcont; nc++)
                {
                  Real fcont = me->min + df*(nc + 0.5);
                  
                  for(BoxIterator bit(curbox); bit.ok(); bit.next())
                    {
                      IntVect ivint = bit();
                      bool coveredCell = false;
                      //i am carefully nesting these so as to not
                      //kill performance in the case of no eb
                      if(me->ebflag)
                        {
                          if(me->coveredCells[ilev].contains(ivint))
                            coveredCell = true;
                        }
                      Vector<Tuple<Real, 2> > vecpoints;
                      Tuple<Real,5> xlocs;
                      Tuple<Real,5> ylocs;
                      getPSCorners(me,ivint,domain,
                                   dx, dy, rimin,rjmin,
                                   xlocs,ylocs);
                      if(!coveredCell)
                        {
                          Tuple<Real,5> flocs;
                          getDataLocs(me,ivint,domain,bigfab,
                                      flocs);
                          for(int iloc = 0; iloc <= 3; iloc++)
                            {
                              Real fint = flocs[iloc];
                              Real xint = xlocs[iloc];
                              Real yint = ylocs[iloc];
                              Real fout = flocs[iloc+1];
                              Real xout = xlocs[iloc+1];
                              Real yout = ylocs[iloc+1];
                              if(
                                 ((fint > fcont)&& (fout < fcont))
                                 ||((fint < fcont)&& (fout > fcont))
                                 )
                                {
                                  Real xpt =  spainterp(fint, fout, fcont, xint, xout);
                                  Real ypt =  spainterp(fint, fout, fcont, yint, yout);

                                  Tuple<Real, 2> tupp;
                                  tupp[0] = xpt;
                                  tupp[1] = ypt;
                                  vecpoints.push_back(tupp);
                                }
                            }
                          if(vecpoints.size() > 0)
                            {
                              Tuple<Real, 2> pt1 = vecpoints[0];
                              ps_strm << pt1[0] << " " << pt1[1] << " mt"
                                      << endl;
                              for(int ivec=0;ivec < int(vecpoints.size()); ivec++)
                                {
                                  Tuple<Real, 2> pt2 = vecpoints[ivec];
                                  ps_strm << pt2[0] << " " << pt2[1] << " lt"
                                          << endl;
                                }
                              ps_strm << "st" << endl;
                            }
                        } //if(!covered)
                      else
                        {
                          //draw black box if cell is covered
                          ps_strm << "0.0 setgray" << endl;

                          Real riloloc = xlocs[0];
                          Real rjloloc = ylocs[0];

                          ps_strm << riloloc    << "   "  <<  rjloloc << "   "
                                  << dx << "   "  << dy << endl;
                          ps_strm << "rf"  << endl;
                          //ps_strm << "st" << endl;
                        }
                    } // end loop over boxes
                } // end loop over contours

              //sets color to black
              ps_strm << "0.0 setgray" << endl;

            }// if the current box is not empty

        } //end loop over boxes in the level
      idepth *= me->vect_ratio[ilev];
      dx /= me->vect_ratio[ilev];
      dy /= me->vect_ratio[ilev];
      domain.refine(me->vect_ratio[ilev]);
    }

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
}
/***********************************************/
void dumpPseudoColor(PostScriptData* me)
{
  ofstream ps_strm(me->filename.c_str());

  /* draw stuff */
  int papersize = 8*72;
  int imin = 0;
  int jmin = 0;
  int imax = papersize;
  int jmax = papersize;
  Real rimin = imin;
  Real rjmin = jmin;
  Real rimax = imax;
  Real rwidth  = rimax - rimin;

  /* obligatory postscript begining */
  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin << " " << imax << " " << jmax << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;
  ps_strm << "/st {stroke} def" << endl;
  ps_strm << "/cp {closepath} def" << endl;

  /* draw box around plot */
  if(!me->doingsubregion)
    {
      ps_strm << "3 setlinewidth" << endl;
      ps_strm << imin << " " << jmin << " mt" << endl;
      ps_strm << imin << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmax << " lt" << endl;
      ps_strm << imax << " " << jmin << " lt" << endl;
      ps_strm << "cp" << endl << "st" << endl;
    }

  /* plot */
  /* upper left corner is origin---ergo weirdness */

  getMaxMinMag(me);
  getAxes(me);
  int idepth = me->idepth;

  Box domain = me->vect_box[0];
  if(me->doingsubregion) domain = me->lev0subbox;
  Real dx = rwidth/domain.size(0);
  Real dy = rwidth/domain.size(1);
  if(dx > dy) dx = dy;
  if(dy > dx) dy = dx;
                
  for(int  ilev = 0; ilev < int(me->vect_mf.size()); ilev++)
    {


      LevelData<FArrayBox>& mfcur = *me->vect_mf[ilev];
      DisjointBoxLayout dbl = mfcur.disjointBoxLayout();
      DataIterator dit = mfcur.dataIterator();
      for(dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& bigfab = mfcur[dit()];
          Box bigbox = bigfab.box();
          Box curbox = dbl.get(dit());

#if(CH_SPACEDIM == 3)
          Box threedbox = curbox;
          IntVect loiv = threedbox.smallEnd();
          IntVect hiiv = threedbox.bigEnd();
          loiv[me->inormal]= idepth;
          hiiv[me->inormal]= idepth;
          curbox= Box(loiv,hiiv);
          curbox &= threedbox;
#endif
          curbox &= domain;
          if(!curbox.isEmpty())
            {
              //stuff is not normalized to one
              for(BoxIterator bit(curbox); bit.ok(); bit.next())
                {
                  IntVect ivint = bit();
                  bool coveredCell = false;
                  //i am carefully nesting these so as to not
                  //kill performance in the case of no eb
                  if(me->ebflag)
                    {
                      if(me->coveredCells[ilev].contains(ivint))
                        coveredCell = true;
                    }
                  Vector<Tuple<Real, 2> > vecpoints;
                  Tuple<Real,5> xlocs;
                  Tuple<Real,5> ylocs;
                  getPSCorners(me,ivint,domain,
                               dx, dy,rimin,rjmin,
                               xlocs,ylocs);

                  Real fval = bigfab(ivint, me->cur_var);
                  Real fmag = me->mag;
                  Real fmin = me->min;
                  Real fscale;
                  if(Abs(fmag) > 1.0e-10) 
                    fscale = (fval-fmin)/fmag;
                  else 
                    fscale = 0;

                  int ir, ig, ib;
                  if(coveredCell)
                    {
                      ir = 0;
                      ig = 0;
                      ib = 0;
                    }
                  else
                    {
                      int icolor = (int)(253.*fscale);
                      icolor =abs(icolor)+2;
                      if((icolor>256)||(icolor<0)) icolor = 0;
                      ir = me->rcol[icolor];
                      ig = me->gcol[icolor];
                      ib = me->bcol[icolor];
                    }
                  Real rcolor = ((Real) ir)/256.;
                  Real gcolor = ((Real) ig)/256.;
                  Real bcolor = ((Real) ib)/256.;
                  ps_strm << setprecision(3)
                          << rcolor << "   "
                          << gcolor << "   "
                          << bcolor << "   "
                          << "rgb" << endl;

                  Real rilo = xlocs[0];
                  //strangeness with coordinate system means
                  // high is low,  war is peace
                  //freedom is slavery and ignorance is strength
                  Real rjlo = ylocs[0];

                  ps_strm << rilo    << "   "  <<  rjlo << "   "
                          << dx << "   "  <<  dy << endl;
                  ps_strm << "rf"  << endl;

                }
              //sets color to black
              if(me->drawboxes)
                {
                  ps_strm << "0.0 setgray" << endl;
                  IntVect ivlo = curbox.smallEnd();
                  IntVect ivhi = curbox.bigEnd();
                  ivlo -= domain.smallEnd();
                  ivhi -= domain.smallEnd();
                  ivhi += IntVect::Unit;
                  
                  Real rilo = rimin+ ivlo[me->axisdir[0]]*dx;
                  Real rjlo = rjmin+ ivlo[me->axisdir[1]]*dy;
                  Real rihi = rimin+ ivhi[me->axisdir[0]]*dx;
                  Real rjhi = rjmin+ ivhi[me->axisdir[1]]*dx;
                  /* draw box */
                  ps_strm << "0.5 setlinewidth" << endl;
                  ps_strm << rilo << " " << rjlo << " mt" << endl;
                  ps_strm << rilo << " " << rjhi << " lt" << endl;
                  ps_strm << rihi << " " << rjhi << " lt" << endl;
                  ps_strm << rihi << " " << rjlo << " lt" << endl;
                  ps_strm << "cp" << endl
                          << "st" << endl;
                }
            } //if curbox is not empty
        } //loop over boxes

      idepth *= me->vect_ratio[ilev];
      dx /= me->vect_ratio[ilev];
      dy /= me->vect_ratio[ilev];
      domain.refine(me->vect_ratio[ilev]);
    }//loop over levels

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
}
/***********************************************/
void dumpColorMap(PostScriptData* me)
{
  ofstream ps_strm(me->palettename.c_str());

  /* draw stuff */
  /* obligatory postscript begining */
  Real papersize = 8*72;
  Real rxdraw = g_paletteb;
  Real rydraw = me->ydraw;
  Real rimin = 0;
  Real rjmin = 0;
  Real rimax = papersize-36;
  Real rjmax = papersize-36;
  Real maxside = Max(rxdraw,rydraw);
  //Real xfac = .5 * rxdraw/maxside;
  Real xfac = 0.25;
  Real yfac = .5 * rydraw/maxside;
  rimax *= xfac;
  rjmax *= yfac;
  /*
    rimax += 36;
    rjmax += 36;
    rimin += 36;
    rjmin += 36;
  */

  int imax = int(rimax);
  int jmax = int(rjmax);
  int imin = int(rimin);
  int jmin = int(rjmin);

  ps_strm << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  // below, 12 is the font size
  ps_strm << "%%BoundingBox: "
          << imin << " " << jmin-12 << " " << imax << " " << jmax+12 << endl;
  ps_strm << "%%EndComments" << endl;
  ps_strm << "%%EndProlog" << endl;

  // abbreviations to make the file smaller
  ps_strm << "/Helvetica 12 selectfont" << endl;
  ps_strm << "/lt {lineto} def" << endl;
  ps_strm << "/mt {moveto} def" << endl;
  ps_strm << "/rf {rectfill} def" << endl;
  ps_strm << "/rgb {setrgbcolor} def" << endl;

  getMaxMinMag(me);
  Real fmin = me->min;
  Real fmag = me->mag;
  int ioffset = 20;

  Real  rbarheight = rjmax - rjmin;

  Real rdycol = rbarheight/256.0;
  Real rdxcol = ioffset;
  Real rilo = rimin;
  Real rjlo = rjmin;
  // the following loop starts at 1, not 0, to avoid displaying color 0,
  // which is always black regardless of the colormap.
  for(int ivec = 1; ivec < 256 ; ivec++)
    {
      Real rcolor = ((Real) me->rcol[ivec])/256.;
      Real gcolor = ((Real) me->gcol[ivec])/256.;
      Real bcolor = ((Real) me->bcol[ivec])/256.;
      ps_strm << rcolor << "   "
              << gcolor << "   "
              << bcolor << "   "
              << "rgb" << endl;
      ps_strm << rilo << " " << rjlo << " mt" << endl;
      ps_strm << rilo    << "   "  <<  rjlo << "   "
              << rdxcol << "   "  << rdycol << endl;
      ps_strm << "rf"  << endl;
      rjlo += rdycol;
    }

  int num_lines = 8;
  Real line_spacing = rbarheight/num_lines;
  Real dval = fmag/num_lines;
  Real value = fmin;
  rilo = rimin + 2*rdxcol;
  // 5 is about half the height of the characters
  rjlo = rjmin - 5;

  ps_strm << "0 setgray" << endl;

  for (int ivec = 0; ivec <= num_lines; ivec++)
    {
      ps_strm << rilo << " " << rjlo << " mt" << endl;
      ps_strm << "( " << value << " ) show" << endl;
      value +=  dval;
      rjlo  += line_spacing;
    }

  /* get out */
  ps_strm << "showpage" << endl;
  ps_strm.close();
}



void  getPSCorners(PostScriptData* me,IntVect& iv, Box& domain, Real dx, Real dy,
                   Real rimin, Real rjmin,Tuple<Real,5>& xlocs,Tuple<Real,5>& ylocs)
{
  //start at lower left, go counter-clockwise
  int im = iv[me->axisdir[0]];
  int jm = iv[me->axisdir[1]];
  IntVect ivdomlo = domain.smallEnd();
  int ivd = ivdomlo[me->axisdir[0]];
  int jvd = ivdomlo[me->axisdir[1]];
  Real yl = rjmin+(jm-jvd)*dy;
  Real yh = rjmin+(jm-jvd+1)*dy;
  Real xl = rimin+(im-ivd)*dx;
  Real xh = rimin+(im-ivd+1)*dx;

  xlocs[0] = xl;
  ylocs[0] = yl;
  xlocs[1] = xl;
  ylocs[1] = yh;
  xlocs[2] = xh;
  ylocs[2] = yh;
  xlocs[3] = xh;
  ylocs[3] = yl;
  xlocs[4] = xlocs[0];
  ylocs[4] = ylocs[0];
}
void getDataLocs(PostScriptData* me,IntVect& iv,Box& domain,
                 FArrayBox& bigstate,Tuple<Real,5>& flocs)
{
  int ivar = me->cur_var;
  SideIterator sit;
  int idir0 = me->axisdir[0];
  int idir1 = me->axisdir[1];
  IntVect iv0 = BASISV(idir0);
  IntVect iv1 = BASISV(idir1);
  Real fmid = bigstate(iv, ivar);
  int nf = 1;
  //lower left corner of cell
  Tuple<IntVect, 3> outiv;

  Real fave = fmid;
  outiv[0] = iv - iv0 - iv1;
  outiv[1] = iv - iv0;
  outiv[2] = iv - iv1;
  for(int i = 0; i < 3; i++)
    {
      if(domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[0] = fave/nf;

  //upper left corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv - iv0 + iv1;
  outiv[1] = iv - iv0;
  outiv[2] = iv + iv1;
  for(int i = 0; i < 3; i++)
    {
      if(domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[1] = fave/nf;

  //upper right corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv + iv0 + iv1;
  outiv[1] = iv + iv0;
  outiv[2] = iv + iv1;
  for(int i = 0; i < 3; i++)
    {
      if(domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[2] = fave/nf;

  //lower right corner of cell
  fave = fmid;
  nf = 1;
  outiv[0] = iv + iv0 - iv1;
  outiv[1] = iv + iv0;
  outiv[2] = iv - iv1;
  for(int i = 0; i < 3; i++)
    {
      if(domain.contains(outiv[i]))
        {
          fave += bigstate(outiv[i], ivar);
          nf++;
        }
    }
  flocs[3] = fave/nf;
  flocs[4] = flocs[0];
}
