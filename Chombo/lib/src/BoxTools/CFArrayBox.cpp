#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include "SPACE.H"
using std::cin;
using std::cout;
using std::cerr;
using std::setw;
using std::setprecision;
using std::ios;
using std::pow;
using std::sqrt;

#include "Misc.H"
#include "CFArrayBox.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

CFArrayBox::CFArrayBox()
  :
  BaseFab<Complex>()
{
}

CFArrayBox::CFArrayBox(const Box& a_box,
                       int        a_n,
                       Complex*   a_alias)
  :
  BaseFab<Complex>(a_box,a_n,a_alias)
{
  // Note: All work is done in BaseFab<Complex> constructor
}

CFArrayBox::~CFArrayBox()
{
}

Real CFArrayBox::norm(const Box& a_subbox,
                      int        a_p,
                      int        a_comp,
                      int        a_numcomp) const
{
  CH_assert(a_comp >= 0 && a_comp + a_numcomp <= nComp());
  CH_assert(a_numcomp > 0);
  CH_assert(a_p >= 0);
  CH_assert(m_domain.contains(a_subbox));

  Real* tmp = 0;
  int tmplen = 0;
  Real nrm = 0;
#ifdef CH_COUNT_FLOPS
  unsigned long long int np=a_subbox.numPts()*a_numcomp;
#endif
  if (a_p == 0)
    {
      // here begins a normed fab function piece
      ForAllThisCPencil(Complex,a_subbox,a_comp,a_numcomp)
        {
          const Complex* row = &thisR;
          if (tmp == 0)
            {
              tmp = new Real[thisLen];
              tmplen = thisLen;
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] = std::abs(row[i]);
                }
            }
          else
            {
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] = Max(tmp[i],std::abs(row[i]));
                }
            }
        } EndForPencil

      nrm = tmp[0];
      for (int i = 1; i < tmplen; i++)
        {
          nrm = Max(nrm, tmp[i]);
        }
#ifdef CH_COUNT_FLOPS
      // 4 flops for std::abs(a+b*i) = sqrt(a*a + b*b)
      ch_flops()+=4*np*a_numcomp;
#endif
      // here it ends
    }
  else if (a_p == 1)
    {
      // here begins a normed fab function piece
      ForAllThisCPencil(Complex,a_subbox,a_comp,a_numcomp)
        {
          const Complex* row = &thisR;
          if (tmp == 0)
            {
              tmp = new Real[thisLen];
              tmplen = thisLen;
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] =std::abs(row[i]);
                }
            }
          else
            {
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] += std::abs(row[i]);
                }
            }
        } EndForPencil

      nrm = tmp[0];
      for (int i = 1; i < tmplen; i++)
        {
          nrm += tmp[i];
        }
      // here it ends
#ifdef CH_COUNT_FLOPS
      // 4 flops for std::abs(a+b*i) = sqrt(a*a + b*b)
      ch_flops()+=4*np*a_numcomp;
#endif
    }
  else if (a_p == 2)
    {
      nrm = std::abs(sqrt(sumPow(a_subbox, 2, a_comp, a_numcomp)));
    }
  else
    {
      // so standard norms weren't good enough for you?
      Real invpwr = 1.0/a_p;
      nrm = std::abs(pow(sumPow(a_subbox, a_p, a_comp, a_numcomp),invpwr));
    }

  delete [] tmp;

  return nrm;
}

Real CFArrayBox::norm(int a_p,
                      int a_comp,
                      int a_numcomp) const
{
  return norm(m_domain,a_p,a_comp,a_numcomp);
}

// utility function used in norms and such only works for a_p >= 2
Complex CFArrayBox::sumPow(const Box& a_subbox,
                           int        a_p,
                           int        a_comp,
                           int        a_numcomp) const
{
  CH_assert(a_p >= 2);
  CH_assert(a_numcomp > 0);
  CH_assert(m_domain.contains(a_subbox));

  Complex sum = Complex::Zero();
  Complex* tmp = NULL;
  int tmplen = 0;

  if (a_p == 2)
  {
    ForAllThisCPencil(Complex,a_subbox,a_comp,a_numcomp)
    {
      const Complex* row = &thisR;
      if (tmp == 0)
      {
        tmp = new Complex[thisLen];
        tmplen = thisLen;
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] = row[i]*row[i];
        }
      }
      else
      {
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] += row[i]*row[i];
        }
      }
    } EndForPencil

    sum = tmp[0];
    for (int i = 1; i < tmplen; i++)
    {
      sum += tmp[i];
    }
  }
  else
  {
    // so standard norms weren't good enough for you?
    // Complex pwr = Complex(a_p);

    ForAllThisCPencil(Complex,a_subbox,a_comp,a_numcomp)
    {
      const Complex* row = &thisR;
      if (tmp == 0)
      {
        tmp = new Complex[thisLen];
        tmplen = thisLen;
        for (int i = 0; i < thisLen; i++)
        {
          // tmp[i] = pow(row[i],pwr);
          std::complex<double> t = pow(row[i],a_p);
          tmp[i].re()=t.real();
          tmp[i].im()=t.imag();
        }
      }
      else
      {
        for (int i = 0; i < thisLen; i++)
        {
          // tmp[i] += pow(row[i],pwr);
          std::complex<double> t = pow(row[i],a_p);
          tmp[i].re() += t.real();
          tmp[i].im() += t.imag();
        }
      }
    } EndForPencil

    sum = tmp[0];
    for (int i = 1; i < tmplen; i++)
    {
      sum += tmp[i];
    }
  }

  delete [] tmp;
#ifdef CH_COUNT_FLOPS
  // 6 flops per Complex multplication
  ch_flops()+=a_subbox.numPts()*a_p*6*a_numcomp;
#endif

  return sum;
}

// Take the dot product of "this" and "a_fab2" over their common box and
// all components.
Complex CFArrayBox::dotProduct(const CFArrayBox& a_fab2) const
{
  const CFArrayBox& fab1 = *this;
  Box commonBox = fab1.box() & a_fab2.box();
  return dotProduct(a_fab2, commonBox);
}

Complex CFArrayBox::dotProduct(const CFArrayBox& a_fab2, const Box& a_box) const
{
  Complex dot = Complex::Zero();
  const CFArrayBox& fab1 = *this;

  int startcomp = 0;
  int endcomp = fab1.nComp()-1;
  int numcomp = endcomp+1;

  CH_assert(fab1.nComp() == a_fab2.nComp());

  ForAllThisCBNNXC(Complex, a_box, startcomp, numcomp, a_fab2, startcomp)
    {
      dot += thisR * a_fab2R;
    } EndForTX
#ifdef CH_COUNT_FLOPS
        // 2 flops per Complex addition, 6 per Complex multiplication
        ch_flops()+=a_box.numPts()*8*numcomp;
#endif

  return dot;
}


// Complex CFArrayBox::min(int a_comp) const
// {
//   Complex *_min_row = 0;
//   int _X_len = 0;

//   ForAllThisCPencil(Complex,m_domain,a_comp,1)
//   {
//     const Complex* _row = &thisR;
//     if (_min_row == 0)
//     {
//       _min_row = new Complex[thisLen];
//       _X_len = thisLen;
//       for (int i = 0; i < thisLen; i++)
//       {
//         _min_row[i] = _row[i];
//       }
//     }
//     else
//     {
//       for (int i = 0; i < thisLen; i++)
//       {
//         _min_row[i] = Min(_row[i],_min_row[i]);
//       }
//     }
//   } EndForPencil;

//   Complex _min = _min_row[0];
//   for (int i = 1; i < _X_len; i++)
//   {
//     _min = Min(_min,_min_row[i]);
//   }

//   delete [] _min_row;

//   return _min;
// }

// Complex CFArrayBox::min(const Box& a_subbox,
//                     int        a_comp) const
// {
//   CH_assert(m_domain.contains(a_subbox));

//   Complex *_min_row = 0;
//   int _X_len = 0;

//   ForAllThisCPencil(Complex,a_subbox,a_comp,1)
//   {
//     const Complex* _row = &thisR;
//     if (_min_row == 0)
//     {
//       _min_row = new Complex[thisLen];
//       _X_len = thisLen;
//       for (int i = 0; i < thisLen; i++)
//       {
//         _min_row[i] = _row[i];
//       }
//     }
//     else
//     {
//       for (int i = 0; i < thisLen; i++)
//       {
//         _min_row[i] = Min(_row[i],_min_row[i]);
//       }
//     }
//   } EndForPencil;

//   Complex _min = _min_row[0];
//   for (int i = 1; i < _X_len; i++)
//   {
//     _min = Min(_min,_min_row[i]);
//   }

//   delete [] _min_row;

//   return _min;
// }

// Complex CFArrayBox::max(int a_comp) const
// {
//   Complex *_max_row = 0;
//   int _X_len = 0;

//   ForAllThisCPencil(Complex,m_domain,a_comp,1)
//   {
//     const Complex* _row = &thisR;
//     if (_max_row== 0)
//     {
//       _max_row = new Complex[thisLen];
//       _X_len = thisLen;
//       for (int i = 0; i < thisLen; i++)
//       {
//         _max_row[i] = _row[i];
//       }
//     }
//     else
//     {
//       for (int i = 0; i < thisLen; i++)
//       {
//         _max_row[i] = Max(_row[i],_max_row[i]);
//       }
//     }
//   } EndForPencil;

//   Complex _max = _max_row[0];
//   for (int i = 1; i < _X_len; i++)
//   {
//     _max = Max(_max,_max_row[i]);
//   }

//   delete [] _max_row;

//   return _max;
// }

// Complex CFArrayBox::max(const Box& a_subbox,
//                     int        a_comp) const
// {
//   CH_assert(m_domain.contains(a_subbox));

//   Complex *_max_row = 0;
//   int _X_len = 0;

//   ForAllThisCPencil(Complex,a_subbox,a_comp,1)
//   {
//     const Complex* _row = &thisR;
//     if (_max_row == 0)
//     {
//       _max_row = new Complex[thisLen];
//       _X_len = thisLen;
//       for (int i = 0; i < thisLen; i++)
//       {
//         _max_row[i] = _row[i];
//       }
//     }
//     else
//     {
//       for (int i = 0; i < thisLen; i++)
//       {
//         _max_row[i] = Max(_row[i],_max_row[i]);
//       }
//     }
//   } EndForPencil;

//   Complex _max = _max_row[0];
//   for (int i = 1; i < _X_len; i++)
//   {
//     _max = Max(_max,_max_row[i]);
//   }

//   delete [] _max_row;

//   return _max;
// }

// IntVect CFArrayBox::minIndex(int a_comp) const
// {
//   IntVect _min_loc(m_domain.smallEnd());
//   Complex _min_val = (*this).operator()(_min_loc,a_comp);

//   ForAllThisCBNN(Complex,m_domain,a_comp,1)
//   {
//     if (thisR < _min_val)
//     {
//       _min_val = thisR;
//       D_EXPR6(_min_loc[0] = iR,
//               _min_loc[1] = jR,
//               _min_loc[2] = kR,
//               _min_loc[3] = _iv[3],
//               _min_loc[4] = _iv[4],
//               _min_loc[5] = _iv[5]);
//     }
//   } EndFor

//   return _min_loc;
// }

// IntVect CFArrayBox::minIndex(const Box& a_subbox,
//                             int        a_comp) const
// {
//   CH_assert(m_domain.contains(a_subbox));

//   IntVect _min_loc(a_subbox.smallEnd());
//   Complex _min_val = (*this).operator()(_min_loc,a_comp);

//   ForAllThisCBNN(Complex,a_subbox,a_comp,1)
//   {
//     if (thisR < _min_val)
//     {
//       _min_val = thisR;
//       D_EXPR6(_min_loc[0] = iR,
//               _min_loc[1] = jR,
//               _min_loc[2] = kR,
//               _min_loc[3] = _iv[3],
//               _min_loc[4] = _iv[4],
//               _min_loc[5] = _iv[5]);
//     }
//   } EndFor

//   return _min_loc;
// }

// IntVect CFArrayBox::maxIndex(int a_comp) const
// {
//   IntVect _max_loc(m_domain.smallEnd());
//   Complex _max_val = (*this).operator()(_max_loc,a_comp);

//   ForAllThisCBNN(Complex,m_domain,a_comp,1)
//   {
//     if (thisR > _max_val)
//     {
//       _max_val = thisR;
//       D_EXPR6(_max_loc[0] = iR,
//               _max_loc[1] = jR,
//               _max_loc[2] = kR,
//               _max_loc[3] = _iv[3],
//               _max_loc[4] = _iv[4],
//               _max_loc[5] = _iv[5]);
//     }
//   } EndFor

//   return _max_loc;
// }

// IntVect CFArrayBox::maxIndex(const Box& a_subbox,
//                             int        a_comp) const
// {
//   CH_assert(m_domain.contains(a_subbox));

//   IntVect _max_loc(a_subbox.smallEnd());
//   Complex _max_val = (*this).operator()(_max_loc,a_comp);

//   ForAllThisCBNN(Complex,a_subbox,a_comp,1)
//   {
//     if (thisR > _max_val)
//     {
//       _max_val = thisR;
//       D_EXPR6(_max_loc[0] = iR,
//               _max_loc[1] = jR,
//               _max_loc[2] = kR,
//               _max_loc[3] = _iv[3],
//               _max_loc[4] = _iv[4],
//               _max_loc[5] = _iv[5]);
//     }
//   } EndFor

//   return _max_loc;
// }



void CFArrayBox::abs()
{
  ForAllThis(Complex)
  {
    thisR = Complex(std::abs(thisR));
  } EndFor
}

void CFArrayBox::abs(int a_comp,
                     int a_numcomp)
{
  ForAllThisNN(Complex,a_comp,a_numcomp)
  {
    thisR = Complex(std::abs(thisR));
  } EndFor
}

void CFArrayBox::abs(const Box& a_subbox,
                     int        a_comp,
                     int        a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNN(Complex,a_subbox,a_comp,a_numcomp)
  {
    thisR = Complex(std::abs(thisR));
  } EndFor
}

//here is where the normedfab functions end

//here is where the arithfab functions begin

Complex CFArrayBox::sum(int a_comp,
                        int a_numcomp) const
{
  Complex *_sum_row = 0;
  int _sum_len = 0;

  ForAllThisCPencil(Complex,m_domain,a_comp,a_numcomp)
  {
    const Complex* _row = &thisR;
    if (_sum_row == 0)
    {
      _sum_row = new Complex[thisLen];
      _sum_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] += _row[i];
      }
    }
  } EndForPencil;

  Complex _sum = _sum_row[0];
  for (int i = 1; i < _sum_len; i++)
  {
    _sum += _sum_row[i];
  }

  delete [] _sum_row;
#ifdef CH_COUNT_FLOPS
  // 2 flops per Complex addition
  ch_flops()+=2*m_domain.numPts()*a_numcomp;
#endif
  return _sum;
}

Complex CFArrayBox::sum(const Box& a_subbox,
                        int        a_comp,
                        int        a_numcomp) const
{
  CH_assert(m_domain.contains(a_subbox));

  Complex *_sum_row = 0;
  int _sum_len = 0;

  ForAllThisCPencil(Complex,a_subbox,a_comp,a_numcomp)
  {
    const Complex* _row = &thisR;
    if (_sum_row == 0)
    {
      _sum_row = new Complex[thisLen];
      _sum_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] = _row[i];
      }
    }
    else
      {
        for (int i = 0; i < thisLen; i++)
        {
          _sum_row[i] += _row[i];
        }
      }
  } EndForPencil;

  Complex _sum = _sum_row[0];
  for (int i = 1; i < _sum_len; i++)
  {
    _sum += _sum_row[i];
  }

  delete [] _sum_row;
#ifdef CH_COUNT_FLOPS
  // 2 flops per Complex addition
  ch_flops()+=2*a_subbox.numPts()*a_numcomp;
#endif
  return _sum;
}

CFArrayBox& CFArrayBox::invert(Complex a_r)
{
  ForAllThis(Complex)
  {
    thisR = a_r / thisR;
  } EndFor

#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex division
      ch_flops()+=11*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::invert(Complex a_r,
                               int  a_comp,
                               int  a_numcomp)
{
  ForAllThisNN(Complex,a_comp,a_numcomp)
  {
    thisR = a_r / thisR;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex division
      ch_flops()+=11*m_domain.numPts()*a_numcomp;
#endif

  return *this;
}

CFArrayBox& CFArrayBox::invert(Complex       a_r,
                               const Box& a_subbox,
                               int        a_comp,
                               int        a_numcomp)
{
  ForAllThisBNN(Complex,a_subbox,a_comp,a_numcomp)
  {
    thisR = a_r / thisR;
  } EndFor

#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex division
      ch_flops()+=11*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::negate(const Box& a_subbox,
                               int        a_comp,
                               int        a_numcomp)
{
  ForAllThisBNN(Complex,a_subbox,a_comp,a_numcomp)
  {
    thisR = - thisR;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex negation
      ch_flops()+=2*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::negate(int a_comp,
                               int a_numcomp)
{
  ForAllThisNN(Complex,a_comp,a_numcomp)
  {
    thisR = - thisR;
  } EndFor

#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex negation
      ch_flops()+=2*m_domain.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::negate()
{
  ForAllThis(Complex)
  {
    thisR = - thisR;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex negation
      ch_flops()+=2*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(Complex       a_r,
                             const Box& a_subbox,
                             int        a_comp,
                             int        a_numcomp)
{
  ForAllThisBNN(Complex,a_subbox,a_comp,a_numcomp)
  {
    thisR += a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=2*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(Complex a_r,
                             int  a_comp,
                             int  a_numcomp)
{
  ForAllThisNN(Complex,a_comp,a_numcomp)
  {
    thisR += a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=2*m_domain.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator += (Complex a_r)
{
  ForAllThis(Complex)
  {
    thisR += a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=2*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator += (const CFArrayBox& a_x)
{
  ForAllThisXC(Complex,a_x)
  {
    thisR += a_xR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=2*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(Complex a_r)
{
  return operator += (a_r);
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_x)
{
  return operator += (a_x);
}

// added bvs Tue May 18, PDT 1999
CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             const Complex&    a_scale)
{
  ForAllThisBNNXC(Complex, a_src.box() & m_domain, 0, m_nvar, a_src, 0)
  {
    thisR += a_srcR * a_scale;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition, 6 flops per Complex multiplication
      ch_flops()+=m_domain.numPts()*nComp()*8;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             const Complex&    a_scale,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXC(Complex, a_src.box() & m_domain, a_destcomp, a_numcomp, a_src, a_srccomp)
  {
    thisR += a_srcR * a_scale;
  } EndForTX

#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition, 6 flops per Complex multiplication
      ch_flops()+=m_domain.numPts()*nComp()*8;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXC(Complex,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=m_domain.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             const Box&       a_subbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Complex,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=a_subbox.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             const Box&       a_srcbox,
                             const Box&       a_destbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXCBN(Complex,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition
      ch_flops()+=a_srcbox.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::plus(const CFArrayBox& a_src,
                             const Box&       a_srcbox,
                             const Box&       a_destbox,
                             const Complex&   a_scale,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXCBN(Complex,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR += a_srcR * a_scale;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex addition, 6 flops per Complex multiplication
      ch_flops()+=a_srcbox.numPts()*a_numcomp*8;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator -= (Complex a_r)
{
  return operator += (-a_r);
}

CFArrayBox& CFArrayBox::operator -= (const CFArrayBox& a_x)
{
  ForAllThisXC(Complex,a_x)
  {
    thisR -= a_xR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex subtraction
      ch_flops()+=m_domain.numPts()*nComp()*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::minus(const CFArrayBox& a_x)
{
  return operator -= (a_x);
}

CFArrayBox& CFArrayBox::minus(const CFArrayBox& a_src,
                              int               a_srccomp,
                              int               a_destcomp,
                              int               a_numcomp)
{
  ForAllThisBNNXC(Complex,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex subtraction
      ch_flops()+=m_domain.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::minus(const CFArrayBox& a_src,
                              const Box&       a_subbox,
                              int              a_srccomp,
                              int              a_destcomp,
                              int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Complex,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex subtraction
      ch_flops()+=a_subbox.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::minus(const CFArrayBox& a_src,
                              const Box&       a_srcbox,
                              const Box&       a_destbox,
                              int              a_srccomp,
                              int              a_destcomp,
                              int              a_numcomp)
{
  ForAllThisBNNXCBN(Complex,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 2 flops per Complex subtraction
      ch_flops()+=a_srcbox.numPts()*a_numcomp*2;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator *= (Complex a_r)
{
  ForAllThis(Complex)
  {
    thisR *= a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
     // 6 flops per Complex Complex multiply
     ch_flops()+=6*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::mult(Complex a_r)
{
  return operator *= (a_r);
}

CFArrayBox& CFArrayBox::mult(Complex a_r,
                             int  a_comp,
                             int  a_numcomp)
{
  ForAllThisNN(Complex,a_comp,a_numcomp)
  {
    thisR *= a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*m_domain.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::mult(Complex    a_r,
                             const Box& a_subbox,
                             int        a_comp,
                             int        a_numcomp)
{
  ForAllThisBNN(Complex,a_subbox,a_comp,a_numcomp)
  {
    thisR *= a_r;
  } EndFor
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator *= (const CFArrayBox &a_x)
{
  ForAllThisXC(Complex,a_x)
  {
    thisR *= a_xR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*m_domain.numPts()*nComp();
#endif
  return *this;
}





template<unsigned char D, typename T1, typename T2, typename F>
struct BFMImpl
{
  static void f(T1* __restrict t1, const T2* __restrict t2, int* lo1, int* lo2, int* stride1, int* stride2, int* counting,F f)
  {
    t1+=lo1[D]*stride1[D-1];
    t2+=lo2[D]*stride2[D-1];
    for(int i=0; i<counting[D]; ++i, t1+=stride1[D-1], t2+=stride2[D-1]) 
      BFMImpl<D-1,T1,T2,F>::f(t1, t2, lo1, lo2, stride1, stride2, counting, f);
  }
};
 
template<typename T1, typename T2, typename F>
struct BFMImpl<0,T1,T2,F>
{
  static void f(T1* __restrict__ t1, const T2* __restrict__ t2, int* lo1, int* lo2, int* stride1, int* stride2, int* counting,F f)
  {
    int trips=counting[0];
    int l1=lo1[0];
    int l2=lo2[0];
    for(int i=0; i!=trips; ++i) f(t1+i+l1,t2[i+l2]);
  }
};

template<unsigned char R> inline void rcount(int* c, const Box& r)
{
  rcount<R-1>(c,r);
  c[R] = r.bigEnd()[R]-r.smallEnd()[R]+1;
}
template<> inline void rcount<0>(int* c, const Box& r)
{
  c[0]=r.bigEnd()[0]-r.smallEnd()[0]+1;
}
//helper template functions for Grid
template<unsigned char R> inline void gridStride(int* stride, int* lo, const Box& b, const Box& r)
{
  gridStride<R-1>(stride, lo, b, r);
  stride[R] = (b.bigEnd()[R]-b.smallEnd()[R]+1)*stride[R-1];
  lo[R] = r.smallEnd()[R]-b.smallEnd()[R];
}

template<> void inline gridStride<0>(int* stride, int* lo, const Box& b, const Box& r)
{
  stride[0] = b.bigEnd()[0]-b.smallEnd()[0]+1;
  lo[0] = r.smallEnd()[0]-b.smallEnd()[0];
}
template<unsigned char D, typename T1, typename T2, typename F>
void BFM(T1* __restrict t1, const T2* __restrict t2, const Box& b1, const Box& b2, int c1, 
         int c2, int ncomp, const Box& r, F f)
{
  CH_assert(b1.contains(r));
  CH_assert(b2.contains(r));
  int lo1[D+1], lo2[D+1]; lo1[D]=c1; lo2[D]=c2;
  int stride1[D], stride2[D];
  int counting[D+1]; counting[D]=ncomp;
  rcount<D-1>(counting,r);
  gridStride<D-1>(stride1, lo1, b1, r);
  gridStride<D-1>(stride2, lo2, b2, r);

  BFMImpl<D,T1,T2,F>::f(t1,t2,lo1,lo2,stride1,stride2,counting,f);
}


CFArrayBox& CFArrayBox::mult (const BaseFab<Real> &a_x)
{
  CH_assert(nComp()==a_x.nComp());
  if(box() == a_x.box())
    {
      const int count=m_nvar*m_numpts;
      const Real* p=a_x.dataPtr(0);
      for(size_t i=0; i!=count; i++)
        m_dptr[i]*=p[i];
#ifdef CH_COUNT_FLOPS
      // 2 flops per Real Complex multiply
      ch_flops()+=2*m_nvar*m_numpts;
#endif
    }
  else
    {
      Box b = box()&a_x.box();
      mult(a_x, b);
    }
  return *this;
}


inline void bFuncCR(Complex* __restrict__ a, Real b){(*a)*=b;} // we can replace this with a lambda in C++11 (bvs, Dec 2015)

CFArrayBox& CFArrayBox::mult(const BaseFab<Real>& a_x,
                             const Box&           a_box)
{
  CH_assert(nComp()==a_x.nComp());
  Complex* __restrict c=dataPtr();
  const Real* __restrict x = a_x.dataPtr();
  const int nvar = m_nvar;
  //with C++11 I can make D default to CH_SPACEDIM and we can pare this down to just BFM (bvs. Dec 2015)
  BFM<CH_SPACEDIM>(c, x,m_domain, a_x.box(), 0, 0, nvar, a_box,
                   [](Complex* a, Real b){(*a)*=b;});
#ifdef CH_COUNT_FLOPS
  // 2 flops per Real Complex multiply
  ch_flops()+=2*m_nvar*a_box.numPts();
#endif
  return *this;
} 

CFArrayBox& CFArrayBox::mult(const CFArrayBox& a_x)
{
  return operator *= (a_x);
}

CFArrayBox& CFArrayBox::mult(const CFArrayBox& a_src,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXC(Complex,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*m_domain.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::mult(const CFArrayBox& a_src,
                             const Box&       a_subbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Complex,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::mult(const CFArrayBox& a_src,
                             const Box&       a_srcbox,
                             const Box&       a_destbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXCBN(Complex,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 flops per Complex Complex multiply
      ch_flops()+=6*a_srcbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::operator /= (Complex a_r)
{
  // Multiply by reciprocal, because Complex Complex
  // multiplications require fewer operations than divisions.
  Complex reciprocal = Complex(1.) / a_r;
  return operator *= (reciprocal);
}

CFArrayBox& CFArrayBox::divide(Complex a_r)
{
  return operator /= (a_r);
}

CFArrayBox& CFArrayBox::divide(Complex a_r,
                               int  a_comp,
                               int  a_numcomp)
{
  // Multiply by reciprocal, because Complex Complex
  // multiplications require fewer operations than divisions.
  Complex reciprocal = Complex(1.) / a_r;
  return mult(reciprocal, a_comp, a_numcomp);
}

CFArrayBox& CFArrayBox::divide(Complex       a_r,
                               const Box& a_subbox,
                               int        a_comp,
                               int        a_numcomp)
{
  // Multiply by reciprocal, because Complex Complex
  // multiplications require fewer operations than divisions.
  Complex reciprocal = Complex(1.) / a_r;
  return mult(reciprocal, a_subbox, a_comp, a_numcomp);
}

CFArrayBox& CFArrayBox::operator /= (const CFArrayBox &a_x)
{
  ForAllThisXC(Complex,a_x)
  {
    thisR /= a_xR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex Complex divide
      ch_flops()+=11*m_domain.numPts()*nComp();
#endif
  return *this;
}

CFArrayBox& CFArrayBox::divide(const CFArrayBox& a_x)
{
  return operator /= (a_x);
}

CFArrayBox& CFArrayBox::divide(const CFArrayBox& a_src,
                               int              a_srccomp,
                               int              a_destcomp,
                               int              a_numcomp)
{
  ForAllThisBNNXC(Complex,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex Complex divide
      ch_flops()+=11*m_domain.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::divide(const CFArrayBox& a_src,
                               const Box&       a_subbox,
                               int              a_srccomp,
                               int              a_destcomp,
                               int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));
  ForAllThisBNNXC(Complex,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex Complex divide
      ch_flops()+=11*a_subbox.numPts()*a_numcomp;
#endif
  return *this;
}

CFArrayBox& CFArrayBox::divide(const CFArrayBox& a_src,
                               const Box&       a_srcbox,
                               const Box&       a_destbox,
                               int              a_srccomp,
                               int              a_destcomp,
                               int              a_numcomp)
{
  ForAllThisBNNXCBN(Complex,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 11 flops per Complex Complex divide
      ch_flops()+=11*a_srcbox.numPts()*a_numcomp;
#endif
  return *this;
}

//-----------------------------------------------------------------------
CFArrayBox&
CFArrayBox::
axby(const CFArrayBox& a_X, const CFArrayBox& a_Y,
     Complex a_A, Complex a_B)
{
  ForAllThisBNNXCBNYCBN(Complex, m_domain, 0, nComp(), a_X, a_X.m_domain, \
                        0, a_Y, a_Y.m_domain, 0)
  {
    thisR = a_A * a_XR + a_B * a_YR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 + 6 + 2 = 14 flops each
      ch_flops()+=14*m_domain.numPts()*nComp();
#endif
  return *this;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
CFArrayBox&
CFArrayBox::
axby(const CFArrayBox& a_X, const CFArrayBox& a_Y,
     Complex a_A, Complex a_B,
     int a_destComp, int a_xComp, int a_yComp)
{
  ForAllThisBNNXCBNYCBN(Complex, m_domain, a_destComp, 1, a_X, a_X.m_domain, \
                        a_xComp, a_Y, a_Y.m_domain, a_yComp)
  {
    thisR = a_A * a_XR + a_B * a_YR;
  } EndForTX
#ifdef CH_COUNT_FLOPS
      // 6 + 6 + 2 = 14 flops each
      ch_flops()+=14*m_domain.numPts();
#endif
  return *this;
}
//-----------------------------------------------------------------------

void CFArrayBox::performCopy(const BaseFab<Complex>& a_src,
                             const Box&           a_srcbox,
                             int                  a_srccomp,
                             const Box&           a_destbox,
                             int                  a_destcomp,
                             int                  a_numcomp)
{
  //BaseFab<Complex>::performCopy(a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp, a_numcomp);
  Box r = a_srcbox&a_destbox;
  BFM<CH_SPACEDIM>(dataPtr(),a_src.dataPtr(), box(), a_src.box(), a_destcomp, a_srccomp,a_numcomp,r,
                   [](Complex* __restrict__ d, Complex s){(*d)=s;});
}
#include "NamespaceFooter.H"
