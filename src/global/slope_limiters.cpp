/*!
  @file   slope_limiters.cpp
  @brief  Implementation of slope_limiters.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "slope_limiters.H"
#include "data_ops.H"


Real slope_limiters::koren(const Real a_slopeL, const Real a_slopeR){
  Real slope = 0.0;


  const Real aa = a_slopeL*a_slopeL;
  const Real ab = a_slopeL*a_slopeR;

  if(ab <= 0.0){
    slope = 0.0;
  }
  else if(aa <= 0.25*ab){
    slope = a_slopeL;
  }
  else if(aa <= 2.5*ab){
    slope = (2*a_slopeL + a_slopeR)/6.0;
  }
  else {
    slope = a_slopeR;
  }
  
  return slope;
}

Real slope_limiters::maxmod(const Real a_slopeL, const Real a_slopeR){
  Real slope = 0.0;

  if (a_slopeL*a_slopeR > 0.0){
    if(Abs(a_slopeL) > Abs(a_slopeR)){
      slope = a_slopeL;
    }
    else {
      slope = a_slopeR;
    }
  }
  
  return slope;
}

Real slope_limiters::minmod(const Real a_slopeL, const Real a_slopeR){
  Real slope = 0.0;

  if (a_slopeL*a_slopeR > 0.0){
    if(Abs(a_slopeL) < Abs(a_slopeR)){
      slope = a_slopeL;
    }
    else {
      slope = a_slopeR;
    }
  }
  
  return slope;
}

Real slope_limiters::superbee(const Real a_slopeL, const Real a_slopeR){
  Real slope = 0.0;

  if(a_slopeL*a_slopeR > 0.0){
    const Real slope1 = Min(2*a_slopeL, 1*a_slopeR);
    const Real slope2 = Min(1*a_slopeL, 2*a_slopeR);

    slope = slope_limiters::maxmod(slope1, slope2);
  }
  return slope;
}

Real slope_limiters::van_leer(const Real a_slopeL, const Real a_slopeR){
  Real slope = 0.0;

  if(a_slopeL*a_slopeR > 0.0){
    const Real slopeC = 0.5*(a_slopeL + a_slopeR);

    slope = data_ops::sgn(slopeC)*Min(2*a_slopeL, Min(2*a_slopeR, slopeC));
  }
  return slope;
}
