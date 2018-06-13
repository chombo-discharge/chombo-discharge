/*!
  @file   rod_coarsen.cpp
  @brief  Implementation of rod_coarsen.H
  @author Robert Marskar
  @date   June 2018
*/

#include "rod_coarsen.H"

#include <ParmParse.H>

rod_coarsen::rod_coarsen(){

  bool enable_coarsen = false;
  
  {
    ParmParse pp("rod_coarsen");
    std::string str;

    pp.query("enable_coarsening", str);
    if(str == "true"){
      enable_coarsen = true;
    }
    else{
      enable_coarsen = false;
    }

    if(enable_coarsen){
      int coarsen_levels = -1;
      Vector<Real> lo1(SpaceDim, 0.0);
      Vector<Real> hi1(SpaceDim, 0.0);

      pp.queryarr("box_lo_corner", lo1, 0, SpaceDim);
      pp.queryarr("box_hi_corner", hi1, 0, SpaceDim);
      pp.query("finest_level", coarsen_levels);

      real_box box1 = real_box(RealVect(D_DECL(lo1[0], lo1[1], lo1[2])), RealVect(D_DECL(hi1[0], hi1[1], hi1[2])));
      
      m_coarsen_boxes.push_back(box1);
      m_coarsen_levels.resize(1, coarsen_levels);

    }
  }
}

rod_coarsen::~rod_coarsen(){

}

