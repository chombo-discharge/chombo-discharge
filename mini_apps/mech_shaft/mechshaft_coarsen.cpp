/*!
  @file   mechshaft_coarsen.cpp
  @brief  Implementation of mechshaft_coarsen.H
  @author Robert Marskar
  @date   June 2018
*/

#include "mechshaft_coarsen.H"

#include <ParmParse.H>

mechshaft_coarsen::mechshaft_coarsen(){

  bool enable_coarsen = false;
  
  {
    ParmParse pp("mechshaft_coarsen");
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
      Vector<Real> lo2(SpaceDim, 0.0);
      Vector<Real> hi2(SpaceDim, 0.0);

      pp.queryarr("box1_lo_corner", lo1, 0, SpaceDim);
      pp.queryarr("box1_hi_corner", hi1, 0, SpaceDim);
      pp.queryarr("box2_lo_corner", lo2, 0, SpaceDim);
      pp.queryarr("box2_hi_corner", hi2, 0, SpaceDim);
      pp.query("finest_level", coarsen_levels);

      real_box box1 = real_box(RealVect(D_DECL(lo1[0], lo1[1], lo1[2])), RealVect(D_DECL(hi1[0], hi1[1], hi1[2])));
      real_box box2 = real_box(RealVect(D_DECL(lo2[0], lo2[1], lo2[2])), RealVect(D_DECL(hi2[0], hi2[1], hi2[2])));
      
      m_coarsen_boxes.push_back(box1);
      m_coarsen_boxes.push_back(box2);
      m_coarsen_levels.resize(2, coarsen_levels);

    }
  }
}

mechshaft_coarsen::~mechshaft_coarsen(){

}

