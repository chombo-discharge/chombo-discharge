/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_GeoCoarsener.cpp
  @brief  Implementation of CD_GeoCoarsener.H
  @author Robert Marskar
*/

// Std includes
#include <sstream>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_GeoCoarsener.H>
#include <CD_NamespaceHeader.H>

GeoCoarsener::GeoCoarsener(){
  m_coarsen_boxes.resize(0);
  m_coarsen_levels.resize(0);

  
  { // Info from input script
    ParmParse pp("GeoCoarsener");

    int num_boxes = 0;
    pp.query("num_boxes", num_boxes);

    if(num_boxes > 0){
      m_coarsen_boxes.resize(num_boxes);
      m_coarsen_levels.resize(num_boxes);
      m_inverse.resize(num_boxes, 0);
      
      const int ndigits = (int) log10((double) num_boxes) + 1;
      
      for (int ibox = 0; ibox < num_boxes; ibox++){
	char* cstr = new char[ndigits];
	sprintf(cstr, "%d", 1+ibox);

	std::string str = "false";
	std::string str1 = "box" + std::string(cstr) + "_lo";
	std::string str2 = "box" + std::string(cstr) + "_hi";
	std::string str3 = "box" + std::string(cstr) + "_lvl";
	std::string str4 = "box" + std::string(cstr) + "_inv";

	Vector<Real> corner_lo(SpaceDim);
	Vector<Real> corner_hi(SpaceDim);
	int finest_box_lvl;

	pp.getarr(str1.c_str(), corner_lo, 0, SpaceDim);
	pp.getarr(str2.c_str(), corner_hi, 0, SpaceDim);
	pp.get(str3.c_str(),    finest_box_lvl);
	pp.query(str4.c_str(),    str);

	const RealVect c1 = RealVect(D_DECL(corner_lo[0], corner_lo[1], corner_lo[2]));
	const RealVect c2 = RealVect(D_DECL(corner_hi[0], corner_hi[1], corner_hi[2]));

	m_coarsen_boxes[ibox]  = real_box(c1,c2);
	m_coarsen_levels[ibox] = finest_box_lvl;
	m_inverse[ibox]        = (str == "true") ? 1 : 0;
	
	delete cstr;
      }
    }
  }
}

GeoCoarsener::~GeoCoarsener(){
}

void GeoCoarsener::coarsenTags(Vector<IntVectSet>& a_tags, const Vector<Real>& a_dx, const RealVect& a_origin) const {
  CH_TIME("GeoCoarsener::coarsenTags");
  if(!(m_coarsen_boxes.size() == m_coarsen_levels.size())){
    pout() << "GeoCoarsener::coarsenTags - m_geoCoarsen is not well defined. Skipping the coarsening step" << endl;
  }
  else{
    if(m_coarsen_boxes.size() > 0){
      for (int lvl = 0; lvl < a_tags.size(); lvl++){
	const int num_coarsen = m_coarsen_boxes.size();
	const IntVectSet tmp = a_tags[lvl];
	for (IVSIterator it(tmp); it.ok(); ++it){
	  const IntVect iv   = it();
	  const RealVect pos = a_origin + RealVect(iv)*a_dx[lvl];

	  bool remove_tag = false;

	  bool inside_coarsen_box  = false;
	  bool inside_inverse_box  = false;
	  bool outside_inverse_box = false;

	  // Check the regular box
	  for (int ibox = 0; ibox < num_coarsen; ibox++){
	    const RealVect lo = m_coarsen_boxes[ibox].get_lo();
	    const RealVect hi = m_coarsen_boxes[ibox].get_hi();

	    const bool inverse     = (m_inverse[ibox] != 0) ? true : false;
	    const bool inside_box  = pos > lo && pos < hi;
	    const bool coarsen_lvl = lvl >= m_coarsen_levels[ibox];

	    if(inside_box && !inverse && coarsen_lvl){
	      inside_coarsen_box = true;
	    }
	  }

	  // Check the inverse boxes
	  for (int ibox = 0; ibox < num_coarsen; ibox++){
	    const RealVect lo = m_coarsen_boxes[ibox].get_lo();
	    const RealVect hi = m_coarsen_boxes[ibox].get_hi();

	    const bool inverse     = (m_inverse[ibox] != 0) ? true : false;
	    const bool inside_box  = pos > lo && pos < hi;
	    const bool coarsen_lvl = lvl >= m_coarsen_levels[ibox];

	    if(inside_box && inverse && coarsen_lvl){ // Protect tag
	      inside_inverse_box = true;
	    }
	    else if(!inside_box && inverse && coarsen_lvl){ // Remove tag
	      outside_inverse_box = true;
	    }
	  }

	  if(inside_coarsen_box && !inside_inverse_box){
	    a_tags[lvl] -= iv;
	  }
	  if(!inside_inverse_box && outside_inverse_box){
	    a_tags[lvl] -= iv;
	  }
	}
      }
    }
  }
}

Vector<real_box> GeoCoarsener::getCoarsenBoxes(){
  return m_coarsen_boxes;
}


Vector<int> GeoCoarsener::getCoarsenLevels(){
  return m_coarsen_levels;
}

#include <CD_NamespaceFooter.H>
