/*!
  @file   cell_tagger.cpp
  @brief  Implementation of cell_tagger.H
  @author Robert marskar
  @date   May. 2018
*/

#include "cell_tagger.H"
#include "data_ops.H"

#include <EBArith.H>
#include <ParmParse.H>

cell_tagger::cell_tagger(){
  CH_TIME("cell_tagger::cell_tagger");
  m_verbosity = 10;
  if(m_verbosity > 5){
    pout() << "cell_tagger::cell_tagger" << endl;
  }

  m_name  = "cell_tagger";
}

cell_tagger::~cell_tagger(){

}

int cell_tagger::get_num_plot_vars(){
  return 0;
}

int cell_tagger::get_buffer(){
  return m_buffer;
}

void cell_tagger::parse_boxes(){
  
  ParmParse pp(m_name.c_str());

  int num_boxes = 0;
  pp.get("num_boxes", num_boxes);

  m_tagboxes.resize(0);
  if(num_boxes > 0){
    m_tagboxes.resize(num_boxes);

    const int ndigits = (int) log10((double) num_boxes) + 1;
      
    for (int ibox = 0; ibox < num_boxes; ibox++){
      char* cstr = new char[ndigits];
      sprintf(cstr, "%d", 1+ibox);

      std::string str1 = "box" + std::string(cstr) + "_lo";
      std::string str2 = "box" + std::string(cstr) + "_hi";

      Vector<Real> corner_lo(SpaceDim);
      Vector<Real> corner_hi(SpaceDim);

      pp.getarr(str1.c_str(), corner_lo, 0, SpaceDim);
      pp.getarr(str2.c_str(), corner_hi, 0, SpaceDim);

      const RealVect c1 = RealVect(D_DECL(corner_lo[0], corner_lo[1], corner_lo[2]));
      const RealVect c2 = RealVect(D_DECL(corner_hi[0], corner_hi[1], corner_hi[2]));

      m_tagboxes[ibox] = real_box(c1,c2);
	
      delete cstr;
    }
  }
}

void cell_tagger::parse_buffer(){

  ParmParse pp(m_name.c_str());
  pp.get("buffer", m_buffer);
  m_buffer = Max(0, m_buffer);
}

void cell_tagger::parse_verbosity(){
  CH_TIME("cell_tagger::parse_verbosity");

  ParmParse pp(m_name.c_str());
  pp.get("verbosity", m_verbosity);
}

void cell_tagger::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp){
  CH_TIME("cell_tagger::write_plot_data");
  if(m_verbosity > 3){
    pout() << "cell_tagger::write_plot_data" << endl;
  }
}

bool cell_tagger::inside_tag_box(const RealVect a_pos){
  bool do_this_refine = (m_tagboxes.size() > 0) ? false : true;
  for (int ibox = 0; ibox < m_tagboxes.size(); ibox++){
    const RealVect lo = m_tagboxes[ibox].get_lo();
    const RealVect hi = m_tagboxes[ibox].get_hi();

    if(a_pos >= lo && a_pos <= hi){
      do_this_refine = true;
    }
  }

  return do_this_refine;
}
