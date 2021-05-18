#include "air_eed.H"
#include "air_eed_species.H"

#include <chrono>
#include <random>

air_eed::eed::eed(){
  m_name = "eed";
  m_chargeNumber = 0;
  ParmParse pp("air_eed");
  std::string str;
  Vector<Real> vec(SpaceDim);

  pp.get("init_energy", m_init_energy);
  pp.get("mobile_electrons", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive = (str == "true") ? true : false;
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air_eed::electron::electron(){
  m_name = "electron";
  m_chargeNumber = -1;
  ParmParse pp("air_eed");
  std::string str;
  Vector<Real> vec(SpaceDim);
  
  pp.get("mobile_electrons", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_electrons", str); m_isDiffusive = (str == "true") ? true : false;
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air_eed::M_plus::M_plus(){
  m_name = "M_plus";
  m_chargeNumber = 1;
  ParmParse pp("air_eed");
  std::string str;
  Vector<Real> vec(SpaceDim);
    
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
  pp.get("uniform_density", m_uniform_density);
  pp.get("seed_density", m_seed_density);
  pp.get("seed_radius",   m_seed_rad);
  pp.getarr("seed_position", vec, 0, SpaceDim); m_seed_pos=RealVect(D_DECL(vec[0], vec[1], vec[2]));
}

air_eed::M_minus::M_minus(){
  m_name = "M_minus";
  m_chargeNumber = -1;
  ParmParse pp("air_eed");
  std::string str;
  
  pp.get("mobile_ions", str);    m_isMobile    = (str == "true") ? true : false;
  pp.get("diffusive_ions", str); m_isDiffusive = (str == "true") ? true : false;
}

air_eed::agg_photon::agg_photon(){
  m_name      = "agg_photon";

  Real p;
  ParmParse pp("air_eed");

  pp.get("pressure", p);
  pp.get("c4v0_X1v0_beer", m_kappa);
  m_kappa = 1./(m_kappa*p);


  // Get excitation effiencies and photoionization efficiencies
  pp.get("c4v0_exc_eff", m_c4v0_exc_eff);
  pp.get("c4v1_exc_eff", m_c4v1_exc_eff);
  pp.get("b1v1_exc_eff", m_b1v1_exc_eff);

  pp.get("c4v0_X1v0_photoi_eff", m_c4v0_X1v0_photoi_eff);
  pp.get("c4v0_X1v1_photoi_eff", m_c4v0_X1v1_photoi_eff);
  pp.get("c4v1_X1v0_photoi_eff", m_c4v1_X1v0_photoi_eff);
  pp.get("c4v1_X1v1_photoi_eff", m_c4v1_X1v1_photoi_eff);
  pp.get("c4v1_X1v2_photoi_eff", m_c4v1_X1v2_photoi_eff);
  pp.get("c4v1_X1v3_photoi_eff", m_c4v1_X1v3_photoi_eff);
  pp.get("b1v1_X1v0_photoi_eff", m_b1v1_X1v0_photoi_eff);
  pp.get("b1v1_X1v1_photoi_eff", m_b1v1_X1v1_photoi_eff);

  
  // Make probabilities and sum them to one
  m_prob_c4v0_X1v0 = m_c4v0_exc_eff*m_c4v0_X1v0_photoi_eff;
  m_prob_c4v0_X1v1 = m_c4v0_exc_eff*m_c4v0_X1v1_photoi_eff;
  m_prob_c4v1_X1v0 = m_c4v1_exc_eff*m_c4v1_X1v0_photoi_eff;
  m_prob_c4v1_X1v1 = m_c4v1_exc_eff*m_c4v1_X1v1_photoi_eff;
  m_prob_c4v1_X1v2 = m_c4v1_exc_eff*m_c4v1_X1v2_photoi_eff;
  m_prob_c4v1_X1v3 = m_c4v1_exc_eff*m_c4v1_X1v3_photoi_eff;
  m_prob_b1v1_X1v0 = m_b1v1_exc_eff*m_b1v1_X1v0_photoi_eff;
  m_prob_b1v1_X1v1 = m_b1v1_exc_eff*m_b1v1_X1v1_photoi_eff;

  //
  Real sum = 0.0;
  sum += m_prob_c4v0_X1v0;
  sum += m_prob_c4v0_X1v1;
  sum += m_prob_c4v1_X1v0;
  sum += m_prob_c4v1_X1v1;
  sum += m_prob_c4v1_X1v2;
  sum += m_prob_c4v1_X1v3;
  sum += m_prob_b1v1_X1v0;
  sum += m_prob_b1v1_X1v1;

  m_prob_c4v0_X1v0 *= 1./sum;
  m_prob_c4v0_X1v1 *= 1./sum;
  m_prob_c4v1_X1v0 *= 1./sum;
  m_prob_c4v1_X1v1 *= 1./sum;
  m_prob_c4v1_X1v2 *= 1./sum;
  m_prob_c4v1_X1v3 *= 1./sum;
  m_prob_b1v1_X1v0 *= 1./sum;
  m_prob_b1v1_X1v1 *= 1./sum;


  // Make the kappa vector
  m_kappas.resize(8);
  pp.get("c4v0_X1v0_beer", m_kappas[0]);
  pp.get("c4v0_X1v1_beer", m_kappas[1]);
  pp.get("c4v1_X1v0_beer", m_kappas[2]);
  pp.get("c4v1_X1v1_beer", m_kappas[3]);
  pp.get("c4v1_X1v2_beer", m_kappas[4]);
  pp.get("c4v1_X1v3_beer", m_kappas[5]);
  pp.get("b1v1_X1v0_beer", m_kappas[6]);
  pp.get("b1v1_X1v1_beer", m_kappas[7]);

  for (int i = 0; i < m_kappas.size(); i++){
    m_kappas[i] = 1./(p*m_kappas[i]);
  }

  // Make the probability vector
  m_probs.resize(9);
  m_probs[0] = 0.0;
  m_probs[1] = m_probs[0] + m_prob_c4v0_X1v0;
  m_probs[2] = m_probs[1] + m_prob_c4v0_X1v1;
  m_probs[3] = m_probs[2] + m_prob_c4v1_X1v0;
  m_probs[4] = m_probs[3] + m_prob_c4v1_X1v1;
  m_probs[5] = m_probs[4] + m_prob_c4v1_X1v2;
  m_probs[6] = m_probs[5] + m_prob_c4v1_X1v3;
  m_probs[7] = m_probs[6] + m_prob_b1v1_X1v0;
  m_probs[8] = m_probs[7] + m_prob_b1v1_X1v1;


  // INIT RNG
  pp.get("rng_seed", m_rng_seed);
  if(m_rng_seed < 0) {
    m_rng_seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  m_rng   = new std::mt19937_64(m_rng_seed);
  m_udist = new std::uniform_real_distribution<Real> (0.0, 1.0);
}

int air_eed::agg_photon::draw_photon_type() const {
  Real r = (*m_udist)(*m_rng);

  int ret = -1;
  for (int t = 0; t < m_probs.size()-1; t++){
    if(r >= m_probs[t] && r <= m_probs[t+1]){
      ret = t;
      break;
    }
  }

#if 1 // Debug
  if(ret < 0){
    MayDay::Abort("stop");
  }
#endif
  
  return ret;
}

Real air_eed::agg_photon::get_random_kappa() const{
  const int type = draw_photon_type();

  return m_kappas[type];
}

Real air_eed::eed::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();
  return m_init_energy*(m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad)));
}

Real air_eed::electron::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
}

Real air_eed::M_plus::initialData(const RealVect a_pos, const Real a_time) const{
  const Real factor = (a_pos - m_seed_pos).vectorLength();

  return m_uniform_density + m_seed_density*exp(-factor*factor/(m_seed_rad*m_seed_rad));
#include "CD_NamespaceFooter.H"
