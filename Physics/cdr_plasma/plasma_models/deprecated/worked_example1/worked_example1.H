#include "plasma_kinetics.H"

class worked_example1 : public plasma_kinetics {

public:

  worked_example1(){
    m_num_species = 1;
    m_num_Photons = 0;
    m_species[0] = RefCountedPtr<species> (new worked_example1::my_species());
  }

  ~worked_example1(){}

  Vector<Real> advection_kinetics::compute_cdr_diffusion_coefficients(const Real&         a_time,
								      const RealVect&     a_pos,
								      const RealVect&     a_E,
								      const Vector<Real>& a_cdr_densities) const {
    return Vector<Real>(m_num_species, 0.0);
  }

  Vector<RealVect> advection_kinetics::compute_cdr_velocities(const Real&         a_time,
							      const RealVect&     a_pos,
							      const RealVect&     a_E,
							      const Vector<Real>& a_cdr_densities) const {
    return Vector<RealVect> velo(m_num_species, a_E);
  }

  Vector<Real> advection_kinetics::compute_cdr_source_terms(const Real              a_time,
							    const RealVect&         a_pos,
							    const RealVect&         a_E,
							    const RealVect&         a_gradE,
							    const Vector<Real>&     a_cdr_densities,
							    const Vector<Real>&     a_rte_densities,
							    const Vector<RealVect>& a_grad_cdr) const {
    return Vector<Real>(m_num_species, 0.0);
  }

  Vector<Real> advection_kinetics::compute_cdr_electrode_fluxes(const Real&         a_time,
								const RealVect&     a_pos,
								const RealVect&     a_normal,
								const RealVect&     a_E,
								const Vector<Real>& a_cdr_densities,
								const Vector<Real>& a_cdr_velocities,
								const Vector<Real>& a_cdr_gradients,
								const Vector<Real>& a_rte_fluxes,
								const Vector<Real>& a_extrap_cdr_fluxes) const {
    Vector<Real> flux(m_num_species, 0.0);
    flux[0] = Max(0.0, a_extrap_cdr_fluxes[0]);
    return flux;
  }

  Vector<Real> advection_kinetics::compute_cdr_dielectric_fluxes(const Real&         a_time,
								 const RealVect&     a_pos,
								 const RealVect&     a_normal,
								 const RealVect&     a_E,
								 const Vector<Real>& a_cdr_densities,
								 const Vector<Real>& a_cdr_velocities,
								 const Vector<Real>& a_cdr_gradients,
								 const Vector<Real>& a_rte_fluxes,
								 const Vector<Real>& a_extrap_cdr_fluxes) const {
    Vector<Real> flux(m_num_species, 0.0);
    flux[0] = Max(0.0, a_extrap_cdr_fluxes[0]);
    return flux;
  }

  Vector<Real> advection_kinetics::compute_rte_source_terms(const Real&         a_time,
							    const RealVect&     a_pos,
							    const RealVect&     a_E,
							    const Vector<Real>& a_cdr_densities) const {
    return Vector<Real>(m_num_Photons, 0.0);
  }

  Real advection_kinetics::initial_sigma(const Real a_time, const RealVect& a_pos) const {
    return 0.0;
  }

  class my_species : public species {
  public:
    my_species(){
      m_name      = "my_species";
      m_isDiffusive = false;
      m_isMobile    = true;
      m_chargeNumber    = 1;

      {
	ParmParse pp("my_kinetics");

	Vector<Real> arr;
	pp.get("width", m_L);
	pp.getarr("center", arr, 0, SpaceDim);
	m_center = RealVect(D_DECL(arr[0], arr[1], arr[2]));
      }
    }

    ~my_species(){

    }

    Real initialData(const RealVect a_pos, const Real a_time) const{
      const Real diff  = (a_pos - m_center).vectorLength();;
      const Real diff2 = diff*diff;
      return exp(-diff2/(2*m_L*m_L));
    }

  protected:

    RealVect m_center;
    Real m_L;
  };
#include "CD_NamespaceFooter.H"
