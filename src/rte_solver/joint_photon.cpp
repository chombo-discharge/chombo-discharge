/*!
  @file   joint_photon.cpp
  @brief  Implementation of joint_photon.H
  @author Robert Marskar
  @date   May 2019
*/

#include "photon.H"
#include "joint_photon.H"


joint_photon::joint_photon() : BinItem() {

}

joint_photon::joint_photon(const Real& a_mass, const RealVect& a_position, const size_t& a_num_photons){
  this->define(a_mass, a_position, a_num_photons);
}

joint_photon::~joint_photon(){

}

void joint_photon::define(const Real& a_mass, const RealVect& a_position, const size_t& a_num_photons){
  m_mass        = a_mass;
  m_position    = a_position;
  m_num_photons = a_num_photons;
}

void joint_photon::add_photon(const photon* const a_photon){
  m_position = m_position*m_mass + a_photon->position()*a_photon->mass();
  m_mass += a_photon->mass();

  m_position /= m_mass;

  m_num_photons++;
}

void joint_photon::clear(){
  m_mass = 0.0;
  m_position = RealVect::Zero;
  m_num_photons = 0;
}


size_t&   joint_photon::num_photons() {return m_num_photons;}
Real&     joint_photon::mass()        {return m_mass;       }
RealVect& joint_photon::position()    {return m_position;   }

const size_t&   joint_photon::num_photons() const {return m_num_photons;}
const Real&     joint_photon::mass()        const {return m_mass;       }
const RealVect& joint_photon::position()    const {return m_position;   }


bool joint_photon::operator == (const joint_photon* a_p) const{ return (*this == *a_p); }
bool joint_photon::operator != (const joint_photon& a_p) const{ return !(*this == a_p); }
bool joint_photon::operator == (const joint_photon& a_p) const{ return ( m_position== a_p.m_position && m_mass == a_p.m_mass);
}

int  joint_photon::size() const { return ( BinItem::size() + sizeof(m_mass) + sizeof(m_position) + sizeof(m_num_photons)); }

void joint_photon::linearOut(void* buf) const {
  Real* buffer = (Real*)buf;
  D_TERM6( *buffer++ = m_position[0];,
	   *buffer++ = m_position[1];,
	   *buffer++ = m_position[2];,
	   *buffer++ = m_position[3];,
	   *buffer++ = m_position[4];,
	   *buffer++ = m_position[5];);
  *buffer++ = m_mass;
  
  size_t* intbuf = (size_t*) buf;
  *intbuf = m_num_photons;
}

void joint_photon::linearIn(void* buf) {
  Real* buffer = (Real*)buf;
  D_TERM6( m_position[0] = *buffer++;,
	   m_position[1] = *buffer++;,
	   m_position[2] = *buffer++;,
	   m_position[3] = *buffer++;,
	   m_position[4] = *buffer++;,
	   m_position[5] = *buffer++;);
  m_mass = *buffer++;
  
  size_t* intbuf = (size_t*) buf;
  m_num_photons = *intbuf;
}
