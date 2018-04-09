/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "RealVect.H"
#include "RefCountedPtr.H"

class dcel_vert;
class dcel_poly;

/*!
  @brief Edge class for handling surface tesselations with DCEL meshes
*/
class dcel_edge {
public:

  /*!
    @brief Constructor
  */
  dcel_edge();

  /*!
    @brief Destructor
  */
  ~dcel_edge();

  /*!
    @brief Define function
  */
  virtual void define(const RefCountedPtr<dcel_vert>& a_vert,
		      const RefCountedPtr<dcel_edge>& a_pair,
		      const RefCountedPtr<dcel_edge>& a_next,
		      const RefCountedPtr<dcel_edge>& a_prev,
		      const RealVect                  a_normal = RealVect::Zero);


    /*!
    @brief Set polygon
  */
  virtual void set_poly(const RefCountedPtr<dcel_poly>& a_poly);

  /*!
    @brief Set vertex
  */
  virtual void set_vert(const RefCountedPtr<dcel_vert>& a_vert);

  /*!
    @brief Set pair
  */
  virtual void set_pair(const RefCountedPtr<dcel_edge>& a_pair);

  /*!
    @brief Set next edge
  */
  virtual void set_next(const RefCountedPtr<dcel_edge>& a_next);

  /*!
    @brief Set previous
  */
  virtual void set_prev(const RefCountedPtr<dcel_edge>& a_prev);
  
  /*!
    @brief Set normal vector
  */
  virtual void set_normal(const RealVect a_normal);
  
  /*!
    @brief Get vertex
  */
  virtual const RefCountedPtr<dcel_vert>& get_vert() const;

    /*!
    @brief Get vertex
  */
  virtual  RefCountedPtr<dcel_vert>& get_vert();

  /*!
    @brief Get other vertex
  */
  virtual const RefCountedPtr<dcel_vert>& get_other_vert() const;

  /*!
    @brief Get other vertex
  */
  virtual RefCountedPtr<dcel_vert>& get_other_vert();

  /*!
    @brief Get the other half edge
  */
  virtual const RefCountedPtr<dcel_edge>& get_pair() const;

  /*!
    @brief Get the other half edge
  */
  virtual RefCountedPtr<dcel_edge>& get_pair();

  /*!
    @brief Get the previous half edge
  */
  virtual const RefCountedPtr<dcel_edge>& get_prev() const;

    /*!
    @brief Get the previous half edge
  */
  virtual RefCountedPtr<dcel_edge>& get_prev();

  /*!
    @brief Get the next half edge
  */
  virtual const RefCountedPtr<dcel_edge>& get_next() const;

  /*!
    @brief Get the next half edge
  */
  virtual RefCountedPtr<dcel_edge>& get_next();

  /*!
    @brief Get polygon
  */
  virtual const RefCountedPtr<dcel_poly>& get_poly() const;

  /*!
    @brief Get polygon
  */
  virtual RefCountedPtr<dcel_poly>& get_poly();
  
  /*!
    @brief Get pseudo-normal
  */
  virtual RealVect get_normal() const;

  /*!
    @brief Get the signed distance between this edge and a point
  */
  virtual Real signed_distance(const RealVect a_x0) const;

protected:

  /*!
    @brief Edge pseudo-normal
  */
  RealVect m_normal;

  /*!
    @brief Vertex at end of edge
  */
  RefCountedPtr<dcel_vert> m_vert; 

  /*!
    @brief Edge pair
  */
  RefCountedPtr<dcel_edge> m_pair;

  /*!
    @brief Previous edge
  */
  RefCountedPtr<dcel_edge> m_prev; 

  /*!
    @brief Next edge
  */
  RefCountedPtr<dcel_edge> m_next;

  /*!
    @brief Reference to polygon
  */
  RefCountedPtr<dcel_poly> m_poly;
};

#endif
