/*!
  @file   tags.H
  @brief  Declaration of a class that holds tags that can be communicated (unlike IntVectSet)
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Segregate implementation and header
*/

#ifndef _tags_
#define _tags_

#include <Interval.H>
#include <Box.H>
#include <DenseIntVectSet.H>

/*!
  @brief Class that holds an IntVectSet that can be communicated
*/
class tags {
public:

  /*!
    @brief Null constructor
  */
  tags();
  
  /*!
    @brief Constructor
  */
  tags(const Box& a_box);

  /*!
    @brief Destructor
  */
  virtual ~tags();

  /*!
    @brief Get IVS
  */
  DenseIntVectSet& get_ivs();

  /*!
    @brief Get IVS
  */
  const DenseIntVectSet& get_ivs() const;

  /*!
    @brief Copy function
  */
  virtual void copy(const Box&        a_region_from,
		    const Interval&   a_dst_interval,
		    const Box&        a_region_to,
		    const tags&       a_src,
		    const Interval&   a_src_interval);

  /*!
    @brief preAllocatable memory or not
  */
  const static int preAllocatable();

  /*!
    @brief Return size in bytes
  */
  int size(const Box& a_box, const Interval& a_comps) const;

  /*!
    @brief Return size of object
  */
  int linearSize() const;

  /*!
    @brief Linearization function
  */
  void linearIn(const void* const a_in_buffer, const Box& R, const Interval& a_comps);

  /*!
    @brief Linearization function
  */
  void linearOut(void* const a_out_buffer, const Box& R, const Interval& a_comps) const;

protected:

  DenseIntVectSet m_ivs;
};
#endif

