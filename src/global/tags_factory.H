/*!
  @file   tags_factory.H
  @brief  Declaration of a class for generating tags objects. 
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Segregate implementation and header
*/

#ifndef _TAGS_FACTORY_
#define _TAGS_FACTORY_

#include "tags.H"

#include <DataIndex.H>

/*!
@brief Class for generating tags object in a LevelData setting
*/
class tags_factory : public DataFactory<tags> {
public:

  /*!
    @brief Constructor
  */
  tags_factory(){
  }

  /*!
    @brief Destructor
  */
  virtual ~tags_factory(){

  }

  /*!
    @brief Create new object
  */
  virtual tags* create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const {
    return new tags(a_box);
  }
  
};

#endif
