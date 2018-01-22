/*!
  @file cdr_gdnv.H
  @brief Godunov solver for CDR equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _CDR_GDNV_
#define _CDR_GDNV_

#include "cdr_solver.H"

/*!
  @brief Abstract class for cdr-equations. This is (almost) a whole solver. 
*/
class cdr_gdnv : public cdr_solver {
public:

  /*!
    @brief Constructor
  */
  cdr_gdnv();

  /*!
    @brief Destructor
  */
  virtual ~cdr_gdnv();

  /*!
    @brief Query number of ghost cells
  */
  virtual int query_ghost() const;
  
protected:

  /*!
    @brief Godunov face extrapolation method
  */
  virtual void extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state);
};

#endif
