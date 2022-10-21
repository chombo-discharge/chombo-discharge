/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Average.H
  @brief  Declaration of averaging methods.
  @author Robert Marskar
*/

#ifndef CD_Average_H
#define CD_Average_H

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Various averaging methods. 
  @details 'Arithmetic'/'Harmonic' mean what you think, while 'Conservative' also incorporate volume, face, 
  and EB face fractions such that the total quantity is conserved. 
*/
enum class Average
{
  Arithmetic,
  Harmonic,
  Conservative,
  Geometric
};

#include <CD_NamespaceFooter.H>

#endif
