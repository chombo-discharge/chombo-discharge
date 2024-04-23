/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBMultigridInterpolator.cpp
  @brief  Implementation of CD_EBMultigridInterpolator.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

EBMultigridInterpolator::EBMultigridInterpolator()
{
  CH_TIME("EBMultigridInterpolator::EBMultigridInterpolator");
}

EBMultigridInterpolator::~EBMultigridInterpolator()
{
  CH_TIME("EBMultigridInterpolator::~EBMultigridInterpolator");
}

#include <CD_NamespaceFooter.H>
