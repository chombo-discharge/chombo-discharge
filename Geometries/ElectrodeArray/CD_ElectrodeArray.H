/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ElectrodeArray.H
  @brief  Electrode array geometry.
  @author Robert Marskar
*/

#ifndef CD_ELECTRODEARRAY_H
#define CD_ELECTRODEARRAY_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that defines an array of electrodes.
*/
class ElectrodeArray : public ComputationalGeometry
{
public:
  ElectrodeArray();
  ~ElectrodeArray();
};

#include <CD_NamespaceFooter.H>

#endif
