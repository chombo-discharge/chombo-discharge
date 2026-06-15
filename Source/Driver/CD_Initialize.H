/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Initialize.H
  @brief  Declaration of the initialization functions for chombo-discharge
  @author Robert Marskar
*/

#ifndef CD_INITIALIZE_H
#define CD_INITIALIZE_H

#if defined(_OPENMP) && defined(CH_USE_MEMORY_TRACKING)
#error "OpenMP and memory tracking are not compatible"
#endif

#if defined(CH_USE_PETSC) && !defined(CH_MPI)
#error "USE_PETSC=TRUE requires USE_MPI=TRUE"
#endif

#if defined(CH_USE_PETSC)
#include <petsc.h>
#endif

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Global chombo-discharge ParmParse input parameter object. This must have global lifetime.
*/
extern ParmParse* dischargeParser;

/**
  @brief Global chombo-discharge input file name
*/
extern std::string dischargeInputFile;

/**
  @brief chombo-discharge main initialization function.
  @details This initializes MPI and possibly also PETSc, and builds the ParmParse table used for fetching argument. 
  @param[in] argc The usual number of arguments on the command line.
  @param[in] argv[] Usual arguments.
*/
#if defined(CH_USE_PETSC)
PetscErrorCode
initialize(int argc, char* argv[]);
#else
void
initialize(int argc, char* argv[]);
#endif

/**
  @brief chombo-discharge finalization functions. Exits MPI and writes timer report if called for.
  @return Exit code (0 on success).
*/
#if defined(CH_USE_PETSC)
PetscErrorCode
finalize();
#else
int
finalize();
#endif

#include <CD_NamespaceFooter.H>

#endif
