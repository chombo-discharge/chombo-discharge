#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERUTIL_H_
#define _FOURTHORDERUTIL_H_

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "MayDay.H"

#include "NamespaceHeader.H"

/// utility functions for 4th-order computations

///
/** convert from cell centers to face centers:
    phi[i+1/2] = 9/16*(phi[i] + phi[i+1]) - 1/16*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFaceCenters(LevelData<FluxBox>& a_faceData,
                                  const LevelData<FArrayBox>& a_cellData);


///
/** convert from cell centers to face centers:
    phi[i+1/2] = 9/16*(phi[i] + phi[i+1]) - 1/16*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFaceCenters(FluxBox& a_faceData,
                                  const FArrayBox& a_cellData);


///
/** convert from cell centers to face centers:
    phi[i+1/2] = 9/16*(phi[i] + phi[i+1]) - 1/16*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFaceCenters(FArrayBox& a_faceData,
                                  const FArrayBox& a_cellData,
                                  const Box& a_faceBox,
                                  int a_dir);

///
/** convert from cell averages to face averages:
    phi[i+1/2] = 7/12*(phi[i] + phi[i+1]) - 1/12*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFace(LevelData<FluxBox>& a_faceData,
                           const LevelData<FArrayBox>& a_cellData);


///
/** convert from cell averages to face averages:
    phi[i+1/2] = 7/12*(phi[i] + phi[i+1]) - 1/12*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFace(FluxBox& a_faceData,
                           const FArrayBox& a_cellData);


///
/** convert from cell averages to face averages:
    phi[i+1/2] = 7/12*(phi[i] + phi[i+1]) - 1/12*(phi[i-1] + phi[i+2])
 */
void fourthOrderCellToFace(FArrayBox& a_faceData,
                           const FArrayBox& a_cellData,
                           const Box& a_faceBox,
                           int a_dir);

///
/**  given face-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order face averages of phi
     (or decrement, a_sgn = -1, to compute face-centered from face-averaged)
*/
void fourthOrderAverage(LevelData<FluxBox>& a_phi,
                        int a_sgn = 1);

///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
*/
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        int a_sgn = 1);


///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
     using one-sided corrections at domain boundaries
*/
void fourthOrderAverage(LevelData<FArrayBox>& a_phi,
                        const ProblemDomain&  a_domain,
                        int a_sgn = 1);


///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
*/
void fourthOrderAverageCell(FArrayBox& a_phi,
                            int a_sgn = 1);

///
/**  given cell-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order cell averages of phi
     (or decrement, a_sgn = -1, to compute cell-centered from cell-averaged)
     using one-sided corrections at domain boundaries
*/
void fourthOrderAverageCell(FArrayBox& a_phi,
                            const ProblemDomain&  a_domain,
                            const Box&            a_bx,
                            int a_sgn = 1);

///
/**  given face-centered phi, increment with Laplacian(phi) correction
     in place to compute 4th-order face averages of phi
     (or decrement, a_sgn = -1, to compute face-centered from face-averaged)
*/
void fourthOrderAverageFace(FluxBox& a_phi,
                            int a_sgn = 1);

///
/**  replace a_phiAvg by
     a_phiCen + a_sgn * D^2_(a_dir)(a_phiAvg)
 */
void fourthOrderAverageCenterFace(FArrayBox& a_phiAvg,
                                  const FArrayBox& a_phiCen,
                                  int a_dir,
                                  int a_sgn = 1);

/// compute fourth-order approximation to u*v
/**
   u and v are assumed to already be fourth-order face-averages
 */
void
fourthOrderMultFace(FluxBox& a_uTimesV,
                    const FluxBox& a_u,
                    const FluxBox& a_v,
                    int a_sgn = 1);

/// compute fourth-order approximation to u*v
/**
   u and v are assumed to already be fourth-order face-averages
 */
void
fourthOrderMultFace(FArrayBox& a_uTimesV,
                    const FArrayBox& a_u,
                    const FArrayBox& a_uForDx,
                    const FArrayBox& a_v,
                    const FArrayBox& a_vForDx,
                    int a_faceDir,
                    int a_sgn = 1);

/// compute fourth-order approximation to u*v
/**
   u and v are assumed to already be fourth-order face-averages
 */
void
fourthOrderMult(LevelData<FluxBox>& a_uTimesV,
                const LevelData<FluxBox>& a_u,
                const LevelData<FluxBox>& a_v,
                int a_sgn = 1);

/// compute fourth-order approximation to u*v
/**
   u and v are assumed to already be fourth-order cell-averages
 */
void
fourthOrderCellProd(LevelData<FArrayBox>& a_uTimesV,
                    const LevelData<FArrayBox>& a_u,
                    const LevelData<FArrayBox>& a_v);

/// Computes <FG> from <F> and <G> using 1-sided derivatives at boundaries
void fourthOrderCellProd(FArrayBox&           a_FG,
                         const FArrayBox&     a_F,
                         const FArrayBox&     a_G,
                         const Box&           a_box,
                         const ProblemDomain& a_domain,
                         const bool           a_oneSidedDG = false,
                         Interval             a_fixedDims = Interval() );

/// this is just a basic wrapper around the INCREMENTGRADPROD fortran
void
incrementGradProduct(FArrayBox& a_gradProduct,
                     const FArrayBox& a_u,
                     const FArrayBox& a_v,
                     const Box& a_region,
                     Real a_dx,
                     Real a_factor,
                     int a_dir);

/// Fourth-order extrapolate cell-averaged interior data to two ghost cell layers at the domain boundary
/**
 * At the domain boundary, perform a fourth-order extrapolation of cell-averaged interior data
 * to fill two layers of ghost cells (we assume an exchange will be performed elsewhere
 * to fill interior ghost cells).  Among other possible uses, this enables the
 * fourthOrderCellToFace utility to then be called without worrying about boundaries.
 */
void
fourthOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const ProblemDomain& a_domain,
                                  const Box&           a_interiorbox);

/// Fill interior ghost cells and boundary ghost cell via extrapolation of cell-averaged data
/**
 * Fill interior ghost cells by copying and boundary ghost cells by extrapolation of cell-averaged data
 */
void
fourthOrderCellGhostFill(LevelData<FArrayBox>& a_data,
                         const ProblemDomain&  a_domain);

// Fill width-two boxes on low and high sides of a_srcbox (in all directions)
// with cell-averaged values extrapolated to fourth-order from values in a_srcbox.
/**
 * Perform a fourth-order extrapolation of cell-averaged a_data values in a_srcbox to set
 * the values of a_data in the adjacent width-two boxes on the low and high
 * sides of a_srcbox.
 */
void
fourthOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox);

// Fill width-one boxes on low and high sides of a_srcbox (in all directions)
// with values extrapolated to second-order from values in a_srcbox.
/**
 * Perform a second-order extrapolation of a_data values in a_srcbox to set
 * the values of a_data in the adjacent width-one boxes on the low and high
 * sides of a_srcbox in all dimensions.
 */
void
secondOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox);

// Fill width-one boxes on low and high sides of a_srcbox in dimension a_dir
// with values extrapolated to second-order from values in a_srcbox.
/**
 * Perform a second-order extrapolation of a_data values in a_srcbox to set
 * the values of a_data in the adjacent width-one boxes on the low and high
 * sides of a_srcbox in dimension a_dir.
 */
void
secondOrderCellExtrap(FArrayBox& a_data,
                      const Box& a_srcbox,
                      int a_dir);

/// Second-order extrapolate at domain boundary (in all directions)
/**
 * Perform a second-order extrapolation of a_data values in a_srcbox to set
 * the values of a_data in the adjacent width-one boxes on the low and high
 * sides of a_srcbox.  The extrapolation is only performed if the adjacent
 * width-one boxes lie outside the domain box (i.e., are physical
 * boundary ghost cells).
 */
void
secondOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const Box&           a_interiorbox,
                                  const ProblemDomain& a_domain);

/// Second-order extrapolate data in transverse directions
/**
 * Perform a second-order extrapolation of the input a_data from
 * the a_srcbox box to adjacent destination boxes of width one on
 * the low and high sides of srcbox in the directions transverse
 * to a_dir.  This can be used, for example, to fill ghost cells
 * prior to computing the transverse gradients (relative to the
 * a_dir direction) on cell faces via the product formula (for
 * computing the fourth-order average of a product in terms of
 * the averages of each factor).
 */
void
secondOrderTransExtrap(FArrayBox& a_data,
                       const int  a_dir,
                       const Box& a_srcbox);

/// Second-order extrapolate data in transverse directions at domain boundary
/**
 * Perform a second-order extrapolation of the input a_data from
 * the a_interiorbox box to adjacent destination boxes of width one on
 * the low and high sides of a_interiorbox in the directions transverse
 * to dir.  The extrapolation is only performed if the adjacent
 * width-one boxes lie outside the domain box (i.e., are physical
 * boundary ghost cells).  This can be used to fill ghost cells
 * prior to computing the transverse gradients (relative to the
 * a_dir direction) on cell faces via the product formula (for
 * computing the fourth-order average of a product in terms of
 * the averages of each factor).
 */
void
secondOrderTransExtrapAtDomainBdry(FArrayBox& a_data,
                                   const int  a_dir,
                                   const Box& a_interiorbox,
                                   const ProblemDomain& a_domain);

/// Second-order extrapolate data to faces in transverse directions at domain boundary
/**
 * Perform a second-order extrapolation of the input a_data from
 * faces of a_interiorbox to adjacent destination boxes of width one on
 * the low and high sides of a_interiorbox in the directions transverse
 * to dir.  The extrapolation is only performed if the adjacent
 * width-one boxes lie outside the domain box (i.e., are physical
 * boundary ghost faces).
 */
void
secondOrderTransExtrapFacesAtDomainBdry(FArrayBox& a_data,
                                        const int  a_dir,
                                        const Box& a_interiorbox,
                                        const ProblemDomain& a_domain);

/// Zeroth-order extrapolate at domain boundary (in all directions)
/**
 * Perform a zeroth-order extrapolation of a_data values in a_srcbox to set
 * the values of a_data in the adjacent width-one boxes on the low and high
 * sides of a_srcbox.  The extrapolation is only performed if the adjacent
 * width-one boxes lie outside the domain box (i.e., are physical
 * boundary ghost cells).
 */
void
zerothOrderCellExtrapAtDomainBdry(FArrayBox&           a_data,
                                  const Box&           a_interiorbox,
                                  const ProblemDomain& a_domain);

/// Determine transverse components of N by averaging from nearby faces
void averageTransverseN(FArrayBox&       a_Nface,
                        const FArrayBox& a_NtransverseFace,
                        const Box&       a_box,
                        const int        a_faceDir,
                        const int        a_transverseDir);

/// Convert from <FG> and <G> to <F>.
/**
 */
void cellFGToCellF(FArrayBox& a_F,
                   const FArrayBox& a_FG, // may have multiple components
                   const FArrayBox& a_G, // 1 component
                   const Box& a_box); // where we want to compute F

/// Computes \<F\> from \<FG\> and \<G\> using 1-sided derivatives at boundaries
void cellFGToCellF(FArrayBox&           a_F,
                   const FArrayBox&     a_FG,
                   const FArrayBox&     a_G,
                   const Box&           a_box,
                   const ProblemDomain& a_domain,
                   const bool           a_oneSidedDG = false,
                   Interval             a_fixedDims = Interval() );

#include "NamespaceFooter.H"

#endif
