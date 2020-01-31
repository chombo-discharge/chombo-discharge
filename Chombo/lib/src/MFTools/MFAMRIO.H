#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DFMartin, Wed, April 27, 2005

#ifndef _MFAMRIO_H_
#define _MFAMRIO_H_

#include "LevelData.H"
#include "REAL.H"
#include "Vector.H"
#include "MFCellFAB.H"
#include "IntVectSet.H"
#include <string>
#include "NamespaceHeader.H"

using std::string;
using std::istream;

#ifdef CH_USE_HDF5

///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a single phase in a MFCellFAB.  This is useful for debugging.
    *a_dataPtr is written to a file named fab.hdf5

*/
void
writeMFFAB(const MFCellFAB* a_dataPtr, int a_phase);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single phase of a MFCellFAB.  This is useful for debugging.
*/
void
viewMFFAB(const MFCellFAB* a_dataPtr, int a_phase);


///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running browse with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single phase of a MFCellFAB.  This is useful for debugging.
*/
void
browseMFFAB(const MFCellFAB* a_dataPtr, int a_phase);

///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a MFCellFAB.  This is useful for debugging.  *a_dataPtr is written
    to the file given by a_filename.
*/
void
writeMFFABname(const MFCellFAB* a_dataPtr,
               int a_phase,
               const char*      a_filename);

void
writeMFLevel(const LevelData<MFCellFAB>* a_dataPtr,
             int a_phase);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single phase of a LevelData<MFCellFAB>.  This is useful for debugging.
*/
void
viewMFLevel(const LevelData<MFCellFAB>* a_dataPtr, int a_phase);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single phase of a LevelData<MFCellFAB>.  This is useful for debugging.
*/
void
browseMFLevel(const LevelData<MFCellFAB>* a_dataPtr, int a_phase);

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single phase of a LevelData<MFCellFAB>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeMFLevelname(const LevelData<MFCellFAB>* a_dataPtr,
                 int a_phase,
                 const char*                 a_filename);

///
/** Write a plotfile using writeEBHDF5.  Useful for debugging.  *a_dataPtr is
    written to a file named MFAMR.hdf5.
*/
void
writeMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr, int a_phase);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default.  This is useful for debugging.
*/
void
viewMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr, int a_phase);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default.  This is useful for debugging.
*/
void
browseMFAMR(const Vector<LevelData<MFCellFAB>* >* a_dataPtr, int a_phase);

///
/** Write a plotfile using writeEBHDF5. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeMFAMRname(const Vector<LevelData<MFCellFAB>* >* a_dataPtr,
               int a_phase,
               const char*                           a_filename);


#endif  //CH_USE_HDF5

#include "NamespaceFooter.H"
#endif
