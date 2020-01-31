#!/bin/sh

#
#
#      _______              __
#     / ___/ /  ___  __ _  / /  ___
#    / /__/ _ \/ _ \/  V \/ _ \/ _ \
#    \___/_//_/\___/_/_/_/_.__/\___/
#    Please refer to Copyright.txt, in Chombo's root directory.
#

#
# Collects the source files we need from the CCSE directory and drops them
# all into the CCSE_all directory.  Adds Chombo-style build infrastructure.
#
CCSE_ALL=../CCSE_all
rm -rf $CCSE_ALL
mkdir -p $CCSE_ALL

STRIPCOMMENTS=../../util/chfpp/ccse/stripcomments.py

CPP_FILES="\
 ./amrlib/BCRec.cpp\
 ./BoxLib/FabArray.cpp\
 ./BoxLib/MultiFab.cpp\
 ./BoxLib/BaseFab.cpp\
 ./BoxLib/FArrayBox.cpp\
 ./BoxLib/BoxDomain.cpp\
 ./BoxLib/BoxArray.cpp\
 ./BoxLib/BoxList.cpp\
 ./BoxLib/Orientation.cpp\
 ./BoxLib/IndexType.cpp\
 ./BoxLib/IntVect.cpp\
 ./BoxLib/Box.cpp\
 ./BoxLib/FPC.cpp\
 ./BoxLib/FabConv.cpp\
 ./BoxLib/BLWorkQueue.cpp\
 ./BoxLib/BLThread.cpp\
 ./BoxLib/CArena.cpp\
 ./BoxLib/BArena.cpp\
 ./BoxLib/Arena.cpp\
 ./BoxLib/VisMF.cpp\
 ./BoxLib/ParallelDescriptor.cpp\
 ./BoxLib/DistributionMapping.cpp\
 ./BoxLib/UseCount.cpp\
 ./BoxLib/Utility.cpp\
 ./BoxLib/ParmParse.cpp\
 ./BoxLib/BoxLib.cpp\
 ./mglib/Laplacian.cpp\
 ./mglib/MultiGrid.cpp\
 ./mglib/Mask.cpp\
 ./mglib/LinOp.cpp\
 ./mglib/InterpBndryData.cpp\
 ./mglib/CGSolver.cpp\
 ./mglib/BndryData.cpp\
 ./mglib/ABecLaplacian.cpp\
 ./bndrylib/RealBox.cpp\
 ./bndrylib/CoordSys.cpp\
 ./bndrylib/Geometry.cpp\
 ./bndrylib/BndryRegister.cpp\
 ./bndrylib/FabSet.cpp\
 ./iamrlib/MacBndry.cpp\
 ./iamrlib/MacOperator.cpp"


for cpp in $CPP_FILES; do
    cat $cpp | grep -v MPI_Finalize | sed -e 's/BL_USE_MPI/CH_MPI/g' -e 's/BL_/CH_/g' > $CCSE_ALL/`basename $cpp`
done

for H in `find . -name "*.H"`; do
    cat $H | sed -e 's/BL_USE_MPI/CH_MPI/g' -e 's/BL_/CH_/g' | python $STRIPCOMMENTS > $CCSE_ALL/`basename $H`
done

# These are the ones compiled by mglib/Test.  Don't try compiling every .F file
# in sight, because some (e.g. PROB_2D.F) depend on nonexistent #include files.
F_FILES="\
 ./BoxLib/SPECIALIZE_2D.F\
 ./BoxLib/SPECIALIZE_3D.F\
 ./mglib/MG_2D.F\
 ./mglib/MG_3D.F\
 ./mglib/LP_2D.F\
 ./mglib/LP_3D.F\
 ./mglib/LO_UTIL.F\
 ./mglib/LO_2D.F\
 ./mglib/LO_3D.F\
 ./mglib/INTERPBNDRYDATA_2D.F\
 ./mglib/INTERPBNDRYDATA_3D.F\
 ./mglib/CG_2D.F\
 ./mglib/CG_3D.F\
 ./mglib/ABec_UTIL.F\
 ./mglib/ABec_2D.F\
 ./mglib/ABec_3D.F\
 ./bndrylib/COORDSYS_2D.F\
 ./bndrylib/COORDSYS_3D.F\
 ./iamrlib/MACOPERATOR_2D.F\
 ./iamrlib/MACOPERATOR_3D.F\
 ./iamrlib/MACPROJ_2D.F\
 ./iamrlib/MACPROJ_3D.F"

f_FILES="\
 ./BoxLib/BLBoxLib_F.f\
 ./BoxLib/BLParmParse_F.f\
 ./BoxLib/BLutil_F.f"

CH_TIMER_FILES="\
 ../BoxTools/CH_Timer.H\
 ../BoxTools/ClockTicks.H\
 ../BoxTools/List.H\
 ../BoxTools/ListImplem.H\
 ../BoxTools/MayDay.H\
 ../BoxTools/NamespaceFooter.H\
 ../BoxTools/NamespaceHeader.H\
 ../BoxTools/Pool.H\
 ../BoxTools/Vector.H"


# CCSE quasi-fortran that needs translation to .cpre and then to .f format:
for F in $F_FILES ; do
    BLFname=`basename $F | cut -d. -f1`.BL_F
    cat $F | sed 's/BL_/CH_/g' | tr "\t" "@" | sed 's/@/            /g' | python $STRIPCOMMENTS > $CCSE_ALL/$BLFname
done
# Tabs screw up the .F-->.cpre translation phase, and result in fortran files
# with code in the forbidden first seven columns of your punch card ;-)
# The stripcomments.py is needed because our $(CPP) is called with -C, to preserve
# "//" which is a legitimate Fortran operator.  Unfortunately, CCSE's .F files often
# have /**/-style comments in them.

# CCSE straight Fortran files:
for f in $f_FILES ; do
    cat $f | sed 's/BL_/CH_/g' > $CCSE_ALL/`basename $f`
done

# CHOMBO include files (for timer instrumentation)
for f in $CH_TIMER_FILES ; do
    ln -s $f $CCSE_ALL
done

cp GNUmakefile.ccse $CCSE_ALL/GNUmakefile

cp -Rf facade $CCSE_ALL
