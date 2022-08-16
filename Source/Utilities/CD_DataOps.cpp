/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataOps.cpp
  @brief  Implementation of CD_DataOps.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBArith.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <CCProjectorF_F.H>
#include <EBLevelDataOpsF_F.H>
#include <PolyGeom.H>

// Our includes
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_ParallelOps.H>
#include <CD_MultifluidAlias.H>
#include <CD_NamespaceHeader.H>

void
DataOps::averageCellVectorToFaceScalar(EBAMRFluxData&               a_faceData,
                                       const EBAMRCellData&         a_cellData,
                                       const Vector<ProblemDomain>& a_domains)
{
  CH_TIME("DataOps::averageCellVectorToFaceScalar(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_faceData.size(); lvl++) {

    CH_assert(a_faceData[lvl]->nComp() == 1);
    CH_assert(a_cellData[lvl]->nComp() == SpaceDim);

    DataOps::averageCellVectorToFaceScalar(*a_faceData[lvl], *a_cellData[lvl], a_domains[lvl]);
  }
}

void
DataOps::averageCellVectorToFaceScalar(LevelData<EBFluxFAB>&       a_faceData,
                                       const LevelData<EBCellFAB>& a_cellData,
                                       const ProblemDomain&        a_domain)
{
  CH_TIME("DataOps::averageCellVectorToFaceScalar");

  CH_assert(a_faceData.nComp() == 1);
  CH_assert(a_cellData.nComp() == SpaceDim);

  // TLDR: This computes the face-centered data from averages of the cell-centered data. Note that a single layer of
  //       ghost faces are filled (the ones that are "tangential" to the box).

  for (DataIterator dit = a_faceData.dataIterator(); dit.ok(); ++dit) {
    EBFluxFAB&       fluxData = a_faceData[dit()];
    const EBCellFAB& cellData = a_cellData[dit()];
    const EBISBox&   ebisbox  = cellData.getEBISBox();
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const Box&       box      = a_cellData.disjointBoxLayout().get(dit());

    constexpr int numTanGhost = 1;

    for (int dir = 0; dir < SpaceDim; dir++) {
      const int faceDir  = dir;
      const int cellComp = dir;
      const int faceComp = 0;

      EBLevelDataOps::averageCellToFace(fluxData[dir],
                                        cellData,
                                        ebgraph,
                                        box,
                                        numTanGhost,
                                        faceDir,
                                        a_domain,
                                        cellComp,
                                        faceComp);
    }
  }
}

void
DataOps::averageCellScalarToFaceScalar(EBAMRFluxData&               a_faceData,
                                       const EBAMRCellData&         a_cellData,
                                       const Vector<ProblemDomain>& a_domains)
{
  CH_TIME("DataOps::averageCellScalarToFaceScalar(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_faceData.size(); lvl++) {

    CH_assert(a_faceData[lvl]->nComp() == 1);
    CH_assert(a_cellData[lvl]->nComp() == 1);

    DataOps::averageCellScalarToFaceScalar(*a_faceData[lvl], *a_cellData[lvl], a_domains[lvl]);
  }
}

void
DataOps::averageCellScalarToFaceScalar(LevelData<EBFluxFAB>&       a_faceData,
                                       const LevelData<EBCellFAB>& a_cellData,
                                       const ProblemDomain&        a_domain)
{
  CH_TIME("DataOps::averageCellScalarToFaceScalar");

  CH_assert(a_faceData.nComp() == 1);
  CH_assert(a_cellData.nComp() == 1);

  // TLDR: This computes the face-centered data from averages of the cell-centered data. Note that a single layer of
  //       ghost faces are filled (the ones that are "tangential" to the box).

  for (DataIterator dit = a_faceData.dataIterator(); dit.ok(); ++dit) {
    EBFluxFAB&       fluxData = a_faceData[dit()];
    const EBCellFAB& cellData = a_cellData[dit()];
    const EBISBox&   ebisbox  = cellData.getEBISBox();
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const Box&       box      = a_cellData.disjointBoxLayout().get(dit());

    constexpr int numTanGhost = 1;

    for (int dir = 0; dir < SpaceDim; dir++) {
      const int faceDir  = dir;
      const int cellComp = 0;
      const int faceComp = 0;

      EBLevelDataOps::averageCellToFace(fluxData[dir],
                                        cellData,
                                        ebgraph,
                                        box,
                                        numTanGhost,
                                        faceDir,
                                        a_domain,
                                        cellComp,
                                        faceComp);
    }
  }
}

void
DataOps::averageCellToFace(EBAMRFluxData&               a_faceData,
                           const EBAMRCellData&         a_cellData,
                           const Vector<ProblemDomain>& a_domains)
{
  CH_TIME("DataOps::averageCellToFace(EBAMRFluxData");

  for (int lvl = 0; lvl < a_faceData.size(); lvl++) {
    DataOps::averageCellToFace(*a_faceData[lvl], *a_cellData[lvl], a_domains[lvl]);
  }
}

void
DataOps::averageCellToFace(LevelData<EBFluxFAB>&       a_faceData,
                           const LevelData<EBCellFAB>& a_cellData,
                           const ProblemDomain&        a_domain)
{
  CH_TIME("DataOps::averageCellToFace(LD<EBFluxFAB)");

  CH_assert(a_faceData.nComp() == a_cellData.nComp());

  // TLDR: This routine computes face-centered data as an arithmetic average of cell-centered data, including cut-cell faces. The kernels are simple arithmetic averages,
  //       only modified on domain boundaries. If a domain face is covered, we just set the covered value to the value in the covered cell just inside the domain. Otherwise,
  //       we iterate through the domain faces and extrapolate (to first order) the cell-centered data to the domain face. If there are not enough cells for extrapolation,
  //       we set the value on the face equal to the value in the vof just inside the domain.

  const int numComp = a_faceData.nComp();

  const DisjointBoxLayout& dbl = a_faceData.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const EBCellFAB& cellData = a_cellData[dit()];
    const EBISBox&   ebisbox  = cellData.getEBISBox();
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const Box        cellBox  = dbl[dit()];

    const BaseFab<Real>& regCellData = cellData.getSingleValuedFAB();

    const IntVectSet ivsIrreg = ebisbox.getIrregIVS(cellBox);

    // Go through faces oriented in direction dir.
    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB&     faceData    = a_faceData[dit()][dir];
      BaseFab<Real>& regFaceData = faceData.getSingleValuedFAB();

      // Face-centered box and irregular region. Note that the irregular kernel only does interior faces. We patch up those later.
      const Box    faceBox = surroundingNodes(cellBox, dir);
      FaceIterator faceit(ivsIrreg, ebgraph, dir, FaceStop::SurroundingNoBoundary);

      // Do each component. Regular cells are done in Fortran while cut-cell faces
      for (int comp = 0; comp < numComp; comp++) {

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          regFaceData(iv, comp) = 0.5 * (regCellData(iv, comp) + regCellData(iv - BASISV(dir), comp));
        };

        // Irregular kernel.
        auto irregularKernel = [&](const FaceIndex& face) -> void {
          const VolIndex& vofLo = face.getVoF(Side::Lo);
          const VolIndex& vofHi = face.getVoF(Side::Hi);

          faceData(face, comp) = 0.5 * (cellData(vofLo, comp) + cellData(vofHi, comp));
        };

        // Run the kernels.
        BoxLoops::loop(faceBox, regularKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }

      // Now fix up domain faces
      for (SideIterator sit; sit.ok(); ++sit) {

        // This box is on the outside of cellBox. I.e. it can also end up being outside of the domain.
        const Box outsideBox = adjCellBox(cellBox, dir, sit(), 1);

        if (!a_domain.contains(outsideBox)) {

          // Define kernel regions. This is a box which is the first strip of cells on the inside of the domain.
          const Box    insideBox = adjCellBox(cellBox, dir, sit(), -1);
          FaceIterator faceit(IntVectSet(insideBox), ebgraph, dir, FaceStop::AllBoundaryOnly);

          // Go through cells in that box and set the covered faces to something reasonable.
          auto regularKernel = [&](const IntVect& iv) -> void {
            // Find face index corresponding to cell index.
            IntVect ivFace = iv;
            if (sit() == Side::Hi) {
              ivFace.shift(BASISV(dir));
            }

            // Put reasonable values in face at domain boundaries.
            for (int comp = 0; comp < numComp; comp++) {
              regFaceData(ivFace, comp) = regCellData(iv, comp);
            }
          };

          // Irregular kernel.
          auto irregularKernel = [&](const FaceIndex& face) -> void {
            // Get the vof just inside the domain on this face.
            const VolIndex& vof = face.getVoF(flip(sit()));

            const Vector<FaceIndex>& furtherFaces = ebisbox.getFaces(vof, dir, flip(sit()));

            for (int comp = 0; comp < numComp; comp++) {

              // If we have more faces available, extrapolate from inside the domain. Otherwise, set the face value
              // to the vof just insid ethe domain.

              if (furtherFaces.size() == 0) { // Set the value to be the value immediately inside the face.
                faceData(face, comp) = cellData(vof, comp);
              }
              else { // Extrapolate.
                const Real valClose = cellData(vof, comp);

                Real valFar = 0.0;
                for (int i = 0; i < furtherFaces.size(); i++) {
                  const VolIndex& farVoF = furtherFaces[i].getVoF(flip(sit()));
                  valFar += cellData(farVoF, comp);
                }
                valFar /= furtherFaces.size();

                faceData(face, comp) = 0.5 * (3.0 * valClose - valFar);
              }
            }
          };

          // Execute the kernels.
          BoxLoops::loop(insideBox, regularKernel);
          BoxLoops::loop(faceit, irregularKernel);
        }
      }
    }
  }
}

void
DataOps::averageFaceToCell(EBAMRCellData&               a_cellData,
                           const EBAMRFluxData&         a_faceData,
                           const Vector<ProblemDomain>& a_domains)
{
  CH_TIME("DataOps::averageFaceToCell(EBAMRCellData, EBAMRFluxData, Vector<ProblemDomain>)");

  for (int lvl = 0; lvl < a_faceData.size(); lvl++) {
    DataOps::averageFaceToCell(*a_cellData[lvl], *a_faceData[lvl], a_domains[lvl]);
  }
}

void
DataOps::averageFaceToCell(LevelData<EBCellFAB>&       a_cellData,
                           const LevelData<EBFluxFAB>& a_fluxData,
                           const ProblemDomain&        a_domain)
{
  CH_TIME("DataOps::averageFaceToCell(LD<EBCellFAB>, LD<EBCellFAB>, ProblemDomain)");

  CH_assert(a_fluxData.nComp() == a_cellData.nComp());

  const int numComp = a_cellData.nComp();

  for (DataIterator dit = a_cellData.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       cellData = a_cellData[dit()];
    const EBFluxFAB& fluxData = a_fluxData[dit()];
    const EBISBox&   ebisbox  = cellData.getEBISBox();
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();

    // Regions for the kernels.
    const Box&        cellBox = a_cellData.disjointBoxLayout()[dit()];
    const IntVectSet& ivs     = ebisbox.getIrregIVS(cellBox);
    VoFIterator       vofit(ivs, ebgraph);

    // Component loop.
    for (int comp = 0; comp < numComp; comp++) {

      // Hooks for single-valued data.
      BaseFab<Real>& cellreg = cellData.getSingleValuedFAB();

      // Regular kernel. We compute phi(cell) = sum(phi(face))/sum(faces).
      auto regularKernel = [&](const IntVect& iv) -> void {
        constexpr Real factor = 1. / (2 * SpaceDim);

        cellreg(iv, comp) = 0.0;
        for (int dir = 0; dir < SpaceDim; dir++) {
          const BaseFab<Real>& facereg = fluxData[dir].getSingleValuedFAB();

          cellreg(iv, comp) += factor * (facereg(iv + BASISV(dir), comp) + facereg(iv, comp));
        }
      };

      // Irregular kernel. Same as the above except that we need to explicitly get the faces. Note that this is an arithmetic average,
      // i.e. not weighted by the face area fractions.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        cellData(vof, comp) = 0.0;

        int numFaces = 0;
        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const std::vector<FaceIndex>& faces = ebisbox.getFaces(vof, dir, sit()).stdVector();

            for (const auto& face : faces) {
              cellData(vof, comp) += fluxData[dir](face, comp);
            }

            numFaces += faces.size();
          }
        }

        cellData(vof, comp) *= 1. / numFaces;
      };

      // Run the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::dotProduct(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2)
{
  CH_TIME("DataOps::dotProduct(mfamrcell)");

  for (int lvl = 0; lvl < a_result.size(); lvl++) {
    DataOps::dotProduct(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void
DataOps::dotProduct(LevelData<MFCellFAB>&       a_result,
                    const LevelData<MFCellFAB>& a_data1,
                    const LevelData<MFCellFAB>& a_data2)
{
  CH_TIME("DataOps::dotProduct(LD<MFCellFAB>)");

  CH_assert(a_data2.nComp() == a_data1.nComp());
  CH_assert(a_result.nComp() == 1);

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit) {
    MFCellFAB&       result = a_result[dit()];
    const MFCellFAB& data1  = a_data1[dit()];
    const MFCellFAB& data2  = a_data2[dit()];
    const Box&       box    = a_result.disjointBoxLayout().get(dit());

    for (int i = 0; i < result.numPhases(); i++) {
      EBCellFAB&       resultPhase = result.getPhase(i);
      const EBCellFAB& data1Phase  = data1.getPhase(i);
      const EBCellFAB& data2Phase  = data2.getPhase(i);

      DataOps::dotProduct(resultPhase, data1Phase, data2Phase, box);
    }
  }
}

void
DataOps::dotProduct(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2)
{
  CH_TIME("DataOps::dotProduct(EBAMRCellData)");

  CH_assert(a_data1[0]->nComp() == a_data2[0]->nComp());
  CH_assert(a_result[0]->nComp() == 1);

  for (int lvl = 0; lvl < a_result.size(); lvl++) {
    DataOps::dotProduct(*a_result[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void
DataOps::dotProduct(LevelData<EBCellFAB>&       a_result,
                    const LevelData<EBCellFAB>& a_data1,
                    const LevelData<EBCellFAB>& a_data2)
{
  CH_TIME("DataOps::dotProduct(LD<EBCellFAB>)");

  CH_assert(a_data1.nComp() == a_data2.nComp());
  CH_assert(a_result.nComp() == 1);

  for (DataIterator dit = a_result.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       result = a_result[dit()];
    const EBCellFAB& data1  = a_data1[dit()];
    const EBCellFAB& data2  = a_data2[dit()];
    const Box&       box    = a_result.disjointBoxLayout().get(dit());

    DataOps::dotProduct(result, data1, data2, box);
  }
}

void
DataOps::dotProduct(EBCellFAB& a_result, const EBCellFAB& a_data1, const EBCellFAB& a_data2, const Box& a_box)
{
  CH_TIME("DataOps::dotProduct(EBCellFAB)");

  CH_assert(a_data1.nComp() == a_data2.nComp());
  CH_assert(a_result.nComp() == 1);

  const int dstComp = 0;
  const int numComp = a_data1.nComp();

  // Regular cells
  BaseFab<Real>&       resultReg = a_result.getSingleValuedFAB();
  const BaseFab<Real>& data1Reg  = a_data1.getSingleValuedFAB();
  const BaseFab<Real>& data2Reg  = a_data2.getSingleValuedFAB();

  // Kernel region.
  const EBISBox&    ebisbox = a_result.getEBISBox();
  const EBGraph&    ebgraph = ebisbox.getEBGraph();
  const IntVectSet& ivs     = ebisbox.getIrregIVS(a_box);
  VoFIterator       vofit(ivs, ebgraph);

  // Regular dot product kernel.
  auto regularKernel = [&](const IntVect& iv) -> void {
    resultReg(iv, dstComp) = 0.0;
    for (int comp = 0; comp < numComp; comp++) {
      resultReg(iv, dstComp) += data1Reg(iv, comp) * data2Reg(iv, comp);
    }
  };

  // Cut-cell kernel. Same as the above but the data access is obviously different (because of multi-valued cells).
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_result(vof, dstComp) = 0.0;

    for (int comp = 0; comp < numComp; comp++) {
      a_result(vof, dstComp) += a_data1(vof, comp) * a_data2(vof, comp);
    }
  };

  // Run the kernels.
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

void
DataOps::incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(MFAMRCellData)");

  CH_assert(a_lhs[0]->nComp() == a_rhs[0]->nComp());

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(LD<MFCellFAB)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  MFLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void
DataOps::incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(EBAMRCellData)");

  CH_assert(a_lhs[0]->nComp() == a_rhs[0]->nComp());

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(LD<EBCellFAB)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void
DataOps::plus(EBAMRCellData&       a_lhs,
              const EBAMRCellData& a_rhs,
              const int            a_srcComp,
              const int            a_dstComp,
              const int            a_numComp)
{
  CH_TIME("DataOps::plus(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::plus(*a_lhs[lvl], *a_rhs[lvl], a_srcComp, a_dstComp, a_numComp);
  }
}

void
DataOps::plus(LevelData<EBCellFAB>&       a_lhs,
              const LevelData<EBCellFAB>& a_rhs,
              const int                   a_srcComp,
              const int                   a_dstComp,
              const int                   a_numComp)
{
  CH_TIME("DataOps::plus(LD<EBCellFAB>)");

  for (DataIterator dit = a_lhs.disjointBoxLayout().dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       lhs = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    lhs.plus(rhs, a_srcComp, a_dstComp, a_numComp);
  }
}

void
DataOps::incr(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(LD<EBFluxData>)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    DataOps::incr(a_lhs[dit()], a_rhs[dit()], a_scale);
  }
}

void
DataOps::incr(EBFluxFAB& a_lhs, const EBFluxFAB& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(EBFluxData)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  EBFluxFAB rhsClone;
  rhsClone.clone(a_rhs);

  for (int dir = 0; dir < SpaceDim; dir++) {
    EBFaceFAB& lhs = a_lhs[dir];

    rhsClone[dir] *= a_scale;
    lhs += rhsClone[dir];
  }
}

void
DataOps::incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(LD<BaseIVFAB>)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int numComp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>&       lhs = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs = a_rhs[dit()];

    // Iteration space
    VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph());

    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < numComp; comp++) {
        lhs(vof, comp) += rhs(vof, comp) * a_scale;
      }
    };

    // Run kernel

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(EBAMRIFData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale)
{
  CH_TIME("DataOps::incr(LD<DomainFluxIFFAB>)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int numComp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    DomainFluxIFFAB&       lhs = a_lhs[dit()];
    const DomainFluxIFFAB& rhs = a_rhs[dit()];

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        BaseIFFAB<Real>&       curLHS = lhs(dir, sit());
        const BaseIFFAB<Real>& curRHS = rhs(dir, sit());

        const IntVectSet& ivs     = curLHS.getIVS();
        const EBGraph&    ebgraph = curLHS.getEBGraph();

        FaceIterator faceit(ivs, ebgraph, dir, FaceStop::SurroundingWithBoundary);

        auto kernel = [&](const FaceIndex& face) -> void {
          for (int comp = 0; comp < numComp; comp++) {
            curLHS(face, comp) += curRHS(face, comp) * a_scale;
          }
        };

        BoxLoops::loop(faceit, kernel);
      }
    }
  }
}

void
DataOps::incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(EBAMRCellData, EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(LD<EBCellFAB>, LD<BaseIVFAB>)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int numComp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&             lhs     = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs     = a_rhs[dit()];
    const EBGraph&         ebgraph = rhs.getEBGraph();
    const IntVectSet&      ivs     = rhs.getIVS();

    // Kernel space
    VoFIterator vofit(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < numComp; comp++) {
        lhs(vof, comp) += rhs(vof, comp) * a_scale;
      }
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(EBAMRIVData, EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::incr(*a_lhs[lvl], *a_rhs[lvl], a_scale);
  }
}

void
DataOps::incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale)
{
  CH_TIME("DataOps::incr(LD<BaseIVFAB>, LD<EBCellFAB>)");

  CH_assert(a_lhs.nComp() == a_rhs.nComp());

  const int numComp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>& lhs = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    // Kernel space
    VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph());

    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < numComp; comp++) {
        lhs(vof, comp) += rhs(vof, comp) * a_scale;
      }
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src)
{
  CH_TIME("DataOps::copy(MFAMRCellData)");

  for (int lvl = 0; lvl < a_dst.size(); lvl++) {
    CH_assert(a_dst[lvl]->nComp() == a_src[lvl]->nComp());

    if (a_src[lvl] != nullptr && a_dst[lvl] != nullptr) {
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void
DataOps::copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src)
{
  CH_TIME("DataOps::copy(EBAMRCellData)");

  for (int lvl = 0; lvl < a_dst.size(); lvl++) {
    CH_assert(a_dst[lvl]->nComp() == a_src[lvl]->nComp());

    if (a_src[lvl] != nullptr && a_dst[lvl] != nullptr) {
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void
DataOps::copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src)
{
  CH_TIME("DataOps::copy(EBAMRIVData)");

  for (int lvl = 0; lvl < a_dst.size(); lvl++) {
    CH_assert(a_dst[lvl]->nComp() == a_src[lvl]->nComp());

    if (a_src[lvl] != nullptr && a_dst[lvl] != nullptr) {
      a_src[lvl]->localCopyTo(*a_dst[lvl]);
    }
  }
}

void
DataOps::divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lhsComp, const int a_rhsComp)
{
  CH_TIME("DataOps::divide(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::divide(*a_lhs[lvl], *a_rhs[lvl], a_lhsComp, a_rhsComp);
  }
}

void
DataOps::divide(LevelData<EBCellFAB>&       a_lhs,
                const LevelData<EBCellFAB>& a_rhs,
                const int                   a_lhsComp,
                const int                   a_rhsComp)
{
  CH_TIME("DataOps::divide(EBAMRCellData)");

  CH_assert(a_lhs.nComp() > a_lhsComp);
  CH_assert(a_rhs.nComp() > a_rhsComp);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       lhs = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    lhs.divide(rhs, a_rhsComp, a_lhsComp, 1);
  }
}

void
DataOps::divideByScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs)
{
  CH_TIME("DataOps::divideByScalar(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::divideByScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::divideByScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("DataOps::divideByScalar(LD<EBCellFAB>)");

  const int numCompsLHS = a_lhs.nComp();

  CH_assert(a_rhs.nComp() == 1);
  CH_assert(a_lhs.nComp() >= 1);

  constexpr int rhsComp = 0;

  for (int lhsComp = 0; lhsComp < numCompsLHS; lhsComp++) {
    DataOps::divide(a_lhs, a_rhs, lhsComp, rhsComp);
  }
}

void
DataOps::divideFallback(EBAMRCellData& a_numerator, const EBAMRCellData& a_denominator, const EBAMRCellData& a_fallback)
{
  CH_TIME("DataOps::divideFallback(EBAMRCellData)");

  for (int lvl = 0; lvl < a_numerator.size(); lvl++) {
    DataOps::divideFallback(*a_numerator[lvl], *a_denominator[lvl], *a_fallback[lvl]);
  }
}

void
DataOps::divideFallback(LevelData<EBCellFAB>&       a_numerator,
                        const LevelData<EBCellFAB>& a_denominator,
                        const LevelData<EBCellFAB>& a_fallback)
{
  CH_TIME("DataOps::LD<EBCellFAB>)");

  CH_assert(a_numerator.nComp() == a_denominator.nComp());
  CH_assert(a_numerator.nComp() == a_fallback.nComp());

  for (DataIterator dit = a_numerator.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       numerator   = a_numerator[dit()];
    const EBCellFAB& denominator = a_denominator[dit()];
    const EBCellFAB& fallback    = a_fallback[dit()];

    // Hooks to single-valued data.
    BaseFab<Real>&       regNumerator   = numerator.getSingleValuedFAB();
    const BaseFab<Real>& regDenominator = denominator.getSingleValuedFAB();
    const BaseFab<Real>& regFallback    = fallback.getSingleValuedFAB();

    // Kernel regions
    const Box        cellBox  = a_numerator.disjointBoxLayout()[dit()];
    const EBISBox&   ebisBox  = numerator.getEBISBox();
    const EBGraph&   ebGraph  = ebisBox.getEBGraph();
    const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);
    VoFIterator      vofit(irregIVS, ebGraph);

    // I need to clone the input data because the regular kernel will screw with it, breaking our cut-cell kernels.
    EBCellFAB cloneNumerator;
    cloneNumerator.clone(numerator);

    // Regular cells
    for (int comp = 0; comp < a_numerator.nComp(); comp++) {

      // Regular kernel.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (std::abs(regDenominator(iv, comp)) > 0.0) {
          regNumerator(iv, comp) /= regDenominator(iv, comp);
        }
        else {
          regNumerator(iv, comp) = regFallback(iv, comp);
        }
      };

      // Irregular kernel. Same as the above but recall that we may have multi-valued cells.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const Real& denom = denominator(vof, comp);

        if (std::abs(denom) > 0.0) {
          numerator(vof, comp) = cloneNumerator(vof, comp) / denom;
        }
        else {
          numerator(vof, comp) = fallback(vof, comp);
        }
      };

      // Run the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::divideFallback(EBAMRCellData& a_numerator, const EBAMRCellData& a_denominator, const Real a_fallback)
{
  CH_TIME("DataOps::divideFallback(EBAMRCellData)");

  for (int lvl = 0; lvl < a_numerator.size(); lvl++) {
    DataOps::divideFallback(*a_numerator[lvl], *a_denominator[lvl], a_fallback);
  }
}

void
DataOps::divideFallback(LevelData<EBCellFAB>&       a_numerator,
                        const LevelData<EBCellFAB>& a_denominator,
                        const Real                  a_fallback)
{
  CH_TIME("DataOps::divdeFallback");

  CH_assert(a_numerator.nComp() == a_denominator.nComp());

  for (DataIterator dit = a_numerator.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       numerator   = a_numerator[dit()];
    const EBCellFAB& denominator = a_denominator[dit()];

    // Hook to single-valued data.
    BaseFab<Real>&       regNumerator   = numerator.getSingleValuedFAB();
    const BaseFab<Real>& regDenominator = denominator.getSingleValuedFAB();

    // Kernel regions
    const Box        cellBox  = a_numerator.disjointBoxLayout()[dit()];
    const EBISBox&   ebisBox  = numerator.getEBISBox();
    const EBGraph&   ebGraph  = ebisBox.getEBGraph();
    const IntVectSet irregIVS = ebisBox.getIrregIVS(cellBox);
    VoFIterator      vofit(irregIVS, ebGraph);

    // I need to clone the input data because the regular kernel will screw with it.
    EBCellFAB cloneNumerator;
    cloneNumerator.clone(numerator);

    // Regular cells
    for (int comp = 0; comp < a_numerator.nComp(); comp++) {

      // Regular kernel.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (std::abs(regDenominator(iv, comp)) > 0.0) {
          regNumerator(iv, comp) /= regDenominator(iv, comp);
        }
        else {
          regNumerator(iv, comp) = a_fallback;
        }
      };

      // Irregular kernel. Same as the above but recall that we may have multi-valued cells.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const Real& denom = denominator(vof, comp);

        if (std::abs(denom) > 0.0) {
          numerator(vof, comp) = cloneNumerator(vof, comp) / denom;
        }
        else {
          numerator(vof, comp) = a_fallback;
        }
      };

      // Execute the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::floor(EBAMRCellData& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::floor(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::floor(*a_lhs[lvl], a_value);
  }
}

void
DataOps::floor(LevelData<EBCellFAB>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::floor(LD<EBCelLFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&     lhs     = a_lhs[dit()];
    const Box      box     = lhs.getRegion(); // Note:: All cells are floored (ghosts also)
    const EBISBox& ebisbox = lhs.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const int      numComp = a_lhs.nComp();

    // Hook to single-valued data.
    BaseFab<Real>& lhs_reg = lhs.getSingleValuedFAB();

    // Irregular kernel region
    VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph);

    // Components loop -- we do all components.
    for (int comp = 0; comp < numComp; comp++) {

      // Regular kernel.
      auto regularKernel = [&](const IntVect& iv) -> void { lhs_reg(iv, comp) = std::max(a_value, lhs_reg(iv, comp)); };

      // Irregular kernel
      auto irregularKernel = [&](const VolIndex& vof) -> void { lhs(vof, comp) = std::max(a_value, lhs(vof, comp)); };

      // Execute the kernels.
      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::floor(EBAMRIVData& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::floor(EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::floor(*a_lhs[lvl], a_value);
  }
}

void
DataOps::floor(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::floor(LD<BaseIVFAB<Real> >)");

  const int numComp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>& lhs = a_lhs[dit()];

    // Irregular kernel region.
    VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph());

    // Do all components.
    for (int comp = 0; comp < numComp; comp++) {

      // Irregular kernel.
      auto irregularKernel = [&](const VolIndex& vof) -> void { lhs(vof, comp) = std::max(a_value, lhs(vof, comp)); };

      // Run kernel
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::max(EBAMRCellData& a_data, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2)
{
  CH_TIME("DataOps::max(EBAMRCellData x3)");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    DataOps::max(*a_data[lvl], *a_data1[lvl], *a_data2[lvl]);
  }
}

void
DataOps::max(LevelData<EBCellFAB>& a_data, const LevelData<EBCellFAB>& a_data1, const LevelData<EBCellFAB>& a_data2)
{
  CH_TIME("DataOps::max(LD<EBCellFAB> x3)");

  const int numComp = a_data.nComp();

  CH_assert(a_data1.nComp() == numComp);
  CH_assert(a_data2.nComp() == numComp);

  const DisjointBoxLayout& dbl = a_data.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {

    EBCellFAB&       data  = a_data[dit()];
    const EBCellFAB& data1 = a_data1[dit()];
    const EBCellFAB& data2 = a_data2[dit()];

    FArrayBox&       dataReg  = data.getFArrayBox();
    const FArrayBox& data1Reg = data1.getFArrayBox();
    const FArrayBox& data2Reg = data2.getFArrayBox();

    auto regularKernel = [&](const IntVect& iv) -> void {
      for (int comp = 0; comp < numComp; comp++) {
        dataReg(iv, comp) = std::max(data1Reg(iv, comp), data2Reg(iv, comp));
      }
    };

    auto irregularKernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < numComp; comp++) {
        data(vof, comp) = std::max(data1(vof, comp), data2(vof, comp));
      }
    };

    // Regions for the kernels.
    const Box&        cellBox = dbl[dit()];
    const EBISBox&    ebisbox = data.getEBISBox();
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs     = ebisbox.getIrregIVS(cellBox);
    VoFIterator       vofit(ivs, ebgraph);

    // Run kernels
    BoxLoops::loop(cellBox, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

void
DataOps::getMaxMin(Real& a_max, Real& a_min, EBAMRCellData& a_data, const int a_comp)
{
  CH_TIME("DataOps::getMaxMin(EBAMRCellData)");

  a_max = -std::numeric_limits<Real>::max();
  a_min = std::numeric_limits<Real>::max();

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    Real lvlMax = -std::numeric_limits<Real>::max();
    Real lvlMin = std::numeric_limits<Real>::max();

    DataOps::getMaxMin(lvlMax, lvlMin, *a_data[lvl], a_comp);

    a_max = std::max(a_max, lvlMax);
    a_min = std::min(a_min, lvlMin);
  }
}

void
DataOps::getMaxMin(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_data, const int a_comp)
{
  CH_TIME("DataOps::getMaxMin(LD<EBCellFAB>)");

  EBLevelDataOps::getMaxMin(a_max, a_min, a_data, a_comp);
}

void
DataOps::getMaxMin(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data)
{
  CH_TIME("DataOps::getMaxMin(Vector<EBAMRCellData>)");

  a_max.resize(a_data.size(), -std::numeric_limits<Real>::max());
  a_min.resize(a_data.size(), std::numeric_limits<Real>::max());

  constexpr int comp    = 0;
  constexpr int numComp = 1;

  for (int i = 0; i < a_data.size(); i++) {
    CH_assert(a_data[i][0]->nComp() == numComp);

    DataOps::getMaxMin(a_max[i], a_min[i], a_data[i], comp);
  }
}

void
DataOps::getMaxMinNorm(Real& a_max, Real& a_min, EBAMRCellData& a_data)
{
  CH_TIME("DataOps::getMaxMinNorm(EBAMRCellData)");

  a_max = -std::numeric_limits<Real>::max();
  a_min = std::numeric_limits<Real>::max();

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    Real lvlMax;
    Real lvlMin;

    DataOps::getMaxMinNorm(lvlMax, lvlMin, *a_data[lvl]);

    a_max = std::max(a_max, lvlMax);
    a_min = std::min(a_min, lvlMin);
  }
}

void
DataOps::getMaxMinNorm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& a_data)
{
  CH_TIME("DataOps::LD<EBCelLFAB>");

  a_max = -std::numeric_limits<Real>::max();
  a_min = std::numeric_limits<Real>::max();

  const int maskComp = 0;
  const int numComp  = a_data.nComp();

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit) {
    const Box&       box     = a_data.disjointBoxLayout().get(dit());
    const EBCellFAB& data    = a_data[dit()];
    const EBISBox&   ebisbox = data.getEBISBox();
    const EBGraph&   ebgraph = ebisbox.getEBGraph();

    // Use a mask for ignoring covered data (which can be bogus).
    EBCellFAB coveredMask(ebisbox, box, 1);
    coveredMask.setVal(1.0);
    coveredMask.setCoveredCellVal(-1.0, 0);

    // Hooks to single-valued data
    const BaseFab<Real>& mask    = coveredMask.getSingleValuedFAB();
    const BaseFab<Real>& dataReg = data.getSingleValuedFAB();

    // Irregular kernel region.
    VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph);

    // Regular kernel. This computes the max/min values from sqrt(x1*x1 + x2*x2 + x3*x3 + ...).
    auto regularKernel = [&](const IntVect& iv) -> void {
      if (mask(iv, maskComp) > 0.0) {
        Real curValue = 0.0;
        for (int comp = 0; comp < numComp; comp++) {
          curValue += dataReg(iv, comp) * dataReg(iv, comp);
        }
        curValue = sqrt(curValue);

        a_max = std::max(a_max, curValue);
        a_min = std::min(a_min, curValue);
      }
    };

    // Irregular kernel. This computes the max/min values from sqrt(x1*x1 + x2*x2 + x3*x3 + ...).
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      Real cur = 0.0;
      for (int comp = 0; comp < numComp; comp++) {
        cur += data(vof, comp) * data(vof, comp);
      }
      cur = sqrt(cur);

      a_max = std::max(a_max, cur);
      a_min = std::min(a_min, cur);
    };

    // Run the kernels.
    BoxLoops::loop(box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }

  // If running with MPI then we need to reduce the result.
  a_max = ParallelOps::max(a_max);
  a_min = ParallelOps::min(a_min);
}

void
DataOps::getMaxMinNorm(Real& a_max, Real& a_min, EBAMRIVData& a_data)
{
  CH_TIME("DataOps::getMaxMinNorm(EBAMRIVData)");

  a_max = -std::numeric_limits<Real>::max();
  a_min = std::numeric_limits<Real>::min();

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    Real lvlMax;
    Real lvlMin;

    DataOps::getMaxMinNorm(lvlMax, lvlMin, *a_data[lvl]);

    a_max = std::max(a_max, lvlMax);
    a_min = std::min(a_min, lvlMin);
  }
}

void
DataOps::getMaxMinNorm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real>>& a_data)
{
  CH_TIME("DataOps::getMaxMinNorm(LD<BaseIVFAB>)");

  a_max = -std::numeric_limits<Real>::max();
  a_min = std::numeric_limits<Real>::max();

  const int numComp = a_data.nComp();

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit) {
    const BaseIVFAB<Real>& data = a_data[dit()];

    // Iteration space
    VoFIterator vofit(data.getIVS(), data.getEBGraph());

    auto kernel = [&](const VolIndex& vof) -> void {
      Real cur = 0.0;
      for (int comp = 0; comp < numComp; comp++) {
        cur += data(vof, comp) * data(vof, comp);
      }
      cur = sqrt(cur);

      a_max = std::max(a_max, cur);
      a_min = std::min(a_min, cur);
    };

    BoxLoops::loop(vofit, kernel);
  }

  // If running with MPI then we need to reduce the result.
  a_max = ParallelOps::max(a_max);
  a_min = ParallelOps::min(a_min);
}

void
DataOps::invert(EBAMRFluxData& a_data)
{
  CH_TIME("DataOps::invert(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    DataOps::invert(*a_data[lvl]);
  }
}

void
DataOps::invert(LevelData<EBFluxFAB>& a_data)
{
  CH_TIME("DataOps::invert(LD<EBFluxFAB>)");

  const int numComp = a_data.nComp();

  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit) {
    EBFluxFAB&       data    = a_data[dit()];
    const EBISBox&   ebisbox = data.getEBISBox();
    const EBGraph&   ebgraph = ebisbox.getEBGraph();
    const Box        box     = a_data.disjointBoxLayout().get(dit());
    const IntVectSet irreg   = ebisbox.getIrregIVS(box);

    // Get faces oriented in direction dir
    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB& data = a_data[dit()][dir];

      // Kernel regions.
      const Box    facebox = surroundingNodes(box, dir);
      FaceIterator faceit(irreg, ebgraph, dir, FaceStop::SurroundingWithBoundary);

      // Need a copy because regular kernel inverts face data also.
      EBFaceFAB cpy(ebisbox, box, dir, numComp);
      cpy.setVal(0.0);
      cpy += data;

      // Hook to single-valued data.
      BaseFab<Real>& dataReg = data.getSingleValuedFAB();

      // Component loop -- we do all components.
      for (int comp = 0; comp < numComp; comp++) {

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void { dataReg(iv, comp) = 1. / dataReg(iv, comp); };

        // Irregular kernel.
        auto irregularKernel = [&](const FaceIndex& face) -> void { data(face, comp) = 1. / data(face, comp); };

        // Run the kernels.
        BoxLoops::loop(facebox, regularKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }
    }
  }
}

void
DataOps::kappaSum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs)
{
  CH_TIME("DataOps::kappaSum");

  CH_assert(a_lhs.nComp() == 1);

  Real mass = 0.;

  constexpr int comp = 0;

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    const Box        box     = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& lhs     = a_lhs[dit()];
    const EBISBox&   ebisbox = lhs.getEBISBox();

    // This is crazy -- we should have a Fortran loop here. This could become a bottleneck for people.
    const EBGraph&   ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    VoFIterator vofit(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void { mass += ebisbox.volFrac(vof) * lhs(vof, comp); };

    BoxLoops::loop(vofit, kernel);
  }

  a_mass = EBLevelDataOps::parallelSum(mass);
}

void
DataOps::kappaScale(EBAMRCellData& a_data)
{
  CH_TIME("DataOps::kappaScale(EBAMRCellData");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    DataOps::kappaScale(*a_data[lvl]);
  }
}

void
DataOps::kappaScale(LevelData<EBCellFAB>& a_data)
{
  CH_TIME("DataOps::kappaScale(LD<EBCellFAB>");

  EBLevelDataOps::kappaWeight(a_data);
}

void
DataOps::kappaScale(MFAMRCellData& a_data)
{
  CH_TIME("DataOps::kappaScale(MFAMRCellData");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    DataOps::kappaScale(*a_data[lvl]);
  }
}

void
DataOps::kappaScale(LevelData<MFCellFAB>& a_data)
{
  CH_TIME("DataOps::kappaScale(LD<MFCellFAB>");

  MFLevelDataOps::kappaWeight(a_data);
}

void
DataOps::volumeScale(EBAMRCellData& a_data, const Vector<Real>& a_dx)
{
  CH_TIME("DataOps::volumeScale(EBAMRCellData, Vector<Real>");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    DataOps::scale(*a_data[lvl], std::pow(a_dx[lvl], SpaceDim));
  }
}

void
DataOps::multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs)
{
  CH_TIME("DataOps::multiply(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("DataOps::multiply(LD<EBCellFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void
DataOps::multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs)
{
  CH_TIME("DataOps::multiply(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs)
{
  CH_TIME("DataOps::multiply(LD<EBFluxFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()] *= a_rhs[dit()];
  }
}

void
DataOps::multiply(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs)
{
  CH_TIME("DataOps::multiply(EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::multiply(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::multiply(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs)
{
  CH_TIME("DataOps::multiply(LD<BaseIVFAB>)");

  const int ncomp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>&       lhs     = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs     = a_rhs[dit()];
    const EBGraph&         ebgraph = lhs.getEBGraph();
    const IntVectSet&      ivs     = lhs.getIVS() & rhs.getIVS();

    // Kernel space
    VoFIterator vofit(ivs, ebgraph);

    // Kernel
    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < ncomp; comp++) {
        lhs(vof, comp) *= rhs(vof, comp);
      }
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::multiplyScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs)
{
  CH_TIME("DataOps::multiplyScalar(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::multiplyScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::multiplyScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("DataOps::multiplyScalar(LD<EBCellFAB>)");

  CH_assert(a_rhs.nComp() == 1);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       lhs = a_lhs[dit()];
    const EBCellFAB& rhs = a_rhs[dit()];

    for (int comp = 0; comp < lhs.nComp(); comp++) {
      lhs.mult(rhs, 0, comp, 1);
    }
  }
}

void
DataOps::multiplyScalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs)
{
  CH_TIME("DataOps::multiplyScalar(EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::multiplyScalar(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::multiplyScalar(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs)
{
  CH_TIME("DataOps::multiplyScalar(LD<BaseIVFAB>)");

  CH_assert(a_rhs.nComp() == 1);

  const int ncomp = a_lhs.nComp();

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>&       lhs     = a_lhs[dit()];
    const BaseIVFAB<Real>& rhs     = a_rhs[dit()];
    const EBGraph&         ebgraph = lhs.getEBGraph();
    const IntVectSet&      ivs     = lhs.getIVS();

    VoFIterator vofit(ivs, ebgraph);

    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < ncomp; comp++) {
        lhs(vof, comp) *= rhs(vof, 0);
      }
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p)
{
  CH_TIME("DataOps::norm");

  Real volume;

  a_norm = EBLevelDataOps::kappaNorm(volume, a_data, EBLEVELDATAOPS_ALLVOFS, a_domain, a_p);
}

Real
DataOps::norm(const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p, const bool a_regOnly)
{
  CH_TIME("DataOps::norm");

  Real volume;
  Real a_norm;

  if (a_regOnly) {
    a_norm = EBLevelDataOps::kappaNorm(volume, a_data, EBLEVELDATAOPS_REGULARVOFS, a_domain, a_p);
  }
  else {
    a_norm = EBLevelDataOps::kappaNorm(volume, a_data, EBLEVELDATAOPS_ALLVOFS, a_domain, a_p);
  }

  return a_norm;
}

void
DataOps::scale(MFAMRCellData& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(MFAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void
DataOps::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(LD<MFCellFAB>)");

  MFLevelDataOps::scale(a_lhs, a_scale);
}

void
DataOps::scale(MFAMRFluxData& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(MFAMRFluxData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void
DataOps::scale(LevelData<MFFluxFAB>& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(LD<MFFluxFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    for (int iphase = 0; iphase < a_lhs[dit()].numPhases(); iphase++) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB& lhs = a_lhs[dit()].getPhase(iphase)[dir];

        lhs *= a_scale;
      }
    }
  }
}

void
DataOps::scale(EBAMRIVData& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(EBAMRIVData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void
DataOps::scale(EBAMRCellData& a_lhs, const Real a_scale)
{
  CH_TIME("DataOps::scale(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void
DataOps::scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale)
{
  CH_TIME("DataOps::scale(LD<EBCellFAB>)");

  EBLevelDataOps::scale(a_lhs, a_scale);
}

void
DataOps::scale(EBAMRFluxData& a_lhs, const Real a_scale)
{
  CH_TIME("DataOps::scale(LD<EBFluxFAB>)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    EBLevelDataOps::scale(*a_lhs[lvl], a_scale);
  }
}

void
DataOps::scale(LevelData<BaseIVFAB<Real>>& a_lhs, const Real& a_scale)
{
  CH_TIME("DataOps::scale(LD<BaseIVFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    BaseIVFAB<Real>& lhs = a_lhs[dit()];

    VoFIterator vofit(lhs.getIVS(), lhs.getEBGraph());

    auto kernel = [&](const VolIndex& vof) -> void {
      for (int comp = 0; comp < a_lhs.nComp(); comp++) {
        lhs(vofit(), comp) *= a_scale;
      }
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::setCoveredValue(EBAMRCellData& a_lhs, const int a_comp, const Real a_value)
{
  CH_TIME("DataOps::setCoveredValue(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setCoveredValue(*a_lhs[lvl], a_comp, a_value);
  }
}

void
DataOps::setCoveredValue(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value)
{
  CH_TIME("DataOps::setCoveredValue(LD<EBCellFAB>)");

  EBLevelDataOps::setCoveredVal(a_lhs, a_comp, a_value);
}

void
DataOps::setCoveredValue(EBAMRCellData& a_lhs, const Real a_value)
{
  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setCoveredValue(*a_lhs[lvl], a_value);
  }
}

void
DataOps::setCoveredValue(LevelData<EBCellFAB>& a_lhs, const Real a_value)
{
  for (int comp = 0; comp < a_lhs.nComp(); comp++) {
    EBLevelDataOps::setCoveredVal(a_lhs, comp, a_value);
  }
}

void
DataOps::setInvalidValue(EBAMRCellData& a_lhs, const Vector<int>& a_refRat, const Real a_value)
{
  CH_TIME("DataOps::setInvalidValud(EBAMRCellData, Vector<int>, Real)");

  const int finestLevel = a_lhs.size() - 1;

  for (int lvl = finestLevel; lvl > 0; lvl--) {
    LevelData<EBCellFAB>&       coarLD = *a_lhs[lvl - 1];
    const LevelData<EBCellFAB>& fineLD = *a_lhs[lvl];

    const DisjointBoxLayout& dblCoar = coarLD.disjointBoxLayout();
    const DisjointBoxLayout& dblFine = fineLD.disjointBoxLayout();

    const int nComp = coarLD.nComp();

    for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
      EBCellFAB& coarData = coarLD[dit()];
      const Box  cellBox  = dblCoar[dit()];

      const EBISBox& ebisBox = coarData.getEBISBox();
      const EBGraph& ebGraph = ebisBox.getEBGraph();

      for (LayoutIterator lit = dblFine.layoutIterator(); lit.ok(); ++lit) {
        const Box fineBox = dblFine[lit()];
        const Box coFiBox = coarsen(fineBox, a_refRat[lvl - 1]);

        const Box overlapBox = cellBox & coFiBox;

        if (!overlapBox.isEmpty()) {

          // Irregular kernel region.
          VoFIterator vofit(IntVectSet(overlapBox), ebGraph);

          // Regular and covered cells' kernel.
          BaseFab<Real>& coarFAB = coarData.getSingleValuedFAB();

          // Do all components.
          for (int comp = 0; comp < nComp; comp++) {

            // Kernel called for all cells.
            auto regularKernel = [&](const IntVect& iv) -> void { coarFAB(iv, comp) = a_value; };

            // Kernel called for irregular cells.
            auto irregularKernel = [&](const VolIndex& vof) -> void {
              for (int comp = 0; comp < nComp; comp++) {
                coarData(vof, comp) = a_value;
              }
            };

            // Execute kernels.
            BoxLoops::loop(overlapBox, regularKernel);
            BoxLoops::loop(vofit, irregularKernel);
          }
        }
      }
    }
  }
}

void
DataOps::setValue(MFAMRCellData&                             a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Vector<Real>&                        a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(MFAMRCellData, std::function)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void
DataOps::setValue(LevelData<MFCellFAB>&                      a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Real                                 a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(LD<MFCellFAB>, std::function)");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {

    MFCellFAB& lhs = a_lhs[dit()];

    for (int i = 0; i < lhs.numPhases(); i++) {
      EBCellFAB&     phaseData    = lhs.getPhase(i);
      BaseFab<Real>& phaseDataFAB = phaseData.getSingleValuedFAB();

      // Kernel regions
      const Box         box     = phaseData.box();
      const EBISBox&    ebisbox = phaseData.getEBISBox();
      const EBGraph&    ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg   = ebisbox.getIrregIVS(box);
      VoFIterator       vofit(irreg, ebgraph);

      // Regular cells
      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

        phaseDataFAB(iv, a_comp) = a_function(pos);
      };

      // Cut-cells.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect& iv = vof.gridIndex();

        const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

        phaseData(vof, a_comp) = a_function(pos);
      };

      // Run kernels.
      BoxLoops::loop(box, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
DataOps::setValue(EBAMRCellData&                             a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Vector<Real>&                        a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(EBAMRCellData, std::function)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void
DataOps::setValue(LevelData<EBCellFAB>&                      a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Real                                 a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(LD<EBCellFAB>, std::function)");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    EBCellFAB&     lhs    = a_lhs[dit()];
    BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();

    // Kernel regions.
    const Box         box     = lhs.box();
    const EBISBox&    ebisbox = lhs.getEBISBox();
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& irreg   = ebisbox.getIrregIVS(box);
    VoFIterator       vofit(irreg, ebgraph);

    // Regular cells
    auto regularKernel = [&](const IntVect& iv) -> void {
      const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

      lhsFAB(iv, a_comp) = a_function(pos);
    };

    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect& iv = vof.gridIndex();

      const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx; // + ebisbox.centroid(vof)*a_dx;

      lhs(vof, a_comp) = a_function(pos);
    };

    // Run kernels.
    BoxLoops::loop(box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

void
DataOps::setValue(EBAMRFluxData&                             a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Vector<Real>&                        a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(EBAMRFluxData, std::function)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void
DataOps::setValue(LevelData<EBFluxFAB>&                      a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Real                                 a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(LD<EBFluxFAB>, std::function)");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB&     lhs    = a_lhs[dit()][dir];
      BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();

      // Kernel regions.
      const Box         box     = lhs.getCellRegion();
      const EBISBox&    ebisbox = lhs.getEBISBox();
      const EBGraph&    ebgraph = ebisbox.getEBGraph();
      const IntVectSet& irreg   = ebisbox.getIrregIVS(box);
      const Box         facebox = surroundingNodes(box, dir);
      FaceIterator      faceit(irreg, ebgraph, dir, FaceStop::SurroundingWithBoundary);

      // Regular cells
      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = a_probLo + RealVect(iv) * a_dx;

        lhsFAB(iv, a_comp) = a_function(pos);
      };

      auto irregularKernel = [&](const FaceIndex& face) -> void {
        const RealVect pos = EBArith::getFaceLocation(face, a_dx, a_probLo);

        lhs(face, a_comp) = a_function(pos);
      };

      // Run the kernels.
      BoxLoops::loop(facebox, regularKernel);
      BoxLoops::loop(faceit, irregularKernel);
    }
  }
}

void
DataOps::setValue(EBAMRIVData&                               a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Vector<Real>&                        a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(EBAMRIVData, std::function)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl], a_comp);
  }
}

void
DataOps::setValue(LevelData<BaseIVFAB<Real>>&                a_lhs,
                  const std::function<Real(const RealVect)>& a_function,
                  const RealVect                             a_probLo,
                  const Real                                 a_dx,
                  const int                                  a_comp)
{
  CH_TIME("DataOps::setValue(LD<BaseIVFAB>, std::function)");

  // As we don't specify where the function should be evaluated, this routine sets a_lhs to be evaluated at the cell center.
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    BaseIVFAB<Real>& lhs = a_lhs[dit()];

    // Kernel region.
    const EBGraph&    ebgraph = lhs.getEBGraph();
    const IntVectSet& irreg   = lhs.getIVS();
    VoFIterator       vofit(irreg, ebgraph);

    // Irregular kernel.
    auto kernel = [&](const VolIndex& vof) {
      const IntVect  iv  = vof.gridIndex();
      const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;

      lhs(vof, a_comp) = a_function(pos);
    };

    BoxLoops::loop(vofit, kernel);
  }
}

void
DataOps::setValue(EBAMRCellData&                                 a_lhs,
                  const std::function<RealVect(const RealVect)>& a_function,
                  const RealVect                                 a_probLo,
                  const Vector<Real>&                            a_dx)
{
  CH_TIME("DataOps::setValue(EBAMRCellData, std::function<RealVect>)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_function, a_probLo, a_dx[lvl]);
  }
}

void
DataOps::setValue(LevelData<EBCellFAB>&                          a_lhs,
                  const std::function<RealVect(const RealVect)>& a_function,
                  const RealVect                                 a_probLo,
                  const Real                                     a_dx)
{
  CH_TIME("DataOps::setValue(LD<EBCellFAB>, std::function<RealVect>)");

  CH_assert(a_lhs.nComp() == SpaceDim);

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    EBCellFAB&     lhs    = a_lhs[dit()];
    BaseFab<Real>& lhsFAB = lhs.getSingleValuedFAB();

    // Kernel regions.
    const Box         box     = lhs.box();
    const EBISBox&    ebisbox = lhs.getEBISBox();
    const EBGraph&    ebgraph = ebisbox.getEBGraph();
    const IntVectSet& irreg   = ebisbox.getIrregIVS(box);
    VoFIterator       vofit(irreg, ebgraph);

    // Regular kernel.
    auto regularKernel = [&](const IntVect& iv) -> void {
      const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;
      const RealVect val = a_function(pos);

      for (int comp = 0; comp < SpaceDim; comp++) {
        lhsFAB(iv, comp) = val[comp];
      }
    };

    // Irregular cells
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect& iv = vof.gridIndex();

      const RealVect pos = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx; // + ebisbox.centroid(vof)*a_dx;
      const RealVect val = a_function(pos);

      for (int comp = 0; comp < SpaceDim; comp++) {
        lhs(vof, comp) = val[comp];
      }
    };

    // Run the kernels.
    BoxLoops::loop(box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

void
DataOps::setValue(EBAMRCellData& a_data, const Real& a_value)
{
  CH_TIME("DataOps::setValue(EBAMRCellData, Real)");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void
DataOps::setValue(EBAMRCellData& a_lhs, const Real a_value, const int a_comp)
{
  CH_TIME("DataOps::setValue(EBAMRCellData, Real, int)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_value, a_comp);
  }
}

void
DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp)
{
  CH_TIME("DataOps::setValue(LD<EBCellFAB>, Real, int)");

  EBLevelDataOps::setVal(a_lhs, a_value, a_comp);
}

void
DataOps::setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::setValue(LD<EBCellFAB>, Real)");

  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
DataOps::setValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::setValue(LD<EBFluxFAB>, Real)");

  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
DataOps::setValue(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::setValue(LD<BaseIVFAB>, Real)");

  EBLevelDataOps::setVal(a_lhs, a_value);
}

void
DataOps::setValue(EBAMRFluxData& a_data, const Real& a_value)
{
  CH_TIME("DataOps::setValue(EBAMRFluxData, Real)");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void
DataOps::setValue(EBAMRIVData& a_data, const Real& a_value)
{
  CH_TIME("DataOps::setValue(EBAMRIVData, Real)");

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    EBLevelDataOps::setVal(*a_data[lvl], a_value);
  }
}

void
DataOps::setValue(MFAMRCellData& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(MFAMRCellData, Real)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void
DataOps::setValue(LevelData<MFCellFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(LD<MFCellFAB>, Real)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    MFCellFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void
DataOps::setValue(MFAMRFluxData& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(MFAMRFluxFAB, Real)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void
DataOps::setValue(LevelData<MFFluxFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(LD<MFFluxFAB>, Real)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    MFFluxFAB& lhs = a_lhs[dit()];
    lhs.setVal(a_value);
  }
}

void
DataOps::setValue(MFAMRIVData& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(MFAMRIVData, Real)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void
DataOps::setValue(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value)
{
  CH_TIME("DataOps::setValue(LD<MFBaseIVFAB>, Real)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].setVal(a_value);
  }
}

void
DataOps::setValue(EBAMRIFData& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::setValue(EBAMRIFData, Real)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::setValue(*a_lhs[lvl], a_value);
  }
}

void
DataOps::setValue(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value)
{
  CH_TIME("DataOps::setValue(LD<DomainFluxIFFAB>, Real)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    DomainFluxIFFAB& lhs = a_lhs[dit()];
    for (SideIterator sit; sit.ok(); ++sit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        BaseIFFAB<Real>& fab = lhs(dir, sit());
        fab.setVal(a_value);
      }
    }
  }
}

void
DataOps::sum(Real& a_value)
{
  CH_TIME("DataOps::sum");

#ifdef CH_MPI
  Real cur    = a_value;
  Real tmp    = a_value;
  int  result = MPI_Allreduce(&cur, &tmp, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("DataOps::sum - communication error on Allreduce");
  }

  a_value = tmp;
#endif
}

void
DataOps::squareRoot(EBAMRFluxData& a_lhs)
{
  CH_TIME("DataOps::squareRoot(EBAMRFluxData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::squareRoot(*a_lhs[lvl]);
  }
}

void
DataOps::squareRoot(LevelData<EBFluxFAB>& a_lhs)
{
  CH_TIME("DataOps::squareRoot(LD<EBFluxFAB>)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    const Box& box = a_lhs.disjointBoxLayout().get(dit());

    for (int dir = 0; dir < SpaceDim; dir++) {
      EBFaceFAB& lhs = a_lhs[dit()][dir];

      // Kernel regions.
      const Box         facebox = surroundingNodes(box, dir);
      const EBISBox&    ebisbox = lhs.getEBISBox();
      const EBGraph&    ebgraph = ebisbox.getEBGraph();
      const IntVectSet& ivs     = ebisbox.getIrregIVS(box);
      FaceIterator      faceit(ivs, ebgraph, dir, FaceStop::SurroundingWithBoundary);

      // Need a copy because regular kernel will invert cells also.
      EBFaceFAB cpy(ebisbox, box, dir, lhs.nComp());
      cpy.setVal(0.0);
      cpy += lhs;

      // Hook to single-valued data.
      BaseFab<Real>& lhs_reg = lhs.getSingleValuedFAB();

      // All comps
      for (int comp = 0; comp < lhs.nComp(); comp++) {

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void { lhs_reg(iv, comp) = sqrt(lhs_reg(iv, comp)); };

        // Irregular kernel. Reaches into the cpy data holder because the regular kernel will have messed
        // with the data.
        auto irregularKernel = [&](const FaceIndex& face) -> void { lhs(face, comp) = sqrt(cpy(face, comp)); };

        // Execute the kernels.
        BoxLoops::loop(facebox, regularKernel);
        BoxLoops::loop(faceit, irregularKernel);
      }
    }
  }
}

void
DataOps::vectorLength(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs)
{
  CH_TIME("DataOps::vectorLength(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::vectorLength(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::vectorLength(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("DataOps::vectorLength(LD<EBCellFAB>");

  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       lhs = a_lhs[dit()];
    const Box&       box = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs = a_rhs[dit()];

    DataOps::vectorLength(lhs, rhs, box);
  }
}

void
DataOps::vectorLength(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box)
{
  CH_TIME("DataOps::vectorLength(EBCellFAB)");

  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  // Component in a_lhs that put the data into.
  constexpr int comp = 0;

  // Mask for skipping computation on covered cells
  const EBISBox& ebisbox = a_lhs.getEBISBox();
  EBCellFAB      coveredMask(ebisbox, a_box, 1);
  coveredMask.setVal(1.0);
  coveredMask.setCoveredCellVal(-1.0, 0);

  // Kernel regions
  VoFIterator vofit(ebisbox.getIrregIVS(a_box), ebisbox.getEBGraph());

  // Hooks to single-valued data.
  BaseFab<Real>&       lhsReg = a_lhs.getSingleValuedFAB();
  const BaseFab<Real>& rhsReg = a_rhs.getSingleValuedFAB();
  const BaseFab<Real>& mask   = coveredMask.getSingleValuedFAB();

  auto regularKernel = [&](const IntVect& iv) -> void {
    lhsReg(iv, comp) = 0.0;

    if (mask(iv, comp) > 0.0) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        lhsReg(iv, comp) += rhsReg(iv, dir) * rhsReg(iv, dir);
      }

      lhsReg(iv, comp) = sqrt(lhsReg(iv, comp));
    }
  };

  // Irregular kernel. Same as the above.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_lhs(vof, comp) = 0.;

    for (int i = 0; i < SpaceDim; i++) {
      a_lhs(vof, comp) += a_rhs(vof, i) * a_rhs(vof, i);
    }

    a_lhs(vof, comp) = sqrt(a_lhs(vof, comp));
  };

  // Run the kernels
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

void
DataOps::vectorLength2(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs)
{
  CH_TIME("DataOps::vectorLength2(EBAMRCellData)");

  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    DataOps::vectorLength2(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

void
DataOps::vectorLength2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs)
{
  CH_TIME("DataOps::vectorLength2(LD<EBCellFAB>)");

  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    EBCellFAB&       lhs = a_lhs[dit()];
    const Box&       box = a_lhs.disjointBoxLayout().get(dit());
    const EBCellFAB& rhs = a_rhs[dit()];

    DataOps::vectorLength2(lhs, rhs, box);
  }
}

void
DataOps::vectorLength2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& a_box)
{
  CH_TIME("DataOps::vectorLength2(EBCellFAB)");

  CH_assert(a_lhs.nComp() == 1);
  CH_assert(a_rhs.nComp() == SpaceDim);

  // Component in a_lhs that put the data into.
  constexpr int comp = 0;

  // Mask for skipping computation on covered cells
  const EBISBox& ebisbox = a_lhs.getEBISBox();
  EBCellFAB      coveredMask(ebisbox, a_box, 1);
  coveredMask.setVal(1.0);
  coveredMask.setCoveredCellVal(-1.0, 0);

  // Kernel regions
  VoFIterator vofit(ebisbox.getIrregIVS(a_box), ebisbox.getEBGraph());

  // Hooks to single-valued data.
  BaseFab<Real>&       lhsReg = a_lhs.getSingleValuedFAB();
  const BaseFab<Real>& rhsReg = a_rhs.getSingleValuedFAB();
  const BaseFab<Real>& mask   = coveredMask.getSingleValuedFAB();

  auto regularKernel = [&](const IntVect& iv) -> void {
    lhsReg(iv, comp) = 0.0;

    if (mask(iv, comp) > 0.0) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        lhsReg(iv, comp) += rhsReg(iv, dir) * rhsReg(iv, dir);
      }
    }
  };

  // Irregular kernel. Same as the above.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    a_lhs(vof, comp) = 0.;

    for (int i = 0; i < SpaceDim; i++) {
      a_lhs(vof, comp) += a_rhs(vof, i) * a_rhs(vof, i);
    }
  };

  // Run the kernels
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

void
DataOps::computeMinValidBox(RealVect& a_lo, RealVect& a_hi, const RealVect a_normal, const RealVect a_centroid)
{
  CH_TIME("DataOps::computeMinValidBox");

  // TLDR: A stupid routine. The below code shifts planes until we "hit" the EB, but we
  //       should be able to compute up with a way better method of doing this.

  const int num_segments = 10;

  // Default values
  a_lo = -0.5 * RealVect::Unit;
  a_hi = 0.5 * RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const RealVect plane_normal = (sit() == Side::Lo) ? BASISREALV(dir) : -BASISREALV(dir);
      const RealVect plane_point  = RealVect::Zero - 0.5 * plane_normal; // Center point on plane
      const RealVect base_shift   = plane_normal / (1.0 * num_segments);

#if CH_SPACEDIM == 2
      Vector<RealVect> corners(2);
      const int        otherDir = (dir + 1) % SpaceDim;
      corners[0]                = plane_point - 0.5 * RealVect(BASISV(otherDir));
      corners[1]                = plane_point + 0.5 * RealVect(BASISV(otherDir));
#elif CH_SPACEDIM == 3
      Vector<RealVect> corners(4);
      const int        otherDir1 = (dir + 1) % SpaceDim;
      const int        otherDir2 = (dir + 2) % SpaceDim;
      corners[0]                 = plane_point - 0.5 * RealVect(BASISV(otherDir1)) - 0.5 * RealVect(BASISV(otherDir2));
      corners[1]                 = plane_point - 0.5 * RealVect(BASISV(otherDir1)) + 0.5 * RealVect(BASISV(otherDir2));
      corners[2]                 = plane_point + 0.5 * RealVect(BASISV(otherDir1)) - 0.5 * RealVect(BASISV(otherDir2));
      corners[3]                 = plane_point + 0.5 * RealVect(BASISV(otherDir1)) + 0.5 * RealVect(BASISV(otherDir2));
#endif

      // Shift corners in direction plane_normal with length base_shift. Keep track of the total
      // displacement of the plane.
      RealVect shift_vector = RealVect::Zero;
      bool     allInside    = DataOps::allCornersInsideEb(corners, a_normal, a_centroid);

      while (allInside) {

        // Shift the corners
        DataOps::shiftCorners(corners, base_shift);
        shift_vector += base_shift;

        // Check if shifted corners are inside EB
        allInside = DataOps::allCornersInsideEb(corners, a_normal, a_centroid);

        // If they are, we can change some components of a_lo
        if (allInside) {
          if (sit() == Side::Lo) {
            a_lo[dir] = -0.5 + shift_vector[dir];
          }
          else if (sit() == Side::Hi) {
            a_hi[dir] = 0.5 + shift_vector[dir];
          }
        }
      }
    }
  }
}

bool
DataOps::allCornersInsideEb(const Vector<RealVect>& a_corners, const RealVect a_normal, const RealVect a_centroid)
{
  CH_TIME("DataOps::allCornersInsideEb");

  bool ret = true;

  // If any point it outside the EB, i.e. inside the domain boundary, return false.
  for (int i = 0; i < a_corners.size(); i++) {
    if (PolyGeom::dot((a_corners[i] - a_centroid), a_normal) > 0.0) {
      ret = false;
    }
  }

  return ret;
}

void
DataOps::shiftCorners(Vector<RealVect>& a_corners, const RealVect& a_distance)
{
  CH_TIME("DataOps::shiftCorners");

  for (int i = 0; i < a_corners.size(); i++) {
    a_corners[i] += a_distance;
  }
}

#include <CD_NamespaceFooter.H>
