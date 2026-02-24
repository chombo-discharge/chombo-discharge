#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_TriangleLookup.H>
#include <CD_Triangle.H>
#include <CD_FieldStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

using T         = Real;
using Meta      = std::array<Real, 3>;
using BV        = EBGeometry::BoundingVolumes::AABBT<T>;
constexpr int K = 4;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);
#warning "Debug code enabled. Remember to fix before completing PR"
#if 0
  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new DiskProfiledPlane());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<FieldStepper<FieldSolverGMG>>(new FieldStepper<FieldSolverGMG>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();
#endif
  const auto triSDF = EBGeometry::Parser::readIntoTriangleBVH<T, Meta, BV, K>("x_low.vtk");

  Triangle tri;
  ChomboDischarge::finalize();
}
