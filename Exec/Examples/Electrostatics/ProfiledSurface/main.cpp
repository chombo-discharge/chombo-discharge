#include <CD_Driver.H>
#include <CD_FieldSolverGMG.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_TriangleCollection.H>
#include <CD_Triangle.H>
#include <CD_FieldStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Electrostatics;

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
  const auto x                  = EBGeometry::Vec3T<Real>::zero();
  const auto triangles          = DataParser::readTriangles("test.vtk", "potential");
  const auto triangleCollection = std::make_shared<TriangleCollection>(triangles);
  const auto closestTriangles   = triangleCollection->getClosestTriangles(x);

  if (procID() == 0) {
    const auto& tri = *(closestTriangles.front().first);
    const auto  y   = tri.projectToTrianglePlane(x);
    const auto  z   = tri.interpolate(y);

    std::cout << "v0 = " << tri.getVertexPositions()[0] << "\n";
    std::cout << "v1 = " << tri.getVertexPositions()[1] << "\n";
    std::cout << "v2 = " << tri.getVertexPositions()[2] << "\n";
    
    std::cout << "num triangles = " << triangles.size() << "\n";
    std::cout << "closest triangles size = " << closestTriangles.size() << "\n";
    std::cout << "x = " << x << "\n";
    std::cout << "y = " << y << "\n";    
    std::cout << "inside(x) = " << tri.isInside(x) << "\n";
    std::cout << "inside(y) = " << tri.isInside(y) << "\n";
    std::cout << "dist(x) = " << tri.signedDistance(x) << "\n";
    std::cout << "dist(y) = " << tri.signedDistance(y) << "\n";
    std::cout << "interp = " << z << "\n";
  }

  ChomboDischarge::finalize();
}
