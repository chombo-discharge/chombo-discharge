// Scratch: an infinite plate, rotated 45 degrees about z, thinner than a cell.  A cell straddling it has
// its two corners along the (1,1) diagonal in the gas and the other two inside the plate, which is the
// alternating configuration -- the fluid is two lumps and the surface enters as two sheets.
#include <CD_Driver.H>
#include <CD_GeometryStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Geometry;

class SlabIF : public BaseIF
{
public:
  SlabIF(const RealVect& a_normal, const RealVect& a_point, const Real a_thickness)
    : m_normal(a_normal / a_normal.vectorLength()), m_point(a_point), m_thickness(a_thickness)
  {}

  Real
  value(const RealVect& a_point) const override
  {
    // positive inside the plate, negative in the gas: an electrode's function is negative outside it
    return 0.5 * m_thickness - std::abs(m_normal.dotProduct(a_point - m_point));
  }

  BaseIF*
  newImplicitFunction() const override
  {
    return static_cast<BaseIF*>(new SlabIF(m_normal, m_point, m_thickness));
  }

protected:
  RealVect m_normal;
  RealVect m_point;
  Real     m_thickness;
};

class ThinPlate : public ComputationalGeometry
{
public:
  ThinPlate()
  {
    ParmParse pp("ThinPlate");

    Real         thickness;
    Real         offset;
    Vector<Real> n(SpaceDim);

    pp.get("thickness", thickness);
    pp.get("offset", offset);
    pp.getarr("normal", n, 0, SpaceDim);

    const RealVect normal = RealVect(D_DECL(n[0], n[1], n[2]));
    const RealVect point  = offset * normal / normal.vectorLength();

    m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new SlabIF(normal, point, thickness)), true));
  }
};

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new ThinPlate());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<GeometryStepper>(new GeometryStepper());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
