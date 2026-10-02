// The plane z = 0 lying on a node plane, with a checkerboard of height a at its nodes. Every face in the node plane
// carries corners b+a, b-a, b+a, b-a, so its signs alternate, and the function at its centre is b. With b < 0 the
// face joins its fluid corners, so the cell above holds two solid tips of height about a in a cell of fluid, and
// the cell below holds a strip of fluid. Dialling a down to round-off reproduces the NoisePlane configuration;
// dialling it up makes the tips real geometry. Either way the cell above is kept: its fluid is one piece and
// fills nearly all of it.
#include <CD_Driver.H>
#include <CD_GeometryStepper.H>

using namespace ChomboDischarge;
using namespace Physics::Geometry;

class SaddleIF : public BaseIF
{
public:
  SaddleIF(const Real a_spacing, const Real a_amplitude, const Real a_offset, const bool a_isolated)
    : m_spacing(a_spacing), m_amplitude(a_amplitude), m_offset(a_offset), m_isolated(a_isolated)
  {}

  Real
  value(const RealVect& a_point) const override
  {
    // positive below the plane (the electrode), negative above it
    Real checker = std::cos(M_PI * a_point[0] / m_spacing) * std::cos(M_PI * a_point[1] / m_spacing);

    // isolated: only the four nodes of the face [0,h]^2 carry the checkerboard
    if (m_isolated) {
      const Real slack = 0.25 * m_spacing;

      for (int dir = 0; dir < SpaceDim - 1; dir++) {
        if (a_point[dir] < -slack || a_point[dir] > m_spacing + slack) {
          checker = 0.0;
        }
      }
    }

    return -a_point[SpaceDim - 1] + m_amplitude * checker + m_offset;
  }

  BaseIF*
  newImplicitFunction() const override
  {
    return static_cast<BaseIF*>(new SaddleIF(m_spacing, m_amplitude, m_offset, m_isolated));
  }

protected:
  Real m_spacing;
  Real m_amplitude;
  Real m_offset;
  bool m_isolated;
};

class Saddle : public ComputationalGeometry
{
public:
  Saddle()
  {
    ParmParse pp("Saddle");

    Real spacing;
    Real amplitude;
    Real offset;
    bool isolated;

    pp.get("spacing", spacing);
    pp.get("amplitude", amplitude);
    pp.get("offset", offset);
    pp.get("isolated", isolated);

    m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new SaddleIF(spacing, amplitude, offset, isolated)), true));
  }
};

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new Saddle());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<GeometryStepper>(new GeometryStepper());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  ChomboDischarge::finalize();
}
