#include <CD_Driver.H>
#include <CD_LookupTable.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  ParmParse pp;

  std::string inputFile;
  int         numPoints;

  pp.get("input_file", inputFile);
  pp.get("num_points", numPoints);

  // Create a standard table and regularize it
  LookupTable1D<Real, 4> table;
  for (int i = 0; i < 10; i++) {
    table.addData(1.0 + i, 1.0 + 2 * i, 1.0 + 3.0 * i, 1.0 + 4.0 * i, 1.0 + 5.0 * i);
  }

  table.prepareTable(0, numPoints, LookupTable::Spacing::Uniform);
  table.prepareTable(0, numPoints, LookupTable::Spacing::Exponential);

  // Slice the table.
  LookupTable1D<Real, 3> slicedTable = table.slice<3>(0, 2, 1, 4);

  const auto& rawData    = table.getRawData();
  const auto& slicedData = slicedTable.getRawData();

  CH_assert(rawData.size() == slicedData.size());

  slicedTable.prepareTable(1, numPoints, LookupTable::Spacing::Uniform);
  slicedTable.prepareTable(1, numPoints, LookupTable::Spacing::Exponential);

  const auto interpRow = slicedTable.interpolate(std::numeric_limits<Real>::max());

  ChomboDischarge::finalize();
}
