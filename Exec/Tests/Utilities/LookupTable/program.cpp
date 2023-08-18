#include <CD_Driver.H>
#include <CD_LookupTable.H>
#include <ParmParse.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Build class options from input script and command line options
  const std::string input_file = argv[1];
  ParmParse         pp(argc - 2, argv + 2, NULL, input_file.c_str());

  // Get the input file name
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

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return 0;
}
