#include <CD_LookupTable.H>

int
main(int argc, char* argv[])
{
  constexpr size_t N = 1;

  LookupTable1D<double, N> table;

  table.addData(1.0,1.0);
  table.addData(2.0,2.0);
  table.addData(3.0,3.0);

  table.prepareTable(0,10,CoordinateSystem::Exponential);

  std::cout << "\n" << std::endl;
  table.outputRawData();
  std::cout << "\n";
  table.outputStructuredData();  

  table.writeRawData("raw.dat");
  table.writeStructuredData("structured.dat");  
  //  table.dumpRawData();
}



