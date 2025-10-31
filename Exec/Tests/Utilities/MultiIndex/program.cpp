#include <ParmParse.H>
#include <CD_MultiIndex.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  constexpr int a_order = 3;

  MultiIndex mi(a_order);

  if (procID() == 0) {

    for (MultiIndex mi(a_order); mi.ok(); ++mi) {
      std::cout << mi.getCurrentIndex() << "\t" << mi.factorial() << std::endl;
    }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
