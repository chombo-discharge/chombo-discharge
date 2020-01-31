
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

#include "graph3.gnu.template.H"

std::string ReplaceString(std::string subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}

int main(int argc, char* argv[])
{
  std::fstream files[200];
  //std::vector<size_t>       workingSet;
  //std::vector<float>        bandwidth;
  if(argc == 1)
    {
      std::cerr<<"usage:  > streamMerge stream.dat.[0..N]"<<std::endl;
      return 0;
    }
  for(int i=0; i<argc-1; ++i)
    {
      files[i].open(argv[i+1], std::ios_base::in);
    }
  unsigned int totalWorkingSet=0;
  while(files[0]>>totalWorkingSet)
    {
      double bandwidth=0;
      double totalBandwidth=0;
      unsigned int workingSet=0;
      files[0]>>totalBandwidth;
      for(int i=1; i<argc-1; ++i)
        {
          files[i]>>workingSet;
          files[i]>>bandwidth;
          totalBandwidth+=bandwidth;
          totalWorkingSet+=workingSet;
        }
      std::cout<<totalWorkingSet<<" "<<totalBandwidth<<"\n";
    }
  
  std::string gplot(graph3_gnu_template);
  
  return 0;
}
      
