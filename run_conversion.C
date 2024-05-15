//#include "LoadLibrary.h"

void run_conversion(int runNumber=42, bool updateSource=true)
{
  //if (LoadLibrary(updateSource)==false) return;

  auto ana = new AnaCC();
  ana -> RunConversion(42,"/home/daquser/data/LiCD2Irrad/SortSi/test_data/");
}
