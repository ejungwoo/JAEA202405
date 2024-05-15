#include "LoadLibrary.C"

void run(int runNumber=42)
{
  //gSystem -> Exec(Form("root -l SortSi.C+\\\(%d\\\)",runNumber));

  if (LoadLibrary(true)==false)
    return false;

  /*
  bool updateSource = true;
  if (updateSource)
  {
    // option: "k" to keep the shared library, "-" to build directly to "build"
    //
    if (gSystem -> CompileMacro("ChannelData.cpp","k-","","build")==false) return;
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so");

    if (gSystem -> CompileMacro("AnaCC.cpp","k-","","build")==false) return;
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so");
  }
  else {
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so");
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so");
  }
  */

  auto ana = new AnaCC();
  ana -> RunConversion(42,"/home/daquser/data/LiCD2Irrad/SortSi/test_data/");
}
