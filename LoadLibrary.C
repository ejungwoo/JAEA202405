bool LoadLibrary(bool updateSource=true)
{
  if (updateSource)
  {
    // option: "k" to keep the shared library, "-" to build directly to "build"
    if (gSystem -> CompileMacro("ChannelData.cpp","k-","","build")==false) return false;
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so");
    if (gSystem -> CompileMacro("AnaCC.cpp","k-","","build")==false) return false;
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so");
  }
  else {
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so");
    gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so");
  }

  return true;
}
