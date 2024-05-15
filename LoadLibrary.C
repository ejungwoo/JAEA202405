bool LoadLibrary(bool updateSource=true)
{
  // option: "k" to keep the shared library, "-" to build directly to "build"
  if (updateSource)
  {
    gSystem -> SetFlagsOpt("-std=c++11");
    if (gSystem -> CompileMacro("ChannelData.cpp","k-","","build")==false) {
      cout << "Failed to create library for ChannelData!" << endl;
      return false;
    }
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so")<0) {
      cout << "Failed to load library for ChannelData!" << endl;
      return false;
    }
    if (gSystem -> CompileMacro("AnaCC.cpp","k-","","build")==false) {
      cout << "Failed to create library for AnaCC!" << endl;
      return false;
    }
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so")<0) {
      cout << "Failed to load library for AnaCC!" << endl;
      return false;
    }
  }
  else {
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so")<0) {
      cout << "Failed to load library for ChannelData!" << endl;
      return false;
    }
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/AnaCC_cpp.so")<0) {
      cout << "Failed to load library for AnaCC!" << endl;
      return false;
    }
  }

  return true;
}
