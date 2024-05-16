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
    if (gSystem -> CompileMacro("Analysis.cpp","k-","","build")==false) {
      cout << "Failed to create library for Analysis!" << endl;
      return false;
    }
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/Analysis_cpp.so")<0) {
      cout << "Failed to load library for Analysis!" << endl;
      return false;
    }
  }
  else {
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/ChannelData_cpp.so")<0) {
      cout << "Failed to load library for ChannelData!" << endl;
      return false;
    }
    if (gSystem -> Load("/home/daquser/data/LiCD2Irrad/SortSi/build/Analysis_cpp.so")<0) {
      cout << "Failed to load library for Analysis!" << endl;
      return false;
    }
  }

  //gSystem -> Load("run_conversion.C");

  return true;
}
