bool CompileAndLoad(const char* name, bool recompile, bool ismacro=false);

bool LoadLibrary(bool recompile=true)
{
    gSystem -> SetFlagsOpt("-std=c++11");
    CompileAndLoad("DetectorSetting", recompile);
    CompileAndLoad("ChannelData", recompile);
    CompileAndLoad("Analysis", recompile);
    //CompileAndLoad("draw_des1", recompile, true);

    return true;
}

bool CompileAndLoad(const char* name, bool recompile, bool ismacro)
{
    // option "k" keep the shared library, option "-" build lib directly to "build"
    TString pathToBuild = "/home/daquser/data/LiCD2Irrad/SortSi/build/";
    if (recompile) {
        if (gSystem -> CompileMacro(Form("%s.%s",name,(ismacro?"C":"cpp")),"k-","","build")==false) {
            cout << "!! Failed to create library for " << name << endl;
            return false;
        }
    }
    if (gSystem -> Load(pathToBuild+Form("%s_%s.so",name,(ismacro?"C":"cpp")))<0) {
        cout << "!! Failed to load library for " << name << endl;
        return false;
    }
    return true;
}

