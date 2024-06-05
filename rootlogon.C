void rootlogon()
{
    int load = gSystem -> Load("build/libJAEAExp.so");
    if (load>=0)
        cout << "JAEAExp" << endl;
}
