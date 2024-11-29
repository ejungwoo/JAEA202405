void rootlogon()
{
    int load = gSystem -> Load("build/libJAEAExp");
    if (load>=0)
        cout << "JAEAExp" << endl;
}
