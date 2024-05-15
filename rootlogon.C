#include "LoadLibrary.C"

void rootlogon()
{
  if (LoadLibrary(true))
    cout << "(rootlogon) Succesfully loaded library" << endl;
  else
    cout << "(rootlogon) Error in LoadLibrary" << endl;
}
