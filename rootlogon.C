#include "LoadLibrary.C"

void rootlogon()
{
  if (LoadLibrary(true))
    cout << "(rootlogon) Succesfully loaded library for JAEA Tandem Experiment Analysis" << endl;
  else
    cout << "(rootlogon) Error in LoadLibrary" << endl;
}
