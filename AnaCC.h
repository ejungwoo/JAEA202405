#ifndef ANA_CRPTCR_CD2_HH
#define ANA_CRPTCR_CD2_HH

#include "TString.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"

#include "ChannelData.cpp"

class AnaCC
{
  public:
    AnaCC() {}
    ~AnaCC() {}

    void Run(int runNo, TString pathToInputFile);

    void ConfigureDateTime();
    void MakeOutputFile();
    void ReadDataFile();
    void FillHistogram();
    void EndOfAnalysis();

  private:
    int fRunNo = -1;
    // somehow below initialization do not work
    TString fPathToInput = "/home/daquser/data/LiCD2Irrad/SortSi/test_data/";
    TString fPathToOutput = "/home/daquser/data/LiCD2Irrad/SortSi/out/";
    TString fDateTime;
    TString fDivisionMax;
    vector<TString> fListOfInputFiles;

    Long64_t fCountEvents;
    TString fFileNameOut;
    TFile* fFileOut;
    TTree* fTreeOut;

    Long64_t fTimeStampPrev;
    int fNumModules;
    int fNumChannels;

    TClonesArray *fChannelArray;
    Int_t fCountChannels = 0;
};

#endif
