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

    void RunConversion(int runNo, TString pathToInputFile);

    void InitializeConversion();
    void ConfigureDateTime();
    void MakeOutputFile();
    void ReadDataFile();
    void EndOfConversion();

    void Draw() {}
    void FillHistogram() {}

  private:
    int fRunNo;
    TString fPathToInput;
    TString fPathToOutput;
    TString fDateTime;
    TString fDivisionMax;
    vector<TString> fListOfInputFiles;

    Long64_t fCountEvents;
    TString fFileNameOut;
    TFile* fFileOut;
    TTree* fTreeOut;

    Long64_t fTimeStampPrev;
    Long64_t fTimeStampSaved;
    int fNumModules;
    int fNumChannels;

    TClonesArray *fChannelArray;
    Long64_t bTimeStampFull;
    Int_t fCountChannels;
    Int_t fCountAllChannels;
    Int_t fCountTimeStampDecrease;
};

#endif
