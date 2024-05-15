#ifndef ANA_CRPTCR_CD2_HH
#define ANA_CRPTCR_CD2_HH

#include "TString.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <fstream>
#include <iostream>
using namespace std;

#include "ChannelData.cpp"

class AnaCC
{
  public:
    AnaCC() { coutx.open("out/dummy_stream"); }
    ~AnaCC() {}

    void RunConversion(int runNo, TString pathToInputFile);
    void RunConversionAndDraw(int runNo, TString pathToInputFile);

    void InitializeConversion();
    void ConfigureDateTime();
    void MakeOutputFile();
    void ReadDataFile();
    void EndOfConversion();

    void InitializeDrawing();
    void Draw() {}
    void FillDrawing() {}

  private:
    std::ofstream coutx;

  // reading
  private:
    int fRunNo;
    TString fPathToInput = "/home/daquser/data/LiCD2Irrad/SortSi/test_data/";
    TString fPathToOutput = "/home/daquser/data/LiCD2Irrad/SortSi/out/";
    TString fDateTime;
    TString fDivisionMax;

    Long64_t fCountEvents = 0;
    TString fFileNameOut;
    TFile* fFileOut;
    TTree* fTreeOut;

    bool fReturnIfNoFile = true;
    Long64_t fTimeStampDecreased = -1; ///< decreased point of time stamp
    Long64_t fTimeStampPrev = -1; ///< previous time stamp
    Long64_t bTimeStamp = -1; ///< branch value for time stamp
    Long64_t bTimeStampDist = -1; ///< branch value for distance do previous time stamp
    const int fNumModules = 7;
    const int fNumChannels = 16;

    TClonesArray *fChannelArray = nullptr;
    Int_t fCountCoincidence[5];
    Int_t fCountChannels = 0;
    Int_t fCountAllChannels = 0;
    Int_t fCountTimeStampDecrease = 0;

  // drawing for data checking
  private:
    bool fFlagDrawing = false;
    TH1D* fHistDTS = nullptr; ///< difference between TS to check TS rate
    TH1D* fHistEnergy = nullptr; ///< energy
    TH2D* fHistEVSCh = nullptr; ///< energy vs channel-id

    int fNumDTS = 500;
    int fMaxDTS = 10000;

    int fNumE = 500;
    int fMaxE = 5000;

    int fNumCh = 112;
    int fMaxCh = 112;

  // drawing for analysis
  private:
    TH1D* fHistChDist = nullptr; ///< single channel distribution
    TH2D* fHistdEVSS1 = nullptr; ///< dE vs dE + S1-energy
    TH1D* fHistTritonE = nullptr; ///< gated triton energy
    TH1D* fHistCrEx = nullptr; ///< 48Cr excitation energy
    TH2D* fHistS3VSEx = nullptr; ///< coincidence S3 vs excitation energy
};

#endif
