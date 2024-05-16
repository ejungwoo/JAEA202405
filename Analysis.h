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
#include "TMath.h"
#include "TList.h"
#include "TSystemDirectory.h"
#include "TApplication.h"

#include <fstream>
#include <iostream>
using namespace std;

#include "ChannelData.cpp"

class Analysis
{
  public:
    Analysis() { coutx.open("out/dummy_stream"); }
    ~Analysis() {}

    void RunConversion(int runNo, TString pathToInputFile);
    void RunConversionOnline(int runNo, TString pathToInputFile);

    void ReadSummaryFile(TString fileName);
    void AnalyzeChannelEnergy(int module, int channelID, bool drawHist=false);
    void AnalyzeChannelEnergy(int gid=-1, bool drawHist=false);

    void InitializeDrawing();

  private:
    void InitializeConversion();
    void ConfigureDateTime();
    void MakeOutputFile();
    void ReadDataFile();
    bool FillDataTree();
    void EndOfConversion();

    void UpdateCvsOnline();
    void SetAttribute(TH1* hist, TPad* pad=(TPad*)nullptr);

  private:
    void Convert2APParameters(double* par, double &mean1, double &sigma1, double &amplitude1, double &mean2, double &sigma2, double &amplitude2);

  public:
    void SetReturnIfNoFile  (bool value) { fReturnIfNoFile = value; }
    void SetIgnoreFileUpdate(bool value) { fIgnoreFileUpdate = value; }
    int  GetGlobalID(UShort_t module, UShort_t channel);
    void GetModCh(int globalID, UShort_t &module, UShort_t &channel);

  private:
    std::ofstream coutx;

  // reading
  private:
    int fRunNo;
    TString fPathToInput = "/home/daquser/data/LiCD2Irrad/analysis/input/";
    TString fPathToOutput = "/home/daquser/data/LiCD2Irrad/analysis/out/";
    TString fDateTime;
    TString fDivisionMax;

    Long64_t fCountEvents = 0;
    TString fFileNameOut;
    TFile* fFileOut = nullptr;
    TTree* fTreeOut = nullptr;

    bool fReturnIfNoFile = false;
    bool fIgnoreFileUpdate = false;
    Long64_t fTimeStampDecreased = -1; ///< decreased point of time stamp
    Long64_t fTimeStampPrev = -1; ///< previous time stamp
    Long64_t bTimeStamp = -1; ///< branch value for time stamp
    Long64_t bTimeStampDist = -1; ///< branch value for distance do previous time stamp
    const int fNumModules = 7;
    const int fNumChannels = 16;

    TClonesArray *fChannelArray = nullptr;
    Int_t fCountMultHit[5];
    Int_t fCountChannels = 0;
    Int_t fCountAllChannels = 0;
    Int_t fCountTimeStampDecrease = 0;

    bool fExitAnalysis = false;
    bool fExitRoot = false;

  private:
    TFile* fFileSummary = nullptr;
    TTree* fTreeSummary = nullptr;
    Long64_t fNumEventsSummary = -1;

    TCanvas* fCvsEnergy = nullptr;

  // mapping
  private:
    /// detector type
    const int kS1Junction   = 0;
    const int kS1Ohmic      = 1;
    const int kS3Junction   = 2;
    const int kS3Ohmic      = 3;
    const int kdEDetector   = 4;
    const int kScintillator = 5;
    const int kFaradayCup   = 6;

    const int fNumS1Junction   = 32;
    const int fNumS1Ohmic      = 0; // 16
    const int fNumS3Junction   = 24;
    const int fNumS3Ohmic      = 32;
    const int fNumdEDetector   = 11;
    const int fNumScintillator = 1;
    const int fNumFaradayCup   = 1;

  // drawing for data checking
  private:
    Long64_t fUpdateDrawingEveryNEvent = -1;
    Long64_t fCountEventsForUpdate = 0;
    TCanvas* fCvsOnline = nullptr;
    TH1D* fHistDTS = nullptr; ///< difference between TS to check TS rate
    TH1D* fHistChCount = nullptr; ///< single channel event count
    TH1D* fHistEnergy = nullptr; ///< energy
    TH2D* fHistEVSCh = nullptr; ///< energy vs channel-id

  private:
    int fNumDTS = 500;
    int fMaxDTS = 10000;
    int fNumE = 8200;
    int fMaxE = 8200;
    int fNumCh = 112;
    int fMaxCh = 112;

  // drawing for analysis
  private:
    TH2D* fHistdEVSS1 = nullptr; ///< dE vs dE + S1-energy
    TH1D* fHistTritonE = nullptr; ///< gated triton energy
    TH1D* fHistCrEx = nullptr; ///< 48Cr excitation energy
    TH2D* fHistS3VSEx = nullptr; ///< coincidence S3 vs excitation energy
};

#endif
