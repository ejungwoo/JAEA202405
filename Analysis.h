#ifndef ANA_CRPTCR_CD2_HH
#define ANA_CRPTCR_CD2_HH

#define _NUMBER_OF_MODULES_ 7
#define _NUMBER_OF_CHANNELS_ 16

#include "TClonesArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TMath.h"
#include "TText.h"
#include "TList.h"
#include "TLine.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TParameter.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TSystemDirectory.h"

#include <fstream>
#include <iostream>
using namespace std;

#include "ChannelData.cpp"

class Analysis
{
  public:
    Analysis();
    ~Analysis() {}

  // general
  public:
    int  GetFakeModule(int module); ///< Get fake module index from real module number 
    int  GetRealModule(int fake); ///< Get real module number from fake module index
    int  GetGlobalID(UShort_t module, UShort_t channel);
    void GetModCh(int globalID, UShort_t &module, UShort_t &channel);
    int  GetDetectorType(int module, int channelID=0);

  // conversion
  public:
    void RunConversion(int runNo, TString pathToInputFile);
    void AddEnergyCalibration(TString name);
    void SetDrawOnline(Long64_t everyNEvents=10000) { fUpdateDrawingEveryNEvent = everyNEvents; }
    void SetSkipTSError(bool ignore=true) { fSkipTSError = ignore; }
    void SetStopAtTSError(bool stop) { fStopAtTSError = stop; }
    void SetFileNumberRange(int d1, int d2) { fFileNumberRange1 = d1, fFileNumberRange2 = d2; }
    void SetEventCountLimit(Long64_t limit) { fEventCountLimit = limit; }
    void SetReturnIfNoFile(bool value) { fReturnIfNoFile = value; }
    void SetIgnoreFileUpdate(bool value) { fIgnoreFileUpdate = value; }
    void SetADCThreshold(UShort_t value) { fADCThreshold = value; }
    void SetOutputFileName(TString name) { fFileNameOut = name; }
    void SetAutoUpdateDrawing(bool value) { fAutoUpdateDrawing = value; }
    void SetShowEnergyConversion(bool value) { fShowEnergyConversion = value; }
    void SetCoincidenceTSRange(int value) { fCoincidenceTSRange = value; }
    void SetCoincidenceMult(int value) { fCoincidenceMult = value; }

  // alpha energy calibration
  public:
    void ReadSummaryFile(TString fileName);
    void AnalyzeAlphaTestModule(int module, bool drawHist=false, TString fileName="");
    bool AnalyzeAlphaTest(int module, int channelID, bool drawHist=false, TVirtualPad* cvs=(TVirtualPad*)nullptr);

  private:
    void InitializeConversion();
    void ConfigureDateTime();
    void MakeOutputFile();
    void ReadDataFile();
    bool FillDataTree();
    void EndOfConversion();

    void WriteRunParameters(TFile* file, int option);
    void UpdateCvsOnline(bool firstDraw=false);
    void SetAttribute(TH1* hist, TVirtualPad* pad, int npad=1, bool td=false);
    void SetAttribute(TH1* hist, int npad=1, bool td=false) { SetAttribute(hist, (TVirtualPad*)nullptr, npad, td); }

    void InitializeMapping();
    void InitializeDrawing();
    void InitializeAlphaAnalysis(TString fileName="");

  private:
    Double_t FxTwoAlpha(Double_t *xy, Double_t *par);
    void Convert2APParameters(double* par, double &mean1, double &sigma1, double &amplitude1, double &mean2, double &sigma2, double &amplitude2);
    TF1* CreateFxTwoAlpha(TString name, double min, double max);

  private:
    std::ofstream coutx;

  // reading
  private:
    int fRunNo;
    TString fRunName;
    TString fPathToInput = "/home/daquser/data/LiCD2Irrad/analysis/input/";
    TString fPathToOutput = "/home/daquser/data/LiCD2Irrad/analysis/out/";
    TString fMapFileName = "/home/daquser/data/LiCD2Irrad/analysis/ModCh.in";
    TString fDateTime;
    int fFileNumberMax;
    int fFileNumberRange1 = 0;
    int fFileNumberRange2 = 10000;

    Long64_t fCountEvents = 0;
    Long64_t fEventCountLimit = -1;
    TString fFileNameOut;
    TFile* fFileOut = nullptr;
    TTree* fTreeOut = nullptr;

    TString fFileNameAlpha;
    TFile* fFileAlpha = nullptr;

    bool fFirstFileOpened = false;
    bool fReturnIfNoFile = false;
    bool fIgnoreFileUpdate = false;
    Long64_t fTimeStampPrevTrue = -1; ///< decreased point of time stamp
    Long64_t fTimeStampPrev = -1; ///< previous time stamp
    Long64_t bTimeStamp = -1; ///< branch value for time stamp
    Long64_t bTimeStampDist = -1; ///< branch value for distance do previous time stamp
    UShort_t bNumChannels = 0;
    UShort_t fADCThreshold = 0;

    const int fNumModules = _NUMBER_OF_MODULES_;
    const int fNumChannels = _NUMBER_OF_CHANNELS_;

    TClonesArray *fChannelArray = nullptr;
    Int_t fCoincidenceCount[5];
    Int_t fCoincidenceTSRange = 0;
    Int_t fCoincidenceMult = -1;
    Int_t fCountChannels = 0;
    Int_t fCountAllChannels = 0;
    Int_t fCountTSError = 0;
    bool fSkipTSError = false;
    bool fStopAtTSError = false;
    bool fExitAnalysis = false;
    bool fExitRoot = false;

    bool fEnergyConversionSet = false;
    TF1* fFxEnergyConversion[_NUMBER_OF_CHANNELS_][_NUMBER_OF_CHANNELS_];

    int kSameEvent = 0;
    int kNextEvent = 1;
    int kTSError = 2;

  private:
    TFile* fFileSummary = nullptr;
    TTree* fTreeSummary = nullptr;
    Long64_t fNumEventsSummary = -1;

    Double_t fAlphaEnergy1 = 5.486;
    Double_t fAlphaEnergy2 = 5.443;
    Double_t fAEBR = 6.65625; ///< alpha energy branching ratio = 85.2 / 12.8 = 6.65625; 

    TF1* fFitAlpha  = nullptr;
    TF1* fFitAlpha1 = nullptr;
    TF1* fFitAlpha2 = nullptr;

  // mapping
  private:
    /// detector type
    TString fDetectorName[_NUMBER_OF_MODULES_+1];
    const int kDummyDetector = 0;
    const int kS1Junction    = 1;
    const int kS1Ohmic       = 2;
    const int kS3Junction    = 3;
    const int kS3Ohmic       = 4;
    const int kdEDetector    = 5;
    const int kScintillator  = 6;
    const int kFaradayCup    = 7;

    const int fNumS1Junction   = 32;
    const int fNumS1Ohmic      = 0; // 16
    const int fNumS3Junction   = 24;
    const int fNumS3Ohmic      = 32;
    const int fNumdEDetector   = 11;
    const int fNumScintillator = 1;
    const int fNumFaradayCup   = 1;

    //int fMapDetectorType[_NUMBER_OF_MODULES_][_NUMBER_OF_CHANNELS_];
    int **fMapDetectorType;
    int **fMapDetectorChannel;
    bool **fMapDetectorReplaced;
    int **fMapDetectorRFMod;
    int **fMapDetectorRFMCh;

  // drawing for data checking
  private:
    bool fAutoUpdateDrawing = false;
    Long64_t fUpdateDrawingEveryNEvent = 0;
    Long64_t fCountEventsForUpdate = 0;
    TCanvas* fCvsOnline = nullptr;
    TVirtualPad* fPadTSDist1 = nullptr;
    TVirtualPad* fPadTSDist2 = nullptr;
    TH1D* fHistTSDist1 = nullptr;
    TH1D* fHistTSDist2 = nullptr;
    TH1D* fHistChCount = nullptr;
    TH1D* fHistADC = nullptr; ///< ADC
    TH1D* fHistE = nullptr; ///< energy
    TH2D* fHistAVSCh = nullptr; ///< ADC vs channel-id
    TH2D* fHistEVSCh = nullptr; ///< energy vs channel-id
    bool fShowEnergyConversion = true;

  private:
    int    fNumADC = 8200;
    int    fMaxADC = 8200;
    int    fNumCh = _NUMBER_OF_MODULES_*_NUMBER_OF_CHANNELS_;
    int    fMaxCh = _NUMBER_OF_MODULES_*_NUMBER_OF_CHANNELS_;
    int    fNumE = 8200;
    double fMaxE = 25;

  // drawing for analysis
  private:
    TH2D* fHistdEVSS1 = nullptr; ///< dE vs dE + S1-energy
    TH1D* fHistTritonE = nullptr; ///< gated triton energy
    TH1D* fHistCrEx = nullptr; ///< 48Cr excitation energy
    TH2D* fHistS3VSEx = nullptr; ///< coincidence S3 vs excitation energy
};

#endif
