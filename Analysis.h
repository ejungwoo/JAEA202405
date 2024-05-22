#ifndef ANA_CRPTCR_CD2_HH
#define ANA_CRPTCR_CD2_HH

#define NUMBER_OF_MODULES 7
#define NUMBER_OF_CHANNELS 16

#define NUMBER_OF_DETECTORS 7
#define NUM_MAX_DCH 32
const int kX   = 0;
const int kdE  = 1;
const int kS1J = 2;
const int kS3J = 3;
const int kS3O = 4;
const int kSC  = 5;
const int kFC  = 6;

#include "TRandom.h"
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
#include "TCutG.h"
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

//#include <chrono>
#include <fstream>
#include <iostream>
//#include <ios>
using namespace std;

#include "ChannelData.cpp"

class Analysis
{
  public:
    Analysis();
    ~Analysis() {}

  // conversion
  public:
    void RunConversion(int runNo);
    void SetBeamEnergy(double beamEnergy) { fBeamEnergy = beamEnergy; }
    void SetConversionFile(TString fileName="");
    void AddAlphaCalibrationFile(TString name);
    void SetDrawOnline(Long64_t everyNEvents=100000) { fUpdateDrawingEveryNEvent = everyNEvents; }
    void SetSkipTSError(bool skip) { fSkipTSError = skip; }
    void SetStopAtTSError(bool stop) { fStopAtTSError = stop; }
    void SetFileNumberRange(int d1, int d2) { fFileNumberRange1 = d1, fFileNumberRange2 = d2; }
    void SetEventCountLimit(Long64_t limit) { fEventCountLimit = limit; }
    void SetLineCountLimit(Long64_t limit) { fLineCountLimit = limit; }
    void SetReturnIfNoFile(bool value) { fReturnIfNoFile = value; }
    void SetIgnoreFileUpdate(bool value) { fIgnoreFileUpdate = value; }
    void SetADCThreshold(UShort_t value) { fADCThreshold = value; }
    void SetAutoUpdateRun(bool value) { fAutoUpdateRun = value; }
    void SetAutoUpdateDrawing(bool value) { fAutoUpdateDrawing = value; }
    void SetShowEnergyConversion(bool value) { fShowEnergyConversion = value; }
    void SetCoincidenceTSRange(int value) { fCoincidenceTSRange = value; }
    void SetCoincidenceMultRange(int r1, int r2) {
      fCoincidenceMultRange1 = r1;
      fCoincidenceMultRange2 = r2;
    }
    void SetfExcludedES1Coincidence(bool value) { fExcludedES1Coincidence = value; }
    void SetdES1Coincidence(bool value) {
      fCoincidenceMultRange1 = -1;
      fCoincidenceMultRange2 = -1;
      fdES1CoincidenceMode = value;
    } 
    void SetdES1S3Coincidence(bool value) {
      fCoincidenceMultRange1 = -1;
      fCoincidenceMultRange2 = -1;
      fdES1S3CoincidenceMode = value;
    } 
    void SetNumADC(int n) { fNumADC = n; }
    void SetNumE(int n) { fNumE = n; }

    void SetLocalDetectorChannelCut(int det, int dch) { fChosenDet = det; fChosenDCh = dch; }
    void SetLocalEnergyRange(double r1, double r2) { fEnergyRange1 = r1; fEnergyRange2 = r2; }

    void AddLocalS1StripHist(int strip) { fS1ChosenS1Strips.push_back(strip); }


  private:
    void ConfigureDateTime();
    void ReadDataFile();
    bool FillDataTree();
    void EndOfConversion();
    void PrintConversionSummary();

    void ResetTriggerParameters();
    void ResetEventParameters();

    bool CheckOpenFileStatus2();
    bool CheckOpenFileStatus1();
    bool CheckDataLineCondition(double adc, int eventStatus, Long64_t timeStamp);
    bool CheckEventCondition(double de, double ee);
    void FillLocalHistograms();
    bool AskContinueRun(TString message="");
    int  fCountAskContinueRun = 0;
    bool UpdateDrawing();
    bool AskUpdateDrawing(TString message="");
    Long64_t GetFileSize(TString fileName);
    TString fInputFileName;

    void WriteRunParameters(TFile* file, int option);
    void UpdateCvsOnline(bool firstDraw=false, bool write=false);
    bool CheckStopFile();
    void SetAttribute(TH1* hist, TVirtualPad* pad, int npad=1, bool is2D=false);
    void SetAttribute(TH1* hist, int npad=1, bool is2D=false) { SetAttribute(hist, (TVirtualPad*)nullptr, npad, is2D); }

    void InitializeAnalysis();
    void InitializeDrawing();

  // drawing for data checking
  private:
    double fBeamEnergy = 0;
    bool fAutoUpdateRun = false;
    bool fAutoUpdateDrawing = false;
    int fUpdateAfterXSec = 5;
    Long64_t fUpdateDrawingEveryNEvent = 0;
    Long64_t fCountEventsForUpdate = 0;
    TCanvas* fCvsOnline = nullptr;
    TCanvas* fCvsOnline2 = nullptr;
    TVirtualPad* fVPadEx      = nullptr;
    TVirtualPad* fVPadChCount = nullptr;
    TVirtualPad* fVPadADC     = nullptr;
    TVirtualPad* fVPadEVSCh   = nullptr;
    TVirtualPad* fVPadEVSS1Strip = nullptr;
    TVirtualPad* fVPadEVSS3Strip = nullptr;
    TVirtualPad* fVPadEVSAngle   = nullptr;
    TVirtualPad* fVPaddEVSE   = nullptr;
    TVirtualPad* fVPadTriggerRate = nullptr;
    TVirtualPad* fVPadEventRate = nullptr;
    TVirtualPad* fVPadBeamCountInTime = nullptr;
    TVirtualPad* fVPadEventCountInTime = nullptr;
    TVirtualPad* fVPadLocalCountInTime = nullptr;
    TVirtualPad* fVPadStripCountInTime = nullptr;
    TVirtualPad* fVPadRatioCountInTime = nullptr;
    TVirtualPad* fVPadTSDist1 = nullptr;
    TVirtualPad* fVPadTSDist2 = nullptr;
    //TVirtualPad* fVPadEVSStrip = nullptr;

    TVirtualPad* fVPadProtonCountInTime = nullptr;
    TVirtualPad* fVPadDeuteronCountInTime = nullptr;
    TVirtualPad* fVPadTritonCountInTime = nullptr;

    TH1D* fHistBeamCountInTime = nullptr;
    TH1D* fHistEventCountInTime = nullptr;
    TH1D* fHistLocalCountInTime = nullptr;
    TH1D* fHistStripCountInTime[17];
    TH1D* fHistRatioCountInTime[17];
    TH1D* fHistProtonCountInTime = nullptr;
    TH1D* fHistDeuteronCountInTime = nullptr;
    TH1D* fHistTritonCountInTime = nullptr;

    TH1D* fHistTSDist1 = nullptr;
    TH1D* fHistTSDist2 = nullptr;
    TH1D* fHistChCount = nullptr;
    TH1D* fHistADC = nullptr; ///< ADC
    TH1D* fHistE = nullptr; ///< energy

    TH2D* fHistAVSCh = nullptr; ///< ADC vs channel-id
    TH2D* fHistEVSCh = nullptr; ///< energy vs channel-id
    TH2D* fHistEVSS1Strip = nullptr;
    TH2D* fHistEVSS3Strip = nullptr;
    TH2D* fHistEVSAngle   = nullptr;

    TH2D* fHistdAVSA = nullptr;
    TH2D* fHistdEVSE = nullptr;
    //TH2D* fHistEVSStrip = nullptr;

    TH1D* fHistTriggerRate = nullptr;
    TH1D* fHistTriggerRateError = nullptr;
    TH1D* fHistEventRate = nullptr;
    TH1D* fHistEventRateError = nullptr;
    TH1D* fHistEx = nullptr;

    bool fShowEnergyConversion = false;

  // general
  public:
    TString GetOutputPath() const { return fPathToOutput; }
    int  GetModuleIndex(int module); ///< Get module index from real module number 
    int  GetGlobalID(UShort_t module, UShort_t channel);
    void GetModCh(int globalID, UShort_t &module, UShort_t &channel);
    void DetectorToModule(int det, int dch, int &midx, int &mch);
    //int  GetDetectorType(int midx, int chid)    const { return fMapDetectorType[midx][chid]; }
    //int  GetDetectorChannel(int midx, int chid) const { return fMapDetectorChannel[midx][chid]; }
    //int  GetModuleGroup(int midx, int chid)     const { return fMapDetectorGroup[midx][chid]; }
    //bool IsModuleReplaced(int midx, int chid)   const { return fMapDetectorReplaced[midx][chid]; }
    //int  GetRFMod(int midx, int chid)           const { return fMapDetectorRFMod[midx][chid]; }
    //int  GetRFMCh(int midx, int chid)           const { return fMapDetectorRFMCh[midx][chid]; }
    //int  GetModuleIndex(int module)             const { return fMapFEToModuleIndex[module]; }

    TString GetDetectorTitle(int midx, int mch=0, bool addChannel=false);
    double GetCalibratedEnergy(int module, int mch, int adc);
    bool fECalErrorWasSentOut = false;

  private:
    int          fNumADC = 8200;
    const int    fMaxADC = 8200;
    const int    fNumCh = NUMBER_OF_MODULES*NUMBER_OF_CHANNELS;
    const int    fMaxCh = NUMBER_OF_MODULES*NUMBER_OF_CHANNELS;

    int          fNumE = 400;
    const double fMaxE = 25;
    const int    fNumdE = 400;
    const double fMaxdE = 12;
    int          fNumEE = 400;
    const double fMaxEE = 35;

    const int    fNumStrips = 16;
    const int    fMaxStrips = 17;

    const int    fNumRate = 200;
    const int    fMaxRate = 5000;

    const int    fNumEx = 500;
    const int    fMaxEx = 15;

    const int    fMaxTime = 60*6;

    const double fSecondPerTS = 200.*1.e-9;
    const double fTSPerSecton = 1./(200.*1.e-9);

  // mapping
  private:
    /// detector type
    const int fMaxDCh = NUM_MAX_DCH;
    TString fDetectorName[NUMBER_OF_DETECTORS];
    const int fNumDetectors  = NUMBER_OF_DETECTORS;
    const int kDummyDetector = kX;
    const int kdEDetector    = kdE;
    const int kS1Junction    = kS1J;
    const int kS3Junction    = kS3J;
    const int kS3Ohmic       = kS3O;
    const int kScintillator  = kSC;
    const int kFaradayCup    = kFC;

    int **fMapDetectorType;
    int **fMapDetectorChannel;
    bool **fMapDetectorReplaced;
    int **fMapDetectorGroup;
    int **fMapDetectorRFMod;
    int **fMapDetectorRFMCh;
    int fMapFEToModuleIndex[20];
    int** fMapDetectorToGlobalID;
    int* fMapGlobalIDToModuleIndex;
    int* fMapGlobalIDToMCh;
    double fMapS1ChToAngle1[33];
    double fMapS1ChToAngle2[33];
    double fMapS1ChToAngle[33];
    double fMapS3ChToAngle1[33];
    double fMapS3ChToAngle2[33];
    double fMapS3ChToAngle[33];
    int fMapS1ChToStrip[33];
    int fMapS3ChToStrip[33];

  // alpha energy calibration
  public:
    void ReadSummaryFile(TString fileName);
    void SetAlphaTestFile(TString fileName="");
    void AnalyzeAlphaTestModule(int module, bool drawHist=false, TString fileName="");
    bool AnalyzeAlphaTest(int module, int mch, bool drawHist=false, TVirtualPad* cvs=(TVirtualPad*)nullptr);
    Double_t FxTwoAlpha(Double_t *xy, Double_t *par);
    void Convert2APParameters(double* par, double &mean1, double &sigma1, double &amplitude1, double &mean2, double &sigma2, double &amplitude2);
    TF1* CreateFxTwoAlpha(TString fileName, double min, double max);

  private:
    TString fFileNameAlpha;
    TFile* fFileAlpha = nullptr;

    Double_t fAlphaEnergy1 = 5.486;
    Double_t fAlphaEnergy2 = 5.443;
    Double_t fAEBR = 6.65625; ///< alpha energy branching ratio = 85.2 / 12.8 = 6.65625; 

    TF1* fFitAlpha  = nullptr;
    TF1* fFitAlpha1 = nullptr;
    TF1* fFitAlpha2 = nullptr;

  public:
    static void MakeCutGFile(int pdt);
    static void CallCutGFile(int pdt);
    //void GetCutGFile(int pdt);
    void SetExcludeDESECutGFile(TString fileName);
    void SetEnergyCutGFile(TString fileName);

  public:
    double EvalEx(double tp, double tt, double theta);

  private:
    bool fTritonCutGIsSet = false;
    TCutG* fProtonCutG = nullptr;
    TCutG* fDeuteronCutG = nullptr;
    TCutG* fTritonCutG = nullptr;
    TCutG* fExcludeDESECutG = nullptr;
    TCutG* fEnergyCutG = nullptr;

   public:
    void SetRunNo(int val) { fRunNo = val; }
    void SetPathToInput (TString val) { fPathToInput = val; }
    void SetPathToOutput(TString val) { fPathToOutput = val;}
    void SetMapFileName (TString val) { fMapFileName = val; }
    void SetDetFileName (TString val) { fDetFileName = val; }

  // reading
  private:
    int fRunNo = -1;
    TString fRunName;
    TString fDateTime;
    TString fPathToInput = "/home/daquser/data/LiCD2Irrad/analysis/input/";
    TString fPathToOutput = "/home/daquser/data/LiCD2Irrad/analysis/out/";
    TString fMapFileName = "/home/daquser/data/LiCD2Irrad/analysis/ModCh.in";
    TString fDetFileName = "/home/daquser/data/LiCD2Irrad/analysis/DetectorSetting.in";
    int fFileNumberMax = 100;
    int fFileNumberRange1 = 0;
    int fFileNumberRange2 = 100;

    ifstream fDataFile;
    const streamsize fDataFileSizeReadMax = 500000000; // 500 MB
    streamsize fDataFileSize;
    streamsize fDataFileSizeOld = 0; // Size of opened Raw-Data file
    UShort_t fDataFileNumber;

    Long64_t fCountEvents = 0;
    Long64_t fCountBeam = 0;

    Long64_t fEventCountLimit = -1;
    Long64_t fLineCountLimit = -1;
    Long64_t fLastDataPos;
    Long64_t fLastDataFileSize;

    TFile* fFileSummary = nullptr;
    TTree* fTreeSummary = nullptr;
    Long64_t fNumEventsSummary = -1;

    TString fFileNameOut;
    TFile* fFileOut = nullptr;
    TTree* fTreeOut = nullptr;

    bool fFirstFileOpened = false;
    bool fReturnIfNoFile = false;
    bool fIgnoreFileUpdate = false;
    Long64_t fTimeStampPrevTrue = -1; ///< decreased point of time stamp
    Long64_t fTimeStampPrev = -1; ///< previous time stamp
    Long64_t fTimeStampLastSec = 0;
    Long64_t fTimeStampLastSecError = 0;
    Int_t fMinuiteBin = 0;

    void ResetFired();
    Short_t** fFiredDCh;
    Short_t* fFiredDetector;

    UShort_t fADCThreshold = 0;
    double fEnergyRange1 = -1;
    double fEnergyRange2 = -1;

    Long64_t bTimeStamp = -1; ///< branch value for time stamp
    Long64_t bTimeStampDist = -1; ///< branch value for distance do previous time stamp
    UShort_t bNumChannels = 0;
    double bdE = -1;
    double bESum = -1;
    double bE3 = -1;

    const int fNumModules = NUMBER_OF_MODULES;
    const int fNumChannels = NUMBER_OF_CHANNELS;

    TClonesArray *fChannelArray = nullptr;
    bool  fExcludedES1Coincidence = false;
    bool  fdES1CoincidenceMode = false;
    bool  fdES1S3CoincidenceMode = false;
    bool  fGoodCoincidenceEvent = false;
    Int_t fCoincidenceCount[5];
    Int_t fCoincidenceTSRange = 0;
    Int_t fCoincidenceMultRange1 = -1;
    Int_t fCoincidenceMultRange2 = -1;
    Int_t fCountChannels = 0;
    Long64_t fCountAllLines = 0;
    Long64_t fCountAllChannels = 0;
    Long64_t fCountTSError = 0;
    bool fSkipTSError = false;
    bool fStopAtTSError = false;
    bool fExitAnalysis = false;
    bool fExitRoot = false;

    Short_t fChosenDet = -1;
    Short_t fChosenDCh = -1;
    vector<int> fS1ChosenS1Strips;

    vector<int> fdEArrayIdx;
    vector<int> fS1ArrayIdx;
    vector<int> fS3ArrayIdx;
    double fdEADC = 0;
    double fS1ADC = 0;
    double fS3ADC = 0;

    bool fEnergyConversionIsSet = false;
    TF1* fFxEnergyConversion[NUMBER_OF_MODULES][NUMBER_OF_CHANNELS];

    int kSameEvent = 0;
    int kNextEvent = 1;
    int kTSError = 2;

    int kSameSecond = 3;
    int fNextSecond = 4;
    int kTimeError = 5;

    int fCountTriggerPerSec = 0;
    int fCountTriggerPerSecError = 0;
    int fCountEventsPerSec = 0;
    int fCountEventsPerSecError = 0;

  // drawing for analysis
  private:
    //TH2D* fHistdEVSS1 = nullptr; ///< dE vs dE + S1-energy
    //TH1D* fHistTritonE = nullptr; ///< gated triton energy
    //TH1D* fHistCrEx = nullptr; ///< 48Cr excitation energy
    //TH2D* fHistS3VSEx = nullptr; ///< coincidence S3 vs excitation energy

  public:
    static Analysis* GetAnalysis();

  private:
    static Analysis* fInstance;
    std::ofstream coutx;
};

#endif
