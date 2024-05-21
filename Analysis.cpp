#include "Analysis.h"

// Replace define value to "coutx" to omit printing
#define coutd cout<<"+\033[0;36m"<<__LINE__<<" "<<__FILE__<<" #\033[0m "
#define coutn cout<<"\033[0;36mNT\033[0m "
#define couti cout<<"\033[0;32m==\033[0m "
#define coutw cout<<"\033[0;33mWR\033[0m "
#define coute cout<<"\033[0;31mER\033[0m "
//#define coutt cout<<"\033[0;33mTSERROR\033[0m "
#define coutt coutx

//#define DEBUG_DATA_LINE
//#define DEBUG_DATA_LINE_CONDITION
//#define DEBUG_EVENT_LINE_CONDITION
//#define DEBUG_EXIT_ANALYSIS

//////////////////////////////////////////////////////////////////////////////
Analysis *getAna() { return Analysis::GetAnalysis(); }
void makeTCut(TString fileName="") { Analysis::MakeTritonCutFile(fileName); }
void callTCut(TString fileName="") { Analysis::CallTritonCutFile(fileName); }
//////////////////////////////////////////////////////////////////////////////

Analysis* Analysis::fInstance = nullptr;
Analysis* Analysis::GetAnalysis() {
  if (fInstance !=nullptr)
    return fInstance;
  return new Analysis();
} 

Analysis::Analysis()
{
  fInstance = this;
  coutx.open("out/dummy_stream");
  InitializeAnalysis();
}

void Analysis::InitializeAnalysis()
{
  std::ifstream mapFile(fMapFileName);
  if (mapFile.fail()) {
    coute << "Cannot open mapping file: " << fMapFileName << endl;
    return;
  }

  fMapDetectorType = new int*[fNumModules];
  fMapDetectorChannel = new int*[fNumModules];
  fMapDetectorReplaced = new bool*[fNumModules];
  fMapDetectorGroup = new int*[fNumModules];
  fMapDetectorRFMod = new int*[fNumModules];
  fMapDetectorRFMCh = new int*[fNumModules];
  for (int midx=0; midx<fNumModules; ++midx) {
    fMapDetectorType[midx] = new int[fNumChannels];
    fMapDetectorChannel[midx] = new int[fNumChannels];
    fMapDetectorReplaced[midx] = new bool[fNumChannels];
    fMapDetectorGroup[midx] = new int[fNumChannels];
    fMapDetectorRFMod[midx] = new int[fNumChannels];
    fMapDetectorRFMCh[midx] = new int[fNumChannels];
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
      fMapDetectorType[midx][iChannel] = kDummyDetector;
      fMapDetectorChannel[midx][iChannel] = -1;
      fMapDetectorReplaced[midx][iChannel] = false;
      fMapDetectorGroup[midx][iChannel]  = 0;
      fMapDetectorRFMod[midx][iChannel] = -1;
      fMapDetectorRFMCh[midx][iChannel] = -1;
    }
  }
  for (int iModule=0; iModule<20; ++iModule)
    fMapFEToModuleIndex[iModule] = -1;


  fMapDetectorToGlobalID = new int*[10];
  for (int det=0; det<10; ++det) {
    fMapDetectorToGlobalID[det] = new int[32];
    for (int dch=0; dch<10; ++dch) {
      fMapDetectorToGlobalID[det][dch] = -1;
    }
  }

  fMapGlobalIDToModuleIndex = new int[fNumCh];
  fMapGlobalIDToMCh = new int[fNumCh];
  for (int gid=0; gid<fNumCh; ++gid) {
    fMapGlobalIDToModuleIndex[gid] = -1;
    fMapGlobalIDToMCh[gid] = -1;
  }

  fDetectorName[kDummyDetector] = "X";
  fDetectorName[kdEDetector   ] = "dEDetector";
  fDetectorName[kS1Junction   ] = "S1Junction";
  fDetectorName[kS3Junction   ] = "S3Junction";
  fDetectorName[kS3Ohmic      ] = "S3Ohmic";
  fDetectorName[kScintillator ] = "Scintillator";
  fDetectorName[kFaradayCup   ] = "FaradayCup";

  TString name;
  string buffer;
  bool replaced;
  UShort_t module0, channelID0, FENumber, mch, dch, dES1Group, midx;

  getline(mapFile, buffer); // skip header
  while (mapFile >> FENumber >> midx >> mch >> name >> dch >> dES1Group >> replaced)
  {
    int det = kDummyDetector;
    for (det=0; det<fNumDetectors; ++det)
      if (fDetectorName[det]==name)
        break;
    if (det==kDummyDetector) {
      coute << "Dummy detector module! " << endl;
      coute << "  " << FENumber << " " << midx << " " << mch << " " << name << " " << dch << endl;
      continue;
    }
    fMapDetectorType[midx][mch] = det;
    fMapDetectorChannel[midx][mch] = dch;
    fMapDetectorReplaced[midx][mch] = replaced;
    fMapDetectorGroup[midx][mch] = dES1Group;
    if (replaced) {
      mapFile >> module0 >> channelID0;
      fMapDetectorRFMod[midx][mch] = module0;
      fMapDetectorRFMCh[midx][mch] = channelID0;
    }
    if (fMapFEToModuleIndex[FENumber]==-1)
      fMapFEToModuleIndex[FENumber] = midx;
    auto gid = GetGlobalID(midx, mch);
    fMapDetectorToGlobalID[det][dch] = gid;
    fMapGlobalIDToModuleIndex[gid] = midx;
    fMapGlobalIDToMCh[gid] = mch;
  }
  mapFile.close();

  for (int midx=0; midx<fNumModules; ++midx)
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      fFxEnergyConversion[midx][iChannel] = nullptr;

  fCountEvents = 0;
  fChannelArray = nullptr;
  fCountChannels = 0;
  fCountAllLines = 0;
  fCountAllChannels = 0;
  fCountTSError = 0;
  fCoincidenceCount[0] = 0;
  fCoincidenceCount[1] = 0;
  fCoincidenceCount[2] = 0;
  fCoincidenceCount[3] = 0;
  fCoincidenceCount[4] = 0;

  if (fSkipTSError )
  {
    coutw << "The program will ignore time-stamp error!" << endl;
    coutw << "The program will ignore time-stamp error!" << endl;
    coutw << "The program will ignore time-stamp error!" << endl;
    coutw << "The program will ignore time-stamp error!" << endl;
  }
}

void Analysis::RunConversion(int runNo, TString pathIn)
{
  if (fEnergyConversionIsSet==false&&fShowEnergyConversion==true) {
    coutw << "Energy conversion is not set! Cannot proceed energy conversion." << endl;
    fShowEnergyConversion = false;
  }

  fRunNo = runNo; 
  fRunName = Form("[RUN%03d]",fRunNo);
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/analysis/input/") : pathIn);

  ConfigureDateTime();
  SetConversionFile();
  InitializeDrawing();
  ReadDataFile();
  EndOfConversion();
}

void Analysis::AddAlphaCalibrationFile(TString name)
{
  couti << "Set energy calibration functions from " << name << endl;

  TFile* file = new TFile(name,"read");
  //file -> ls();

  for (int midx=0; midx<fNumModules; ++midx) {
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
      auto f1 = (TF1*) file -> Get(Form("AToE_%d_%d",midx,iChannel));
      if (f1!=nullptr)
        fFxEnergyConversion[midx][iChannel] = f1;
    }
  }

  for (int midx=0; midx<fNumModules; ++midx) {
    bool missing = false;
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      if (fFxEnergyConversion[midx][iChannel]==nullptr) {
        missing = true;
        break;
      }
    if (!missing) break;
    coute << "Missing ECal: ";
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      if (fFxEnergyConversion[midx][iChannel]==nullptr)
        cout << "(" << midx << "," << iChannel << ") ";
    cout << endl;
  }

  fEnergyConversionIsSet = true;
  //file -> Close();
}

void Analysis::ReadSummaryFile(TString fileName)
{
  fFileSummary = new TFile(fileName,"read");
  if (fFileSummary -> IsZombie()) {
    coute << "File " << fileName << " is zombie!" << endl;
    fFileSummary = nullptr;
    return;
  }
  fTreeSummary = (TTree*) fFileSummary -> Get("event");
  if (fTreeSummary==nullptr) {
    coute << "Tree from " << fileName << " is zombie!" << endl;
    return;
  }
  fTreeSummary -> SetBranchAddress("ts",&bTimeStamp);
  fTreeSummary -> SetBranchAddress("nch",&bNumChannels);
  fTreeSummary -> SetBranchAddress("de",&bdE);
  fTreeSummary -> SetBranchAddress("ee",&bESum);
  fTreeSummary -> SetBranchAddress("s3",&bE3);
  fTreeSummary -> SetBranchAddress("ch",&fChannelArray);
  fNumEventsSummary = fTreeSummary -> GetEntries();

  //fFileSummary -> ls();
  //fRunNo = ((TParameter<int>*)fFileSummary->Get("runNo"))->GetVal();
  fRunNo = TString(((TNamed*) fFileSummary->Get("runNo"))->GetTitle()).Atoi();
  couti << fileName << "(" << fRunNo << ") containing " << fNumEventsSummary << " events" << endl;
}

Double_t Analysis::FxTwoAlpha(Double_t *xy, Double_t *par)
{
  // fAlphaEnergy1 = 5.486 MeV (85.2 %)
  // fAlphaEnergy2 = 5.443 MeV (12.8 %)

  double x = xy[0];

  double ADCOffset = par[0];
  double energyResolution = par[1];

  double mean1 = par[2]; // ADC mean of 5.486 MeV peak
  double meanP1 = mean1 - ADCOffset; // pure ADC without ADC-offset
  double sigma1 = energyResolution * meanP1;
  double amplitude1 = par[3];

  double ADCEnergyRatio = meanP1 / fAlphaEnergy1;

  double mean2 = fAlphaEnergy2 * ADCEnergyRatio + ADCOffset;
  double meanP2 = mean2 - ADCOffset; // pure ADC without ADC-offset
  double sigma2 = energyResolution * meanP2;
  double amplitude2 = amplitude1 * sigma1 / fAEBR / sigma2;

  double value1 = amplitude1*exp(-0.5*((x-mean1)*(x-mean1)/sigma1/sigma1));
  double value2 = amplitude2*exp(-0.5*((x-mean2)*(x-mean2)/sigma2/sigma2));
  double value  = value1 + value2;

  return value;
}

void Analysis::Convert2APParameters(double* par, double &mean1, double &sigma1, double &amplitude1, double &mean2, double &sigma2, double &amplitude2)
{
  double ADCOffset = par[0];
  double energyResolution = par[1];

  mean1 = par[2]; // ADC mean of 5.486 MeV peak
  double meanP1 = mean1 - ADCOffset; // pure ADC without ADC-offset
  sigma1 = energyResolution * meanP1;
  amplitude1 = par[3];

  double ADCEnergyRatio = meanP1 / fAlphaEnergy1;

  mean2 = fAlphaEnergy2 * ADCEnergyRatio + ADCOffset;
  double meanP2 = mean2 - ADCOffset; // pure ADC without ADC-offset
  sigma2 = energyResolution * meanP2;
  amplitude2 = amplitude1 * sigma1 / 6.65625 / sigma2; // 85.2 / 12.8 = 6.65625; 
}

TF1* Analysis::CreateFxTwoAlpha(TString name, double min, double max)
{
  auto f1 = new TF1(name,this,&Analysis::FxTwoAlpha,min,max,4,"Analysis","FxTwoAlpha");
  return f1;
}

void Analysis::SetAlphaTestFile(TString fileName)
{
  if (fFileAlpha!=nullptr)
    return;
  if (fileName.IsNull())
    fileName = Form("RUN%03d.alpha.root",fRunNo);
  if (fFileNameAlpha.IsNull())
    fFileNameAlpha = fPathToOutput+fileName;
  cout << endl;
  couti << "Creating output file " << fFileNameAlpha << endl;
  fFileAlpha = new TFile(fFileNameAlpha,"recreate");
  WriteRunParameters(fFileAlpha,1);
}

void Analysis::AnalyzeAlphaTestModule(int midx, bool drawAnalysis, TString fileName)
{
  SetAlphaTestFile(fileName);

  TCanvas* cvs = nullptr;
  if (drawAnalysis) {
    cvs = new TCanvas(Form("cvs%d",midx),GetDetectorTitle(midx) + Form("  M(%d)",midx),1600,900);
    cvs -> Divide(4,4);
  }
  bool moduleDataExist = false;
  for (auto mch=0; mch<fNumChannels; ++mch)
  {
    TVirtualPad* pad = nullptr;
    if (drawAnalysis)
      pad = cvs->cd(mch+1);

    bool channelDataExist = AnalyzeAlphaTest(midx, mch, drawAnalysis, pad);
    moduleDataExist = (moduleDataExist || channelDataExist);
  }

  if (drawAnalysis)
  {
    if (!moduleDataExist) {
      coute << "Module " << midx << " is empty!" << endl;
        delete cvs;
    }
    else if (fFileAlpha!=nullptr) {
      fFileAlpha -> cd();
      cvs -> Write();
    }
  }

  cout << endl;
  if (fFileAlpha!=nullptr)
    couti << fFileAlpha -> GetName() << endl;
}

bool Analysis::AnalyzeAlphaTest(int midx, int mch, bool drawAnalysis, TVirtualPad* cvs)
{
  if (fTreeSummary==nullptr) {
    coute << "Summary tree is nullptr! Run ReadSummaryFile before this method!";
    return false;
  }
  
  bool isSingleDrawing = false;
  if (cvs==nullptr)
    isSingleDrawing = true;

  auto detectorTitle = GetDetectorTitle(midx,mch,1) + Form(" M(%d,%d)",midx,mch);
  couti << detectorTitle << endl;

  auto gid = GetGlobalID(midx, mch);
  if (fFitAlpha==nullptr)
  {
    fFitAlpha = CreateFxTwoAlpha("fxTwoAlpha",0,fMaxADC);
    fFitAlpha1 = new TF1("fxAlpha1","gaus(0)",0,fMaxADC);
    fFitAlpha2 = new TF1("fxAlpha2","gaus(0)",0,fMaxADC);
    fFitAlpha1 -> SetLineColor(kBlue);
    fFitAlpha2 -> SetLineColor(kGreen+1);
  }
  double mean1, sigma1, amplitude1, mean2, sigma2, amplitude2;

  TString nameHist = Form("histE%d",gid);
  TString titleHist = detectorTitle + ";ADC;count";
  auto histE = new TH1D(nameHist,titleHist,fNumADC,0,fMaxADC);
  fTreeSummary -> Project(nameHist, "ch.adc", Form("ch.adc>1000&&ch.gid==%d",gid));
  if (histE->GetEntries()==0) {
    coute << "Data is empty!" << endl;
    if (fFileAlpha!=nullptr) {
      auto fxADCToEnergyDummy = new TF1(Form("AToE_%d_%d",midx,mch),"0.005*(x-0.0000)",0,10000);
      fxADCToEnergyDummy -> Write();
    }
    return false;
  }

  double mean = histE -> GetMean();
  double sig = histE -> GetStdDev();

  double bin = histE->GetMaximumBin();
  double yMax = histE->GetBinContent(bin);
  double xPeak = histE->GetBinCenter(bin);

  auto x1 = xPeak - 80;
  auto x2 = xPeak + 80;
  fFitAlpha -> SetParameters(10,0.005,xPeak,yMax);
  fFitAlpha -> SetParLimits(0,0,100);
  histE -> Fit(fFitAlpha,"Q0N","",x1,x2);
  auto parameters = fFitAlpha -> GetParameters();
  Convert2APParameters(parameters, mean1, sigma1, amplitude1, mean2, sigma2, amplitude2);
  double ADCOffset = parameters[0];
  if (ADCOffset<0.001)
    ADCOffset = 0;
  double ADCEnergyRatio = (mean1 - ADCOffset) / fAlphaEnergy1;
  double energyADCRatio = fAlphaEnergy1 / (mean1-ADCOffset);
  double energyResolution = parameters[1] * 2.354;
  double energyResolutionKeV = energyResolution * 1000*fAlphaEnergy1;
  cout << "   Energy res. : " << 100*energyResolution << " % (" << Form("%.2f KeV)",energyResolutionKeV) << endl;
  cout << "   Energy/ADC. : " << energyADCRatio << endl;
  cout << "   ADC offset  : " << ADCOffset << endl;
  cout << Form("   AMS(%.3f)",fAlphaEnergy1) <<" : " << amplitude1 << ", " << mean1 << ", " << sigma1 << endl;
  cout << Form("   AMS(%.3f)",fAlphaEnergy2) <<" : " << amplitude2 << ", " << mean2 << ", " << sigma2 << endl;
  fFitAlpha1 -> SetParameters(amplitude1, mean1, sigma1);
  fFitAlpha2 -> SetParameters(amplitude2, mean2, sigma2);

  if (fFileAlpha!=nullptr) {
    auto fxADCToEnergy = new TF1(Form("AToE_%d_%d",midx,mch),Form("%f*(x-%.4f)",energyADCRatio,ADCOffset),0,10000);
    fxADCToEnergy -> Write();
  }

  if (drawAnalysis)
  {
    TString nameCvs = Form("cvsE%d",gid);
    if (cvs==nullptr)
      cvs = new TCanvas(nameCvs,"",800,450);

    if (isSingleDrawing) SetAttribute(histE,cvs,1);
    else SetAttribute(histE,cvs,16);

    histE -> GetXaxis() -> SetRangeUser(mean2-8*sigma2,mean1+8*sigma1);
    histE -> Draw();
    fFitAlpha -> DrawCopy("samel");
    fFitAlpha1 -> DrawCopy("samel");
    fFitAlpha2 -> DrawCopy("samel");

    const char* percent = "%";
    if (1)
    {
      auto legend = new TLegend(0.18,0.30,0.47,0.85);
      legend -> SetBorderSize(0);
      legend -> SetFillStyle(0);
      legend -> SetTextSize(0.07);
      legend -> AddEntry("",Form("ER = %0.2f %s", 100*energyResolution, percent),"");
      legend -> AddEntry("",Form("(%0.2f KeV)", energyResolutionKeV),"");
      if (ADCOffset>0.001) legend -> AddEntry("",Form("ADC offset = %.3f", ADCOffset),"");
      legend -> AddEntry(fFitAlpha1,Form("1) %0.2f", mean1),"l");
      legend -> AddEntry(fFitAlpha2,Form("2) %0.2f", mean2),"l");
      legend -> Draw();
    }
    else
    {
      auto legend = new TLegend(0.15,0.40,0.47,0.85);
      legend -> SetBorderSize(0);
      legend -> SetFillStyle(0);
      legend -> SetTextSize(0.05);
      legend -> AddEntry("",Form("Resolution = %0.2f %s", 100*energyResolution, percent),"");
      if (ADCOffset>0.001) legend -> AddEntry("",Form("ADC offset = %.3f", ADCOffset),"");
      legend -> AddEntry(fFitAlpha1,Form("M(%.3f) = %0.2f", fAlphaEnergy1,mean1),"l");
      legend -> AddEntry(fFitAlpha2,Form("M(%.3f) = %0.2f", fAlphaEnergy2,mean2),"l");
      legend -> Draw();
    }

    if (isSingleDrawing && fFileAlpha!=nullptr) {
      fFileAlpha -> cd();
      cvs -> Write();
    }
  }

  return true;
}

void Analysis::InitializeDrawing()
{
  if (fFileOut!=nullptr)
    fFileOut -> cd();

  if (fUpdateDrawingEveryNEvent>0)
  {
    fCvsOnline1 = new TCanvas("cvsOnline",fRunName+" online update canvas 1",1800,800);
    fCvsOnline1 -> SetMargin(0,0,0,0);
    fCvsOnline1 -> Divide(3,2,0,0);

    fCvsOnline1 -> cd(2) -> Divide(2,2,0,0);
    fPadTSDist1 = fCvsOnline1 -> cd(2) -> cd(1);
    fPadTSDist2 = fCvsOnline1 -> cd(2) -> cd(2);
    fPadTrgRate = fCvsOnline1 -> cd(2) -> cd(3);
    fPadEvtRate = fCvsOnline1 -> cd(2) -> cd(4);
    fPadEx      = fCvsOnline1 -> cd(3);
    fPadChCount = fCvsOnline1 -> cd(1);
    fPadADC     = fCvsOnline1 -> cd(5);
    fPadEVSCh   = fCvsOnline1 -> cd(4);
    fPaddEVSE   = fCvsOnline1 -> cd(6);
    fPadEVSCh   -> SetLogz();
    fPaddEVSE   -> SetLogz();
  }

  fHistTriggerRate      = new TH1D("histTrgRateG",";nch*trigger/s;count",fNumRate,0,fMaxRate);
  fHistTriggerRateError = new TH1D("histTrgRateE","Bad trigger rate;trigger/s;count",fNumRate,0,fMaxRate);
  fHistTriggerRateError -> SetLineStyle(2);
  //fHistTriggerRate      -> SetFillColor(kGray);
  fHistTriggerRate      -> SetLineColor(kBlack);
  fHistTriggerRateError -> SetLineColor(kBlack);

  fHistEventRate      = new TH1D("histEvtRateG",";event/s;count",fNumRate,0,fMaxRate);
  fHistEventRateError = new TH1D("histEvtRateE","Bad event rate;event/s;count",fNumRate,0,fMaxRate);
  fHistEventRateError -> SetLineStyle(2);
  fHistEventRate      -> SetFillColor(29);
  fHistEventRate      -> SetLineColor(kBlue);
  fHistEventRateError -> SetLineColor(kBlue);
  //fHistEventRate      -> SetStats(0);

  if (fPadTrgRate==fPadEvtRate) {
    //fHistTriggerRate -> SetTitle("Trigger*NCh (black) / Event (blue) rate;trigger/s;count");
    fHistEventRate -> SetTitle("Trigger*NCh (black) / Event (blue) rate;trigger/s;count");
    //fHistEventRate -> SetTitle(";Trigger*NCh/s (black) / Event/s (blue);count");
  }

  //fHistTSDist1 = new TH1D("histTSDist1","TS distance (+);TS-dist;count", 20,-1,20);
  fHistTSDist1 = new TH1D("histTSDist1","TS distance;TS-dist;count", 20,0,20);
  fHistTSDist2 = new TH1D("histTSDist2","TS distance (TS-Error);TS-dist;count", 50,-1000000,0);
  fHistTSDist1 -> SetStats(0);
  //fHistTSDist2 -> SetStats(0);

  fHistChCount = new TH1D("histChCount","Channel count;;", fNumCh,0,fNumCh);
  fHistChCount -> SetStats(0);
  fHistChCount -> SetFillColor(29);

  fHistADC = new TH1D("histADC","ADC (all channels);ADC",fNumADC,0,fMaxADC);
  fHistADC -> SetFillColor(29);
  fHistE = new TH1D("histEnergy",";energy (MeV)",fNumE,0,fMaxE);
  fHistE -> SetFillColor(29);

  fHistAVSCh = new TH2D("histAVSCh", "ADC vs module-ch;;ADC", fNumCh,0,fNumCh,fNumADC,0,fMaxADC);
  fHistEVSCh = new TH2D("histEVSCh", ";;energy (MeV)", fNumCh,0,fNumCh,fNumE,0,fMaxE);
  fHistdAVSA = new TH2D("histdAVSA", "dA vs dA + S1;dA + S1A;dA", fNumADC, 0, 2*fMaxADC, fNumADC, 0, fMaxADC);
  fHistdEVSE = new TH2D("histdEVSE", ";dE + S1E (MeV);dE (MeV)", fNumE,0,fMaxE+fMaxdE, fNumdE, 0, fMaxdE);
  fHistAVSCh -> SetStats(0);
  fHistEVSCh -> SetStats(0);
  fHistdAVSA -> SetStats(0);
  fHistdEVSE -> SetStats(0);

  fHistEx = new TH1D("histEx",";Cr excitation energy;count",fNumEx,0,fMaxEx);
  fHistEx -> SetFillColor(29);

  if (fUpdateDrawingEveryNEvent>0) {
    SetAttribute(fHistTriggerRate,fPadTrgRate,16);
    SetAttribute(fHistTriggerRateError,fPadTrgRate,16);
    SetAttribute(fHistEventRate,fPadEvtRate,16);
    SetAttribute(fHistEventRateError,fPadEvtRate,16);
    SetAttribute(fHistTSDist1,fPadTSDist1,16);
    SetAttribute(fHistTSDist2,fPadTSDist2,16);
    SetAttribute(fHistChCount,fPadChCount);
    SetAttribute(fHistEx,fPadEx);
    SetAttribute(fHistADC,fPadADC);
    SetAttribute(fHistE,fPadADC);
    SetAttribute(fHistAVSCh,fPadEVSCh,1,true);
    SetAttribute(fHistEVSCh,fPadEVSCh,1,true);
    SetAttribute(fHistdAVSA,fPaddEVSE,1,true);
    SetAttribute(fHistdEVSE,fPaddEVSE,1,true);
  }

  //auto TakeCareOfLabels = [this](TH1* hist)
  //{
  //  for (int midx=0; midx<fNumModules; ++midx) {
  //    for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
  //      if (fMapDetectorReplaced[midx][iChannel])
  //      {
  //        auto detectorName = fDetectorName[fMapDetectorType[midx][iChannel]];
  //        auto dch = fMapDetectorChannel[midx][iChannel];
  //        auto gid = GetGlobalID(midx,iChannel);
  //        TString title = Form("%s%d",detectorName.Data(),dch);
  //        hist -> GetXaxis() -> SetBinLabel(gid+1,title);
  //      }
  //    }
  //  }
  //};
  //TakeCareOfLabels(fHistChCount);

  UpdateCvsOnline(true);
}

void Analysis::SetAttribute(TH1* hist, TVirtualPad* pad, int npad, bool is2D)
{
  auto ax = hist -> GetXaxis();
  auto ay = hist -> GetYaxis();

  if (npad==2)
  {
    if (pad!=nullptr) {
      pad -> SetLeftMargin(0.15);
      pad -> SetRightMargin(0.10);
      pad -> SetBottomMargin(0.12);
      pad -> SetTopMargin(0.12);
    }
    gStyle -> SetTitleH(0.07);
    ax -> SetLabelSize(0.04);
    ay -> SetLabelSize(0.04);
    ax -> SetTitleSize(0.05);
    ay -> SetTitleSize(0.05);
    ax -> SetTitleOffset(1.0);
    ay -> SetTitleOffset(1.50);
  }
  else if (npad>=16)
  {
    if (pad!=nullptr) {
      pad -> SetLeftMargin(0.12);
      pad -> SetRightMargin(0.05);
      pad -> SetBottomMargin(0.12);
      pad -> SetTopMargin(0.12);
    }
    gStyle -> SetTitleH(0.08);
    ax -> SetLabelSize(0.06);
    ay -> SetLabelSize(0.06);
    ax -> SetTitleSize(0.07);
    ay -> SetTitleSize(0.07);
    ax -> SetTitleOffset(0.8);
    ay -> SetTitleOffset(0.82);
  }
  else
  {
    if (pad!=nullptr) {
      pad -> SetLeftMargin(0.11);
      pad -> SetRightMargin(0.05);
      if (is2D)
        pad -> SetRightMargin(1.2);
      pad -> SetBottomMargin(0.13);
      pad -> SetTopMargin(0.11);
    }
    gStyle -> SetTitleH(0.06);
    ax -> SetLabelSize(0.04);
    ay -> SetLabelSize(0.04);
    ax -> SetTitleSize(0.05);
    ay -> SetTitleSize(0.05);
    ax -> SetTitleOffset(1.1);
    ay -> SetTitleOffset(1.0);
  }
}

void Analysis::ConfigureDateTime()
{
  cout << endl;
  couti << Form("Looking for RUN%03d in %s",fRunNo,fPathToInput.Data()) << " ..." << endl;

  TList *listOfFiles = TSystemDirectory(fPathToInput,fPathToInput).GetListOfFiles();
  TIter next(listOfFiles);
  TSystemFile* fileObject;
  while ((fileObject=(TSystemFile*)next()))
  {
    int idx = 0;
    if (fileObject->IsDirectory()) continue;
    TString fileName = fileObject -> GetName();
    if (fileName.Index(Form("RUN%03d",fRunNo))!=0) continue;
    if (fileName.EndsWith(".dat")==false) continue;
    fDateTime = fileName(7,12);
    int fileNumber = TString(fileName(25,3)).Atoi();
    if (fileNumber>fFileNumberMax)
      fFileNumberMax = fileNumber;
    couti << fileName << endl;
  }
  if (fFileNumberMax>fFileNumberRange2)
    fFileNumberMax=fFileNumberRange2;
}

void Analysis::WriteRunParameters(TFile* file, int option)
{
  file -> cd();
  (new TNamed("runNo", Form("%d",fRunNo))) -> Write();
  if (option==0) {
    (new TNamed("Draw every n-event   :", Form("%d",fUpdateDrawingEveryNEvent))) -> Write();
    (new TNamed("Return if no file    :", Form("%d",fReturnIfNoFile))) -> Write();
    (new TNamed("Ignore file update   :", Form("%d",fIgnoreFileUpdate))) -> Write();
    (new TNamed("Ignore TS decrease   :", Form("%d",fSkipTSError))) -> Write();
    (new TNamed("Event count limit    :", Form("%d",fEventCountLimit))) -> Write();
    (new TNamed("File number range    :", Form("%d %d",fFileNumberRange1,fFileNumberRange2))) -> Write();
    (new TNamed("ADC Threshold        :", Form("%d",fADCThreshold))) -> Write();
    (new TNamed("Coincidence dTS cut  :", Form("%d",fCoincidenceTSRange))) -> Write();
    (new TNamed("Coincidence mult rng :", Form("%d %d",fCoincidenceMultRange1,fCoincidenceMultRange2))) -> Write();
    (new TNamed("dES1 Coincidence mode:", Form("%d",fdES1CoincidenceMode))) -> Write();
    (new TNamed("dE13 Coincidence mode:", Form("%d",fdES1S3CoincidenceMode))) -> Write();
    if (fTritonCutG!=nullptr) {
      auto tritonCutGCopy = (TCutG*) fTritonCutG -> Clone("tritonCutGCopy");
      tritonCutGCopy -> Write();
    }
  }
}

void Analysis::SetConversionFile(TString fileName)
{
  if (fFileOut!=nullptr)
    return;

  if (fileName.IsNull())
    fFileNameOut = fPathToOutput+Form("RUN%03d.summary.root",fRunNo);
  else
    fFileNameOut = fPathToOutput+fileName;

  cout << endl;
  couti << "Creating output file " << fFileNameOut << endl;
  fFileOut = new TFile(fFileNameOut,"recreate");
  fTreeOut = new TTree("event","");
  fChannelArray = new TClonesArray("ChannelData",10);
  fTreeOut -> Branch("ts",&bTimeStamp);
  fTreeOut -> Branch("nch",&bNumChannels);
  fTreeOut -> Branch("de",&bdE);
  fTreeOut -> Branch("ee",&bESum);
  fTreeOut -> Branch("s3",&bE3);
  fTreeOut -> Branch("ch",&fChannelArray);

  WriteRunParameters(fFileOut,0);
}

void Analysis::ReadDataFile()
{
  fDataFileNumber = fFileNumberRange1; // Number of Files read
  UShort_t FENumber, mch, adc, tsGroup;
  Long64_t tsLocal, timeStamp;
  Int_t numData = 0;
  char buffer[256];
  //int countOpenFileTrials = 0;

  ChannelData* channelData = NULL;
  fTimeStampPrev = -1;

  while (true)
  {
    //if (CheckOpenFileStatus1()==false) break;
    if (CheckOpenFileStatus2()==false) break;

    fChannelArray -> Clear("C");
    Long64_t countEventsSingleFile = 0;
    fCountChannels = 0;
    Long64_t countLine = 0;
    fdEArrayIdx.clear();
    fS1ArrayIdx.clear();
    fS3ArrayIdx.clear();
    fdEADC = 0.;
    fS1ADC = 0.;
    fS3ADC = 0.;

    while (fDataFile >> buffer)
    {
      fLastDataPos = fDataFile.tellg();

      countLine++;
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &FENumber, &mch, &adc, &tsLocal, &tsGroup);
      if (numData != 5) {
        coute << Form("Number of data in line is not 5 (%d)!",numData) << endl;
        coute << FENumber<<" "<<mch<<" "<<adc<<" "<<tsLocal<<" "<<tsGroup << endl;
        continue;
      }
#ifdef DEBUG_DATA_LINE
      coutd << fCountAskContinueRun << endl;
      if (fCountAskContinueRun==1) coutd << FENumber<<" "<<mch<<" "<<adc<<" "<<tsLocal<<" "<<tsGroup << endl;
#endif

      timeStamp = tsLocal;
      if (tsGroup>0) timeStamp += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      if (bTimeStamp>=0) {
        int dts = int(timeStamp)- int(fTimeStampPrevTrue);
        if (dts>=0) fHistTSDist1 -> Fill(dts);
        else { fHistTSDist2 -> Fill(dts); }
      }

      int eventStatus = kNextEvent;
      if (!fSkipTSError)
      {
        if (timeStamp-fTimeStampPrev<=fCoincidenceTSRange) eventStatus = kSameEvent;
        else if (timeStamp<fTimeStampPrev) eventStatus = kNextEvent;
        else if (timeStamp>fTimeStampPrev) eventStatus = kNextEvent;
      }
      else
      {
        if (timeStamp-fTimeStampPrev<=fCoincidenceTSRange) eventStatus = kSameEvent;
        else if (timeStamp<fTimeStampPrev) eventStatus = kNextEvent;
        else if (timeStamp>fTimeStampPrev) eventStatus = kTSError;
      }

      int timeStatus = kSameSecond;
      double dSecond = (timeStamp-fTimeStampLastSec)*fSecondPerTS;
      if (dSecond<0) timeStatus = kTimeError;
      else if (dSecond>1) timeStatus = fNextSecond;
      else timeStatus = kSameSecond;

      if (timeStatus==kSameSecond) {
        fCountTriggerPerSec++;
      }
      else if (timeStatus==kTimeError) {
        fCountTriggerPerSecError++;
      }
      else if (timeStatus==fNextSecond)
      {
        fHistTriggerRate -> Fill(fCountTriggerPerSec);
        fHistTriggerRateError -> Fill(fCountTriggerPerSecError);
        fCountTriggerPerSecError = 0;
        fCountTriggerPerSec = 0;
        fCountTriggerPerSec++;
        fTimeStampLastSec = timeStamp;
        fTimeStampLastSecError = timeStamp;
      }

      if (CheckDataLineCondition(adc,eventStatus,timeStamp)==false)
        continue;
      if (fExitAnalysis)
        break;

      if (eventStatus==kSameEvent)
      {
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);
      }
      else if (eventStatus==kNextEvent)
      {
        bool eventIsFilled = true;
        if (countEventsSingleFile>0) {
          eventIsFilled = FillDataTree();
          if (UpdateDrawing()==false)
            break;
        }
        if (fExitAnalysis)
          break;

        if (eventIsFilled)
        {
          fCountEvents++;
          countEventsSingleFile++;

          if (timeStatus==kSameSecond) {
            fCountEventsPerSec++;
          }
          else if (timeStatus==kTimeError) {
            fCountEventsPerSecError++;
          }
          else if (timeStatus==fNextSecond)
          {
            fHistEventRate -> Fill(fCountEventsPerSec);
            fHistEventRateError -> Fill(fCountEventsPerSecError);
            fCountEventsPerSecError = 0;
            fCountEventsPerSec = 0;
            fCountEventsPerSec++;
          }
        }

        fCountChannels = 0;
        fChannelArray -> Clear("C");
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);

        fTimeStampPrev = timeStamp;
        fTimeStampPrevTrue = timeStamp;
      }

      Short_t midx = GetModuleIndex(FENumber);
      int gid = GetGlobalID(midx, mch);
      double energy = 0;
      if (fEnergyConversionIsSet)
          energy = GetCalibratedEnergy(midx,mch,adc);
      auto det = fMapDetectorType[midx][mch];
      auto dch = fMapDetectorChannel[midx][mch];
      channelData -> SetData(midx,mch,det,dch,gid,adc,energy,tsLocal,tsGroup,timeStamp);
      fCountAllChannels++;

      if (fdES1CoincidenceMode || fdES1S3CoincidenceMode)
      {
        if (fMapDetectorType[midx][mch]==kS3Junction) {
          fS3ArrayIdx.push_back(fCountChannels-1);
          fS3ADC = adc;
        }
        if (fMapDetectorType[midx][mch]==kS1Junction) {
          fS1ArrayIdx.push_back(fCountChannels-1);
          fS1ADC = adc;
        }
        if (fMapDetectorType[midx][mch]==kdEDetector) {
          fdEArrayIdx.push_back(fCountChannels-1);
          fS1ADC = adc;
        }
      }

      if (fUpdateDrawingEveryNEvent>=0)
      {
        fHistChCount -> Fill(gid);
        fHistADC -> Fill(adc);
        fHistAVSCh -> Fill(gid, adc);
        if (fEnergyConversionIsSet) {
          fHistEVSCh -> Fill(gid, energy);
          fHistE -> Fill(energy);
        }
      }

      if (fEventCountLimit>0 && fCountEvents>=fEventCountLimit) {
        couti << "Event count limit at " << fCountEvents << "!" << endl;
#ifdef DEBUG_EXIT_ANALYSIS
        coutn << "Exit ana: " << "Event count limit at " << fCountEvents << "!" << endl;
#endif
        fExitAnalysis = true;
      }
    }

    if (!fExitAnalysis && fCountEvents>0) 
      FillDataTree();

    if (fExitAnalysis)
      break;

    couti << "End of file!" << endl;
    couti << "Number of events from current file: " << countEventsSingleFile << endl;

    PrintConversionSummary();

  }
  if (fDataFile.is_open())
    fDataFile.close();

  UpdateCvsOnline();
}

bool Analysis::CheckOpenFileStatus2()
{
  if (!fInputFileName.IsNull() && fDataFile.is_open())
  {
    fLastDataFileSize = GetFileSize(fInputFileName);
    if (!fIgnoreFileUpdate)
       if (AskContinueRun("Waiting for possible update of the data file")==false) // after waiting here
          return false;

    Long64_t currentDataFileSize = GetFileSize(fInputFileName);
    if (currentDataFileSize>fLastDataFileSize) {
      couti << fInputFileName << " size increased by " << currentDataFileSize - fLastDataFileSize << endl;
      fDataFile.clear();
      fDataFile.seekg(fLastDataPos);
      return true;
    }
    else {
      couti << "No update! Closing file" << endl;
      fDataFile.close(); // and continue;
    }
  }

  {
    if (fFileNumberMax >=0 && fDataFileNumber>fFileNumberMax) {
      couti << "File number " << fDataFileNumber << " exceeded maximum " << fFileNumberMax << endl;
      return false;
    }
    fInputFileName = Form("%s/RUN%03d_%s_list_%03d.dat",
        fPathToInput.Data(),fRunNo,fDateTime.Data(),fDataFileNumber);
    fDataFileNumber++;

    cout << endl;
    couti << "Reading " << fInputFileName << endl;

    bool waitedForFileUpdate = false;
    bool firstLoop = true;

    while (true)
    {
      if (firstLoop)
        firstLoop = false;
      else {
        if (AskContinueRun("Wait for the new file?")==false)
          return false;
      }

      if (fExitAnalysis)
        return false;

      fDataFile.open(fInputFileName);
      if (fDataFile.fail())
      {
        fDataFile.close();
        if (waitedForFileUpdate) {
          coute << "Failed to open file!" << endl;
          return false;
        }
        else {
          coutw << "There is no file!" << endl;
          if (fReturnIfNoFile)
            return false;
          waitedForFileUpdate = true;
          continue;
        }
      }

      //couti << "Good!" << endl;
      break;
    }
  }

  return true;
}

bool Analysis::CheckOpenFileStatus1()
{
  if (fFileNumberMax >=0 && fDataFileNumber>fFileNumberMax) {
    couti << "Number of inputs " << fDataFileNumber << " exceeded maximum file number " << fFileNumberMax << endl;
    return false;
  }
  fInputFileName = Form("%s/RUN%03d_%s_list_%03d.dat",
      fPathToInput.Data(),fRunNo,fDateTime.Data(),fDataFileNumber);
  fDataFileNumber++;

  cout << endl;
  couti << "Reading " << fInputFileName << endl;

  int countOpenFileTrials = 0;
  while (true)
  {
    fDataFile.open(fInputFileName);
    if (fDataFile.fail())
    {
      fDataFile.close();
      if (countOpenFileTrials>10) {
        coute << "Failed to open file!" << endl;
        break;
      }
      else {
        countOpenFileTrials++;
        coutw << "There is no file!" << endl;
        if (fReturnIfNoFile)
          break;
        couti << "waiting(" << countOpenFileTrials << ") 3s ..." << endl;
        sleep(3);
        continue;
      }
    }

    //couti << "Good!" << endl;

    fDataFileSize = fDataFile.seekg(0, ios::end).tellg(); // obtain the size of file
    fDataFile.seekg(0, ios::beg); // rewind

    if (fDataFileSize>fDataFileSizeReadMax || (countOpenFileTrials!=0 && fDataFileSize==fDataFileSizeOld))
    {
      // let's read file!
      countOpenFileTrials = 0;
      break;
    }
    else if(countOpenFileTrials==0){ // first try
      countOpenFileTrials++;
      fDataFile.close();
      fDataFileSizeOld = fDataFileSize;
      if (fFirstFileOpened) {
        couti << "waiting(" << countOpenFileTrials << ") 3s ..." << endl;
        sleep(3);
      }
      fFirstFileOpened = true;
    }
    else if(fDataFileSize > fDataFileSizeOld){ // writing the file is still continued ...
      if (fIgnoreFileUpdate)
        break;
      fDataFile.close();
      fDataFileSizeOld = fDataFileSize;
      couti << "waiting(" << countOpenFileTrials << ") 60s ..." << endl;
      sleep(60);
    }
  }

  if (countOpenFileTrials != 0) {
    coutw << fInputFileName << " is not found! exit." << endl;
    return false;
  }

  return true;
}

bool Analysis::CheckDataLineCondition(double adc, int eventStatus, Long64_t timeStamp)
{
  fCountAllLines++;
  if (adc < fADCThreshold) {
#ifdef DEBUG_DATA_LINE_CONDITION
    if (fCountAskContinueRun==1) coutd << "Bad data: " << adc << " < " << fADCThreshold << endl;
#endif
    return false;
  }

  if (eventStatus==kTSError) // time-stamp decreased
  {
    if (timeStamp<fTimeStampPrevTrue) {
      coutt << "Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
      fTimeStampPrevTrue = timeStamp;
      ++fCountTSError;
    }
    if (fStopAtTSError) {
      coute << "Exit due to Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
#ifdef DEBUG_EXIT_ANALYSIS
        coutn << "Exit ana: " << "Exit due to Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
#endif
      fExitAnalysis = true;
    }
#ifdef DEBUG_DATA_LINE_CONDITION
    if (fCountAskContinueRun==1) coutd << "Bad data: " << "TS-Error!" << endl;
#endif
    return false;
  }
#ifdef DEBUG_DATA_LINE_CONDITION
  if (fCountAskContinueRun==1) coutd << "Good data" << endl;
#endif

  return true;
}

bool Analysis::CheckEventCondition(double de, double ee)
{
  if (!fGoodCoincidenceEvent)
    return false;
  fGoodCoincidenceEvent = false;

  if (fTritonCutGIsSet)
  {
    if (fTritonCutG->IsInside(ee,de))
      return true;
    else {
#ifdef DEBUG_EVENT_LINE_CONDITION
      coutd << "Bad event: (" << ee << ", " << de << ") is outside of triton cut" << endl;
#endif
      return false;
    }
  }

  return true;
}

bool Analysis::FillDataTree()
{
  if      (fCountChannels==0) fCoincidenceCount[0]++;
  else if (fCountChannels==1) fCoincidenceCount[1]++;
  else if (fCountChannels==2) fCoincidenceCount[2]++;
  else if (fCountChannels==3) fCoincidenceCount[3]++;
  else if (fCountChannels>=4) fCoincidenceCount[4]++;

  //coutd << fdES1CoincidenceMode << " " << fdES1S3CoincidenceMode << " " << fCoincidenceMultRange1 << endl;
  if (fdES1CoincidenceMode) {
    if (fdEArrayIdx.size()==1 && fS1ArrayIdx.size()==1)
      fGoodCoincidenceEvent = true;
    else {
#ifdef DEBUG_EVENT_LINE_CONDITION
      coutd << "Bad event: " << fdEArrayIdx.size() << " " << fS1ArrayIdx.size() << endl;
#endif
      fGoodCoincidenceEvent = false;
    }
  }
  else if (fdES1S3CoincidenceMode) {
    if (fdEArrayIdx.size()==1 && fS1ArrayIdx.size()==1 && fS3ArrayIdx.size()==1 )
      fGoodCoincidenceEvent = true;
    else {
#ifdef DEBUG_EVENT_LINE_CONDITION
      coutd << "Bad event: " << fdEArrayIdx.size() << " " << fS1ArrayIdx.size() << " " << fS3ArrayIdx.size() << endl;
#endif
      fGoodCoincidenceEvent = false;
    }
  }
  else if (fCoincidenceMultRange1>0) {
    if (fCountChannels>=fCoincidenceMultRange1 && fCountChannels<=fCoincidenceMultRange2)
      fGoodCoincidenceEvent = true;
    else {
#ifdef DEBUG_EVENT_LINE_CONDITION
      coutd << "Bad event: mult " << fCountChannels << " not in range " << fCoincidenceMultRange1 << " - " << fCoincidenceMultRange2 << endl;
#endif
      fGoodCoincidenceEvent = false;
    }
  }
  else
    fGoodCoincidenceEvent = true;

  bdE = -1;
  bESum = -1;
  bE3 = -1;
  if ((fdES1CoincidenceMode || fdES1CoincidenceMode) && fGoodCoincidenceEvent)
  {
    auto chdE = (ChannelData*) fChannelArray -> At(fdEArrayIdx[0]);
    auto chS1 = (ChannelData*) fChannelArray -> At(fS1ArrayIdx[0]);
    int groupdE = fMapDetectorGroup[chdE->midx][chdE->mch];
    int groupS1 = fMapDetectorGroup[chS1->midx][chS1->mch];
    bdE = chdE->energy;
    bESum = chdE->energy + chS1->energy;
    if (groupdE==groupS1) {
      fHistdAVSA -> Fill(chdE->adc,    chdE->adc    + chS1->adc);
      fHistdEVSE -> Fill(chdE->energy, chdE->energy + chS1->energy);
    }
    if (fS3ArrayIdx.size()==1)
    {
      auto chS3 = (ChannelData*) fChannelArray -> At(fS3ArrayIdx[0]);
      int groupS3 = fMapDetectorGroup[chS3->midx][chS3->mch];
      bE3 = chdE->energy + chS3->energy;
    }
    fdEArrayIdx.clear();
    fS1ArrayIdx.clear();
    fS3ArrayIdx.clear();
    fdEADC = 0.;
    fS1ADC = 0.;
    fS3ADC = 0.;
  }

  if (CheckEventCondition(bdE, bESum)==false)
    return false;

  if (fTritonCutGIsSet)
  {
    //auto detectorTheta;
    //auto ex = EvalEx(fBeamEnergy, bESum, detectorTheta);
    //fHistEx -> Fill(ex);
  }

  bTimeStamp = fTimeStampPrev;
  bNumChannels = fCountChannels;
  fTreeOut -> Fill();

  return true;
}

Long64_t Analysis::GetFileSize(TString fileName)
{
  std::ifstream file(fileName,std::ios::binary | std::ios::ate);
  Long64_t fileSize = file.tellg();
  file.close();
  return fileSize;
}

bool Analysis::AskContinueRun(TString message)
{
  fCountAskContinueRun++;

  std::string userInput0;
  if (fAutoUpdateRun) {
    if (CheckStopFile())
      userInput0 = "stop";
    else
      userInput0 = "wait";
  }
  else {
    if (!message.IsNull()) couti << message << endl;
    cout << "\033[0;32m" << Form("== (%d) Enter / stop / auto / wait: ",fCountEvents) << "\033[0m";
    std::getline(std::cin, userInput0);
  }
  TString userInput = userInput0;
  userInput.ToLower();

  if (userInput=="stop" || userInput=="exit" || userInput=="x") {
#ifdef DEBUG_EXIT_ANALYSIS
    coutn << "Exit ana: " << "user input " << userInput << endl;
#endif
    fExitAnalysis = true;
    return false;
  }
  else if (userInput.Index(".qqq")==0 || userInput=="qqq") {
    gApplication -> Terminate();
  }
  else if (userInput.Index(".q")==0 || userInput=="q") {
#ifdef DEBUG_EXIT_ANALYSIS
    coutn << "Exit ana: " << "user input " << userInput << endl;
#endif
    fExitAnalysis = true;
    fExitRoot = true;
    return false;
  }
  else
  {
    if (userInput=="auto")
    {
      fAutoUpdateRun = true;
      couti << "automatically update after waiting " << fUpdateAfterXSec << " s" << endl;
    }
    if (userInput=="wait")
    {
      //fUpdateAfterXSec = 0;
      //couti << "Reading all events" << endl;
    }
    if (userInput.IsDec() && userInput.Atoi()>0)
    {
      fUpdateAfterXSec = userInput.Atoi();
      couti << Form("Waiting next %d seconds",fUpdateAfterXSec) << endl;
    }

    int nn = 1;
    if      (fUpdateAfterXSec>1000) nn = 1000;
    else if (fUpdateAfterXSec>100) nn = 100;
    else if (fUpdateAfterXSec>10) nn = 10;
    else if (fUpdateAfterXSec>0) nn = 1;

    for (int sec=0; sec<fUpdateAfterXSec; ++sec)
    {
      if (sec%nn==0) {
        if (sec>0) {
          cout << "\r";
          cout << "                        ";
          cout << "\r";
        }
        cout << "   Waiting ... " << sec << std::flush;
      }
      sleep(1);
    }
    cout << endl;
  }
  return true;
}

bool Analysis::UpdateDrawing()
{
  if (fUpdateDrawingEveryNEvent>=0)
  {
    fCountEventsForUpdate++;

    // fUpdateDrawingEveryNEvent==0 will read all events and draw
    if (fUpdateDrawingEveryNEvent!=0 && fCountEventsForUpdate>=fUpdateDrawingEveryNEvent)
    {
      UpdateCvsOnline();
      if (AskUpdateDrawing("Continue update drawing?")==false)
        return false;
    }
  }
  return true;
}

bool Analysis::AskUpdateDrawing(TString message)
{
  fCountEventsForUpdate = 0;
  std::string userInput0;
  if (fAutoUpdateDrawing) {
    if (CheckStopFile())
      userInput0 = "stop";
    else {
      userInput0 = "";
      couti << "Auto update drawing ..." << endl;
    }
  }
  else {
    if (!message.IsNull()) couti << message << endl;
    cout << "\033[0;32m" << Form("== (%d) Enter / stop / auto / all: ",fCountEvents) << "\033[0m";
    std::getline(std::cin, userInput0);
  }
  TString userInput = userInput0;
  userInput.ToLower();
  if (userInput=="stop" || userInput=="exit" || userInput=="x") {
#ifdef DEBUG_EXIT_ANALYSIS
    coutn << "Exit ana: " << "user input " << userInput << endl;
#endif
    fExitAnalysis = true;
    return false;
  }
  else if (userInput.Index(".qqq")==0 || userInput=="qqq") {
    gApplication -> Terminate();
  }
  else if (userInput.Index(".q")==0 || userInput=="q") {
#ifdef DEBUG_EXIT_ANALYSIS
    coutn << "Exit ana: " << "user input " << userInput << endl;
#endif
    fExitAnalysis = true;
    fExitRoot = true;
    return false;
  }
  else
  {
    if (userInput=="auto")
    {
      fAutoUpdateDrawing = true;
      couti << Form("automatically update every %d events",fUpdateDrawingEveryNEvent) << endl;
    }
    if (userInput=="all")
    {
      fUpdateDrawingEveryNEvent = 0;
      couti << "Reading all events" << endl;
    }
    if (userInput.IsDec() && userInput.Atoi()>0)
    {
      fUpdateDrawingEveryNEvent = userInput.Atoi();
      couti << Form("Reading next %d events",fUpdateDrawingEveryNEvent) << endl;
    }
  }

  return true;
}

void Analysis::UpdateCvsOnline(bool firstDraw)
{
  if (fUpdateDrawingEveryNEvent<=0)
    return;

  if (fCvsOnline1==nullptr)
    return;

  auto DrawModuleBoundary = [this](double yMax) {
    for (int midx=0; midx<fNumModules; ++midx) {
      if (midx!=0) {
        auto line = new TLine(midx*fNumChannels,0,midx*fNumChannels,yMax);
        line -> SetLineColor(kBlack);
        line -> Draw("samel");
      }
      for (int iDiv : {1,2,3}) {
        auto line = new TLine(midx*fNumChannels+iDiv*4,0,midx*fNumChannels+iDiv*4,yMax*0.5);
        line -> SetLineColor(kGray+1);
        line -> SetLineStyle(2);
        line -> Draw("samel");
      }
      TString detectorName = TString(fDetectorName[fMapDetectorType[midx][0]](0,3));
      auto tt = new TText((midx+0.5)*fNumChannels,-yMax*0.1,detectorName);
      tt -> SetTextAlign(22);
      tt -> Draw("same");
    }
  };

  auto TakeCareOfStatsBox = [this](TH1* hist) {
    TPaveStats* box = dynamic_cast<TPaveStats*>(hist -> FindObject("stats"));
    if (box!=nullptr) {
      box -> Draw();
      //box -> SetBorderSize(0);
    }
  };

  if (fPadTrgRate==fPadEvtRate)
  {
    fPadTrgRate -> cd();
    //double maxTrigger = fHistTriggerRate->GetBinContent(fHistTriggerRate->GetMaximumBin());
    //double maxEvent = fHistEventRate->GetBinContent(fHistEventRate->GetMaximumBin());
    //if (maxTrigger<maxEvent)
    //  fHistTriggerRate -> SetMaximum(maxTrigger*1.05);
    //else
    //  fHistEventRate -> SetMaximum(maxEvent*1.05);

    fHistEventRate -> Draw();
    fHistTriggerRate -> Draw("same");
    fHistTriggerRateError -> Draw("same");
    fHistEventRateError -> Draw("same");
  }
  else {
    fPadTrgRate -> cd();
    fHistTriggerRate -> Draw();
    fHistTriggerRateError -> Draw("same");
    fPadEvtRate -> cd();
    fHistEventRate -> Draw();
    fHistEventRateError -> Draw("same");
  }

  fPadTSDist1 -> cd();
  fHistTSDist1 -> Draw();
  if (fPadTSDist2!=nullptr)
  {
    fPadTSDist2 -> cd();
    fHistTSDist2 -> Draw();
  }

  //fCvsOnline1 -> cd(1);
  fPadChCount -> cd();
  fHistChCount -> Draw();
  DrawModuleBoundary(fHistChCount->GetMaximum()*1.05);
  TakeCareOfStatsBox(fHistChCount);

  //fCvsOnline1 -> cd(2);
  fPadADC -> cd();
  TH1D* histEOnline = fHistADC;
  if (fShowEnergyConversion)
    histEOnline = fHistE;
  histEOnline -> Draw();

  //fCvsOnline1 -> cd(3);
  fPadEVSCh -> cd();
  TH2D* hist2DOnline = fHistAVSCh;
  double yMax = fMaxADC;
  if (fShowEnergyConversion) {
    hist2DOnline = fHistEVSCh;
    yMax = fMaxE;
  }
  hist2DOnline -> Draw("colz");
  DrawModuleBoundary(yMax);
  TakeCareOfStatsBox(hist2DOnline);

  //fCvsOnline1 -> cd(4);
  fPaddEVSE -> cd();
  TH2D* histdEEOnline = fHistdAVSA;
  if (fShowEnergyConversion) {
    histdEEOnline = fHistdEVSE;
  }
  histdEEOnline -> Draw("colz");

  fPadEx -> cd();
  fHistEx -> Draw();

  if (fCvsOnline2!=nullptr) {
    fCvsOnline2 -> Modified();
    fCvsOnline2 -> Update();
  }
  fCvsOnline1 -> Modified();
  fCvsOnline1 -> Update();
}

bool Analysis::CheckStopFile()
{
  if (access("stop",F_OK)!=-1) {
    coutd << "found stop file" << endl;
    coutn << "Please remove stop file!" << endl;
    coutn << "Please remove stop file!" << endl;
    coutn << "Please remove stop file!" << endl;
    coutn << "Please remove stop file!" << endl;
    coutn << "Please remove stop file!" << endl;
    return true;
  }
  return false;
}

void Analysis::EndOfConversion()
{
  fFileOut -> cd();
  cout << endl;
  couti << "End of conversion!" << endl;
  couti << "Writting tree ..." << endl;
  fTreeOut -> Write();

  if (fCvsOnline1!=nullptr)
  {
    UpdateCvsOnline();
    fCvsOnline1 -> Write();
  }
  fHistChCount -> Write();
  fHistADC -> Write();
  fHistAVSCh -> Write();
  fHistEVSCh -> Write();

  PrintConversionSummary();

  couti << "End of conversion!" << endl;
  couti << "Output file name: " << fFileNameOut << endl;

  if (fExitRoot)
    gApplication -> Terminate();
}

void Analysis::PrintConversionSummary()
{
  cout << endl;
  couti << "Number of events: " << fCountEvents << endl;
  couti << "Number of channels: " << fCountAllChannels << " (from " << fCountAllLines << ")" << endl;
  couti << "Number of times TS-error occured: " << fCountTSError << endl;
  couti << "Number of events with 0   coincidence channels: " << fCoincidenceCount[0] << endl;
  couti << "Number of events with 1   coincidence channels: " << fCoincidenceCount[1] << endl;
  couti << "Number of events with 2   coincidence channels: " << fCoincidenceCount[2] << endl;
  couti << "Number of events with 3   coincidence channels: " << fCoincidenceCount[3] << endl;
  couti << "Number of events with =>4 coincidence channels: " << fCoincidenceCount[4] << endl;
}

int Analysis::GetModuleIndex(int FENumber)
{
  return fMapFEToModuleIndex[FENumber];
}

int Analysis::GetGlobalID(UShort_t midx, UShort_t mch)
{
  int globalID = midx*fNumChannels + mch;
  return globalID;
}

void Analysis::GetModCh(int globalID, UShort_t &midx, UShort_t &mch)
{
  midx = globalID / fNumChannels;
  mch = globalID % fNumChannels;
}

void Analysis::DetectorToModule(int det, int dch, int &midx, int &mch)
{
  auto gid = fMapDetectorToGlobalID[det][dch];
  midx = fMapGlobalIDToModuleIndex[gid];
  mch = fMapGlobalIDToMCh[gid];
}

TString Analysis::GetDetectorTitle(int midx, int mch, bool addChannel)
{
  int det = fMapDetectorType[midx][mch];
  int dch = fMapDetectorChannel[midx][mch];
  TString detectorTitle = fDetectorName[det];
  if (addChannel)
    detectorTitle = detectorTitle + Form("-%d",dch);
  return detectorTitle;
}

double Analysis::GetCalibratedEnergy(int midx, int mch, int adc)
{
  if (fFxEnergyConversion[midx][mch]==nullptr) {
    if (!fECalErrorWasSentOut) {
      coute << "Energy conversion (" << midx << " " << mch << ") is nullptr! ..." << endl;
      fECalErrorWasSentOut = true;
    }
    return 0;
  }
  return fFxEnergyConversion[midx][mch] -> Eval(adc);
}

void Analysis::MakeTritonCutFile(TString fileName)
{
  if (fileName.IsNull())
    fileName = "tritonCutG.root";
  fileName = getAna()->GetOutputPath() + fileName;
  cout << "Creating " << fileName << endl;
  auto file = new TFile(fileName,"recreate");
  TCutG* tritonCutG = (TCutG*) gROOT->GetListOfSpecials()->FindObject("CUTG");
  tritonCutG -> SetName("tritonCutG");
  file -> cd();
  tritonCutG -> SetVarX("de");
  tritonCutG -> SetVarY("ee");
  //tritonCutG -> SetVarX("ch.gid");
  //tritonCutG -> SetVarY("ch.energy");
  
  
  tritonCutG -> Write();
}

void Analysis::CallTritonCutFile(TString fileName)
{
  if (fileName.IsNull())
    fileName = "out/tritonCutG.root";
  couti << "Set triton graphic cut from " << fileName << endl;

  TDirectory* currentDirectory;
  if (gDirectory!=nullptr)
    currentDirectory = gDirectory;
  TFile* file = new TFile(fileName,"read");
  TCutG* tritonCutG = (TCutG*) file -> Get("tritonCutG");
  tritonCutG -> SetVarX("de");
  tritonCutG -> SetVarY("ee");
  //tritonCutG -> SetVarX("ch.gid");
  //tritonCutG -> SetVarY("ch.energy");
  if (currentDirectory!=nullptr)
    currentDirectory -> cd();
}

void Analysis::SetTritonCutFile(TString fileName)
{
  couti << "Set triton graphic cut from " << fileName << endl;

  TFile* file = new TFile(fileName,"read");
  if (file->IsZombie()) {
    coute << fileName << " is zombie!!" << endl;
    return;
  }

  fTritonCutG = (TCutG*) file -> Get("tritonCutG");
  fTritonCutG -> SetVarX("de");
  fTritonCutG -> SetVarY("ee");

  fTritonCutGIsSet = true;
}

const double mp = 938.791;
const double m48cr = 44669.3;
const double m50cr = 46524.8;
const double mt = 2809.45;
double Analysis::EvalEx(double tp, double tt, double theta)
{
  double mtPtt = 2809.45+tt;
  double mtMtt = 2809.45*tt;
  double val1 = 2.25559e+09;
  double val2 = 93049.6*mtPtt - 56327.5*mtMtt;
  double val3 = 2*sqrt((mp+tp)*(mp+tp)-881329);
  double val4 = sqrt(mtPtt*mtPtt-7.89301e+06)*cos(theta);
  double ex = sqrt(val1 - val2 + val3 * val4) - m48cr;
  return ex;

  /*
  {
    double val1 = mp*mp + m50cr*m50cr + mt+mt + 2*(mp+tp)*m50cr;
    double val2 = 2*(mt+tt)*m50cr - 2*(mp*tp)*(mt*tt);
    double val3 = 2*sqrt((mp+tp)*(mp+tp)-mp*mp) * sqrt((mt+tt)*(mt+tt)-mt*mt)*cos(theta);
    double ex = sqrt(val1-val2+val3) - m48cr;

    return ex;
  }
  */
}
