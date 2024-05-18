#include "Analysis.h"

// Replace define value to "coutx" to omit printing
#define coutd cout<<"+\033[0;36m"<<__LINE__<<" "<<__FILE__<<" #\033[0m "
#define couti cout<<"\033[0;32m==\033[0m "
#define coutw cout<<"\033[0;33mWR\033[0m "
#define coute cout<<"\033[0;31mER\033[0m "
//#define coutt cout<<"\033[0;33mTSERROR\033[0m "
#define coutt coutx

Analysis::Analysis()
{
  coutx.open("out/dummy_stream");
  InitializeMapping();
}

void Analysis::InitializeMapping()
{
  std::ifstream mapFile(fMapFileName);
  if (mapFile.fail()) {
    coute << "Cannot open mapping file: " << fMapFileName << endl;
    return;
  }

  fMapDetectorType = new int*[fNumModules];
  fMapDetectorChannel = new int*[fNumModules];
  fMapDetectorReplaced = new bool*[fNumModules];
  fMapDetectorRFMod = new int*[fNumModules];
  fMapDetectorRFMCh = new int*[fNumModules];
  for (int iModule=0; iModule<fNumModules; ++iModule) {
    fMapDetectorType[iModule] = new int[fNumChannels];
    fMapDetectorChannel[iModule] = new int[fNumChannels];
    fMapDetectorReplaced[iModule] = new bool[fNumChannels];
    fMapDetectorRFMod[iModule] = new int[fNumChannels];
    fMapDetectorRFMCh[iModule] = new int[fNumChannels];
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
      fMapDetectorType[iModule][iChannel] = kDummyDetector;
      fMapDetectorChannel[iModule][iChannel] = -1;
      fMapDetectorReplaced[iModule][iChannel] = false;
      fMapDetectorRFMod[iModule][iChannel] = -1;
      fMapDetectorRFMCh[iModule][iChannel] = -1;
    }
  }

  TString name;
  string buffer;
  bool replaced;
  UShort_t module0, channelID0;
  UShort_t module, channelID, detectorChannel;
  getline(mapFile, buffer); // skip header
  while (mapFile >> module >> channelID >> name >> detectorChannel >> replaced)
  {
    //couti << module << " " << channelID << " " << name << " " << detectorChannel << endl;
    int detectorType = kDummyDetector;
    if (name=="S1Junction")   detectorType = kS1Junction;
    if (name=="S1Ohmic")      detectorType = kS1Ohmic;
    if (name=="S3Junction")   detectorType = kS3Junction;
    if (name=="S3Ohmic")      detectorType = kS3Ohmic;
    if (name=="dEDetector")   detectorType = kdEDetector;
    if (name=="Scintillator") detectorType = kScintillator;
    if (name=="FaradayCup")   detectorType = kFaradayCup;
    int fake = GetFakeModule(module);
    fMapDetectorType[fake][channelID] = detectorType;
    fMapDetectorChannel[fake][channelID] = detectorChannel;
    fMapDetectorReplaced[fake][channelID] = replaced;
    if (replaced) {
      mapFile >> module0 >> channelID0;
      fMapDetectorRFMod[fake][channelID] = module0;
      fMapDetectorRFMCh[fake][channelID] = channelID0;
    }
  }
  mapFile.close();

  fDetectorName[kDummyDetector] = "X";
  fDetectorName[kS1Junction   ] = "S1J";
  fDetectorName[kS1Ohmic      ] = "S1O";
  fDetectorName[kS3Junction   ] = "S3J";
  fDetectorName[kS3Ohmic      ] = "S3O";
  fDetectorName[kdEDetector   ] = "dE";
  fDetectorName[kScintillator ] = "SC";
  fDetectorName[kFaradayCup   ] = "FC";
}

void Analysis::RunConversion(int runNo, TString pathIn)
{
  fRunNo = runNo; 
  fRunName = Form("RUN%03d",fRunNo);
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/analysis/input/") : pathIn);

  InitializeConversion();
  ConfigureDateTime();
  MakeOutputFile();
  InitializeDrawing();
  ReadDataFile();
  EndOfConversion();
}

void Analysis::AddEnergyCalibration(TString name)
{
  couti << "Set energy calibration functions from " << name << endl;
  for (int fake=0; fake<fNumModules; ++fake)
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      fFxEnergyConversion[fake][iChannel] = nullptr;

  TFile* file = new TFile(name,"read");

  for (int fake=0; fake<fNumModules; ++fake) {
    int module = GetRealModule(fake);
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel) {
      fFxEnergyConversion[fake][iChannel] = (TF1*) file -> Get(Form("ae%02d%02d",module,iChannel));
    }
  }

  for (int fake=0; fake<fNumModules; ++fake) {
    int module = GetRealModule(fake);
    bool missing = false;
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      if (fFxEnergyConversion[fake][iChannel]==nullptr) {
        missing = true;
        break;
      }
    if (!missing) break;
    coute << "Missing ECal: ";
    for (int iChannel=0; iChannel<fNumChannels; ++iChannel)
      if (fFxEnergyConversion[fake][iChannel]==nullptr)
        cout << "(" << module << "," << iChannel << ") ";
    cout << endl;
  }

  fEnergyConversionSet = true;
  file -> Close();
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
  fTreeSummary -> SetBranchAddress("ch",&fChannelArray);
  fNumEventsSummary = fTreeSummary -> GetEntries();

  //fFileSummary -> ls();
  //fRunNo = ((TParameter<int>*)fFileSummary->Get("run"))->GetVal();
  fRunNo = TString(((TNamed*) fFileSummary->Get("run"))->GetTitle()).Atoi();
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

void Analysis::InitializeAlphaAnalysis(TString fileName)
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

void Analysis::AnalyzeAlphaTestModule(int module, bool drawAnalysis, TString fileName)
{
  InitializeAlphaAnalysis(fileName);
  TCanvas* cvs = nullptr;
  if (drawAnalysis) {
    cvs = new TCanvas(Form("cvsEModule%d",module),Form("energy of module %d",module),1600,900);
    cvs -> Divide(4,4);
  }
  bool moduleDataExist = false;
  for (auto channelID=0; channelID<fNumChannels; ++channelID)
  {
    TVirtualPad* pad = nullptr;
    if (drawAnalysis)
      pad = cvs->cd(channelID+1);

    bool channelDataExist = AnalyzeAlphaTest(module, channelID, drawAnalysis, pad);
    moduleDataExist = (moduleDataExist || channelDataExist);
  }

  if (drawAnalysis)
  {
    if (!moduleDataExist) {
      coute << "Module " << module << " is empty!" << endl;
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

bool Analysis::AnalyzeAlphaTest(int module, int channelID, bool drawAnalysis, TVirtualPad* cvs)
{
  if (fTreeSummary==nullptr) {
    coute << "Summary tree is nullptr! Run ReadSummaryFile before this method!";
    return false;
  }
  
  bool isSingleDrawing = false;
  if (cvs==nullptr)
    isSingleDrawing = true;

  couti << Form("channel (%d,%d)",module,channelID) << endl;

  auto gid = GetGlobalID(module, channelID);
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
  TString titleHist = Form("channel (%d,%d) energy distribution;ADC;count",module,channelID);
  auto histE = new TH1D(nameHist,titleHist,fNumADC,0,fMaxADC);
  fTreeSummary -> Project(nameHist, "ch.energy", Form("ch.gid==%d",gid));
  if (histE->GetEntries()==0) {
    coute << "Data is empty!" << endl;
    //0.001810*(x-0.0000)
    if (fFileAlpha!=nullptr) {
      auto fxADCToEnergyDummy = new TF1(Form("ae%02d%02d",module,channelID),"0.001810*(x-0.0000)",0,10000);
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
    //auto fxADCToEnergy = new TF1(Form("ae%02d%02d",module,channelID),"[0]*(x-[1])",0,10000);
    //fxADCToEnergy -> SetParameter(0,energyADCRatio);
    //fxADCToEnergy -> SetParameter(1,ADCOffset);
    auto fxADCToEnergy = new TF1(Form("ae%02d%02d",module,channelID),Form("%f*(x-%.4f)",energyADCRatio,ADCOffset),0,10000);
    fxADCToEnergy -> Write();
    //(new TParameter<double>("ADCOffset",ADCOffset)) -> Write();
    //(new TParameter<double>("ADCEnergyRatio",ADCEnergyRatio)) -> Write();
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

void Analysis::InitializeConversion()
{
  fPathToOutput = "/home/daquser/data/LiCD2Irrad/analysis/out/";
  fCountEvents = 0;
  fChannelArray = nullptr;
  fCountChannels = 0;
  fCountAllChannels = 0;
  fCountTSError = 0;
  fReturnIfNoFile = true;
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

void Analysis::InitializeDrawing()
{
  if (fFileOut!=nullptr)
    fFileOut -> cd();

  if (fUpdateDrawingEveryNEvent>0)
  {
    fCvsOnline = new TCanvas("cvsOnline",fRunName+" online update canvas",1600,900);
    fCvsOnline -> Divide(2,2);
    fCvsOnline -> cd(3) -> SetLogz();
    fCvsOnline -> cd(4) -> Divide(2,1);
    fPadTSDist1 = fCvsOnline -> cd(4) -> cd(1);
    fPadTSDist2 = fCvsOnline -> cd(4) -> cd(2);
  }

  fHistTSDist1 = new TH1D("histTSDist1",fRunName+" TS distance (+);TS-dist;count", 20,0,20);
  fHistTSDist2 = new TH1D("histTSDist2",fRunName+" TS distance (-);TS-dist;count", 50,-1000000,0);
  fHistChCount = new TH1D("histChCount",fRunName+" channel count;global-ch;event count", fNumCh,0,fNumCh);
  fHistChCount -> SetFillColor(29);
  fHistADC = new TH1D("histADC",fRunName+" ADC (all channels);ADC",fNumADC,0,fMaxADC);
  fHistADC -> SetFillColor(29);
  fHistE = new TH1D("histEnergy",fRunName+" energy (all channels);energy (MeV)",fNumE,0,fMaxE);
  fHistE -> SetFillColor(29);
  fHistAVSCh = new TH2D("histAVSCh", fRunName+" ADC vs channel-id;global-ch;ADC", fNumCh,0,fNumCh,fNumADC,0,fMaxADC);
  fHistEVSCh = new TH2D("histEVSCh", fRunName+" energy vs channel-id;global-ch;energy (MeV)", fNumCh,0,fNumCh,fNumE,0,fMaxE);

  if (fUpdateDrawingEveryNEvent>0) {
    SetAttribute(fHistTSDist1,fPadTSDist1,2);
    SetAttribute(fHistTSDist2,fPadTSDist2,2);
    SetAttribute(fHistChCount,fCvsOnline->cd(1));
    SetAttribute(fHistADC,fCvsOnline->cd(2));
    SetAttribute(fHistE,fCvsOnline->cd(2));
    SetAttribute(fHistAVSCh,fCvsOnline->cd(3),1,true);
    SetAttribute(fHistEVSCh,fCvsOnline->cd(3),1,true);
  }

  UpdateCvsOnline(true);
}

void Analysis::SetAttribute(TH1* hist, TVirtualPad* pad, int npad, bool td)
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
      if (td)
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
  (new TNamed("run", Form("%d",fRunNo))) -> Write();
  if (option==0) {
    (new TNamed("Energy threshold   :", Form("%d",fADCThreshold))) -> Write();
    (new TNamed("Draw every n-event :", Form("%d",fUpdateDrawingEveryNEvent))) -> Write();
    (new TNamed("Return if no file  :", Form("%d",fReturnIfNoFile))) -> Write();
    (new TNamed("Ignore file update :", Form("%d",fIgnoreFileUpdate))) -> Write();
    (new TNamed("Ignore TS decrease :", Form("%d",fSkipTSError))) -> Write();
    (new TNamed("Event count limit  :", Form("%d",fEventCountLimit))) -> Write();
    (new TNamed("File number range  :", Form("%d %d",fFileNumberRange1,fFileNumberRange2))) -> Write();
    (new TNamed("ADC Threshold      :", Form("%d",fADCThreshold))) -> Write();
    (new TNamed("Coincidence dTS cut:", Form("%d",fCoincidenceTSRange))) -> Write();
  }
}

void Analysis::MakeOutputFile()
{
  if (fFileNameOut.IsNull())
    fFileNameOut = fPathToOutput+Form("RUN%03d.summary.root",fRunNo);
  cout << endl;
  couti << "Creating output file " << fFileNameOut << endl;
  fFileOut = new TFile(fFileNameOut,"recreate");
  fTreeOut = new TTree("event","");
  fChannelArray = new TClonesArray("ChannelData",10);
  fTreeOut -> Branch("ts",&bTimeStamp);
  fTreeOut -> Branch("nch",&bNumChannels);
  fTreeOut -> Branch("ch",&fChannelArray);

  WriteRunParameters(fFileOut,0);
}

void Analysis::ReadDataFile()
{
  ifstream fileIn;
  streamsize fileSize;
  streamsize fileSizeOld = 0; // Size of opened Raw-Data file
  const streamsize fileSizeMax = 500000000; // 500 MB
  UShort_t countInputs = fFileNumberRange1; // Number of Files read
  UShort_t module, channelID, adc;
  Long64_t timeStampLocal, timeStamp;
  Int_t tsGroup;
  Int_t numData = 0;
  char buffer[256];
  int countOpenFileTrials = 0;

  ChannelData* channelData = NULL;
  fTimeStampPrev = -1;

  while(1)
  {
    if (fFileNumberMax >=0 && countInputs>fFileNumberMax) {
      couti << "Number of inputs " << countInputs << " exceed maximum input number " << fFileNumberMax << endl;
      break;
    }

    TString fileNameInput = TString::Format("%s/RUN%03d_%s_list_%03d.dat", fPathToInput.Data(),fRunNo,fDateTime.Data(),countInputs);
    cout << endl;
    couti << "Reading " << fileNameInput << endl;

    countOpenFileTrials = 0;
    while (1) // endless loop (2) to check the status of the Raw-Data file //
    {
      fileIn.open(fileNameInput);
      if (fileIn.fail())
      {
        fileIn.close();
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

      fileSize = fileIn.seekg(0, ios::end).tellg(); // obtain the size of file
      fileIn.seekg(0, ios::beg); // rewind

      if (fileSize>fileSizeMax || (countOpenFileTrials!=0 && fileSize==fileSizeOld)) // good file or final file
      {
        countOpenFileTrials = 0;
        countInputs++;
        break;
      }
      else if(countOpenFileTrials==0){ // first try
        countOpenFileTrials++;
        fileIn.close();
        fileSizeOld = fileSize;
        if (fFirstFileOpened) {
          couti << "waiting(" << countOpenFileTrials << ") 3s ..." << endl;
          sleep(3);
        }
        fFirstFileOpened = true;
      }
      else if(fileSize > fileSizeOld){ // writing the file is still continued ...
        if (fIgnoreFileUpdate)
          break;
        fileIn.close();
        fileSizeOld = fileSize;
        couti << "waiting(" << countOpenFileTrials << ") 60s ..." << endl;
        sleep(60);
      }
    }

    if (countOpenFileTrials != 0) {
      coutw << fileNameInput << " is not found! exit." << endl;
      break;
    }

    //couti << "File size is " << fileSize << endl;

    fChannelArray -> Clear("C");
    Long64_t countEventsSingleFile = 0;
    fCountChannels = 0;
    Long64_t countLine = 0;
    while (fileIn >> buffer)
    {
      countLine++;
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &module, &channelID, &adc, &timeStampLocal, &tsGroup);
      if (numData != 5) {
        coute << TString::Format("Number of data in line is not 5 (%d) %u %u %u %lld %u", numData, module, channelID, adc, timeStampLocal, tsGroup) << endl;
        continue;
      }

      if (adc < fADCThreshold)
        continue;

      timeStamp = timeStampLocal;
      if (tsGroup>0)
        timeStamp += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      if (bTimeStamp>=0) {
        int dts = int(timeStamp)- int(fTimeStampPrevTrue);
        if (dts>=0) fHistTSDist1 -> Fill(dts);
        else fHistTSDist2 -> Fill(dts);
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

      if (eventStatus==kTSError) // time-stamp decreased
      {
        if (timeStamp<fTimeStampPrevTrue) {
          coutt << "Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          fTimeStampPrevTrue = timeStamp;
          ++fCountTSError;
        }
        if (fStopAtTSError) {
          coute << "Return due to Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          return;
        }
        continue;
      }
      else if (eventStatus==kSameEvent) // same event
      {
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);
      }
      else if (eventStatus==kNextEvent) // next event
      {
        if (countEventsSingleFile>0)
          FillDataTree();
        if (fExitAnalysis)
          break;

        fCountChannels = 0;
        fChannelArray -> Clear("C");
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);

        fTimeStampPrev = timeStamp;
        fTimeStampPrevTrue = timeStamp;
        fCountEvents++;
        countEventsSingleFile++;
      }

      int gid = GetGlobalID(module, channelID);
      channelData -> SetData(module,channelID,gid,adc,timeStampLocal,tsGroup,timeStamp);
      fCountAllChannels++;

      if (fUpdateDrawingEveryNEvent>=0)
      {
        fHistChCount -> Fill(gid);
        fHistADC -> Fill(adc);
        fHistAVSCh -> Fill(gid, adc);
        if (fEnergyConversionSet)
        {
          double energy = fFxEnergyConversion[GetFakeModule(module)][channelID] -> Eval(adc);
          fHistEVSCh -> Fill(gid, energy);
          fHistE -> Fill(energy);
        }
      }

      if (fEventCountLimit>0 && fCountEvents>=fEventCountLimit) {
        couti << "Event count limit at " << fCountEvents << "!" << endl;
        fExitAnalysis = true;
      }
    }

    if (!fExitAnalysis && fCountEvents>0) 
      FillDataTree();

    if (fExitAnalysis)
      break;

    couti << "Number of events from last file: " << countEventsSingleFile << endl;

    couti << "Number of all events: " << fCountEvents << endl;
    couti << "Number of all channels: " << fCountAllChannels << endl;
    couti << "Number of times TS-error occured: " << fCountTSError << endl;
    couti << "Number of events with 0   coincidence channels: " << fCoincidenceCount[0] << endl;
    couti << "Number of events with 1   coincidence channels: " << fCoincidenceCount[1] << endl;
    couti << "Number of events with 2   coincidence channels: " << fCoincidenceCount[2] << endl;
    couti << "Number of events with 3   coincidence channels: " << fCoincidenceCount[3] << endl;
    couti << "Number of events with =>4 coincidence channels: " << fCoincidenceCount[4] << endl;

    fileIn.close();
  }

  UpdateCvsOnline();
}

bool Analysis::FillDataTree()
{
  if      (fCountChannels==0) fCoincidenceCount[0]++;
  else if (fCountChannels==1) fCoincidenceCount[1]++;
  else if (fCountChannels==2) fCoincidenceCount[2]++;
  else if (fCountChannels==3) fCoincidenceCount[3]++;
  else if (fCountChannels>=4) fCoincidenceCount[4]++;

  if (fCoincidenceMult>0)
    if (fCountChannels!=fCoincidenceMult)
      return false;

  bTimeStamp = fTimeStampPrev;
  bNumChannels = fCountChannels;
  fTreeOut -> Fill();

  if (fUpdateDrawingEveryNEvent>=0)
  {
    fCountEventsForUpdate++;

    // fUpdateDrawingEveryNEvent==0 will read all events and draw
    if (fUpdateDrawingEveryNEvent!=0 && fCountEventsForUpdate>=fUpdateDrawingEveryNEvent)
    {
      UpdateCvsOnline();
      fCountEventsForUpdate = 0;
      std::string userInput0;
      if (fAutoUpdateDrawing)
        userInput0 = "";
      else {
        cout << "\033[0;32m" << Form("== (%d) Enter / stop / auto / all: ",fCountEvents,fUpdateDrawingEveryNEvent) << "\033[0m";
        std::getline(std::cin, userInput0);
      }
      TString userInput = userInput0;
      userInput.ToLower();
      if (userInput=="stop" || userInput=="exit" || userInput=="x")
        fExitAnalysis = true;
      else if (userInput.Index(".q")==0 || userInput=="q") {
        fExitAnalysis = true;
        fExitRoot = true;
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
    }
  }
  return true;
}

void Analysis::UpdateCvsOnline(bool firstDraw)
{
  if (fUpdateDrawingEveryNEvent<=0)
    return;

  if (fCvsOnline==nullptr)
    return;

  auto DrawModuleBoundary = [this](double yMax) {
    for (int iModule=0; iModule<fNumModules; ++iModule) {
      if (iModule!=0) {
        auto line = new TLine(iModule*fNumChannels,0,iModule*fNumChannels,yMax);
        line -> SetLineColor(kBlack);
        line -> Draw("samel");
      }
      for (int iDiv : {1,2,3}) {
        auto line = new TLine(iModule*fNumChannels+iDiv*4,0,iModule*fNumChannels+iDiv*4,yMax*0.5);
        line -> SetLineColor(kGray+1);
        line -> SetLineStyle(2);
        line -> Draw("samel");
      }
      //auto tt = new TText((iModule+0.5)*fNumChannels,yMax*0.1,Form("%d",GetRealModule(iModule)));
      auto tt = new TText((iModule+0.5)*fNumChannels,yMax*0.1,fDetectorName[GetDetectorType(iModule)]);
      //tt -> SetTextColor(kGray);
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

  fPadTSDist1 -> cd();
  fHistTSDist1 -> Draw();
  fPadTSDist2 -> cd();
  fHistTSDist2 -> Draw();

  fCvsOnline -> cd(1);
  fHistChCount -> Draw();
  DrawModuleBoundary(fHistChCount->GetMaximum()*1.05);
  TakeCareOfStatsBox(fHistChCount);

  fCvsOnline -> cd(2);
  TH1D* histEOnline = fHistADC;
  if (fShowEnergyConversion)
    histEOnline = fHistE;
  histEOnline -> Draw();

  fCvsOnline -> cd(3);
  TH2D* hist2DOnline = fHistAVSCh;
  double yMax = fMaxADC;
  if (fShowEnergyConversion) {
    hist2DOnline = fHistEVSCh;
    yMax = fMaxE;
  }
  hist2DOnline -> Draw("colz");
  DrawModuleBoundary(yMax);
  TakeCareOfStatsBox(hist2DOnline);
  fCvsOnline -> Modified();
  fCvsOnline -> Update();
}

void Analysis::EndOfConversion()
{
  fFileOut -> cd();
  fTreeOut -> Write();

  if (fCvsOnline!=nullptr)
  {
    UpdateCvsOnline();
    fCvsOnline -> Write();
  }
  fHistChCount -> Write();
  fHistADC -> Write();
  fHistAVSCh -> Write();
  fHistEVSCh -> Write();

  cout << endl;
  couti << "End of conversion!" << endl;
  couti << "Number of events: " << fCountEvents << endl;
  couti << "Number of all channels: " << fCountAllChannels << endl;
  couti << "Number of times TS-error occured: " << fCountTSError << endl;
  couti << "Number of events with 0   coincidence channels: " << fCoincidenceCount[0] << endl;
  couti << "Number of events with 1   coincidence channels: " << fCoincidenceCount[1] << endl;
  couti << "Number of events with 2   coincidence channels: " << fCoincidenceCount[2] << endl;
  couti << "Number of events with 3   coincidence channels: " << fCoincidenceCount[3] << endl;
  couti << "Number of events with =>4 coincidence channels: " << fCoincidenceCount[4] << endl;
  couti << "Output file name: " << fFileNameOut << endl;

  if (fExitRoot)
    gApplication -> Terminate();
}

int Analysis::GetFakeModule(int module)
{
  if (module>2) return (module - 6);
  return module;
}

int Analysis::GetRealModule(int fake)
{
  if (fake>2) return (fake + 6);
  return fake;
}

int Analysis::GetGlobalID(UShort_t module, UShort_t channelID)
{
  if (module>2) module = module - 6;
  int globalID = module*fNumChannels + channelID;
  return globalID;
}

void Analysis::GetModCh(int globalID, UShort_t &module, UShort_t &channelID)
{
  module = globalID / fNumChannels;
  if (module>2) module = module + 6;
  channelID = globalID % fNumChannels;
}

int Analysis::GetDetectorType(int iModule, int channelID)
{
  return fMapDetectorType[iModule][channelID];
}
