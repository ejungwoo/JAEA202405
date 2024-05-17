#include "Analysis.h"

// Replace define value to "coutx" to omit printing
#define coutd cout<<"+\033[0;36m"<<__LINE__<<" "<<__FILE__<<" #\033[0m "
#define couti cout<<"\033[0;32m==\033[0m "
#define coutw cout<<"\033[0;33mWR\033[0m "
#define coute cout<<"\033[0;31mER\033[0m "
//#define coutt cout<<"\033[0;33mTSERROR\033[0m "
#define coutt coutx

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

void Analysis::ReadSummaryFile(TString fileName)
{
  fFileSummary = new TFile(fileName,"read");
  fTreeSummary = (TTree*) fFileSummary -> Get("event");
  fTreeSummary -> SetBranchAddress("ts",&bTimeStamp);
  //fTreeSummary -> SetBranchAddress("tsDist",&bTimeStampDist);
  fTreeSummary -> SetBranchAddress("ch",&fChannelArray);
  fNumEventsSummary = fTreeSummary -> GetEntries();

  //fFileSummary -> ls();
  //fRunNo = ((TParameter<int>*)fFileSummary->Get("run"))->GetVal();
  fRunNo = TString(((TNamed*) fFileSummary->Get("run"))->GetTitle()).Atoi();
  couti << fileName << "(" << fRunNo << ") containing " << fNumEventsSummary << " events" << endl;
}

Double_t Analysis::FxTwoAlpha(Double_t *xy, Double_t *par)
{
  // 5.486 MeV (85.2 %)
  // 5.443 MeV (12.8 %)

  double x = xy[0];

  double ADCOffset = par[0];
  double energyResolution = par[1];

  double mean1 = par[2]; // ADC mean of 5.486 MeV peak
  double meanP1 = mean1 - ADCOffset; // pure ADC without ADC-offset
  double sigma1 = energyResolution * meanP1;
  double amplitude1 = par[3];

  double ADCEnergyRatio = meanP1 / 5.486;

  double mean2 = 5.443 * ADCEnergyRatio + ADCOffset;
  double meanP2 = mean2 - ADCOffset; // pure ADC without ADC-offset
  double sigma2 = energyResolution * meanP2;
  double amplitude2 = amplitude1 * sigma1 / 6.65625 / sigma2; // 85.2 / 12.8 = 6.65625; 

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

  double ADCEnergyRatio = meanP1 / 5.486;

  mean2 = 5.443 * ADCEnergyRatio + ADCOffset;
  double meanP2 = mean2 - ADCOffset; // pure ADC without ADC-offset
  sigma2 = energyResolution * meanP2;
  amplitude2 = amplitude1 * sigma1 / 6.65625 / sigma2; // 85.2 / 12.8 = 6.65625; 
}

TF1* Analysis::CreateFxTwoAlpha(TString name, double min, double max)
{
  auto f1 = new TF1(name,this,&Analysis::FxTwoAlpha,min,max,4,"Analysis","FxTwoAlpha");
  return f1;
}

void Analysis::InitializeAlphaAnalysis()
{
  if (fFileNameAlpha.IsNull())
    fFileNameAlpha = fPathToOutput+Form("RUN%03d.alpha.root",fRunNo);
  cout << endl;
  couti << "Creating output file " << fFileNameAlpha << endl;
  fFileAlpha = new TFile(fFileNameAlpha,"recreate");
  WriteRunParameters(fFileAlpha,1);
}

void Analysis::AnalyzeAlphaTestModule(int module, bool drawAnalysis)
{
  InitializeAlphaAnalysis();
  TCanvas* cvs = nullptr;
  if (drawAnalysis) {
    cvs = new TCanvas(Form("cvsEModule%d",module),Form("energy of module %d",module),1600,900);
    cvs -> Divide(4,4);
  }
  bool moduleDataExist = false;
  for (auto channelID=0; channelID<16; ++channelID)
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

  cout << endl;
  couti << Form("channel (%d,%d)",module,channelID) << endl;

  auto gid = GetGlobalID(module, channelID);
  if (fFitAlpha==nullptr)
  {
    fFitAlpha = CreateFxTwoAlpha("fxTwoAlpha",0,fMaxE);

    fFitAlpha1 = new TF1("fxAlpha1","gaus(0)",0,fMaxE);
    fFitAlpha2 = new TF1("fxAlpha2","gaus(0)",0,fMaxE);
    fFitAlpha1 -> SetLineColor(kBlue);
    fFitAlpha2 -> SetLineColor(kGreen+1);
  }
  double mean1, sigma1, amplitude1, mean2, sigma2, amplitude2;

  TString nameHist = Form("histE%d",gid);
  TString titleHist = Form("channel (%d,%d) energy distribution;energy;count",module,channelID);
  auto histE = new TH1D(nameHist,titleHist,fNumE,0,fMaxE);
  fTreeSummary -> Project(nameHist, "ch.energy", Form("ch.gid==%d",gid));
  if (histE->GetEntries()==0) {
    coute << "Data is empty!" << endl;
    return false;
  }

  double mean = histE -> GetMean();
  double sig = histE -> GetStdDev();

  double bin = histE->GetMaximumBin();
  double yMax = histE->GetBinContent(bin);
  double xPeak = histE->GetBinCenter(bin);

  fFitAlpha -> SetParameters(mean, sig, yMax);
  auto x1 = xPeak - 100;
  auto x2 = xPeak + 100;
  fFitAlpha -> SetParameters(10,0.003,xPeak,yMax);
  fFitAlpha -> SetParLimits(0,0,100);
  histE -> Fit(fFitAlpha,"Q0N","",x1,x2);
  auto parameters = fFitAlpha -> GetParameters();
  Convert2APParameters(parameters, mean1, sigma1, amplitude1, mean2, sigma2, amplitude2);
  double ADCOffset = parameters[0];
  if (ADCOffset<0.001)
    ADCOffset = 0;
  double ADCEnergyRatio = (mean1 - ADCOffset) / 5.486;
  double energyADCRatio = 5.486 / (mean1-ADCOffset);
  double energyResolution = 100*parameters[1] * 2.354;
  double energyResolutionKeV = 0;
  couti << "Energy res. : " << energyResolution << " % (" << energyResolutionKeV << ")" << endl;
  couti << "Energy/ADC. : " << energyADCRatio << endl;
  couti << "ADC offset  : " << ADCOffset << endl;
  couti << "AMS(5.486)  : " << amplitude1 << ", " << mean1 << ", " << sigma1 << endl;
  couti << "AMS(5.443)  : " << amplitude2 << ", " << mean2 << ", " << sigma2 << endl;
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
      legend -> AddEntry("",Form("ER = %0.2f %s", energyResolution, percent),"");
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
      legend -> AddEntry("",Form("Resolution = %0.2f %s", energyResolution, percent),"");
      if (ADCOffset>0.001) legend -> AddEntry("",Form("ADC offset = %.3f", ADCOffset),"");
      legend -> AddEntry(fFitAlpha1,Form("M(5.486) = %0.2f", mean1),"l");
      legend -> AddEntry(fFitAlpha2,Form("M(5.443) = %0.2f", mean2),"l");
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
  fCountMultHit[0] = 0;
  fCountMultHit[1] = 0;
  fCountMultHit[2] = 0;
  fCountMultHit[3] = 0;
  fCountMultHit[4] = 0;

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
  }

  fHistChCount = new TH1D("histChCount",fRunName+" channel count;global-ch;event count", fNumCh,0,fNumCh);
  fHistChCount -> SetFillColor(29);
  fHistEnergy = new TH1D("histEnergy",fRunName+" energy (all channels);energy",fNumE,0,fMaxE);
  fHistEnergy -> SetFillColor(29);
  fHistEVSCh = new TH2D("histEVSCh", fRunName+" energy vs channel-id;global-ch;energy", fNumCh,0,fNumCh,fNumE,0,fMaxE);

  if (fUpdateDrawingEveryNEvent>0) {
    SetAttribute(fHistChCount,fCvsOnline->cd(1));
    SetAttribute(fHistEnergy,fCvsOnline->cd(2));
    SetAttribute(fHistEVSCh,fCvsOnline->cd(3),1,true);
  }

  UpdateCvsOnline(true);
}

void Analysis::SetAttribute(TH1* hist, TVirtualPad* pad, int npad, bool td)
{
  auto ax = hist -> GetXaxis();
  auto ay = hist -> GetYaxis();

  if (npad>=16)
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
  //(new TParameter<int>("run", fRunNo)) -> Write();
  (new TNamed("run", Form("%d",fRunNo))) -> Write();
  if (option==0) {
    //(new TParameter<int>("Energy threshold",   fEnergyThreshold)) -> Write();
    //(new TParameter<int>("Draw every n-event", fUpdateDrawingEveryNEvent)) -> Write();
    //(new TParameter<int>("Return if no file",  fReturnIfNoFile)) -> Write();
    //(new TParameter<int>("Ignore file update", fIgnoreFileUpdate)) -> Write();
    //(new TParameter<int>("Ignore TS decrease", fSkipTSError)) -> Write();
    //(new TParameter<int>("Event count limit",  fEventCountLimit)) -> Write();
    //(new TParameter<int>("File number range",  fFileNumberRange1,fFileNumberRange2)) -> Write();
    (new TNamed("Energy threshold   :", Form("%d",fEnergyThreshold))) -> Write();
    (new TNamed("Draw every n-event :", Form("%d",fUpdateDrawingEveryNEvent))) -> Write();
    (new TNamed("Return if no file  :", Form("%d",fReturnIfNoFile))) -> Write();
    (new TNamed("Ignore file update :", Form("%d",fIgnoreFileUpdate))) -> Write();
    (new TNamed("Ignore TS decrease :", Form("%d",fSkipTSError))) -> Write();
    (new TNamed("Event count limit  :", Form("%d",fEventCountLimit))) -> Write();
    (new TNamed("File number range  :", Form("%d %d",fFileNumberRange1,fFileNumberRange2))) -> Write();
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
  //fTreeOut -> Branch("tsDist",&bTimeStampDist);
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
  UShort_t module, channelID, energy;
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
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &module, &channelID, &energy, &timeStampLocal, &tsGroup);
      if (numData != 5) {
        coute << TString::Format("Number of data in line is not 5 (%d) %u %u %u %lld %u", numData, module, channelID, energy, timeStampLocal, tsGroup) << endl;
        continue;
      }

      if (energy < fEnergyThreshold)
        continue;

      timeStamp = timeStampLocal;
      if (tsGroup>0)
        timeStamp += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      bool tsError = false;
      bool nextEvent = false;
      bool sameEvent = false;
      if (fSkipTSError)
      {
        if (timeStamp<fTimeStampPrev) nextEvent = true;
        if (timeStamp>fTimeStampPrev) nextEvent = true;
        if (timeStamp==fTimeStampPrev) sameEvent = true;
      }
      else
      {
        if (timeStamp<fTimeStampPrev) tsError = true;
        if (timeStamp>fTimeStampPrev) nextEvent = true;
        if (timeStamp==fTimeStampPrev) sameEvent = true;
      }

      if (tsError) // time-stamp decreased
      {
        if (timeStamp<fTimeStampDecreased) {
          coutt << "Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          fTimeStampDecreased = timeStamp;
          ++fCountTSError;
        }
        if (fStopAtTSError) {
          coute << "Return due to Time-stamp error (" << fCountTSError << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          return;
        }
        continue;
      }
      else if (sameEvent) // same event
      {
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);
      }
      else if (nextEvent) // next event
      {
        if (countEventsSingleFile>0)
          FillDataTree();
        if (fExitAnalysis)
          break;

        fCountChannels = 0;
        fChannelArray -> Clear("C");
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);

        fTimeStampPrev = timeStamp;
        fTimeStampDecreased = timeStamp;
        fCountEvents++;
        countEventsSingleFile++;
      }

      int gid = GetGlobalID(module, channelID);
      channelData -> SetData(module,channelID,gid,energy,timeStampLocal,tsGroup,timeStamp);
      fCountAllChannels++;

      if (fUpdateDrawingEveryNEvent>=0)
      {
        fHistChCount -> Fill(gid);
        fHistEnergy -> Fill(energy);
        fHistEVSCh -> Fill(gid, energy);
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
    couti << "Number of events with 0   coincidence channels: " << fCountMultHit[0] << endl;
    couti << "Number of events with 1   coincidence channels: " << fCountMultHit[1] << endl;
    couti << "Number of events with 2   coincidence channels: " << fCountMultHit[2] << endl;
    couti << "Number of events with 3   coincidence channels: " << fCountMultHit[3] << endl;
    couti << "Number of events with =>4 coincidence channels: " << fCountMultHit[4] << endl;

    fileIn.close();
  }

  UpdateCvsOnline();
}

bool Analysis::FillDataTree()
{
  if      (fCountChannels==0) fCountMultHit[0]++;
  else if (fCountChannels==1) fCountMultHit[1]++;
  else if (fCountChannels==2) fCountMultHit[2]++;
  else if (fCountChannels==3) fCountMultHit[3]++;
  else if (fCountChannels>=4) fCountMultHit[4]++;
  if (bTimeStamp<0) bTimeStampDist = -1;
  else  bTimeStampDist = fTimeStampPrev - bTimeStamp;
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
        auto line = new TLine(iModule*16,0,iModule*16,yMax);
        line -> SetLineColor(kBlack);
        line -> Draw("samel");
      }
      for (int iDiv : {1,2,3}) {
        auto line = new TLine(iModule*16+iDiv*4,0,iModule*16+iDiv*4,yMax*0.5);
        line -> SetLineColor(kGray+1);
        line -> SetLineStyle(2);
        line -> Draw("samel");
      }
      auto tt = new TText((iModule+0.5)*16,yMax*0.1,Form("%d",GetRealModuleNumberFromIdx(iModule)));
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

  fCvsOnline -> cd(1);
  fHistChCount -> Draw();
  DrawModuleBoundary(fHistChCount->GetMaximum()*1.05);
  TakeCareOfStatsBox(fHistChCount);

  fCvsOnline -> cd(2);
  fHistEnergy -> Draw();

  fCvsOnline -> cd(3);
  fHistEVSCh -> Draw("colz");
  DrawModuleBoundary(8200);
  TakeCareOfStatsBox(fHistEVSCh);
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
  fHistEnergy -> Write();
  fHistEVSCh -> Write();

  cout << endl;
  couti << "End of conversion!" << endl;
  couti << "Number of events: " << fCountEvents << endl;
  couti << "Number of all channels: " << fCountAllChannels << endl;
  couti << "Number of times TS-error occured: " << fCountTSError << endl;
  couti << "Number of events with 0   coincidence channels: " << fCountMultHit[0] << endl;
  couti << "Number of events with 1   coincidence channels: " << fCountMultHit[1] << endl;
  couti << "Number of events with 2   coincidence channels: " << fCountMultHit[2] << endl;
  couti << "Number of events with 3   coincidence channels: " << fCountMultHit[3] << endl;
  couti << "Number of events with =>4 coincidence channels: " << fCountMultHit[4] << endl;
  couti << "Output file name: " << fFileNameOut << endl;

  if (fExitRoot)
    gApplication -> Terminate();
}

int Analysis::GetRealModuleNumberFromIdx(int iModule)
{
  if (iModule>2) return (iModule + 6);
  return iModule;
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
