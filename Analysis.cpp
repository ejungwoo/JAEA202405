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
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/analysis/input/") : pathIn);

  InitializeConversion();
  ConfigureDateTime();
  MakeOutputFile();
  ReadDataFile();
  EndOfConversion();
}

void Analysis::RunConversionOnline(int runNo, TString pathIn)
{
  fRunNo = runNo; 
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/analysis/input/") : pathIn);

  InitializeConversion();
  ConfigureDateTime();
  MakeOutputFile();
  InitializeDrawing(); //
  ReadDataFile();
  EndOfConversion();
}

void Analysis::ReadSummaryFile(TString fileName)
{
  fFileSummary = new TFile(fileName,"read");
  fTreeSummary = (TTree*) fFileSummary -> Get("event");
  fTreeSummary -> SetBranchAddress("ts",&bTimeStamp);
  fTreeSummary -> SetBranchAddress("tsDist",&bTimeStampDist);
  fTreeSummary -> SetBranchAddress("channel",&fChannelArray);
  fNumEventsSummary = fTreeSummary -> GetEntries();

  couti << fileName << " containing " << fNumEventsSummary << " events" << endl;
}

void Analysis::AnalyzeChannelEnergy(int module, int channelID, bool drawHist)
{
  auto gid = GetGlobalID(module, channelID);
  AnalyzeChannelEnergy(gid);
}

Double_t FxAlphaPeak2(Double_t *xy, Double_t *par)
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

  double ADCToEnergy = meanP1 / 5.486;

  double mean2 = 5.443 * ADCToEnergy + ADCOffset;
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

  double ADCToEnergy = meanP1 / 5.486;

  mean2 = 5.443 * ADCToEnergy + ADCOffset;
  double meanP2 = mean2 - ADCOffset; // pure ADC without ADC-offset
  sigma2 = energyResolution * meanP2;
  amplitude2 = amplitude1 * sigma1 / 6.65625 / sigma2; // 85.2 / 12.8 = 6.65625; 
}

void Analysis::AnalyzeChannelEnergy(int gid, bool drawHist)
{
  //TF1* fitAlpha = new TF1("AlphaSinglePeak","gaus(0)",0,fMaxE);
  TF1* fitAlpha = new TF1("AlphaTwoPeak",FxAlphaPeak2,0,fMaxE,4);
  TF1* fitAlpha1 = new TF1("Alpha1","gaus(0)",0,fMaxE);
  TF1* fitAlpha2 = new TF1("Alpha2","gaus(0)",0,fMaxE);
  fitAlpha1 -> SetLineColor(kBlue);
  fitAlpha2 -> SetLineColor(kGreen+1);
  double mean1, sigma1, amplitude1, mean2, sigma2, amplitude2;

  if (gid>=0)
  {
    TString nameHist = Form("histE%d",gid);
    auto histE = new TH1D(nameHist,Form("channel %d energy distribution;energy;count",gid),fNumE,0,fMaxE);
    fTreeSummary -> Project(nameHist, "channel.energy", Form("channel.gid==%d",gid));
    if (drawHist)
    {
      if (fCvsEnergy==nullptr)
        fCvsEnergy = new TCanvas("cvsEnergy","",800,450);
      double mean = histE -> GetMean();
      double sig = histE -> GetStdDev();

      double bin = histE->GetMaximumBin();
      double yMax = histE->GetBinContent(bin);
      double xPeak = histE->GetBinCenter(bin);

      fitAlpha -> SetParameters(mean, sig, yMax);
      auto x1 = xPeak - 100;
      auto x2 = xPeak + 100;
      fitAlpha -> SetParameters(0,0.003,xPeak,yMax);
      fitAlpha -> SetParLimits(0,0,100);
      histE -> Fit(fitAlpha,"Q0N","",x1,x2);
      auto parameters = fitAlpha -> GetParameters();
      Convert2APParameters(parameters, mean1, sigma1, amplitude1, mean2, sigma2, amplitude2);
      couti << "ADC Offset : " << parameters[0] << endl;
      couti << "Energy resolution : " << 100*parameters[1] * 2.354 << " % (" << parameters[1] << ")" << endl;
      couti << "AMS(5.486) : " << amplitude1 << ", " << mean1 << ", " << sigma1 << endl;
      couti << "AMS(5.443) : " << amplitude2 << ", " << mean2 << ", " << sigma2 << endl;
      fitAlpha1 -> SetParameters(amplitude1, mean1, sigma1);
      fitAlpha2 -> SetParameters(amplitude2, mean2, sigma2);

      histE -> GetXaxis() -> SetRangeUser(mean1-5*sigma1,mean1+5*sigma1);
      histE -> Draw();
      fitAlpha -> Draw("samel");
      fitAlpha1 -> Draw("samel");
      fitAlpha2 -> Draw("samel");
    }
  }
}

void Analysis::InitializeConversion()
{
  fPathToOutput = "/home/daquser/data/LiCD2Irrad/analysis/out/";
  fCountEvents = 0;
  fChannelArray = nullptr;
  fCountChannels = 0;
  fCountAllChannels = 0;
  fCountTimeStampDecrease = 0;
  fReturnIfNoFile = true;
  fCountMultHit[0] = 0;
  fCountMultHit[1] = 0;
  fCountMultHit[2] = 0;
  fCountMultHit[3] = 0;
  fCountMultHit[4] = 0;
}

void Analysis::InitializeDrawing()
{
  if (fFileOut!=nullptr)
    fFileOut -> cd();

  fUpdateDrawingEveryNEvent = 100000;
  fCvsOnline = new TCanvas("cvsOnline","online update canvas",1600,900);
  fCvsOnline -> Divide(2,2);

  fHistChCount = new TH1D("histChCount","channel count;global-ch;event count", fNumCh,0,fNumCh);
  fHistDTS     = new TH1D("histDTS",   "difference between TS;dTS",fNumDTS,0,fMaxDTS);
  fHistEnergy  = new TH1D("histEnergy","energy (all channels);energy",fNumE,0,fMaxE);
  fHistEVSCh   = new TH2D("histEVSCh", "energy vs channel-id;global-ch;energy", fNumCh,0,fNumCh,fNumE,0,fMaxE);

  UpdateCvsOnline();
}

void Analysis::SetAttribute(TH1* hist, TPad* pad)
{
}

void Analysis::ConfigureDateTime()
{
  cout << endl;
  couti << Form("Looking for RUN%03d in %s",fRunNo,fPathToInput.Data()) << " ..." << endl;

  TList *listOfFiles = TSystemDirectory(fPathToInput,fPathToInput).GetListOfFiles();
  TIter next(listOfFiles);
  TSystemFile* fileObject;
  fDivisionMax = 10;
  while ((fileObject=(TSystemFile*)next()))
  {
    int idx = 0;
    if (fileObject->IsDirectory()) continue;
    TString fileName = fileObject -> GetName();
    if (fileName.Index(Form("RUN%03d",fRunNo))!=0) continue;
    if (fileName.EndsWith(".dat")==false) continue;
    fDateTime = fileName(7,12);
    int division = TString(fileName(25,3)).Atoi();
    if (division>fDivisionMax)
      fDivisionMax = division;
    couti << fileName << endl;
  }
}

void Analysis::MakeOutputFile()
{
  fFileNameOut = fPathToOutput+Form("RUN%03d.summary.root",fRunNo);
  cout << endl;
  couti << "Creating output file " << fFileNameOut << endl;
  fFileOut = new TFile(fFileNameOut,"recreate");
  fTreeOut = new TTree("event","");
  fChannelArray = new TClonesArray("ChannelData",10);
  fTreeOut -> Branch("ts",&bTimeStamp);
  fTreeOut -> Branch("tsDist",&bTimeStampDist);
  fTreeOut -> Branch("channel",&fChannelArray);
}

void Analysis::ReadDataFile()
{
  ifstream fileIn;
  streamsize fileSize;
  streamsize fileSizeOld = 0; // Size of opened Raw-Data file
  const streamsize fileSizeMax = 500000000; // 500 MB
  UShort_t countInputs = 0; // Number of Files read
  UShort_t module, channelID, energy, tsGroup;
  Long64_t timeStampLocal, timeStamp;
  Int_t numData = 0;
  char buffer[256];
  int countOpenFileTrials = 0;

  ChannelData* channelData = NULL;
  fTimeStampPrev = -1;

  while(1)
  {
    if (fDivisionMax >=0 && countInputs>fDivisionMax) {
      couti << "Number of inputs " << countInputs << " exceed maximum input number " << fDivisionMax << endl;
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
        couti << "waiting(" << countOpenFileTrials << ") 3s ..." << endl;
        sleep(3);
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

    couti << "File size is " << fileSize << endl;

    fChannelArray -> Clear("C");
    Long64_t countEventsSingleFile = 0;
    fCountChannels = 0;
    while (fileIn >> buffer)
    {
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &module, &channelID, &energy, &timeStampLocal, &tsGroup);
      if (numData != 5) {
        coute << TString::Format("Number of data in line is not 5 (%d) %u %u %u %lld %u", numData, module, channelID, energy, timeStampLocal, tsGroup) << endl;
        continue;
      }

      timeStamp = timeStampLocal;
      if (tsGroup>0)
        timeStamp += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      if (timeStamp>fTimeStampPrev) // next event
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
      else if (timeStamp==fTimeStampPrev) // same event
      {
        channelData = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels++);
      }
      else if (timeStamp<fTimeStampPrev) // time-stamp decreased
      {
        if (timeStamp<fTimeStampDecreased) {
          coutt << "Time-stamp decreased (" << fCountTimeStampDecrease << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          fTimeStampDecreased = timeStamp;
          ++fCountTimeStampDecrease;
        }
        continue;
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
    }

    if (!fExitAnalysis && fCountEvents>0) 
      FillDataTree();
    if (fExitAnalysis)
      break;

    couti << "Number of channels in current file: " << countEventsSingleFile << endl;
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
  fTreeOut -> Fill();

  if (fUpdateDrawingEveryNEvent>=0)
  {
    fHistDTS -> Fill(bTimeStampDist);
    fCountEventsForUpdate++;

    // fUpdateDrawingEveryNEvent==0 will read all events and draw
    if (fUpdateDrawingEveryNEvent!=0 && fCountEventsForUpdate>=fUpdateDrawingEveryNEvent)
    {
      UpdateCvsOnline();
      fCountEventsForUpdate = 0;
      cout << "\033[0;32m" << Form("== (%d) Enter / stop / all: ",fCountEvents,fUpdateDrawingEveryNEvent) << "\033[0m";
      std::string userInput0;
      std::getline(std::cin, userInput0);
      TString userInput = userInput0;
      userInput.ToLower();
      if (userInput=="stop" || userInput=="exit" || userInput=="x")
        fExitAnalysis = true;
      else if (userInput.Index(".q")==0 ) {
        fExitAnalysis = true;
        fExitRoot = true;
      }
      else
      {
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

void Analysis::UpdateCvsOnline()
{
  if (fCvsOnline!=nullptr)
  {
    fCvsOnline -> cd(1);
    fHistDTS -> Draw();
    fCvsOnline -> cd(2);
    fHistChCount -> Draw();
    fCvsOnline -> cd(3);
    fHistEnergy -> Draw();
    fCvsOnline -> cd(4);
    fHistEVSCh -> Draw("colz");
    fCvsOnline -> Modified();
    fCvsOnline -> Update();
  }
}

void Analysis::EndOfConversion()
{
  fFileOut -> cd();
  fTreeOut -> Write();

  if (fCvsOnline!=nullptr)
  {
    UpdateCvsOnline();
    fCvsOnline -> Write();
    fHistDTS -> Write();
    fHistChCount -> Write();
    fHistEnergy -> Write();
    fHistEVSCh -> Write();
  }

  cout << endl;
  couti << "End of conversion!" << endl;
  couti << "Number of events: " << fCountEvents << endl;
  couti << "Number of all channels: " << fCountAllChannels << endl;
  couti << "Number of times time-stamp decreased: " << fCountTimeStampDecrease << endl;
  couti << "Number of events with 0   coincidence channels: " << fCountMultHit[0] << endl;
  couti << "Number of events with 1   coincidence channels: " << fCountMultHit[1] << endl;
  couti << "Number of events with 2   coincidence channels: " << fCountMultHit[2] << endl;
  couti << "Number of events with 3   coincidence channels: " << fCountMultHit[3] << endl;
  couti << "Number of events with =>4 coincidence channels: " << fCountMultHit[4] << endl;
  couti << "Output file name: " << fFileNameOut << endl;

  if (fExitRoot)
    gApplication -> Terminate();
}

int Analysis::GetGlobalID(UShort_t module, UShort_t channelID)
{
  int globalID = module*fNumChannels + channelID;
  return globalID;
}

void Analysis::GetModCh(int globalID, UShort_t &module, UShort_t &channelID)
{
  module = globalID / fNumChannels;
  channelID = globalID % fNumChannels;
}
