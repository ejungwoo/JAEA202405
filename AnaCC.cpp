#define coutd cout<<"+\033[0;36m"<<__LINE__<<" "<<__FILE__<<" #\033[0m "
#define couti cout<<"\033[0;32m==\033[0m "
#define coutw cout<<"\033[0;33mww\033[0m "
#define coute cout<<"\033[0;31mee\033[0m "
//#define coutt cout<<"\033[0;33mts-error\033[0m "
#define coutt coutx

#include "AnaCC.h"
#include "TList.h"
#include "TSystemDirectory.h"

void AnaCC::RunConversion(int runNo, TString pathIn)
{
  fRunNo = runNo; 
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/SortSi/test_data/") : pathIn);

  InitializeConversion();
  ConfigureDateTime();
  MakeOutputFile();
  ReadDataFile();
  EndOfConversion();
}

void AnaCC::RunConversionAndDraw(int runNo, TString pathIn)
{
  fRunNo = runNo; 
  fPathToInput = (pathIn.IsNull() ? TString("/home/daquser/data/LiCD2Irrad/SortSi/test_data/") : pathIn);

  InitializeDrawing(); //
  InitializeConversion();
  ConfigureDateTime();
  MakeOutputFile();
  ReadDataFile();
  EndOfConversion();
  Draw(); //
}

void AnaCC::InitializeConversion()
{
  fPathToOutput = "/home/daquser/data/LiCD2Irrad/SortSi/out/";
  fCountEvents = 0;
  fChannelArray = nullptr;
  fCountChannels = 0;
  fCountAllChannels = 0;
  fCountTimeStampDecrease = 0;
  fReturnIfNoFile = true;
  fCountCoincidence[0] = 0;
  fCountCoincidence[1] = 0;
  fCountCoincidence[2] = 0;
  fCountCoincidence[3] = 0;
  fCountCoincidence[4] = 0;
}

void AnaCC::InitializeDrawing()
{
  fFlagDrawing = true;
  fHistDTS    = new TH1D("fHistDTS",   "difference between TS",fNumDTS,0,fMaxDTS);
  fHistEnergy = new TH1D("fHistEnergy","energy",               fNumE,0,fMaxE);
  fHistEVSCh  = new TH2D("fHistEVSCh", "energy vs channel-id", fNumE,0,fMaxE,fNumCh,0,fNumCh);
}

void AnaCC::ConfigureDateTime()
{
  cout << endl;
  couti << Form("Looking for RUN%03d in %s",fRunNo,fPathToInput.Data()) << " ..." << endl;
  cout << endl;

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

void AnaCC::MakeOutputFile()
{
  fFileNameOut = fPathToOutput+Form("RUN%03d.ana.root",fRunNo);
  cout << endl;
  couti << "Creating output file " << fFileNameOut << endl;
  fFileOut = new TFile(fFileNameOut,"recreate");
  fTreeOut = new TTree("event","");
  fChannelArray = new TClonesArray("ChannelData",10);
  fTreeOut -> Branch("ts",&bTimeStamp);
  fTreeOut -> Branch("tsDist",&bTimeStampDist);
  fTreeOut -> Branch("ch",&fChannelArray);
}

void AnaCC::ReadDataFile()
{
  ifstream fileIn;
  streamsize fileSize;
  streamsize fileSizeOld = 0; // Size of opened Raw-Data file
  const streamsize fileSizeMax = 500000000; // 500 MB
  UShort_t countInputs = 0; // Number of Files read
  UShort_t module, channel, energy, tsGroup;
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
          coute << "There is no file!" << endl;
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
        fileIn.close();
        fileSizeOld = fileSize;
        couti << "waiting(" << countOpenFileTrials << ") 60s ..." << endl;
        sleep(60);
      }
    }

    if (countOpenFileTrials != 0) {
      coute << fileNameInput << " is not found! exit." << endl;
      break;
    }

    cout << endl;
    couti << "File size is " << fileSize << endl;

    fChannelArray -> Clear("C");
    Long64_t countEventsSingleFile = 0;
    fCountChannels = 0;
    while (fileIn >> buffer)
    {
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &module, &channel, &energy, &timeStampLocal, &tsGroup);
      if (numData != 5) {
        coute << TString::Format("Number of data in line is not 5 (%d) %u %u %u %lld %u", numData, module, channel, energy, timeStampLocal, tsGroup) << endl;
        continue;
      }

      timeStamp = timeStampLocal;
      if (tsGroup>0)
        timeStamp += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      if (timeStamp>fTimeStampPrev) // next event
      {
        if (countEventsSingleFile>0)
        {
          if (fCountChannels==0) {
            fCountCoincidence[0]++;
            coutd << fCountEvents << endl;
            coutd << fTimeStampPrev << endl;
            coutd << fCountChannels << endl;
          }
          else if (fCountChannels==1) fCountCoincidence[1]++;
          else if (fCountChannels==2) fCountCoincidence[2]++;
          else if (fCountChannels==3) fCountCoincidence[3]++;
          else if (fCountChannels>=4) fCountCoincidence[4]++;
          if (bTimeStamp<0) bTimeStampDist = -1;
          else  bTimeStampDist = fTimeStampPrev - bTimeStamp;
          bTimeStamp = fTimeStampPrev;
          fTreeOut -> Fill();
        }

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
      else if (timeStamp<fTimeStampPrev) // time stamp decreased
      {
        if (timeStamp<fTimeStampDecreased) {
          coutt << "Time stamp decreased (" << fCountTimeStampDecrease << ") " << fTimeStampPrev << " -> " << timeStamp << endl;
          fTimeStampDecreased = timeStamp;
          ++fCountTimeStampDecrease;
        }
        continue;
      }

      channelData -> SetData(module,channel,energy,timeStampLocal,tsGroup,timeStamp);
      fCountAllChannels++;
    }

    if (fCountEvents>0) {
      if      (fCountChannels==0) fCountCoincidence[0]++;
      else if (fCountChannels==1) fCountCoincidence[1]++;
      else if (fCountChannels==2) fCountCoincidence[2]++;
      else if (fCountChannels==3) fCountCoincidence[3]++;
      else if (fCountChannels>=4) fCountCoincidence[4]++;
      bTimeStamp = fTimeStampPrev;
      fTreeOut -> Fill();
    }

    couti << "Number of channels in current file: " << countEventsSingleFile << endl;
    fileIn.close();
  }
}

void AnaCC::EndOfConversion()
{
  fFileOut -> cd();
  fTreeOut -> Write();

  cout << endl;
  couti << "End of conversion!" << endl;
  couti << "Number of events: " << fCountEvents << endl;
  couti << "Number of all channels: " << fCountAllChannels << endl;
  couti << "Number of times time stamp decreased: " << fCountTimeStampDecrease << endl;
  couti << "Number of events with 0   coincidence channels: " << fCountCoincidence[0] << endl;
  couti << "Number of events with 1   coincidence channels: " << fCountCoincidence[1] << endl;
  couti << "Number of events with 2   coincidence channels: " << fCountCoincidence[2] << endl;
  couti << "Number of events with 3   coincidence channels: " << fCountCoincidence[3] << endl;
  couti << "Number of events with =>4 coincidence channels: " << fCountCoincidence[4] << endl;
  couti << "Output file name: " << fFileNameOut << endl;
}
