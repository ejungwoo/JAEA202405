#include "AnaCC.h"
#include "TList.h"
#include "TSystemDirectory.h"

#include <fstream>
#include <iostream>
using namespace std;

void AnaCC::Run(int runNo, TString pathToInputFile)
{
  fRunNo = runNo; 
  fPathToOutput = "/home/daquser/data/LiCD2Irrad/SortSi/test_data/";
  if (pathToInputFile.IsNull()==false)
    fPathToInput = pathToInputFile;
  fPathToOutput = "/home/daquser/data/LiCD2Irrad/SortSi/out/";
  fNumModules = 7;
  fNumChannels = 16;
  fChannelArray = NULL;
  fCountChannels = 0;

  ConfigureDateTime();
  MakeOutputFile();
  ReadDataFile();
  //FillHistogram();
  //EndOfAnalysis();
}

void AnaCC::ConfigureDateTime()
{
  cout << endl;
  cout << Form("==Looking for RUN%03d in %s",fRunNo,fPathToInput.Data()) << endl;
  cout << "..." << endl;

  TList *listOfFiles = TSystemDirectory(fPathToInput,fPathToInput).GetListOfFiles();
  TIter next(listOfFiles);
  TSystemFile* fileObject;
  int fDivisionMax = 0;
  while ((fileObject=(TSystemFile*)next()))
  {
    int idx = 0;
    if (fileObject->IsDirectory()) continue;
    TString fileName = fileObject -> GetName();
    if (fileName.EndsWith(".dat")==false) continue;
    fDateTime = fileName(7,12);
    int division = TString(fileName(25,3)).Atoi();
    if (division>fDivisionMax)
      fDivisionMax = division;
    cout << fileName << endl;
    fListOfInputFiles.push_back(fileName);
  }
}

void AnaCC::MakeOutputFile()
{
  fFileNameOut = fPathToOutput+Form("RUN%03d.ana.root",fRunNo);
  cout << endl;
  cout << "== Creating output file " << fFileNameOut << endl;
  fFileOut = new TFile(fFileNameOut,"recreate");
  fTreeOut = new TTree("event","");
  fChannelArray = new TClonesArray("ChannelData",10);
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
  Long64_t timeStamp, timeStampFull;
  Int_t numData = 0;
  char buffer[256];
  int countOpenFileTrials = 0;

  //for (int i=0; i<fListOfInputFiles.size(); i++)
  //{
  //  TString fileNameInput = fListOfInputFiles.at(i);
  //  cout << "Reading " << fileNameInput << endl;
  //}

  while(1)
  {
    if (fDivisionMax >=0 && countInputs>fDivisionMax) {
      cout << "Number of inputs " << countInputs << " exceed maximum input number " << fDivisionMax << endl;
      break;
    }

    TString fileNameInput = TString::Format("%s/RUN%03d_%s_list_%03d.dat", fPathToInput.Data(),fRunNo,fDateTime.Data(),countInputs);
    cout << endl;
    cout << "== Reading " << fileNameInput << endl;

    while (1) // endless loop (2) to check the status of the Raw-Data file //
    {
      countOpenFileTrials = 0;
      fileIn.open(fileNameInput);
      if (fileIn.fail())
      {
        fileIn.close();
        if(countOpenFileTrials > 10) {
          cout << "Failed to open file!" << endl;
          break;
        }
        else{
          countOpenFileTrials++;
          cout << "waiting(" << countOpenFileTrials << ") 10s ..." << endl;
          sleep(10);
          continue;
        }
      }

      fileSize = fileIn.seekg(0, ios::end).tellg(); // obtain the size of file
      fileIn.seekg(0, ios::beg); // rewind

      if (fileSize>fileSizeMax || (countOpenFileTrials!=0 && fileSize==fileSizeOld))
      // good file or final file
      {
        countOpenFileTrials = 0;
        countInputs++;
        break;
      }
      else if(countOpenFileTrials==0){ // first try
        countOpenFileTrials++;
        fileIn.close();
        fileSizeOld = fileSize;
        cout << "waiting(" << countOpenFileTrials << ") 10s ..." << endl;
        sleep(10);
      }
      else if(fileSize > fileSizeOld){ // writing the file is still continued ...
        fileIn.close();
        fileSizeOld = fileSize;
        cout << "waiting(" << countOpenFileTrials << ") 60s ..." << endl;
        sleep(60);
      }
    }

    if (countOpenFileTrials != 0) {
      cout << fileNameInput << " is Not Found!" << endl;
      break;
    }

    cout << endl;
    if (countInputs == 1) cout << "== Start reading first file!" << endl;
    else cout << "== Start reading" << endl;

    fChannelArray -> Clear("C");
    fCountChannels = 0;
    while (FileIn >> buffer)
    {
      numData = (Int_t) sscanf(buffer, "%hu,%hu,%hu,%lld,%hu", &module, &channel, &energy, &timeStamp, &tsGroup);
      if (numData != 5) {
        cout << TString::Format(" Warning: # of data in line is not 5 (%d) %u %u %u %lld %u", numData, module, channel, energy, timeStamp, tsGroup) << endl;
         break;
      }

      timeStampFull = timeStamp;
      if (tsGroup>0)
        timeStampFull += (Long64_t)(tsGroup) * 0xFFFFFFFFFF; 

      if (timeStampFull<fTimeStampPrev)
      {
        cout << "!!!!!!!! TIME STAMP DECREASED !!!!!!!!" << endl;
        cout << "Exiting from current input file" << endl;
        continue;
      }
      else if (timeStampFull>fTimeStampPrev) // next event
      {
        if (fCountEvents>0)
          fTreeOut -> Fill();
        fChannelArray -> Clear("C");
        fCountChannels = 0;
        fTimeStampPrev = timeStampFull;
        fCountEvents++;
      }

      ChannelData* data = (ChannelData*) fChannelArray -> ConstructedAt(fCountChannels);
      data -> SetData(module,channel,energy,timeStamp,tsGroup,timeStampFull);
      ++fCountChannels;
    }
  }
}

void AnaCC::FillHistogram()
{
}

void AnaCC::EndOfAnalysis()
{
}
