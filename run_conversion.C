void run_conversion(int runNo=-1, double energy=-1)
{
  if (runNo<0 && gApplication->Argc()>=5 && TString(gApplication->Argv()[4]).IsDec())
    runNo = TString(gApplication->Argv()[4]).Atoi();
  if (runNo<0 && gApplication->Argc()>=6 && TString(gApplication->Argv()[5]).IsFloat())
    energy = TString(gApplication->Argv()[4]).Atof();
  if (runNo<0) { cout << "runNo is not given!" << endl; return; }

  auto ana = new Analysis();
  ana -> SetPathToInput("/home/daquser/data/LiCD2Irrad/");
  ana -> SetPathToOutput("/home/daquser/data/LiCD2Irrad/analysis/out/");
  ana -> SetMapFileName("/home/daquser/data/LiCD2Irrad/analysis/ModCh.in");
  ana -> SetDetFileName("/home/daquser/data/LiCD2Irrad/analysis/DetectorSetting.in");

  // general
  ana -> SetReturnIfNoFile(true);
  ana -> SetIgnoreFileUpdate(true);
  ana -> SetAutoUpdateDrawing(true);
  ana -> SetAutoUpdateRun(true);
  ana -> SetDrawOnline(400000);
  ana -> SetSkipTSError(false);
  ana -> SetStopAtTSError(false);
  ana -> SetADCThreshold(50);
  //ana -> SetEventCountLimit(1000000);
  //ana -> SetLineCountLimit(125);
  //ana -> SetFileNumberRange(0,10);
  //ana -> SetDrawOnline(1000000);

  // alpha calibration
  ana -> AddAlphaCalibrationFile("out/RUN011.alpha.root");
  ana -> AddAlphaCalibrationFile("out/RUN015.alpha.root");
  ana -> SetShowEnergyConversion(true);

  ana -> SetCoincidenceTSRange(15);
  //ana -> SetTritonCutFile("out/tritonCutG.root");
  //ana -> SetCoincidenceMultRange(2,3);
  //ana -> SetdES1Coincidence(true);
  //ana -> SetdES1S3Coincidence(true);

  // conditions for local histograms (which will not affect event collection)
  //ana -> SetLocalDetectorChannelCut(kS1J,10);
  ana -> SetLocalEnergyRange(4,6);
  ana -> AddLocalS1StripHist(1);
  ana -> AddLocalS1StripHist(4);
  ana -> AddLocalS1StripHist(8);
  //ana -> AddLocalS1StripHist(12);

  //ant -> SetExcludeDESECutGFile();
  ana -> SetEnergyCutGFile("out/cutG5.root"); // CD2
  //ana -> SetEnergyCutGFile("out/cutG7.root"); // CH2

  //ana -> SetConversionFile("dummy.root");
  //ana -> SetBeamEnergy(energy);
  ana -> RunConversion(runNo);
}
