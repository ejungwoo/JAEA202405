void run_conversion(int runNo=-1)
{
  if (runNo<0 && gApplication->Argc()>=5 && TString(gApplication->Argv()[4]).IsDec())
    runNo = TString(gApplication->Argv()[4]).Atoi();
  if (runNo<0) { cout << "runNo is not given!" << endl; return; }

  auto ana = new Analysis();

  ana -> SetReturnIfNoFile(true);
  ana -> SetIgnoreFileUpdate(false);
  ana -> SetAutoUpdateDrawing(true);
  ana -> SetAutoUpdateRun(false);
  ana -> SetDrawOnline(50000);
  //ana -> SetDrawOnline(1000000);

  ana -> AddAlphaCalibrationFile("out/RUN011.alpha.root");
  ana -> AddAlphaCalibrationFile("out/RUN015.alpha.root");
  ana -> SetShowEnergyConversion(false);

  ana -> SetSkipTSError(false);
  ana -> SetStopAtTSError(false);
  ana -> SetADCThreshold(500);

  //ana -> SetCoincidenceTSRange(1);
  //ana -> SetTritonCutFile("out/tritonCutG.root");
  //ana -> SetCoincidenceMultRange(2,3);
  //ana -> SetdES1Coincidence(true);
  //ana -> SetdES1S3Coincidence(true);

  //ana -> SetEventCountLimit(100);
  //ana -> SetFileNumberRange(0,10);

  //ana -> SetConversionFile("dummy.root");
  ana -> RunConversion(runNo,"/home/daquser/data/LiCD2Irrad/");
}
