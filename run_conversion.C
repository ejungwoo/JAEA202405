void run_conversion(int runNo=11)
{
  auto ana = new Analysis();
  ana -> SetReturnIfNoFile(true); // If true, Do not wait for data file to apear
  ana -> SetIgnoreFileUpdate(false); // If true, do not wait for data file to finish

  //ana -> SetEventCountLimit(1000);
  ana -> SetSkipTSError(true);
  ana -> SetStopAtTSError(true);
  //ana -> SetFileNumberRange(0,10);
  ana -> SetEnergyThreshold(500);
  ana -> SetDrawOnline(200000);
  ana -> SetAutoUpdateDrawing(true);
  ana -> RunConversion(runNo,"/home/daquser/data/LiCD2Irrad/");
}
