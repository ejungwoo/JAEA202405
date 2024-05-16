void run_conversion(int runNumber)
{
  auto ana = new Analysis();
  ana -> SetReturnIfNoFile(true); // If true, Do not wait for data file to apear
  ana -> SetIgnoreFileUpdate(false); // If true, do not wait for data file to finish

  ana -> RunConversion(runNumber,"/home/daquser/data/LiCD2Irrad/");
  //ana -> RunConversionOnline(3,"/home/daquser/data/LiCD2Irrad/");
}
