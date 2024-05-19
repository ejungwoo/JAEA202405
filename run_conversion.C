void run_conversion(int runNo=14)
{
  auto ana = new Analysis();

  ana -> SetReturnIfNoFile(true); // If true, Do not wait for data file to apear
  ana -> SetIgnoreFileUpdate(false); // If true, do not wait for data file to finish
  ana -> AddEnergyCalibration("out/RUN011.alpha.root"); // (Si) output file from calibrate_energy.C
  ana -> AddEnergyCalibration("out/RUN013.alpha.root"); // (dE) output file from calibrate_energy.C
  ana -> SetShowEnergyConversion(true); // if energy calibration file is good, user can chose this option

  ana -> SetSkipTSError(false); // skip time-stamp error
  ana -> SetStopAtTSError(false); // stop program if time-stamp occur
  ana -> SetADCThreshold(10); // do not write data below ADC threshold
  //ana -> SetDrawOnline(50000); // update online canvas with given event number
  ana -> SetDrawOnline(1000000); // update online canvas with given event number
  //ana -> SetAutoUpdateDrawing(true); // do not stop to ask continue
  //ana -> SetdES1Coincidence(true);
  //ana -> SetCoincidenceTSRange(1); // time-stamp range to reconstruct coincidence
  //ana -> SetCoincidenceMult(2); // only write event with given number of coincidence channels in event

  //ana -> SetEventCountLimit(100); // read event until given event count limit
  //ana -> SetFileNumberRange(0,10); // set range of file number

  ana -> SetConversionFile("dummy.root");
  ana -> RunConversion(runNo,"/home/daquser/data/LiCD2Irrad/");
}
