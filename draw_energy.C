void draw_energy(int runNumber=3)
{
  auto ana = new Analysis();

  //ana -> InitializeDrawing();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNumber));
  ana -> AnalyzeChannelEnergy(1,true);
}
