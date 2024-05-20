void calibrate_alpha(int runNo=-1)
{
  if (runNo<0 && gApplication->Argc()>=5 && TString(gApplication->Argv()[4]).IsDec())
    runNo = TString(gApplication->Argv()[4]).Atoi();
  if (runNo<0) { cout << "runNo is not given!" << endl; return; }

  bool drawAnalysis = true;

  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNo));

  if (runNo==11) {
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
  }
  else if (runNo==13) {
    ana -> SetNumADC(500);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
  else if (runNo<=15) {
    ana -> SetNumADC(500);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
  else {
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
}
