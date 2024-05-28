void calibrate_alpha(int runNo=515)
{
  bool drawAnalysis = true;

  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNo));

  if (runNo==511) {
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
  }
  else if (runNo==513) {
    ana -> SetNumADC(500);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
  else if (runNo<=515) {
    ana -> SetNumADC(500);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
  else {
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    //ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
  }
}
