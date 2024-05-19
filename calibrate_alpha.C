void calibrate_alpha(int runNo=13)
{
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
  else {
    ana -> SetAlphaTestFile("dummy.root");
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
  }
}
