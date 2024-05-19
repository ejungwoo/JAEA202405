void calibrate_alpha(int runNo=11)
{
  bool drawAnalysis = true;
  TString fileName = "dummy.root";
  //TString fileName = ""; // for creating auto-named calibration file

  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNo));

  ana -> AnalyzeAlphaTestModule(0,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(1,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(2,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(3,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(4,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(5,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(6,drawAnalysis,fileName);
}
