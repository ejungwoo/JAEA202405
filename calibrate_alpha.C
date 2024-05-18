void calibrate_alpha(int runNo=11)
{
  bool drawAnalysis = true;
  TString fileName = "dummy.root";
  //TString fileName = ""; // for creating auto-named calibration file

  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNo));

  ana -> AnalyzeAlphaTestModule( 9,drawAnalysis,fileName);
  return;
  ana -> AnalyzeAlphaTestModule( 1,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule( 2,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule( 9,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(10,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(11,drawAnalysis,fileName);
  ana -> AnalyzeAlphaTestModule(12,drawAnalysis,fileName);
}
