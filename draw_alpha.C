void draw_alpha(int runNo=11)
{
  bool drawAnalysis = true;
  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("out/RUN%03d.summary.root",runNo));

  //ana -> AnalyzeAlphaTest(0,14,drawAnalysis);

  ana -> AnalyzeAlphaTestModule( 0,drawAnalysis);
  return;
  ana -> AnalyzeAlphaTestModule( 1,drawAnalysis);
  ana -> AnalyzeAlphaTestModule( 2,drawAnalysis);
  ana -> AnalyzeAlphaTestModule( 9,drawAnalysis);
  ana -> AnalyzeAlphaTestModule(10,drawAnalysis);
  ana -> AnalyzeAlphaTestModule(11,drawAnalysis);
  ana -> AnalyzeAlphaTestModule(12,drawAnalysis);
}
