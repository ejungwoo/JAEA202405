#include "jshort.h"

void run_alpha_calibration(int runNo=15)
{
  bool drawAnalysis = true;

  TString path_to_output = "data_converted";

  auto ana = new Analysis();
  ana -> ReadSummaryFile(Form("%s/RUN%03d.summary.root",path_to_output.Data(),runNo));
  ana -> SetPathToOutput("data_calibrated/");

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
    ana -> AnalyzeAlphaTestModule(0,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(1,drawAnalysis);
    //ana -> AnalyzeAlphaTestModule(2,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(3,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(4,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(5,drawAnalysis);
    ana -> AnalyzeAlphaTestModule(6,drawAnalysis);
  }
}
