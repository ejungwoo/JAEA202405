DetectorSetting *getDet() { return DetectorSetting::GetDetectorSetting(); }

void setRun(int runNo) { Analysis::GetAnalysis()->SetRunNo(runNo); }
Analysis *getAna() { return Analysis::GetAnalysis(); }
void makeCutG(TString name) { Analysis::MakeCutGFile(name); }
void makeCutG(int pdt) { Analysis::MakeCutGFile(pdt); }
void callCuts() {
    Analysis::CallCutGFile(1);
    Analysis::CallCutGFile(2);
    Analysis::CallCutGFile(3);
}
