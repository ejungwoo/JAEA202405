/*************************************************************
 * macro: SortSi.C                                           *
 * written by Jung Woo Lee                                   *
 * modified from below information                           *
 *************************************************************/

/*************************************************************
 * macro: SortSi02.C                                         *
 * written by H. Makii                                       *
 * Time-stamp: <2018-02-28 14:22:49 makii>                   *
 * check Time Stamps,  make histograms,                      *
 * search peak for Si detectors,                             *
 *  - single histograms for Si detectors                     *
 *  - (E-DE)-DE histograms                                   *
 * check gain stability for Si detectors as a function of Ts *
 *                                                           *
 * make TTree containing data for DE, E, and MWPC(Q,X,Y,dT), *
 *                        and Neutron Detectors (PH,PSD,TOF) *
 * Trigger : Si Deleta-E                                     *
 *                                                           *
 * Usage:                                                    *
 * SortSi02("RUN035",       RunName No.                      *
 *          "201503061604", Date and Time                    *
 *          -1)             Last File No.                    *
 *                                                           *
 * Raw-Data files : PathIn/RUN035_201503061604_list_???.dat  *
 *  if LastFile < 0, open all existing files,                *
 *   open 000 - (Last File) otherwise                        *
 *                                                           *
 * output files   : PathOut/RUN035.log                       *
 *                          RUN035_hist.root                 *
 *                          RUN035_tree.root                 *
 *                                                           *
 * This macro was tested by using ROOT 5.34/36               *
 *************************************************************/

//#define eout cout<<"+\033[0;36m"<<__LINE__<<"\033[0m "
#define eout cout

#include "SortSi.h"
#include "TSystemDirectory.h"
#include "TSystemDirectory.h"

void SortSi(int RunNumber=41);
void SortSi(TString RunName, TString DateTime, Int_t FLast);

/////////////////////////////////////////////////////////////
// Specify which E (Ea or Eb) should coincide with Delta-E //
/////////////////////////////////////////////////////////////
void InitUseDE(void){
  for(Int_t m=0; m<16; m++){
    TString stmp;
    stmp = sEaDEUse(m,1);
    EaDEUse[m] = (Bool_t)( stmp.Atoi() );
    stmp = sEbDEUse(m,1);
    EbDEUse[m] = (Bool_t)( stmp.Atoi() );
  }
}

////////////////////////////////////////////
// Dictionary for Ch assignment for MWPCs //
////////////////////////////////////////////
void InitChMWPC(void){
  TString stmp;
  for(Int_t m=0; m<4; m++){
    stmp = sUseChQ(m*2,2);
    ChQ[ stmp.Atoi() ] = m;
    stmp = sUseChPosX(m*2,2);
    ChX[ stmp.Atoi() ] = m;
    stmp = sUseChPosY(m*2,2);
    ChY[ stmp.Atoi() ] = m;
    if(m<3){
      stmp = sUseChDt(m*2,2);
      ChDt[ stmp.Atoi() ] = m;
    }
    stmp = sUseChMCP(m*2,2);
    ChMCP[ stmp.Atoi() ] = m;
  }
}

///////////////////////////////////////////////////////////
// Find the peak in TH2 *hin and Fill PID plot TH2 *hout //
///////////////////////////////////////////////////////////
void PeakPosition2(Bool_t fPID, TH2 *hin, TH2 *hout, Double_t &PeakEE, Double_t &PeakDE)
{
  Int_t     binx,biny,binz;
  Double_t  fmin,fmax;
  TF1      *fTmp;
  hin->GetXaxis()->SetRangeUser(2000,5000);
  hin->GetYaxis()->SetRangeUser(1200,3500);
  hin->GetMaximumBin(binx,biny,binz);
  PeakEE = hin->GetXaxis()->GetBinCenter(binx);
  PeakDE = hin->GetYaxis()->GetBinCenter(biny);
  //
  TH1D *pj_EE = hin->ProjectionX("pj_EE");
  fmin  = pj_EE->GetXaxis()->GetBinCenter(1);
  fmin *= 2.;
  fmax  = pj_EE->GetXaxis()->GetBinCenter( hin->GetNbinsX() );
  fTmp = new TF1("fTmp","gaus",fmin,fmax);
  fTmp->SetParameter(0, hin->GetMaximumBin() );
  fTmp->SetParameter(1, PeakEE);
  fTmp->SetParameter(2, 25.);
  pj_EE->Fit("fTmp","QN","",PeakEE-25.,PeakEE+25.);
  PeakEE = fTmp->GetParameter(1);
  //
  TH1D *pj_DE = hin->ProjectionY("pj_DE",binx-25,binx+25);
  fmin  = pj_DE->GetXaxis()->GetBinCenter(1);
  fmin *= 2.;
  fmax  = pj_DE->GetXaxis()->GetBinCenter( hin->GetNbinsX() );
  fTmp = new TF1("fTmp","gaus",fmin,fmax);
  fTmp->SetParameter(0, hin->GetMaximumBin() );
  fTmp->SetParameter(1, PeakDE);
  fTmp->SetParameter(2, 25.);
  pj_DE->Fit("fTmp","QN","",PeakDE-25.,PeakDE+25.);
  PeakDE = fTmp->GetParameter(1);
  //
  TPolyMarker *pm = (TPolyMarker*)hin->GetListOfFunctions()->FindObject("TPolyMarker");
  if(pm)
    hin->GetListOfFunctions()->Remove(pm); delete pm;
  pm = new TPolyMarker(1);
  hin->GetListOfFunctions()->Add(pm);
  pm ->SetPoint(0,PeakEE,PeakDE);
  pm ->SetMarkerStyle(23);
  pm ->SetMarkerColor(kBlue);
  pm ->SetMarkerSize(1.3);
  //
  delete fTmp;
  delete pj_EE;
  delete pj_DE;
  //
  hin->GetXaxis()->UnZoom();
  hin->GetYaxis()->UnZoom();
  if(fPID){
    Int_t    xnbins = hout->GetNbinsX();
    Double_t xlow   = hout->GetXaxis()->GetBinLowEdge(1);
    Double_t xup    = hout->GetXaxis()->GetBinUpEdge(xnbins);
    Int_t    ynbins = hout->GetNbinsY();
    Double_t ylow   = hout->GetYaxis()->GetBinLowEdge(1);
    Double_t yup    = hout->GetYaxis()->GetBinUpEdge(ynbins);
    Double_t bx,by,newbx,newby;
    vector<Double_t> x,y,w;
    for(Int_t   binx=1; binx<=hin->GetNbinsX(); binx++){
      for(Int_t biny=1; biny<=hin->GetNbinsY(); biny++){
        Double_t content = hin->GetBinContent(binx,biny);
        if(content == 0) continue;
        bx    = hin->GetXaxis()->GetBinCenter(binx);
        newbx = bx / PeakEE;
        if(newbx < xlow || newbx > xup) continue;
        by    = hin->GetYaxis()->GetBinCenter(biny);
        newby = by / PeakDE;
        if(newby < ylow || newby > yup) continue;
        x.push_back(newbx);
        y.push_back(newby);
        w.push_back(content);
      }
    }
    hout->FillN( (Int_t)x.size(), x.data(), y.data(), w.data(), 1 );
    x.clear();
    y.clear();
    w.clear();
  }
  //  else
  //    cout << Form(" %s is not contained in PID\n", hin->GetName() );
}

///////////////////////
// Print # of Errors //
///////////////////////
void PrintErrors(Bool_t fWrite, std::ofstream & fout,
		 Long64_t **err, Long64_t **event, Long64_t **err2)
{
  // fWrite = 1 : end of program, 0 : each Raw-Data file
  ErrR = new Double_t*[NMod];
  for(Int_t l=0; l<NMod; l++){
    ErrR[l] = new Double_t[16];
    for(Int_t m=0; m<16; m++)
      ErrR[l][m] = 0.;
  }
  if(fWrite){
    ///////////////////////////////////////////////////////
    // Error of Time Stamp (same Time Stamp) -> Terminal //
    ///////////////////////////////////////////////////////
    cout << " Error of Time Stamp (same Time Stamp)" << endl;
    for(Int_t l=0; l<NMod; l++){
      cout << TString::Format("  Mod. No. ... %02d\n   ", l);
      for(Int_t m=0; m<16; m++){
        if(err[l][m]>0)
          cout << TString::Format("\033[7m%5lld (%10lld)  \033[0m", err[l][m],event[l][m]);
        else
          cout << TString::Format("%5lld (%10lld)  ", err[l][m],event[l][m]);
        if( (m+1) % 4 == 0)
          cout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-3), m);
      }
      cout << endl;
    }
    ///////////////////////////////////////////////
    // (Error) / (Valid Event) ratio -> Log file //
    ///////////////////////////////////////////////
    cout << endl;
    fout << endl;
    fout << " Error of Time Stamp (same Time Stamp)" << endl;
    for(Int_t l=0; l<NMod; l++){
      fout << TString::Format("  Mod. No. ... %02d\n   ", l);
      for(Int_t m=0; m<16; m++){
        if(event[l][m] != 0){
          ErrR[l][m] = ((Double_t)err[l][m] /  (Double_t)event[l][m]*100.0);
          fout << TString::Format("%.4f  ", ErrR[l][m]);
        }
        else
          fout << TString::Format("*.****  ");
        if( (m+1) % 8 == 0)
          fout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-7), m);
      }
      fout << endl;
    }
    //////////////////////////////////////////////////////////////
    // Error of Time Stamp (decrease of Time Stamp) -> Log file //
    //////////////////////////////////////////////////////////////
    fout << endl << " Error of Time Stamp (decrease of Time Stamp)" << endl;
    for(Int_t l=0; l<NMod; l++){
      fout << TString::Format("  Mod. No. ... %02d\n   ", l);
      for(Int_t m=0; m<16; m++){
        fout << TString::Format("%5lld  ", err2[l][m]);
        if( (m+1) % 8 == 0)
          fout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-7), m);
      }
      fout << endl;
    }
  }
  ///////////////////////////////////////////////
  // (Error) / (Valid Event) ratio -> Terminal //
  ///////////////////////////////////////////////
  cout << " Error of Time Stamp (same Time Stamp)" << endl;
  for(Int_t l=0; l<NMod; l++){
    cout << TString::Format("  Mod. No. ... %02d\n   ", l);
    for(Int_t m=0; m<16; m++){
      if(event[l][m] != 0){
        ErrR[l][m] = ((Double_t)err[l][m] /  (Double_t)event[l][m]*100.0);
        if(err[l][m]>0)
          cout << TString::Format("\033[7m%.4f  \033[0m", ErrR[l][m]);
        else
          cout << TString::Format("%.4f  ", ErrR[l][m]);
        if( (ErrR[l][m] > MinSameErr || err2[l][m] > MinDecErr) && HName[l][m] != "" )
          ErrFlag++;
      }
      else
        cout << TString::Format("*.****  ");
      if( (m+1) % 8 == 0)
        cout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-7), m);
    }
    cout << endl;
  }
  if(fWrite){
    //////////////////////////////////////////////////////////////
    // Error of Time Stamp (decrease of Time Stamp) -> Terminal //
    //////////////////////////////////////////////////////////////
    cout << endl << " Error of Time Stamp (decrease of Time Stamp)" << endl;
    for(Int_t l=0; l<NMod; l++){
      cout << TString::Format("  Mod. No. ... %02d\n   ", l);
      for(Int_t m=0; m<16; m++){
        if(err2[l][m]>0)
          cout << TString::Format("\033[7m%5lld  \033[0m", err2[l][m]);
        else
          cout << TString::Format("%5lld  ", err2[l][m]);
        if( (m+1) % 8 == 0)
          cout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-7), m);
      }
      cout << endl;
    }
  }
  for(Int_t l=0; l<NMod; l++){
    delete [] ErrR[l]; ErrR[l] = 0;
  }
  delete [] ErrR; ErrR = 0;
}

///////////////////////////////////////
// Print # of Events stored in TTree //
///////////////////////////////////////
void PrintEventTree(std::ofstream & fout,
		    Long64_t eventTrig, Long64_t eventStored, Long64_t **eventTree, Long64_t **event)
{
  cout << endl << " Number of Events stored in TTree (Total Event)" << endl;
  cout << TString::Format("  -- Number of Trigger : %lld -- \n", eventTrig);
  fout << endl << " Number of Events stored in TTree (Total Event)" << endl;
  fout << TString::Format("  -- Number of Trigger : %lld -- \n", eventTrig);
  for(Int_t l=0; l<NMod; l++){
    cout << TString::Format("  Mod. No. ... %02d\n   ", l);
    fout << TString::Format("  Mod. No. ... %02d\n   ", l);
    for(Int_t m=0; m<16; m++){
      cout << TString::Format("%8lld (%10lld)  ", eventTree[l][m],event[l][m]);
      fout << TString::Format("%8lld (%10lld)  ", eventTree[l][m],event[l][m]);
      if( (m+1) % 4 == 0){
        cout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-3), m);
        fout << TString::Format("(for Ch = %02d - %02d )\n   ",(m-3), m);
      }
    }
    cout << endl;
    fout << endl;
  }
  cout << "  Require conincidence between Si Delta-E and Si E\n"
       << TString::Format("   -> %lld events were stored (%.3f of total) --\n\n",
        eventStored, (Double_t)eventStored/(Double_t)eventTrig*100.0);
  fout << "  Require conincidence between Si Delta-E and Si E\n"
       << TString::Format("   -> %lld events were stored (%.3f of total) --\n\n",
        eventStored, (Double_t)eventStored/(Double_t)eventTrig*100.0);
}

void ConfigureDateTime(TString RunName, TString &DateTime, Int_t &FLast)
{
  eout << endl;
  TList *listOfFiles = TSystemDirectory("test_data","test_data").GetListOfFiles();
  /*
  TIter next(listOfFiles);
  TSystemFile* fileObject;
  int divisionMax = 0;
  while ((fileObject=(TSystemFile*)next()))
  {
    if (fileObject->IsDirectory()) continue;
    TString fileName = fileObject -> GetName();
    if (fileName.EndsWith(".dat")==false) continue;
    DateTime = fileName(7,12);
    int division = TString(fileName(25,3)).Atoi();
    if (division>divisionMax)
      divisionMax = division;
    cout << "Found file " << fileName << " : " <<  DateTime << " / " << division << endl;
  }
  cout << divisionMax << endl;
  FLast = divisionMax;
  */
}

void SortSi(int RunNumber)
{
  SortSi(TString(Form("RUN%03d",RunNumber)),"",-1);
}

void SortSi(TString RunName, TString DateTime, Int_t FLast)
{
  /*
  if (DateTime.IsNull())
    ConfigureDateTime(RunName,DateTime,FLast);
    eout << DateTime << endl;
    eout << FLast << endl;
  return;

  ////////////////
  // Start time //
  ////////////////
  timer.Start();
  //////////////
  // Log file //
  //////////////
  FLog.open( Form("%s/%s.log",PathOut.Data(), RunName.Data() ) );
  eout << Form("Logger name is %s/%s.log",PathOut.Data(), RunName.Data() )  << endl;
  if( FLog.fail() ){
    cout << TString::Format(" Can't open output file: %s/%s.log \n\n",PathOut.Data(),RunName.Data() );
    FLog.close();
    exit(1);
  }
  FLog << TString::Format(" SortSi(%s, %s, %d)\n",RunName.Data(), DateTime.Data(), FLast ) << endl;
  ///////////////////////////////////////////////////////////
  // Definition of Mod and Ch -> name & title of histogram //
  ///////////////////////////////////////////////////////////
  DefFile.open( DefFileName.Data() );
  eout << "Mod -> Ch map is " << DefFileName.Data() << endl;
  if( DefFile.fail() ){
    cout << TString::Format(" Can't open iput file: %s \n\n", DefFileName.Data() );
    FLog.close();
    DefFile.close();
    exit(1);
  }
  FLog << TString::Format("  Mod and Ch definition file: %s", DefFileName.Data() ) << endl;
  HName = new TString*[NMod];
  HTitl = new TString*[NMod];
  for(Int_t l=0; l<NMod; l++){
    HName[l] = new TString[16];
    HTitl[l] = new TString[16];
    for(Int_t m=0; m<16; m++){
      HName[l][m] = "";
      HTitl[l][m] = "";
    }
  }
  getline(DefFile, BufTmp);// skip header
  getline(DefFile, BufTmp);// 
  FLog << "   Mod\tCh\tName\tTitle" << endl;
  while(DefFile >> Mod >> Ch >> HNameBuf){
    // BufTmp (= Title of Histogram) contains space...
    getline(DefFile,BufTmp);
    // Erase (space or tab) in begining of BufTmp
    size_t pos;
    while( (pos = BufTmp.find_first_of(" 　\t") ) == 0){
      BufTmp.erase(BufTmp.begin() );
      if(BufTmp.empty() ) break;
    }
    HName[Mod][Ch] = HNameBuf;
    HTitl[Mod][Ch] = BufTmp;
    FLog << TString::Format("   %02d\t%02d\t%s\t%s",Mod,Ch,
        HName[Mod][Ch].Data(),HTitl[Mod][Ch].Data()) << endl;
  }
  FLog << endl;
  DefFile.close();

  ///////////////////////////////////////////////////////////////////////////////////////
  // Create and open output file, make directory in output file, and define histograms //
  ///////////////////////////////////////////////////////////////////////////////////////
  FHist = new TFile( Form("%s/%s_hist.root",PathOut.Data(), RunName.Data() ), "RECREATE");
  eout << "Histogram file name is " << FHist->GetName() << endl;
  CdMod = new TDirectory*[NMod];
  Hist  = new TH1D**[NMod];
  HistTimeElapsed= new TH1D**[NMod];
  for(Int_t l=0; l<NMod; l++){
    CdMod[l] = FHist->mkdir( Form("Mod%02d",l) );
    CdMod[l] ->cd();
    Hist[ l] = new TH1D*[16];
    HistTimeElapsed[l] = new TH1D*[16];
    for(Int_t m=0; m<16; m++){
      if(HName[l][m] != ""){
        Hist[l][m] = new TH1D(HName[l][m].Data(),HTitl[l][m].Data(),8192, -0.5, 8191.5);
        HNameBuf = TString::Format("%s_TS",HName[l][m].Data() );
        HTitlBuf = TString::Format("The time that elapsed after previous event (%s)",
            HTitl[l][m].Data() );
        HistTimeElapsed[l][m] = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),100, 0.5, 100.5);
      }
    }
  }
  return;

  // Histograms for Ts Difference
  CdTDiff = FHist->mkdir( NameCdTDiff.Data() );
  CdTDiff->cd();
  HTDiffEa = new TH1D("TDiffEa", "Ts diff. betw. #Delta E and Total-E-A",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffEa->GetXaxis()->SetTitle("Ts(Ea) - Ts(Trigger) [ch]");
  HTDiffEb = new TH1D("TDiffEb", "Ts diff. betw. #Delta E and Total-E-B",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffEb->GetXaxis()->SetTitle("Ts(Eb) - Ts(Trigger) [ch]"); HTDiffDE = new TH1D("TDiffDE", "Ts diff. betw. #Delta E and other #Delta E", (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffDE->GetXaxis()->SetTitle("Ts(#Delta E) - Ts(Trigger) [ch]");
  HTDiffQ  = new TH1D("TDiffQ",  "Ts diff. betw. #Delta E and MWPC(Q)",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffQ->GetXaxis()->SetTitle("Ts(Q) - Ts(Trigger) [ch]");
  HTDiffPos= new TH1D("TDiffPos","Ts diff. betw. #Delta E and MWPC(Pos)",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffPos->GetXaxis()->SetTitle("Ts(Pos) - Ts(Trigger) [ch]");
  HTDiffDt = new TH1D("TDiffDt", "Ts diff. betw. #Delta E and MWPC( #Delta T)",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffDt->GetXaxis()->SetTitle("Ts(#Delta T) - Ts(Trigger) [ch]");
  HTDiffMCP= new TH1D("TDiffMCP","Ts diff. betw. #Delta E and #Delta T(MWPC-MCP)",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffMCP->GetXaxis()->SetTitle("Ts(#Delta T[MWPC-MCP]) - Ts(Trigger) [ch]");
  HTDiffPH = new TH1D("TDiffPH", "Ts diff. betw. #Delta E and PH of ND",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffPH->GetXaxis()->SetTitle("Ts(PH of ND) - Ts(Trigger) [ch]");
  HTDiffPSD = new TH1D("TDiffPSD","Ts diff. betw. #Delta E and PSD of ND",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffPSD->GetXaxis()->SetTitle("Ts(PSD of ND) - Ts(Trigger) [ch]");
  HTDiffTOF = new TH1D("TDiffTOF","Ts diff. betw. #Delta E and TOF of ND",
      (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
  HTDiffTOF->GetXaxis()->SetTitle("Ts(TOF of ND) - Ts(Trigger) [ch]");
  for(Int_t m=0; m<12; m++){
    HNameBuf = TString::Format("TDiffEaDE%d",m+1);
    HTitlBuf = TString::Format("Ts diff. betw. #Delta E-%d and Total-E-A",m+1);
    HTDiffEaEach[m]=new TH1D(HNameBuf.Data(),HTitlBuf.Data(),
        (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
    HTDiffEaEach[m]->GetXaxis()->SetTitle("Ts(Ea) - Ts(Trigger) [ch]");
    //
    HNameBuf = TString::Format("TDiffEbDE%d",m+1);
    HTitlBuf = TString::Format("Ts diff. betw. #Delta E-%d and Total-E-B",m+1);
    HTDiffEbEach[m]=new TH1D(HNameBuf.Data(),HTitlBuf.Data(),
        (Int_t)(DTsMax-DTsMin)+1,(Double_t)DTsMin-0.5,(Double_t)DTsMax+0.5);
    HTDiffEbEach[m]->GetXaxis()->SetTitle("Ts(Ea) - Ts(Trigger) [ch]");
  }
  // Histograms for Peak Positions in (E-DE)-DE histograms
  CdPeak= FHist->mkdir( NameCdPeak.Data() );
  CdPeak->cd();
  // for time plotting
  Year          = DateTime( 0,4);
  Month         = DateTime( 4,2);
  Day           = DateTime( 6,2);
  Hour          = DateTime( 8,2);
  Minute        = DateTime(10,2);
  TimeFormat   += TString::Format("%%F%s-%s-%s %s:%s:00 GMT",
      Year.Data(),Month.Data(),Day.Data(),Hour.Data(),Minute.Data() );
  // 
  for(Int_t m=0; m<16; m++){
    // Peak Positions for Ea (X: Ts, Y: NoDE) m = NoEa
    HNameBuf = TString::Format("PeakPosEa%02dDE",                        m+1);
    HTitlBuf = TString::Format("Peak Positions for Ea-%02d and #Delta E",m+1);
    HPeakEaDE[m] = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1., 16,0.5,16.5);
    HPeakEaDE[m]  ->GetXaxis()->SetTimeDisplay(1);
    HPeakEaDE[m]  ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
    HPeakEaDE[m]  ->GetXaxis()->SetLabelOffset(0.025);
    HPeakEaDE[m]  ->GetYaxis()->SetTitle("Number of #Delta E");
    // Peak Positions for Eb (X: Ts, Y: NoDE) m = NoEb
    HNameBuf = TString::Format("PeakPosEb%02dDE",                        m+1);
    HTitlBuf = TString::Format("Peak Positions for Eb-%02d and #Delta E",m+1);
    HPeakEbDE[m] = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1., 16,0.5,16.5);
    HPeakEbDE[m]  ->GetXaxis()->SetTimeDisplay(1);
    HPeakEbDE[m]  ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
    HPeakEbDE[m]  ->GetXaxis()->SetLabelOffset(0.025);
    HPeakEbDE[m]  ->GetYaxis()->SetTitle("Number of #Delta E");
    // Peak Positions for DE (X: Ts, Y: NoEa) m = NoDE
    HNameBuf = TString::Format("PeakPosDE%02dEa",                        m+1);
    HTitlBuf = TString::Format("Peak Positions for #Delta E-%02d and Ea",m+1);
    HPeakDEEa[m] = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1., 16,0.5,16.5);
    HPeakDEEa[m]  ->GetXaxis()->SetTimeDisplay(1);
    HPeakDEEa[m]  ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
    HPeakDEEa[m]  ->GetXaxis()->SetLabelOffset(0.025);
    HPeakDEEa[m]  ->GetYaxis()->SetTitle("Number of Ea");
    // Peak Positions for DE (X: Ts, Y: NoEb) m = NoDE
    HNameBuf = TString::Format("PeakPosDE%02dEb",                        m+1);
    HTitlBuf = TString::Format("Peak Positions for #Delta E-%02d and Eb",m+1);
    HPeakDEEb[m] = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1., 16,0.5,16.5);
    HPeakDEEb[m]  ->GetXaxis()->SetTimeDisplay(1);
    HPeakDEEb[m]  ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
    HPeakDEEb[m]  ->GetXaxis()->SetLabelOffset(0.025);
    HPeakDEEb[m]  ->GetYaxis()->SetTitle("Number of Eb");
  }
  // Histograms for PID
  InitUseDE(); // Specify which E (Ea or Eb) should coincide with Delta-E
  CdPID = FHist->mkdir( NameCdPID.Data() );
  CdPID->cd();
  for(NoDE=0; NoDE<16; NoDE++){
    if( (HName[ModDE][NoDE] != "") && ( EaDEUse[NoDE] || EbDEUse[NoDE] ) ){
      HNameBuf = TString::Format("E_%s_PID",HName[ModDE][NoDE].Data() );
      if( EaDEUse[NoDE] )
        HTitlBuf = TString::Format("(Ea - #Delta E) - #Delta E plot (#Delta E-%d)",NoDE+1);
      else if( EbDEUse[NoDE] )
        HTitlBuf = TString::Format("(Eb - #Delta E) - #Delta E plot (#Delta E-%d)",NoDE+1);
      else if( EaDEUse[NoDE] && EbDEUse[NoDE] )
        HTitlBuf = TString::Format("(E - #Delta E) - #Delta E plot (#Delta E-%d)",NoDE+1);
      HistEDE[NoDE] = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),
          NbinsChX,XLowPIDCh,XUpPIDCh, NbinsChY, YLowPIDCh,YUpPIDCh);
      HistEDE[NoDE]->GetXaxis()->SetTitle("E - #Delta E (arb. unit)");
      HistEDE[NoDE]->GetYaxis()->SetTitle("#Delta E (arb. unit)");
      //	cout << HNameBuf.Data() << endl << endl;
    }
  }
  HNameBuf = "E_DE_PID";
  //  HTitlBuf = Form("sum of (E - #Delta E) - #Delta E PID plot (Except for #Delta E-9) [%s]",RunName.Data() );
  HTitlBuf = Form("sum of (E - #Delta E) - #Delta E PID plot [%s]",RunName.Data() );
  HistEDESum = new TH2D(HNameBuf.Data(),HTitlBuf.Data(),
      NbinsChX,XLowPIDCh,XUpPIDCh, NbinsChY, YLowPIDCh,YUpPIDCh);
  HistEDESum->GetXaxis()->SetTitle("E - #Delta E (arb. unit)");
  HistEDESum->GetYaxis()->SetTitle("#Delta E (arb. unit)");
  FHist->cd();
  // Histograms for Counting Rates
  CdCr= FHist->mkdir( NameCdCr.Data() );
  // Ea
  HNameBuf = "HistCrEa";
  HTitlBuf = "Counting Rate of #Sigma E_{a}";
  HistCrEa = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrEa ->GetXaxis()->SetTimeDisplay(1);
  HistCrEa ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrEa ->GetXaxis()->SetLabelOffset(0.025);
  // Eb
  HNameBuf = "HistCrEb";
  HTitlBuf = "Counting Rate of #Sigma E_{b}";
  HistCrEb = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrEb ->GetXaxis()->SetTimeDisplay(1);
  HistCrEb ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrEb ->GetXaxis()->SetLabelOffset(0.025);
  // DE
  for(Int_t m=0; m<16; m++){
    HNameBuf    = TString::Format("HistCrDE%02d", m+1);
    HTitlBuf    = TString::Format("Counting Rate of #Deleta E-%02d", m+1);
    HistCrDE[m] = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
    HistCrDE[m] ->GetXaxis()->SetTimeDisplay(1);
    HistCrDE[m] ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
    HistCrDE[m] ->GetXaxis()->SetLabelOffset(0.025);
  }
  // FC
  HNameBuf = "HistCrFC";
  HTitlBuf = TString::Format("Counting Rate of FC");
  HistCrFC = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrFC ->GetXaxis()->SetTimeDisplay(1);
  HistCrFC ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrFC ->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-A and MWPC-B)
  HNameBuf    = "HistCrDt_AB";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-A and MWPC-B)"; 
  HistCrDt[0] = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrDt[0] ->GetXaxis()->SetTimeDisplay(1);
  HistCrDt[0] ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrDt[0] ->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-B and MWPC-C)
  HNameBuf    = "HistCrDt_BC";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-B and MWPC-C)"; 
  HistCrDt[1] = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrDt[1] ->GetXaxis()->SetTimeDisplay(1);
  HistCrDt[1] ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrDt[1] ->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-A and MWPC-D)
  HNameBuf    = "HistCrDt_AD";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-A and MWPC-D)"; 
  HistCrDt[2] = new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrDt[2] ->GetXaxis()->SetTimeDisplay(1);
  HistCrDt[2] ->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrDt[2] ->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-A and MCP2)
  HNameBuf    = "HistCrCMP_A2";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-A and MCP2)";
  HistCrMCP[0]= new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrMCP[0]->GetXaxis()->SetTimeDisplay(1);
  HistCrMCP[0]->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrMCP[0]->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-B and MCP1)
  HNameBuf    = "HistCrCMP_B1";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-B and MCP1)";
  HistCrMCP[1]= new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrMCP[1]->GetXaxis()->SetTimeDisplay(1);
  HistCrMCP[1]->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrMCP[1]->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-C and MCP1)
  HNameBuf    = "HistCrCMP_C1";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-C and MCP1)";
  HistCrMCP[2]= new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrMCP[2]->GetXaxis()->SetTimeDisplay(1);
  HistCrMCP[2]->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrMCP[2]->GetXaxis()->SetLabelOffset(0.025);
  // Dt (MWPC-D and MCP2)
  HNameBuf    = "HistCrCMP_D2";
  HTitlBuf    = "Counting Rates of #Delta T (betw. MWPC-D and MCP2)";
  HistCrMCP[3]= new TH1D(HNameBuf.Data(),HTitlBuf.Data(),1,0.,1.);
  HistCrMCP[3]->GetXaxis()->SetTimeDisplay(1);
  HistCrMCP[3]->GetXaxis()->SetTimeFormat( TimeFormat.Data() );
  HistCrMCP[3]->GetXaxis()->SetLabelOffset(0.025);

  ///////////////////////////////////////////
  // Create and open output file for TTree //
  ///////////////////////////////////////////
  // Channel assignment for MWPCs
  InitChMWPC();
  // Output file
  FTree   = new TFile( Form("%s/%s_tree.root",PathOut.Data(), RunName.Data() ), "RECREATE");
  // Make TTree
  OutTree = new TTree( NameOutTree.Data(), RunName.Data() );
  // Branches
  OutTree->Branch("DataT", &DataTT, "DataTT/L"); // L: Long64_t
  OutTree->Branch("NumDE", &NumDE,  "NumDE/s");  // s: UShort_t
  OutTree->Branch("NumEa", &NumEa,  "NumEa/s"); 
  OutTree->Branch("NumEb", &NumEb,  "NumEb/s"); 
  OutTree->Branch("NumDt", &NumDt,  "NumDt/s"); 
  OutTree->Branch("AdcDE",  AdcDE,  "AdcDE[16]/s");
  OutTree->Branch("AdcEa",  AdcEa,  "AdcEa[16]/s");
  OutTree->Branch("AdcEb",  AdcEb,  "AdcEb[16]/s");
  OutTree->Branch("QdcQ",   QdcQ,   "QdcQ[4]/s");
  OutTree->Branch("TdcX",   TdcX,   "TdcX[4]/s");
  OutTree->Branch("TdcY",   TdcY,   "TdcY[4]/s");
  OutTree->Branch("TdcDt",  TdcDt,  "TdcDt[4]/s");
  OutTree->Branch("TdcMCP", TdcMCP, "TdcMCP[4]/s");
  OutTree->Branch("PhNdR",  PhNdR,  "PhNdR[17]/s");
  OutTree->Branch("PhNdL",  PhNdL,  "PhNdL[17]/s");
  OutTree->Branch("PsdNdR", PsdNdR, "PsdNdR[17]/s");
  OutTree->Branch("PsdNdL", PsdNdL, "PsdNdL[17]/s");
  OutTree->Branch("TofNdR", TofNdR, "TofNdR[17]/s");
  OutTree->Branch("TofNdL", TofNdL, "TofNdL[17]/s");
  ////////////
  // Canvas //
  ////////////
  HNameBuf = Form("SortSi [%s]",RunName.Data() );
  Canvas = new TCanvas("Canvas",HNameBuf.Data(),1400,1000);
  Canvas->Divide(2,1);//   0      1     2       3       4     5        6       7
  Color_t lcolor[] = {kBlack,  kGray,  kRed, kGreen,   kYellow,kBlue,kMagenta,kOrange,
    kGreen+3,kCyan+2,kCyan,kYellow+2,kRed+2, kSpring,kAzure,kTeal};
  Canvas->cd(1);
  gPad  ->Divide(1,2);
  Canvas->cd(2);
  gPad  ->Divide(1,4);
  //////////////////////
  // Number of Events //
  //////////////////////
  Event    = new Long64_t*[NMod];
  Err      = new Long64_t*[NMod];
  Err2     = new Long64_t*[NMod];
  EventTree= new Long64_t*[NMod];
  TS       = new Long64_t*[NMod];
  DA       = new UShort_t*[NMod];
  WF       = new UShort_t*[NMod];
  EachEvent= new Double_t*[NMod];
  for(Int_t l=0; l<NMod; l++){
    Event[l]    = new Long64_t[16];
    Err[  l]    = new Long64_t[16];
    Err2[ l]    = new Long64_t[16];
    EventTree[l]= new Long64_t[16];
    TS[   l]    = new Long64_t[16];
    DA[   l]    = new UShort_t[16];
    WF[   l]    = new UShort_t[16];
    EachEvent[l]= new Double_t[16];
    for(Int_t m=0; m<16; m++){
      Event[    l][m] = 0;
      Err[      l][m] = 0;
      Err2[     l][m] = 0;
      EventTree[l][m] = 0;
      TS[       l][m] = 0;
      DA[       l][m] = 0;
      WF[       l][m] = 0;
      EachEvent[l][m] = 0.;
    }
  }
  EventTrig   = 0;
  EventStored = 0;
  //////////////////////////////////////////////////////////
  // Initialize vector for Last Time Stamp of each Module //
  //////////////////////////////////////////////////////////
  for(Int_t l=0; l<NMod; l++)
    TsLastMod.push_back(0xFFFFFFFFFFFF); // 48bit (max. of Time Stamp: 40bit)
  /////////////////////////////////////////////////////////////////////
  // endless loop (1) to read of Raw-Data files written in given RunName //
  /////////////////////////////////////////////////////////////////////
  while(1){
    ////////////////
    // Last File? //
    ////////////////
    if(FLast >=0 && FNum>FLast) break;
    /////////////////////////////////
    // Name of Input Raw-Data File //
    /////////////////////////////////
    FileInName = TString::Format("%s/%s_%s_list_%03d.dat",
        PathIn.Data(),RunName.Data(),DateTime.Data(),FNum);
    /////////////////////////////////////////////////////
    // Open Input Raw-Data File & check the size of it //
    /////////////////////////////////////////////////////
    while(1){ // endless loop (2) to check the status of the Raw-Data file //
      FileIn.open( FileInName.Data() );
      if( FileIn.fail() ){
        FileIn.close();
        if(FOpen > 10)
          break;
        else{
          FOpen++;
          sleep(10);
          continue;
        }
      }
      SSize = FileIn.seekg(0, ios::end).tellg(); // obtain the size of file
      FileIn.seekg(0, ios::beg);                 // rewind
      if(SSize > SSizeMax || (FOpen != 0 && SSize==SSizeOld) ) { // good file or final file
        FOpen = 0;
        FNum++;
        break;
      }
      else if(FOpen == 0){ // first try
        FOpen++;
        FileIn.close();
        SSizeOld = SSize;
        sleep(10);
      }
      else if(SSize > SSizeOld){ // writing the file is still continued ...
        FileIn.close();
        SSizeOld = SSize;
        sleep(60);
      }
    } // end of endless loop (2)
    if(FOpen != 0){ // Next file is not found in the directory (PathIn)
      cout << TString::Format("  %s is Not Found .. \n\n",   FileInName.Data() );
      break;
    }
    if(FNum == 1){ // First File
      FLog << TString::Format(" Raw-Data files are : %s", FileInName.Data() ) << endl;
    }
    else{
      FLog << TString::Format("                      %s", FileInName.Data() ) << endl;
    }
    cout << TString::Format(" Raw-Data file : %s\n", FileInName.Data() );
    ////////////////////////////////////////////
    // Number of Event for Each Raw-Data file //
    ////////////////////////////////////////////
    for(Int_t l=0; l<NMod; l++){
      for(Int_t m=0; m<16; m++){
        EachEvent[l][m] = 0.;
      }
    }
    EachTsFirst = 0;
    ////////////////////////
    // Read Raw-Data file //
    ////////////////////////
    while( FileIn >> Buf ){
      // Buffer -> Data
      NBuf = (Int_t)sscanf(Buf, "%hu,%hu,%hu,%lld,%hu", &Mod, &Ch, &DataA, &DataT, &Evt);
      // Mod,Ch, DataA, Evt: UShort_t, DataT: Long64_t
      // Check # of data in line
      if(NBuf != 5){
        cout << TString::Format(" Warning: # of data in line is not 5 (%d) %u %u %u %lld %u",
            NBuf, Mod, Ch, DataA, DataT, Evt) << endl; break;
      }
      // Check overflow of Time Stamp
      if(Evt>0)
        DataT += (Long64_t)(Evt) * 0xFFFFFFFFFF;
      // Total Number of Events
      Event[Mod][Ch]++;
      // Number of Events for each Raw-Data file
      EachEvent[Mod][Ch]++;
      EachTsLast = DataT;
      if(EachTsFirst == 0)
        EachTsFirst = DataT;
      // check time stamp
      if(TS[Mod][Ch]>0){ // Not first event for Mod,Ch
        if( DataT > TS[Mod][Ch] ){      // Time Stamp increased
          if(HName[Mod][Ch] != ""){
            //////////////////////////////////////////////
            // Fill histograms for almost all the event //
            //////////////////////////////////////////////
            Hist[Mod][Ch]->Fill( (Double_t)(DataA) );
            HistTimeElapsed[Mod][Ch]->Fill( (Double_t)(DataT - TS[Mod][Ch]) );
            ///////////////////////////////////////////////////////
            // Store the data for the previous event in multimap //
            ///////////////////////////////////////////////////////
            // data -> structure
            Dbuf.Mod  = Mod;
            Dbuf.Ch   = Ch;
            Dbuf.DataA= DA[Mod][Ch];
            WF[Mod][Ch]++;
            // structure -> multimap (for all Mod & Ch)
            DmapA.insert( p(TS[Mod][Ch],Dbuf) );
            TsLastMod[Mod] = TS[Mod][Ch];
            // structure -> multimap (for trigger event only)
            if( (Mod==ModDE) && ( EaDEUse[Ch] || EbDEUse[Ch] ) )
              DmapT.insert( p(TS[Mod][Ch],Dbuf) );
          }
        }
        else if( DataT == TS[Mod][Ch]){ // Same time stamp as previous event
          if(Err[Mod][Ch] < 5){
            cout << TString::Format("  same time stamp ( %lld ) for Mod=%u, Ch=%u, DataA=%u\n",
                DataT, Mod, Ch, DataA);
            FLog << "  ********************"
              << TString::Format("* same time stamp ( %lld ) for Mod=%u, Ch=%u, DataA=%u",
                  DataT, Mod, Ch, DataA) << endl;
          }
          Err[Mod][Ch]++;
        }
        else{ // Time Stamp decreased ... [ DataT < TS[Mod][Ch] ]
          if(Err2[Mod][Ch] < 3){
            cout << "\n ====================================================== \n"
              << " Warning: Time Stamp decreased !!\n"
              << TString::Format(" Mod=%u, Ch=%u : ( %lld -> %lld )\n",
                  Mod,Ch,TS[Mod][Ch],DataT)
              << " ====================================================== \n";
            FLog << "                      "
              << TString::Format("  Mod=%u, Ch=%u : ( %lld -> %lld )",
                  Mod,Ch,TS[Mod][Ch],DataT) << endl;
          }
          Err2[Mod][Ch]++;
        }
      }
      ///////////////////////
      // Save current Data //
      ///////////////////////
      if(DataT >= TS[Mod][Ch]){
        TS[Mod][Ch] = DataT;
        DA[Mod][Ch] = DataA;
        WF[Mod][Ch] = 0;
      }
    }
    ///////////////////////////////////////////////////////////////////////
    // End of Each Raw-Data file                                         //
    // Display (Error) / (Valid Event) ratio obtained up to current file //
    ///////////////////////////////////////////////////////////////////////
    PrintErrors(0, FLog, Err, Event, Err2);
    /////////////////////////
    // Close Raw-Data file //
    /////////////////////////
    FileIn.close();
    ////////////////////////////////
    // Define (E-DE)-DE histogams //
    ////////////////////////////////
    for(NoDE=0; NoDE<16; NoDE++){
      for(NoEa=0; NoEa<16; NoEa++){
        if( (HName[ModEa][NoEa] != "") && (HName[ModDE][NoDE] != "") &&
            (Event[ModEa][NoEa] != 0)  && (Event[ModDE][NoDE] != 0)  &&
            (NoEa <= ChLastEa)         && EaDEUse[NoDE] ){
          HNameBuf = TString::Format("%s_%s_%03d",HName[ModEa][NoEa].Data(),HName[ModDE][NoDE].Data(),
              FNum-1);
          HTitlBuf = TString::Format("(E- #Delta E) - #Delta E plot (Ea-%d,  #Delta E-%d)(%03d)",
              NoEa+1,NoDE+1,FNum-1);
          histEaDE[NoEa][NoDE]=new TH2I(HNameBuf.Data(),HTitlBuf.Data(),
              NbinsChX,XLowCh,XUpCh, NbinsChY,YLowCh,YUpCh);
        }
      }
      for(NoEb=0; NoEb<16; NoEb++){
        if( (HName[ModEb][NoEb] != "") && (HName[ModDE][NoDE] != "") &&
            (Event[ModEb][NoEb] != 0)  && (Event[ModDE][NoDE] != 0)  &&
            (NoEb <= ChLastEb)         && EbDEUse[NoDE] ){
          HNameBuf = TString::Format("%s_%s_%03d",HName[ModEb][NoEb].Data(),HName[ModDE][NoDE].Data(),
              FNum-1);
          HTitlBuf = TString::Format("(E - #Delta E) - #Delta E plot (Eb-%d,  #Delta E-%d)(%03d)",
              NoEb+1,NoDE+1,FNum-1);
          histEbDE[NoEb][NoDE]=new TH2I(HNameBuf.Data(),HTitlBuf.Data(),
              NbinsChX,XLowCh,XUpCh, NbinsChY,YLowCh,YUpCh);
          //	  cout << HNameBuf.Data() << endl;
        }
      }
    }
    ////////////////////////////////
    // Fill old event in multimap // 次の event が入力がなく保留となった event の処理
    ////////////////////////////////
    for(Int_t l=0; l<NMod; l++){
      for(Int_t m=0; m<16; m++){
        if( (DataT-TS[l][m]>5000000) && WF[l][m]==0 && TS[l][m]>0 && (HName[l][m] != "") ){
          // Time Stamp of 5000000 correspond 1 sec (if TimeBase = 200).
          // Fill histograms for almost all the event
          Hist[l][m]->Fill( (Double_t)(DataA) );
          HistTimeElapsed[l][m]->Fill( (Double_t)(DataT - TS[l][m]) );
          // data -> structure
          Dbuf.Mod  = (UShort_t)l;
          Dbuf.Ch   = (UShort_t)m;
          Dbuf.DataA= DA[l][m];
          WF[l][m]++;
          // structure -> multimap (for all Mod & Ch)
          DmapA.insert( p(TS[l][m],Dbuf) );
          if( TsLastMod[l]<TS[l][m] )
            TsLastMod[l] = TS[l][m];
          // structure -> multimap (for trigger event only)
          if( (l==ModDE) && (EaDEUse[m] || EbDEUse[m] ) )
            DmapT.insert( p(TS[l][m],Dbuf) );
        }
      }
    }
    ///////////////////////////
    // Sort multimap & fill  //
    ///////////////////////////
    LastTsSafe = *min_element( TsLastMod.begin(),TsLastMod.end() );
    cout << " LastTsSafe ... " << LastTsSafe << endl
      << " DataT ........ " << DataT      << endl;
    //    for(Int_t l=0; l<NMod; l++)
    //      cout << "   " << l << "\t" << TsLastMod[l] << endl;
    itrTrig = DmapT.begin();
    while( itrTrig != DmapT.end() ){
      // Data of Trigger
      DataTT = itrTrig->first;
      ChT    = itrTrig->second.Ch;
      ModT   = itrTrig->second.Mod;
      DataAT = itrTrig->second.DataA;
      if(DataTT > (LastTsSafe-DTsMax) )
        break;
      // setting up Time Window for Prompt Event
      TsWinMin   = DataTT + DTsMin;
      TsWinMax   = DataTT + DTsMax;
      itrAll     = DmapA.lower_bound(TsWinMin);
      if(itrAll != DmapA.begin() )
        --itrAll;
      itrLast    = DmapA.upper_bound(TsWinMax);
      if(itrLast!= DmapA.end()   )
        ++itrLast;
      // Initialize data for TTree
      NumDE  = 1; // for Trigger
      NumEa  = 0;
      NumEb  = 0;
      NumDt  = 0;
      for(Int_t m=0; m<16; m++){
        AdcDE[m] = 0;
        AdcEa[m] = 0;
        AdcEb[m] = 0;
        if(m<4){
          QdcQ[  m]=0;
          TdcX[  m]=0;
          TdcY[  m]=0;
          TdcDt[ m]=0;
          TdcMCP[m]=0;
        }
      }
      for(Int_t m=0; m<17; m++){
        PhNdR[ m]=0;
        PhNdL[ m]=0;
        PsdNdR[m]=0;
        PsdNdL[m]=0;
        TofNdR[m]=0;
        TofNdL[m]=0;
      }
      AdcDE[ChT]=DataAT; // for Trigger
      EventTrig++;
      EventTree[ModT][ChT]++;
      // Sort DmapA for Prompt Event
      while(itrAll != itrLast){
        // multimap -> data
        DataT = itrAll->first;
        Mod   = itrAll->second.Mod;
        Ch    = itrAll->second.Ch;
        DataA = itrAll->second.DataA;
        TDiff = (Int_t)(DataT-DataTT);
        if( (TDiff >= DTsMin) && (TDiff <= DTsMax) ){
          // data -> data for TTree
          if( Mod==ModEa ){                                             // Ea
            HTDiffEa->Fill( (Double_t)(TDiff) );
            HTDiffEaEach[ChT]->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinEa) && (TDiff<TsMaxEa) ){
              EventTree[Mod][Ch]++;
              NumEa++;
              AdcEa[Ch]=DataA;
              if( EaDEUse[ChT] && (Ch<=ChLastEa) )   // Ea - DeltaE histogram
                histEaDE[Ch][ChT]->Fill( (Double_t)(DataA), (Double_t)(DataAT) );
            }
          }
          else if( Mod==ModEb ){                                        // Eb
            HTDiffEb->Fill( (Double_t)(TDiff) );
            HTDiffEbEach[ChT]->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinEb) && (TDiff<TsMaxEb) ){
              EventTree[Mod][Ch]++;
              NumEb++;
              AdcEb[Ch]=DataA;
              if( EbDEUse[ChT] && (Ch<=ChLastEb) )   // Eb - DeltaE histogram
                histEbDE[Ch][ChT]->Fill( (Double_t)(DataA), (Double_t)(DataAT) );
            }
          }
          else if( (Mod==ModDE) && (Ch>=ChFirstDE) && (Ch<=ChLastDE) ){ // DE
            if( Ch != ChT ){ // data of DE other than trigger ch
              HTDiffDE->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinDE) && (TDiff<TsMaxDE) ){
                EventTree[Mod][Ch]++;
                NumDE++;
                AdcDE[Ch]=DataA;
              }
            }
          }
          else if( (Mod==ModQ) && (Ch>=ChFirstQ) && (Ch<=ChLastQ) && 
              (HName[Mod][Ch] != "") ){                            //Q(MWPC)
                HTDiffQ->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinQ) && (TDiff<TsMaxQ) ){
                  // Ch -> No. of Detector
                  itrMWPC = ChQ.find(Ch);
                  if(itrMWPC != ChQ.end() ){
                    EventTree[Mod][Ch]++;
                    ChUse = itrMWPC->second;
                    QdcQ[ChUse]= DataA;
                  }
                  else
                    cout << Form("  Mod=%02d Ch=%02d is not Pulse height of MWPC\n",Mod,Ch);
                }
              }
          else if( (Mod==ModPos) && (Ch>=ChFirstPos) && (Ch<=ChLastPos) && 
              (HName[Mod][Ch] != "") ){                            //Position(MWPC)
                HTDiffPos->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinPos) && (TDiff<TsMaxPos) ){
                  // Ch -> No. of Detector
                  itrMWPC = ChX.find(Ch);
                  if(itrMWPC != ChX.end() ){
                    EventTree[Mod][Ch]++;
                    ChUse = itrMWPC->second;
                    TdcX[ChUse]= DataA;
                  }
                  else{ // itrMWPC == ChX.end()
                    itrMWPC = ChY.find(Ch);
                    if(itrMWPC != ChY.end() ){
                      EventTree[Mod][Ch]++;
                      ChUse = itrMWPC->second;
                      TdcY[ChUse]= DataA;
                    }
                    else
                      cout << Form("  Mod=%02d Ch=%02d is not Position of MWPC\n",Mod,Ch);
                  }
                }
              }
          else if( (Mod==ModDt) && (Ch>=ChFirstDt) && (Ch<=ChLastDt) && 
              (HName[Mod][Ch] != "") ){                            //dT(MWPC)
                HTDiffDt->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinDt) && (TDiff<TsMaxDt) ){
                  // Ch -> No. of Detector
                  itrMWPC = ChDt.find(Ch);
                  if(itrMWPC != ChDt.end() ){
                    EventTree[Mod][Ch]++;
                    NumDt++;
                    ChUse = itrMWPC->second;
                    TdcDt[ChUse]= DataA;
                  }
                  else
                    cout << Form("  Mod=%02d Ch=%02d is not dT of MWPC\n",Mod,Ch);
                }
              }
          else if( (Mod==5) && (Ch==15) && (HName[Mod][Ch] != "") ){  //dT(MWPC)
            HTDiffDt->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinDt) && (TDiff<TsMaxDt) ){
              // Ch -> No. of Detector
              itrMWPC = ChDt.find(Ch);
              if(itrMWPC != ChDt.end() ){
                EventTree[Mod][Ch]++;
                NumDt++;
                ChUse = itrMWPC->second;
                TdcDt[ChUse]= DataA;
              }
              else
                cout << Form("  Mod=%02d Ch=%02d is not dT of MWPC\n",Mod,Ch);
            }
          }
          else if( (Mod==ModMCP) && (Ch>=ChFirstMCP) && (Ch<=ChLastMCP) &&
              (HName[Mod][Ch] != "") ){                            //dT(MWPC-MCP)
                HTDiffMCP->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinMCP) && (TDiff<TsMaxMCP) ){
                  // Ch -> No. of Detector
                  itrMWPC = ChMCP.find(Ch);
                  if(itrMWPC != ChMCP.end() ){
                    EventTree[Mod][Ch]++;
                    ChUse = itrMWPC->second;
                    TdcMCP[ChUse]= DataA;
                  }
                  else
                    cout << Form("  Mod=%02d Ch=%02d is not DeltaT (MWPC-MCP)\n", Mod,Ch);
                }
              }
          else if( Mod==13 && Ch==0 ){                             // ND-R01 PH
            HTDiffPH ->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
              PhNdR[0]=DataA;
            }
          }
          else if( Mod==13 && Ch==1 ){                             // ND-R09 PH
            HTDiffPH ->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
              PhNdR[8]=DataA;
            }
          }
          else if( (Mod==ModNdR) || ( Mod==(ModNdR+1) ) ||
              ( ( Mod==(ModNdR+2) ) && Ch<2 ) ){              // ND-R
            if(Ch % 2 == 0){ // Pulse Height
              HTDiffPH ->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
                // Number of Detector
                ChND = (Mod-ModNdR)*8 + Ch/2;
                // Fill Pulse Height
                PhNdR[ChND]=DataA;
              }
            }
            else{            // PSD
              HTDiffPSD->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinPSD) && (TDiff<TsMaxPSD) ){
                // Number of Detector
                ChND = (Mod-ModNdR)*8 + (Ch-1)/2;
                // Fill PSD
                PsdNdR[ChND]=DataA;
              }
            }
          }   // else if( (Mod==ModNdR) ...)
          else if( ( (Mod==ModNdL)       && (Ch>1) ) ||
              ( Mod==(ModNdL+1) )         ||
              ( ( Mod==(ModNdL+2) ) && (Ch<4) ) ){            // ND-L
            if(Ch % 2 == 0){ // Pulse Height
              HTDiffPH ->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
                // Number of Detector
                ChND = (Mod-ModNdL)*8 + Ch/2-1;
                // Fill Pulse Height
                PhNdL[ChND]=DataA;
              }
            }
            else{            // PSD
              HTDiffPSD->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinPSD) && (TDiff<TsMaxPSD) ){
                // Number of Detector
                ChND = (Mod-ModNdL)*8 + (Ch-1)/2-1;
                // Fill PSD
                PsdNdL[ChND]=DataA;
              }
            }
          }   // else if( (Mod==ModNdL)...)
          else if( ( Mod== ModTofR) ||
              ( Mod==(ModTofR+1) && Ch<1) ){               // TOF(ND-R)
                HTDiffTOF->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinTOF) && (TDiff<TsMaxTOF) ){
                  // Number of Detector
                  ChND = (Mod-ModTofR)*16+Ch;
                  // Fill TOF
                  TofNdR[ChND]=DataA;
                }
              }   // else if( ( Mod== ModTofR) ...)
          else if( ( Mod==(ModTofL)   && Ch>0) ||
              ( Mod==(ModTofL+1) && Ch<1) ){               // TOF(ND-L)
                HTDiffTOF->Fill( (Double_t)(TDiff) );
                if( (TDiff>TsMinTOF) && (TDiff<TsMaxTOF) ){
                  // Number of Detector
                  if(Mod==ModTofL)
                    ChND = Ch-1;
                  else
                    ChND = 16;
                  // Fill TOF
                  TofNdL[ChND]=DataA;
                }
              }   // else if( ( Mod==(ModTofL) ....) )
        }
        // next event
        ++itrAll;
      }
      // Fill TTree
      if( (NumEa > 0) || (NumEb > 0) ){
        OutTree->Fill();
        EventStored++;
      }
      // next trigger
      ++itrTrig;
    }
    ////////////////////////////////////////////
    // Last Time Stamp for Triger in multimap //
    ////////////////////////////////////////////
    VecTS.push_back(DataTT); 
    /////////////////////////////////////////////////////////
    // remove sorted elements from the multimap containers //
    /////////////////////////////////////////////////////////
    itrLast = DmapA.upper_bound(LastTsSafe-DTsMax);
    cout << " Size of multimap for 1 list data file ... " << DmapA.size() << " (All)"     << endl;
    cout << "                                       ... " << DmapT.size() << " (Trigger)" << endl;
    DmapA.erase(DmapA.begin(),itrLast);
    DmapT.erase(DmapT.begin(),itrTrig);
    cout << " Size of multimap remaining .............. " << DmapA.size() << " (All)" << endl;
    cout << "                                       ... " << DmapT.size() << " (Trigger)"
      << endl << endl;
    ////////////////////////////////////////
    // Search Peak in (E-DE)-DE histogams //
    ////////////////////////////////////////
    //    cout << "NoE\tNoDE\tPeakEE\tPeakDE" << endl;
    HistEDESum->Reset();
    for(NoDE=0; NoDE<16; NoDE++){
      Double_t PeakEE, PeakDE;
      Bool_t   fPID;
      for(NoEa=0; NoEa<16; NoEa++){
        if( (HName[ModEa][NoEa] != "") && (HName[ModDE][NoDE] != "") &&
            (Event[ModEa][NoEa] != 0)  && (Event[ModDE][NoDE] != 0)  &&
            (NoEa <= ChLastEa)         && EaDEUse[NoDE] ){
          if(NoEa >= ChFirstEa)
            fPID = kTRUE;
          else
            fPID = kFALSE;
          PeakPosition2(fPID, histEaDE[NoEa][NoDE],HistEDE[NoDE],PeakEE,PeakDE);
          //	  cout << NoEa+1 << "\t" << NoDE+1 << "\t" << PeakEE << "\t" << PeakDE << endl;
          VecEaDETmp[NoEa][NoDE] = PeakEE;
          VecDEEaTmp[NoDE][NoEa] = PeakDE;
        }
      }
      //
      for(NoEb=0; NoEb<16; NoEb++){
        if( (HName[ModEb][NoEb] != "") && (HName[ModDE][NoDE] != "") &&
            (Event[ModEb][NoEb] != 0)  && (Event[ModDE][NoDE] != 0)  &&
            (NoEb <= ChLastEb)         && EbDEUse[NoDE] ){
          if(NoEb >= ChFirstEb)
            fPID = kTRUE;
          else
            fPID = kFALSE;
          PeakPosition2(fPID, histEbDE[NoEb][NoDE],HistEDE[NoDE],PeakEE,PeakDE);
          //	  cout << NoEb+1 << "\t" << NoDE+1 << "\t" << PeakEE << "\t" << PeakDE << endl;
          VecEbDETmp[NoEb][NoDE] = PeakEE;
          VecDEEbTmp[NoDE][NoEb] = PeakDE;
        }
      }
      //      cout << endl;
      if( (HName[ModDE][NoDE] != "") && ( EaDEUse[NoDE] || EbDEUse[NoDE] ) 
          && (Event[ModDE][NoDE] != 0) ) //&& NoDE != 8)
            HistEDESum->Add(HistEDE[NoDE]);
    }
    HNameBuf  = Form("SortSi [%s 000-%03d]",RunName.Data(),FNum-1 );
    Canvas    ->SetTitle(HNameBuf.Data() );
    Canvas    ->cd(1);
    gPad      ->cd(1);
    gPad      ->SetLogz(1);
    HistEDESum->ResetStats();
    HistEDESum->Draw("colz");
    Canvas->Update();
    TPaveStats *st = (TPaveStats*)HistEDESum->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.7);
    st->SetY2NDC(0.9);
    Canvas->Modified();
    VecPeakEaDE.push_back(VecEaDETmp);
    VecPeakEbDE.push_back(VecEbDETmp);
    VecPeakDEEa.push_back(VecDEEaTmp);
    VecPeakDEEb.push_back(VecDEEbTmp);
    /////////////////////////////////////////////////////////////////////
    // Plot Peak positions in (E-DE)-DE histograms as a function of Ts //
    /////////////////////////////////////////////////////////////////////
    NbinsT = (Int_t)VecTS.size();  // = FNum+2
    XbinsT = new Double_t[NbinsT]; // [FNum+2]
    YbinsT = new Double_t[17];
    for(Int_t n=0; n<NbinsT; n++)
      XbinsT[n] = (Double_t)VecTS[n]*TimeBase*1.e-9; // Low edges (sec)
    for(Int_t n=0; n<17; n++)
      YbinsT[n] = (Double_t)n + 0.5; //cout << n << "\t" << n+0.5 << endl;
    // Define histogram and fill it
    for(Int_t m=0; m<16; m++){
      // Peak Positions for Ea (X: Ts, Y: NoDE) m = NoEa
      HPeakEaDE[m]->Reset();
      HPeakEaDE[m]->SetBins(NbinsT-1,XbinsT,16,YbinsT);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        for(Int_t NoDE=0; NoDE<16; NoDE++){
          if( VecPeakEaDE[nx][m][NoDE] > 0){
            Xbin.push_back( HPeakEaDE[m]->GetXaxis()->GetBinCenter(nx+1) );
            Ybin.push_back( Double_t(NoDE+1) );
            Wbin.push_back( VecPeakEaDE[nx][m][NoDE] );
          }
        }
      }
      HPeakEaDE[m]->FillN( (Int_t)Xbin.size(), Xbin.data(), Ybin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Ybin.clear();
      Wbin.clear();
      // Peak Positions for Eb (X: Ts, Y: NoDE) m = NoEb
      HPeakEbDE[m]->Reset();
      HPeakEbDE[m]->SetBins(NbinsT-1,XbinsT,16,YbinsT);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        for(Int_t NoDE=0; NoDE<16; NoDE++){
          if( VecPeakEbDE[nx][m][NoDE] > 0){
            Xbin.push_back( HPeakEbDE[m]->GetXaxis()->GetBinCenter(nx+1) );
            Ybin.push_back( Double_t(NoDE+1) );
            Wbin.push_back( VecPeakEbDE[nx][m][NoDE] );
          }
        }
      }
      HPeakEbDE[m]->FillN( (Int_t)Xbin.size(), Xbin.data(), Ybin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Ybin.clear();
      Wbin.clear();
      // Peak Positions for DE (X: Ts, Y: NoEa) m = NoDE
      HPeakDEEa[m]->Reset();
      HPeakDEEa[m]->SetBins(NbinsT-1,XbinsT,16,YbinsT);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        for(Int_t NoEa=0; NoEa<16; NoEa++){
          if( VecPeakDEEa[nx][m][NoEa] > 0){
            Xbin.push_back( HPeakDEEa[m]->GetXaxis()->GetBinCenter(nx+1) );
            Ybin.push_back( Double_t(NoEa+1) );
            Wbin.push_back( VecPeakDEEa[nx][m][NoEa] );
          }
        }
      }
      HPeakDEEa[m]->FillN( (Int_t)Xbin.size(), Xbin.data(), Ybin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Ybin.clear();
      Wbin.clear();
      if( ( HPeakDEEa[m]->Integral() ) > 0){
        HNameBuf = TString::Format("PjHPeakDE%02d",m+1);
        PjHPeakDE[m] = HPeakDEEa[m]->ProjectionX(HNameBuf.Data(),PjChUseE,PjChUseE);
        PjHPeakDE[m]->SetStats(0);
      }
      // Peak Positions for DE (X: Ts, Y: NoEb) m = NoDE
      HPeakDEEb[m]->Reset();
      HPeakDEEb[m]->SetBins(NbinsT-1,XbinsT,16,YbinsT);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        for(Int_t NoEb=0; NoEb<16; NoEb++){
          if( VecPeakDEEb[nx][m][NoEb] > 0){
            Xbin.push_back( HPeakDEEb[m]->GetXaxis()->GetBinCenter(nx+1) );
            Ybin.push_back( Double_t(NoEb+1) );
            Wbin.push_back( VecPeakDEEb[nx][m][NoEb] );
          }
        }
      }
      HPeakDEEb[m]->FillN( (Int_t)Xbin.size(), Xbin.data(), Ybin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Ybin.clear();
      Wbin.clear();
      if( ( HPeakDEEb[m]->Integral() ) > 0){
        HNameBuf = TString::Format("PjHPeakDE%02d",m+1);
        PjHPeakDE[m] = HPeakDEEb[m]->ProjectionX(HNameBuf.Data(),PjChUseE,PjChUseE);
        PjHPeakDE[m]->SetStats(0);
      }
    }
    delete [] XbinsT;
    delete [] YbinsT;
    Canvas->cd(1);
    gPad  ->cd(2);
    gPad  ->SetGridx(1);
    gPad  ->SetGridy(1);
    Double_t nFactor;
    Bool_t drawFirst = kTRUE;
    if(LegPj) delete LegPj;
    LegPj = new TLegend(0.1, 0.1, 0.4, 0.4);
    LegPj->SetNColumns(2);
    for(NoDE=0; NoDE<16; NoDE++){
      if( ( HPeakDEEa[NoDE]->Integral() ) > 0 || ( HPeakDEEb[NoDE]->Integral() ) > 0){
        PjHPeakDE[NoDE]->SetLineColor(lcolor[NoDE]);
        PjHPeakDE[NoDE]->SetLineWidth(2);
        nFactor = PjHPeakDE[NoDE]->GetBinContent(1);
        if(drawFirst){
          PjHPeakDE[NoDE]->Scale(1. / nFactor);
          PjHPeakDE[NoDE]->GetYaxis()->SetRangeUser(0.85, 1.05);
          HTitlBuf = TString::Format("Ts Dependency of Peak Positions for #Delta E - E-%02d [%s]",
              PjChUseE+1, RunName.Data() );
          PjHPeakDE[NoDE]->SetTitle(HTitlBuf.Data() );
          PjHPeakDE[NoDE]->Draw();
          drawFirst = kFALSE;
        }
        else{
          PjHPeakDE[NoDE]->Scale(1. / nFactor);
          PjHPeakDE[NoDE]->Draw("same");
        }
        LegPj->AddEntry(PjHPeakDE[NoDE],Form("#Delta E-%02d",NoDE+1),"l");
      }
    }
    LegPj ->Draw();
    Canvas->Update();
    Canvas->cd();
    ///////////////////////////////////////////////////////
    // Save (E-DE)-DE Histogamrs for first Raw-Data file //
    ///////////////////////////////////////////////////////
    if(FNum ==1){
      CdEDEFirst = FHist->mkdir( NameCdEDEFirst.Data() );
      CdEDEFirst->cd();
      for(NoDE=0; NoDE<16; NoDE++){
        for(NoEa=0; NoEa<16; NoEa++){
          if( (HName[ModEa][NoEa] != "") && (HName[ModDE][NoDE] != "") &&
              (Event[ModEa][NoEa] != 0)  && (Event[ModDE][NoDE] != 0)  &&
              (NoEa <= ChLastEa)         && EaDEUse[NoDE] )
            histEaDE[NoEa][NoDE]->Write();
        }
        for(NoEb=0; NoEb<16; NoEb++){
          if( (HName[ModEb][NoEb] != "") && (HName[ModDE][NoDE] != "") &&
              (Event[ModEb][NoEb] != 0)  && (Event[ModDE][NoDE] != 0)  &&
              (NoEb <= ChLastEb)         && EbDEUse[NoDE] )
            histEbDE[NoEb][NoDE]->Write();
        }
      }
      // PID plot
      for(NoDE=0; NoDE<16; NoDE++){
        if( (HName[ModDE][NoDE] != "") && ( EaDEUse[NoDE] || EbDEUse[NoDE] ) && 
            (Event[ModDE][NoDE] != 0) ){
          HistEDE[NoDE]->ResetStats();
          HistEDE[NoDE]->Write();
        }
      }
      HistEDESum->ResetStats();
      HistEDESum->Write();
    }
    /////////////////////////////////////
    // Deallocate (E-DE)-DE histograms //
    /////////////////////////////////////
    for(NoDE=0; NoDE<16; NoDE++){
      for(NoEa=0; NoEa<16; NoEa++){
        if( (HName[ModEa][NoEa] != "") && (HName[ModDE][NoDE] != "") ) delete histEaDE[NoEa][NoDE];
      }
      for(NoEb=0; NoEb<16; NoEb++){
        if( (HName[ModEb][NoEb] != "") && (HName[ModDE][NoDE] != "") ) delete histEbDE[NoEb][NoDE];
      }
    }
    //////////////////////////
    // Check Counting Rates //
    //////////////////////////
    DeltaTs = (Double_t)(EachTsLast - EachTsFirst)*TimeBase*1.e-9; // ch->sec
    // Ea
    CountSum=0;
    for(NoEa=0; NoEa<16; NoEa++){
      if( HName[ModEa][NoEa] != "" )
        CountSum+=EachEvent[ModEa][NoEa];
    }
    VecCrEa.push_back(CountSum/DeltaTs);
    // Eb
    CountSum=0;
    for(NoEb=0; NoEb<16; NoEb++){
      if( HName[ModEb][NoEb] != "" )
        CountSum+=EachEvent[ModEb][NoEb];
    }
    VecCrEb.push_back(CountSum/DeltaTs);
    // DE
    for(NoDE=0; NoDE<16; NoDE++){
      if( HName[ModDE][NoDE] != "" )
        VecCrDETmp[NoDE] = EachEvent[ModDE][NoDE]/DeltaTs;
    }
    VecCrDE.push_back(VecCrDETmp);
    // FC
    if( HName[ModFC][ChFC] != "" )
      VecCrFC.push_back(EachEvent[ModFC][ChFC]/DeltaTs);
    // Dt (MWPC)
    for(Int_t m=ChFirstDt; m<=ChLastDt; m++){
      itrMWPC = ChDt.find(m);
      ChUse   = itrMWPC->second;
      if( (HName[ModDt][m] != "") && EachEvent[ModDt][m] != 0 )
        VecCrDtTmp[ChUse] = EachEvent[ModDt][m]/DeltaTs;
      if( (HName[5][m] != "") && EachEvent[5][m] != 0 )
        VecCrDtTmp[ChUse] = EachEvent[5][m]/DeltaTs;
    }
    VecCrDt.push_back(VecCrDtTmp);
    // Dt (MWPC-MCP)
    for(Int_t m=ChFirstMCP; m<=ChLastMCP; m++){
      itrMWPC = ChMCP.find(m);
      ChUse   = itrMWPC->second;
      if( (HName[ModDt][m] != "") && EachEvent[ModDt][m] != 0 )
        VecCrMCPTmp[ChUse] = EachEvent[ModMCP][m]/DeltaTs;
    }
    VecCrMCP.push_back(VecCrMCPTmp);
    /////////////////////////////////////////////
    // Plot Counting Rates as a function of Ts //
    /////////////////////////////////////////////
    XbinsT = new Double_t[NbinsT]; // [FNum+2]
    for(Int_t n=0; n<NbinsT; n++)
      XbinsT[n] = (Double_t)VecTS[n]*TimeBase*1.e-9; // Low edges (sec)
    // Ea
    HistCrEa->Reset();
    HistCrEa->GetYaxis()->UnZoom();
    HistCrEa->SetStats(0);
    HistCrEa->SetBins(NbinsT-1,XbinsT);
    for(Int_t nx=0; nx<NbinsT-1; nx++){
      if( VecCrEa[nx] > 0){
        Xbin.push_back( HistCrEa->GetXaxis()->GetBinCenter(nx+1) );
        Wbin.push_back( VecCrEa[nx] );
      }
    }
    HistCrEa->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
    Xbin.clear();
    Wbin.clear();
    // Eb
    HistCrEb->Reset();
    HistCrEb->GetYaxis()->UnZoom();
    HistCrEb->SetStats(0);
    HistCrEb->SetBins(NbinsT-1,XbinsT);
    for(Int_t nx=0; nx<NbinsT-1; nx++){
      if( VecCrEb[nx] > 0){
        Xbin.push_back( HistCrEb->GetXaxis()->GetBinCenter(nx+1) );
        Wbin.push_back( VecCrEb[nx] );
      }
    }
    HistCrEb->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
    Xbin.clear();
    Wbin.clear();
    MaxEaEb = HistCrEa->GetMaximum();
    MinEaEb = HistCrEa->GetMinimum();
    if( (HistCrEb->GetMaximum() ) > MaxEaEb)
      MaxEaEb = HistCrEb->GetMaximum();
    if( (HistCrEb->GetMinimum() ) < MinEaEb)
      MinEaEb = HistCrEb->GetMinimum();
    // DE
    drawFirst = kTRUE;
    for(NoDE=0; NoDE<16; NoDE++){
      HistCrDE[NoDE]->Reset();
      HistCrDE[NoDE]->GetYaxis()->UnZoom();
      HistCrDE[NoDE]->SetBins(NbinsT-1,XbinsT);
      HistCrDE[NoDE]->SetStats(0);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        if( VecCrDE[nx][NoDE] > 0){
          Xbin.push_back( HistCrDE[NoDE]->GetXaxis()->GetBinCenter(nx+1) );
          Wbin.push_back( VecCrDE[nx][NoDE] );
        }
      }
      HistCrDE[NoDE]->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Wbin.clear();
      if(drawFirst){
        MaxDE=HistCrDE[NoDE]->GetMaximum();
        MinDE=HistCrDE[NoDE]->GetMinimum();
        drawFirst = kFALSE;
      }
      else{
        if( (HistCrDE[NoDE]->GetMaximum() ) > MaxDE )
          MaxDE=HistCrDE[NoDE]->GetMaximum();
        if( (HistCrDE[NoDE]->GetMinimum() ) < MinDE )
          MinDE=HistCrDE[NoDE]->GetMinimum();
      }
    }
    // FC
    HistCrFC->Reset();
    HistCrFC->SetStats(0);
    HistCrFC->SetBins(NbinsT-1,XbinsT);
    for(Int_t nx=0; nx<NbinsT-1; nx++){
      if( VecCrFC[nx] > 0){
        Xbin.push_back( HistCrFC->GetXaxis()->GetBinCenter(nx+1) );
        Wbin.push_back( VecCrFC[nx] );
      }
    }
    HistCrFC->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
    Xbin.clear();
    Wbin.clear();
    // Dt (MWPC - MWPC)
    drawFirst = kTRUE;
    for(ChUse=0; ChUse<3; ChUse++){
      HistCrDt[ChUse]->Reset();
      HistCrDt[ChUse]->GetYaxis()->UnZoom();
      HistCrDt[ChUse]->SetBins(NbinsT-1,XbinsT);
      HistCrDt[ChUse]->SetStats(0);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        if( VecCrDt[nx][ChUse] > 0){
          Xbin.push_back( HistCrDt[ChUse]->GetXaxis()->GetBinCenter(nx+1) );
          Wbin.push_back( VecCrDt[nx][ChUse] );
        }
      }
      HistCrDt[ChUse]->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Wbin.clear();
      if(drawFirst){
        MaxDt=HistCrDt[ChUse]->GetMaximum();
        MinDt=HistCrDt[ChUse]->GetMinimum();
        drawFirst = kFALSE;
      }
      else{
        if( (HistCrDt[ChUse]->GetMaximum() ) > MaxDt )
          MaxDt=HistCrDt[ChUse]->GetMaximum();
        if( (HistCrDt[ChUse]->GetMinimum() ) < MinDt )
          MinDt=HistCrDt[ChUse]->GetMinimum();
      }
    }
    // Dt (MWPC - MCP)
    drawFirst = kTRUE;
    for(ChUse=0; ChUse<4; ChUse++){
      HistCrMCP[ChUse]->Reset();
      HistCrMCP[ChUse]->GetYaxis()->UnZoom();
      HistCrMCP[ChUse]->SetBins(NbinsT-1,XbinsT);
      HistCrMCP[ChUse]->SetStats(0);
      for(Int_t nx=0; nx<NbinsT-1; nx++){
        if( VecCrMCP[nx][ChUse] > 0){
          Xbin.push_back( HistCrMCP[ChUse]->GetXaxis()->GetBinCenter(nx+1) );
          Wbin.push_back( VecCrMCP[nx][ChUse] );
        }
      }
      HistCrMCP[ChUse]->FillN( (Int_t)Xbin.size(), Xbin.data(), Wbin.data(), 1 );
      Xbin.clear();
      Wbin.clear();
      if(drawFirst){
        MaxMCP=HistCrMCP[ChUse]->GetMaximum();
        MinMCP=HistCrMCP[ChUse]->GetMinimum();
        drawFirst = kFALSE;
      }
      else{
        if( (HistCrMCP[ChUse]->GetMaximum() ) > MaxMCP )
          MaxMCP=HistCrMCP[ChUse]->GetMaximum();
        if( (HistCrMCP[ChUse]->GetMinimum() ) < MinMCP )
          MinMCP=HistCrMCP[ChUse]->GetMinimum();
      }
    }
    //
    delete [] XbinsT;
    Canvas  ->cd(2);
    gPad    ->cd(1);
    gPad    ->SetGridx(1);
    gPad    ->SetGridy(1);
    HistCrFC->SetLineColor(kBlack);
    HistCrFC->SetLineWidth(2);
    HistCrFC->Draw();
    //
    Canvas  ->cd(2);
    gPad    ->cd(2);
    gPad    ->SetGridx(1);
    gPad    ->SetGridy(1);
    HistCrEa->SetLineColor(kRed);
    HistCrEb->SetLineColor(kBlue);
    HistCrEa->SetLineWidth(2);
    HistCrEb->SetLineWidth(2);
    HistCrEa->GetYaxis()->SetRangeUser(MinEaEb*0.95,MaxEaEb*1.02);
    HistCrEa->Draw();
    HistCrEb->Draw("same");
    if(LegEaEb) delete LegEaEb;
    LegEaEb = new TLegend(0.10,0.10,0.18,0.25);
    LegEaEb->AddEntry(HistCrEa, "#Sigma E_{a}");
    LegEaEb->AddEntry(HistCrEb, "#Sigma E_{b}");
    LegEaEb->Draw();
    //
    Canvas  ->cd(2);
    gPad    ->cd(3);
    gPad    ->SetGridx(1);
    gPad    ->SetGridy(1);
    drawFirst = kTRUE;
    for(Int_t NoDE=0; NoDE<16; NoDE++){
      HistCrDE[NoDE]->SetLineColor(lcolor[NoDE]);
      HistCrDE[NoDE]->SetLineWidth(2);
      if( HistCrDE[NoDE]->Integral() > 0){
        if(drawFirst){
          HTitlBuf = TString::Format("Ts Dependency of Counting Rate for #Delta E - E-%02d [%s]",
              PjChUseE+1, RunName.Data() );
          HistCrDE[NoDE]->SetTitle(HTitlBuf.Data() );
          HistCrDE[NoDE]->GetYaxis()->SetRangeUser(MinDE*0.95, MaxDE*1.05);
          HistCrDE[NoDE]->Draw();
          drawFirst = kFALSE;
        }
        else
          HistCrDE[NoDE]->Draw("same");
      }
    }
    //
    Canvas  ->cd(2);
    gPad    ->cd(4);
    gPad    ->SetGridx(1);
    gPad    ->SetGridy(1);
    HistCrMCP[0]->SetTitle("Ts Dependency of Counting Rate for #Delta T (MWPC-MCP)");
    HistCrMCP[0]->SetLineColor(kMagenta);
    HistCrMCP[1]->SetLineColor(kRed);
    HistCrMCP[2]->SetLineColor(kBlue);
    HistCrMCP[3]->SetLineColor(kGreen);
    HistCrMCP[0]->SetLineWidth(2);
    HistCrMCP[1]->SetLineWidth(2);
    HistCrMCP[2]->SetLineWidth(2);
    HistCrMCP[3]->SetLineWidth(2);
    HistCrMCP[0]->GetYaxis()->SetRangeUser(0.0, MaxMCP*1.2);
    HistCrMCP[0]->Draw();
    HistCrMCP[1]->Draw("same");
    HistCrMCP[2]->Draw("same");
    HistCrMCP[3]->Draw("same");
    if(LegMCP) delete LegMCP;
    LegMCP = new TLegend(0.10,0.10,0.25,0.40);
    LegMCP->AddEntry(HistCrMCP[0], "MWPC-A & MCP-2");
    LegMCP->AddEntry(HistCrMCP[1], "MWPC-B & MCP-1");
    LegMCP->AddEntry(HistCrMCP[2], "MWPC-C & MCP-1");
    LegMCP->AddEntry(HistCrMCP[3], "MWPC-D & MCP-2");
    LegMCP->Draw();

    //
    Canvas->Update();
    Canvas->cd();
    ///////////////////
    // Check ErrFlag //
    ///////////////////
    if(ErrFlag>0) {cout << " ErrFlag > 0" << endl; break;}
    ////////////////////////////////
    // Compless the Raw-Data file //
    ////////////////////////////////
    //    gSystem->Exec( Form("bzip2 -9 %s", FileInName.Data() ) );
    ////////////////
    // Final File //
    ////////////////
    if(SSize <= SSizeMax) {cout << " SSize <= SSizeMax" << endl; break;}
  } // end of endless loop(1)
  /////////////////////////////////////////////////////
  // Fill the last event for each Mod&Ch to multimap //
  /////////////////////////////////////////////////////
  for(Int_t l=0; l<NMod; l++){ // 次の event が入力がなく保留となった event の処理
    for(Int_t m=0; m<16; m++){
      if( WF[l][m]==0 && TS[l][m]>0 && (HName[l][m] != "") ){
        // Fill histograms for almost all the event
        Hist[l][m]->Fill( (Double_t)(DataA) );
        HistTimeElapsed[l][m]->Fill( (Double_t)(DataT - TS[l][m]) );
        // data -> structure
        Dbuf.Mod  = (UShort_t)l;
        Dbuf.Ch   = (UShort_t)m;
        Dbuf.DataA= DA[l][m];
        // structure -> multimap (for all Mod & Ch)
        DmapA.insert( p(TS[l][m],Dbuf) );
        if( TsLastMod[l]<TS[l][m] )
          TsLastMod[l] = TS[l][m];
        // structure -> multimap (for trigger event only)
        if( (l==ModDE) && (EaDEUse[m] || EbDEUse[m] ) )
          DmapT.insert( p(TS[l][m],Dbuf) );
      }
    }
  }
  ///////////////////////////////////
  // remaining events in multimaps //
  ///////////////////////////////////
  itrTrig = DmapT.begin();
  while( itrTrig != DmapT.end() ){
    // Data of Trigger
    DataTT = itrTrig->first;
    ModT   = itrTrig->second.Mod;
    ChT    = itrTrig->second.Ch;
    DataAT = itrTrig->second.DataA;
    // setting up Time Window;
    TsWinMin   = DataTT + DTsMin;
    TsWinMax   = DataTT + DTsMax;
    itrAll     = DmapA.lower_bound(TsWinMin);
    if(itrAll != DmapA.begin() )
      --itrAll;
    itrLast    = DmapA.upper_bound(TsWinMax);
    if(itrLast!= DmapA.end()   )
      ++itrLast;
    // Initialize data for TTree
    NumDE  = 1; // Trigger
    NumEa  = 0;
    NumEb  = 0;
    NumDt  = 0;
    for(Int_t m=0; m<16; m++){
      AdcDE[m] = 0;
      AdcEa[m] = 0;
      AdcEb[m] = 0;
      if(m<4){
        QdcQ[  m]=0;
        TdcX[  m]=0;
        TdcY[  m]=0;
        TdcDt[ m]=0;
        TdcMCP[m]=0;
      }
    }
    for(Int_t m=0; m<17; m++){
      PhNdR[ m]=0;
      PhNdL[ m]=0;
      PsdNdR[m]=0;
      PsdNdL[m]=0;
      TofNdR[m]=0;
      TofNdL[m]=0;
    }
    AdcDE[ChT]=DataAT; // Trigger
    EventTrig++;
    EventTree[ModT][ChT]++;
    // Sort DmapA
    while(itrAll != itrLast){
      // multimap -> data
      DataT = itrAll->first;
      Mod   = itrAll->second.Mod;
      Ch    = itrAll->second.Ch;
      DataA = itrAll->second.DataA;
      TDiff = (Int_t)(DataT-DataTT);
      if( (TDiff >= DTsMin) && (TDiff <= DTsMax) ){
        // data -> data for TTree
        if( Mod==ModEa ){                                             // Ea
          HTDiffEa->Fill( (Double_t)(TDiff) );
          HTDiffEaEach[ChT]->Fill( (Double_t)(TDiff) );
          if( (TDiff>TsMinEa) && (TDiff<TsMaxEa) ){
            EventTree[Mod][Ch]++;
            NumEa++;
            AdcEa[Ch]=DataA;
          }
        }
        else if( Mod==ModEb ){                                        // Eb
          HTDiffEb->Fill( (Double_t)(TDiff) );
          HTDiffEbEach[ChT]->Fill( (Double_t)(TDiff) );
          if( (TDiff>TsMinEb) && (TDiff<TsMaxEb) ){
            EventTree[Mod][Ch]++;
            NumEb++;
            AdcEb[Ch]=DataA;
          }
        }
        else if( (Mod==ModDE) && (Ch>=ChFirstDE) && (Ch<=ChLastDE) ){ // DE
          if( Ch != ChT ){ // data of DE other than trigger ch
            HTDiffDE->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinDE) && (TDiff<TsMaxDE) ){
              EventTree[Mod][Ch]++;
              NumDE++;
              AdcDE[Ch]=DataA;
            }
          }
        }
        else if( (Mod==ModQ) && (Ch>=ChFirstQ) && (Ch<=ChLastQ) && 
            (HName[Mod][Ch] != "") ){                            //Q(MWPC)
              HTDiffQ->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinQ) && (TDiff<TsMaxQ) ){
                // Ch -> No. of Detector
                itrMWPC = ChQ.find(Ch);
                if(itrMWPC != ChQ.end() ){
                  EventTree[Mod][Ch]++;
                  ChUse = itrMWPC->second;
                  QdcQ[ChUse]= DataA;
                }
                else
                  cout << Form("  Mod=%02d Ch=%02d is not Pulse height of MWPC\n",Mod,Ch);
              }
            }
        else if( (Mod==ModPos) && (Ch>=ChFirstPos) && (Ch<=ChLastPos) && 
            (HName[Mod][Ch] != "") ){                            //Position(MWPC)
              HTDiffPos->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinPos) && (TDiff<TsMaxPos) ){
                // Ch -> No. of Detector
                itrMWPC = ChX.find(Ch);
                if(itrMWPC != ChX.end() ){
                  EventTree[Mod][Ch]++;
                  ChUse = itrMWPC->second;
                  TdcX[ChUse]= DataA;
                }
                else{ // itrMWPC == ChX.end()
                  itrMWPC = ChY.find(Ch);
                  if(itrMWPC != ChY.end() ){
                    EventTree[Mod][Ch]++;
                    ChUse = itrMWPC->second;
                    TdcY[ChUse]= DataA;
                  }
                  else
                    cout << Form("  Mod=%02d Ch=%02d is not Position of MWPC\n",Mod,Ch);
                }
              }
            }
        else if( (Mod==ModDt) && (Ch>=ChFirstDt) && (Ch<=ChLastDt) && 
            (HName[Mod][Ch] != "") ){                            //dT(MWPC)
              HTDiffDt->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinDt) && (TDiff<TsMaxDt) ){
                // Ch -> No. of Detector
                itrMWPC = ChDt.find(Ch);
                if(itrMWPC != ChDt.end() ){
                  EventTree[Mod][Ch]++;
                  NumDt++;
                  ChUse = itrMWPC->second;
                  TdcDt[ChUse]= DataA;
                }
                else
                  cout << Form("  Mod=%02d Ch=%02d is not dT of MWPC\n",Mod,Ch);
              }
            }
        else if( (Mod==ModMCP) && (Ch>=ChFirstMCP) && (Ch<=ChLastMCP) &&
            (HName[Mod][Ch] != "") ){                            //dT(MWPC-MCP)
              HTDiffMCP->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinMCP) && (TDiff<TsMaxMCP) ){
                // Ch -> No. of Detector
                itrMWPC = ChMCP.find(Ch);
                if(itrMWPC != ChMCP.end() ){
                  EventTree[Mod][Ch]++;
                  ChUse = itrMWPC->second;
                  TdcMCP[ChUse]= DataA;
                }
                else
                  cout << Form("  Mod=%02d Ch=%02d is not DeltaT (MWPC-MCP)\n", Mod,Ch);
              }
            }
        else if( (Mod==ModNdR) || ( Mod==(ModNdR+1) ) ||
            ( ( Mod==(ModNdR+2) ) && Ch<2 ) ){              // ND-R
          if(Ch % 2 == 0){ // Pulse Height
            HTDiffPH ->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
              // Number of Detector
              ChND = (Mod-ModNdR)*8 + Ch/2;
              // Fill Pulse Height
              PhNdR[ChND]=DataA;
            }
          }
          else{            // PSD
            HTDiffPSD->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPSD) && (TDiff<TsMaxPSD) ){
              // Number of Detector
              ChND = (Mod-ModNdR)*8 + (Ch-1)/2;
              // Fill PSD
              PsdNdR[ChND]=DataA;
            }
          }
        }   // else if( (Mod==ModNdR) ...)
        else if( ( (Mod==ModNdL)       && (Ch>1) ) ||
            ( Mod==(ModNdL+1) )         ||
            ( ( Mod==(ModNdL+2) ) && (Ch<4) ) ){            // ND-L
          if(Ch % 2 == 0){ // Pulse Height
            HTDiffPH ->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPH) && (TDiff<TsMaxPH) ){
              // Number of Detector
              ChND = (Mod-ModNdL)*8 + Ch/2-1;
              // Fill Pulse Height
              PhNdL[ChND]=DataA;
            }
          }
          else{            // PSD
            HTDiffPSD->Fill( (Double_t)(TDiff) );
            if( (TDiff>TsMinPSD) && (TDiff<TsMaxPSD) ){
              // Number of Detector
              ChND = (Mod-ModNdL)*8 + (Ch-1)/2-1;
              // Fill PSD
              PsdNdL[ChND]=DataA;
            }
          }
        }   // else if( (Mod==ModNdL)...)
        else if( ( Mod== ModTofR) ||
            ( Mod==(ModTofR+1) && Ch<1) ){               // TOF(ND-R)
              HTDiffTOF->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinTOF) && (TDiff<TsMaxTOF) ){
                // Number of Detector
                ChND = (Mod-ModTofR)*16+Ch;
                // Fill TOF
                TofNdR[ChND]=DataA;
              }
            }   // else if( ( Mod== ModTofR) ...)
        else if( ( Mod==(ModTofL)   && Ch>0) ||
            ( Mod==(ModTofL+1) && Ch<1) ){               // TOF(ND-L)
              HTDiffTOF->Fill( (Double_t)(TDiff) );
              if( (TDiff>TsMinTOF) && (TDiff<TsMaxTOF) ){
                // Number of Detector
                if(Mod==ModTofL)
                  ChND = Ch-1;
                else
                  ChND = 16;
                // Fill TOF
                TofNdL[ChND]=DataA;
              }
            }   // else if( ( Mod==(ModTofL) ....) )
      }
      // next event
      ++itrAll;
    }
    // Fill TTree
    if( (NumEa > 0) || (NumEb > 0) ){
      OutTree->Fill();
      EventStored++;
    }
    // next Trigger
    ++itrTrig;
  }
  if( !(DmapA.empty() ) ) DmapA.clear();
  if( !(DmapT.empty() ) ) DmapT.clear();
  ////////////////////////////////////////
  // End of RUN                         //
  // Display number of Time-Stamp error //
  ////////////////////////////////////////
  PrintErrors(1, FLog, Err, Event, Err2);
  //////////////////////
  // write histograms //
  //////////////////////
  for(Int_t l=0; l<NMod; l++){
    CdMod[l] ->cd();
    for(Int_t m=0; m<16; m++){
      if( (HName[l][m] != "") && (Event[l][m] != 0) ){
        // Write Histograms
        Hist[l][m]->Write();
        HistTimeElapsed[l][m]->Write();
      }
    }
  }
  // Histograms for Ts Difference
  CdTDiff  ->cd();
  HTDiffEa ->Write(); delete HTDiffEa;
  HTDiffEb ->Write(); delete HTDiffEb;
  HTDiffDE ->Write(); delete HTDiffDE;
  HTDiffQ  ->Write(); delete HTDiffQ;
  HTDiffPos->Write(); delete HTDiffPos;
  HTDiffDt ->Write(); delete HTDiffDt;
  HTDiffMCP->Write(); delete HTDiffMCP;
  HTDiffPH ->Write(); delete HTDiffPH;
  HTDiffPSD->Write(); delete HTDiffPSD;
  HTDiffTOF->Write(); delete HTDiffTOF;
  for(Int_t m=0; m<12; m++){
    HTDiffEaEach[m]->Write(); delete HTDiffEaEach[m];
    HTDiffEbEach[m]->Write(); delete HTDiffEbEach[m];
  }
  // Peak Positions in (E-DE)-DE histograms
  CdPeak->cd();
  for(Int_t m=0; m<16; m++){
    if( (HPeakEaDE[m]->Integral() ) > 0 )
      HPeakEaDE[m]->Write();
    delete HPeakEaDE[m];
    if( (HPeakEbDE[m]->Integral() ) > 0 )
      HPeakEbDE[m]->Write();
    delete HPeakEbDE[m];
    if( (HPeakDEEa[m]->Integral() ) > 0 )
      HPeakDEEa[m]->Write();
    delete HPeakDEEa[m];
    if( (HPeakDEEb[m]->Integral() ) > 0 )
      HPeakDEEb[m]->Write();
    delete HPeakDEEb[m];
  }
  // PID plot
  CdPID->cd();
  for(NoDE=0; NoDE<16; NoDE++){
    if( (HName[ModDE][NoDE] != "") && ( EaDEUse[NoDE] || EbDEUse[NoDE] ) && (Event[ModDE][NoDE] != 0) ){
      HistEDE[NoDE]->ResetStats();
      HistEDE[NoDE]->Write();
    }      
    delete  HistEDE[NoDE];
  }
  HistEDESum->ResetStats();
  HistEDESum->Write();
  Canvas->Write();
  Canvas->Close();
  delete LegPj;
  delete HistEDESum;
  // Counting Rate
  CdCr->cd();
  if( (HistCrEa->Integral() ) > 0 )
    HistCrEa->Write();
  delete HistCrEa;
  if( (HistCrEb->Integral() ) > 0 )
    HistCrEb->Write();
  delete HistCrEb;
  for(NoDE=0; NoDE<16; NoDE++){
    if( (HistCrDE[NoDE]->Integral() ) > 0 )
      HistCrDE[NoDE]->Write();
    delete HistCrDE[NoDE];
  }    
  ///////////////////////////////////////
  // Print # of Events stored in TTree //
  ///////////////////////////////////////
  PrintEventTree(FLog,EventTrig,EventStored,EventTree,Event);
  ///////////////
  // Stop time //
  ///////////////
  timer.Stop();
  RTime = timer.RealTime();
  CTime = timer.CpuTime();
  cout << TString::Format("\n Real Time = %d sec, CPU Time = %d sec\n",(Int_t)RTime,(Int_t)CTime) 
    << endl;
  FLog << TString::Format("\n Real Time = %d sec, CPU Time = %d sec\n",(Int_t)RTime,(Int_t)CTime)
    << endl;
  ///////////////////
  // Check ErrFlag //
  ///////////////////
  if(ErrFlag>0){
    printf("\033[7m");
    cout << "                                              \n"
      << "                                              \n"
      << "   #######  #######  ######    #####    ###   \n"
      << "   #     #  #     #  #     #  #     #   ###   \n"
      << "   #     #  #     #  #     #  #         ###   \n"
      << "   #     #  #     #  ######    #####     #    \n"
      << "   #     #  #     #  #              #         \n"
      << "   #     #  #     #  #        #     #   ###   \n"
      << "   #######  #######  #         #####    ###   \n"
      << "                                              \n"
      << "                                              \n\n";
    printf("\033[0m");
    FLog << "  *******************************************\n"
      << "  *                                         *\n"
      << "  *  Too much error of Time-Stamp occured!! *\n"
      << "  *                                         *\n"
      << "  *******************************************\n" << endl;
  }
  ////////////////////////
  // Close output files //
  ////////////////////////
  // histograms
  FHist->Close();
  delete FHist;   // Hist[l][m], CdMod[l], and so on are deallocated automatically
  // TTree
  FTree = OutTree->GetCurrentFile();
  //  FTree->Print();
  FTree->Write("",TObject::kOverwrite); // save only the latest version of the tree (by Mark-san)
  FTree->Close();
  delete FTree;
  // Log
  FLog.close();
  //////////////////////////////////////////////////////////////
  // Output some parameters will be used for further analysis //
  //////////////////////////////////////////////////////////////
  FPar.open( Form("%s/%s.para",PathOut.Data(), RunName.Data() ) );
  // Path for output
  FPar << PathOut.Data() << endl;
  // Name of RunName
  FPar << RunName.Data()     << endl;
  // Time Base (ns)
  FPar << TimeBase       << endl;
  // Ea(ChFist - ChLast) -> PID
  FPar << ChFirstEa      << endl;
  FPar << ChLastEa       << endl;
  // Eb(ChFist - ChLast) -> PID
  FPar << ChFirstEb      << endl;
  FPar << ChLastEb       << endl;
  // DE(ChFist - ChLast) -> PID
  FPar << ChFirstDE      << endl;
  FPar << ChLastDE       << endl;
  // combination of DE and E
  FPar << sEaDEUse.Data()<< endl;
  FPar << sEbDEUse.Data()<< endl;
  // close the file
  FPar.close();
  ////////////////////////////////////////
  // Deallocate storage space of arrays //
  ////////////////////////////////////////
  delete [] CdMod;
  delete [] Hist;
  for(Int_t l=0; l<NMod; l++){
    delete [] HName[l]; HName[l]=0;
    delete [] HTitl[l]; HTitl[l]=0;
    delete [] Event[l]; Event[l]=0;
    delete [] Err[  l]; Err[  l]=0;
    delete [] Err2[ l]; Err2[ l]=0;
    delete [] TS[   l]; TS[   l]=0;
    delete [] DA[   l]; DA[   l]=0;
    delete [] WF[   l]; WF[   l]=0;
  }
  delete [] HName;
  delete [] HTitl;
  delete [] Event;
  delete [] Err;
  delete [] Err2;
  delete [] TS;
  delete [] DA;
  delete [] WF;
  //
  exit(0);
  */
}
