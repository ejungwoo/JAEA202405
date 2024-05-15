#ifndef LICD2IRRAD_SORTSI_H
#define LICD2IRRAD_SORTSI_H

/********************************************
 * Header file for SortSi                   *
 ********************************************/
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TH1D.h>
#include <TH2I.h>
#include <TH2D.h>
#include <TFile.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TPolyMarker.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm> // need for std::min_element

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<unsigned short>+;
#endif


////////////////////////////
// Path for Raw-Data file //
////////////////////////////
const TString PathIn = "/home/daquser/data/LiCD2Irrad";
//////////////////////////////////
// Number Of Modules Used (<20) //
//////////////////////////////////
const UShort_t   NMod = 16;
////////////////////////////////
// Max. size of Raw-Data file //
////////////////////////////////
const streamsize SSizeMax=500000000; // 500 MB
////////////////////
// Time Base (ns) //
////////////////////
const Double_t TimeBase = 200.;
///////////
// Timer //
///////////
TStopwatch timer;
Double_t   RTime, CTime;
///////////////////////////
// Path for output files //
///////////////////////////
const TString PathOut= "/home/daquser/data/LiCD2Irrad/SortSi/out";
//////////////
// Log File //
//////////////
ofstream FLog;
ofstream FPar;
//////////////////////////////
// Definition of Mod and Ch //
//////////////////////////////
const TString DefFileName = "/home/daquser/data/LiCD2Irrad/SortSi/ModCh.in";
ifstream      DefFile;
/////////////////////////////////////////
// Raw histograms and output root file //
/////////////////////////////////////////
TH1D       ***Hist;
TString       HNameBuf, HTitlBuf;
TFile        *FHist;
string        BufTmp;
TDirectory  **CdMod;
TString     **HName;
TString     **HTitl;
/////////////////////////////////////////////////////
// the time that has elapsed after previous events //
/////////////////////////////////////////////////////
TH1D ***HistTimeElapsed;
//////////////////////////////////////////
// Number of events, errors, and so on. //
//////////////////////////////////////////
Long64_t **Event, **Err, **Err2, **EventTree, EventTrig, EventStored;
Double_t **ErrR;
////////////////////////////////////////////////////
// Prameters for Operation of input Raw-Data file //
////////////////////////////////////////////////////
UShort_t   FOpen = 0;           // Flag for file open
UShort_t   FNum  = 0;           // Number of Files read
Int_t      ErrFlag = 0;         // Error Flag
TString    FileInName;          // Name of Input Raw-Data File
ifstream   FileIn;              // input stream for Raw-Data file
streamsize SSize, SSizeOld = 0; // Size of opened Raw-Data file
// Parameters to switch the EffFlag on
const Double_t MinSameErr = 1.; // 1 % for Same Ts
const Int_t    MinDecErr  = 10; // 10 times for decrease of Ts
// Data buffer
char     Buf[256];
Int_t    NBuf = 0;
UShort_t Mod, Ch, DataA, Evt, **DA, **WF;
Long64_t DataT;
Long64_t **TS;

///////////////////////////////////////////////
// Max. # of Peaks to be found in histograms //
///////////////////////////////////////////////
const Int_t NPeaks = 10;
/////////////////////////////////////
// Peak positions for Si detectors //
/////////////////////////////////////
// Si Total-Ea
const UShort_t ModEa    =    0;
const UShort_t ChFirstEa=    1; // 0-15->Ea01-Ea16
const UShort_t ChLastEa =   10; // Ea(ChFist - ChLast) -> PID
// Si Total-Eb
const UShort_t ModEb =       1;
const UShort_t ChFirstEb=    2; // 0-15->Eb01-Eb16
const UShort_t ChLastEb =   12; // Eb(ChFist - ChLast) -> PID
// Si Delta-E
const UShort_t ModDE    =    2;
const UShort_t ChFirstDE=    0; // 0-15->DE01-DE16
const UShort_t ChLastDE =   11;
///////////////////////
// for time plotting //
///////////////////////
TString   Year;
TString   Month;
TString   Day;
TString   Hour;
TString   Minute;
TString   TimeFormat("#splitline{%m/%d}{%H:%M}");
Int_t     NbinsT;
Double_t *XbinsT;
Double_t *YbinsT;
vector<Double_t> Xbin,Ybin, Wbin;
////////////////////////////////////////////////////////
// Peak positions for (E-DE)-DE histograms & PID plot //
////////////////////////////////////////////////////////
Int_t NoDE, NoEa, NoEb;
//                        1234567890123456
const TString sEaDEUse = "1111010000000000";
const TString sEbDEUse = "0000001111010000";
Bool_t         EaDEUse[16];
Bool_t         EbDEUse[16];
const Int_t    NbinsChX = 1250;
const Double_t XLowCh   =  502.;  // 4 ch / bin
const Double_t XUpCh    = 5502.;
const Double_t XLowPIDCh= 0.1505;
const Double_t XUpPIDCh = 1.4005;
const Int_t    NbinsChY =  800;
const Double_t YLowCh   =  402.;
const Double_t YUpCh    = 3602.;  // 4 ch / bin
const Double_t YLowPIDCh= 0.201;
const Double_t YUpPIDCh = 1.801;
TCanvas       *Canvas;
const TString  NameCdEDEFirst = "E_DE1stFile";
TDirectory    *CdEDEFirst;
const TString  NameCdPeak     = "PeakPositions";
TDirectory    *CdPeak;
const TString  NameCdPID      = "PID";
TDirectory    *CdPID;
TLegend       *LegPj;
TH2I *histEaDE[ 16][16];// [NoEa][NoDE] : Ea vs DE histograms
TH2I *histEbDE[ 16][16];// [NoEb][NoDE] : Eb vs DE histograms
TH2D *HistEDE[  16];    // [NoDE]       : Ea or Eb vs DE(sum) histogram (PID)
TH2D *HistEDESum;       // sum of HistEDE[0-15]
TH2D *HPeakEaDE[16];    // [NoEa]       : Peak Positions for Ea (X: Ts, Y: NoDE)
TH2D *HPeakEbDE[16];    // [NoEb]       : Peak Positions for Eb (X: Ts, Y: NoDE)
TH2D *HPeakDEEa[16];    // [NoDE]       : Peak Positions for DE (X: Ts, Y: NoEa)
TH2D *HPeakDEEb[16];    // [NoDE]       : Peak Positions for DE (X: Ts, Y: NoEb)
TH1D *PjHPeakDE[16];    // [NoDE]       : Projection of HPeakDEEa or HPeakDEEb
const UShort_t PjChUseE = 2; // 0-15->Ea(b)01-Ea(b)16, which Ea (or Eb) used for HPeakDEE? 
vector<Long64_t> VecTS(1,0); // [FNum]
vector< vector< vector<Double_t> > > VecPeakEaDE; // [FNum-1][NoEa][NoDE]: Peak Pos. for Ea
vector< vector< vector<Double_t> > > VecPeakEbDE; // [FNum-1][NoEb][NoDE]: Peak Pos. for Eb
vector< vector< vector<Double_t> > > VecPeakDEEa; // [Fnum-1][NoDE][NoEa]: Peak Pos. for DE(w/Ea)
vector< vector< vector<Double_t> > > VecPeakDEEb; // [Fnum-1][NoDE][NoEb]: Peak Pos. for DE(w/Eb)
vector< vector<Double_t> > VecEaDETmp( 16, vector<Double_t>(16,0) ); // [NoEa][NoDE]
vector< vector<Double_t> > VecEbDETmp( 16, vector<Double_t>(16,0) ); // [NoEb][NoDE]
vector< vector<Double_t> > VecDEEaTmp( 16, vector<Double_t>(16,0) ); // [NoDE][NoEa]
vector< vector<Double_t> > VecDEEbTmp( 16, vector<Double_t>(16,0) ); // [NoDE][NoEa]
/////////////////////////////////
// sort of data using multimap //
/////////////////////////////////
struct ModChData                       {UShort_t Mod, Ch, DataA;};
struct ModChData                       Dbuf;
multimap<Long64_t,ModChData>           DmapA, DmapT;
multimap<Long64_t,ModChData>::iterator itrAll,itrTrig, itrLast;
typedef pair<Long64_t,ModChData>       p;
vector< Long64_t >                     TsLastMod;
Long64_t                               LastTsSafe, DataTT, TsWinMin, TsWinMax;
UShort_t                               ModT, ChT, DataAT;
const Long64_t                         DTsMin = -1000; // (* TimeBase ns)
const Long64_t                         DTsMax = +1000;
Int_t                                  TDiff;
/////////////////////////////////////////////////////////////
// Histograms for Ts Difference between Trigger and others //
/////////////////////////////////////////////////////////////
const TString  NameCdTDiff = "TDiff";
TDirectory    *CdTDiff;
TH1D          *HTDiffEa; // Ea - DE
TH1D          *HTDiffEaEach[16];
TH1D          *HTDiffEb; // Eb - DE
TH1D          *HTDiffEbEach[16];
TH1D          *HTDiffDE; // DE - DE
TH1D          *HTDiffQ;  // Q(MWPC)       - DE
TH1D          *HTDiffPos;// Position(MWPC)- DE
TH1D          *HTDiffDt; // Dt(MWPC)      - DE
TH1D          *HTDiffMCP;// dT(MCP)       - DE
TH1D          *HTDiffPH; // PH(ND)        - DE
TH1D          *HTDiffPSD;// PSD(ND)       - DE
TH1D          *HTDiffTOF;// TOF(ND)       - DE
//////////////////////////////////
// Channel assignment for MWPCs //
//////////////////////////////////
UShort_t         ChUse;
const UShort_t   ModQ      =  3;
const UShort_t   ChFirstQ  =  0; // 0-15
const UShort_t   ChLastQ   =  4;
const TString   sUseChQ    = "00010304"; // Channels of [Q1=00][Q2=01][Q3=03][Q4=04]
map<Int_t,Int_t> ChQ;
//
const UShort_t   ModPos    =  4;
const UShort_t   ChFirstPos=  1; // 0-15
const UShort_t   ChLastPos = 11;
const TString   sUseChPosX = "01030609"; // Channels of [X1][X2][X3][X4]
const TString   sUseChPosY = "02050811"; // Channels of [Y1][Y2][Y3][Y4]
map<Int_t,Int_t> ChX,ChY;
//
const UShort_t   ModDt     =  4;
const UShort_t   ChFirstDt = 13; // 0-15
const UShort_t   ChLastDt  = 15;
const TString   sUseChDt   = "131415";   // Channels of [A-B][B-C][A-D]
map<Int_t,Int_t> ChDt;
//
const UShort_t   ModMCP    =  8;
const UShort_t   ChFirstMCP= 12; // 0-15
const UShort_t   ChLastMCP = 15;
const TString   sUseChMCP  = "12131415"; // Channels of [MWPC1-MCP2][MWPC2-MCP1][MWPC3-MCP1][MWPC4-MCP2]
map<Int_t,Int_t> ChMCP;
map<Int_t,Int_t>::iterator itrMWPC;
////////////////////////////////////////////
// Channel assiment for Neutron Detectors //
////////////////////////////////////////////
const UShort_t  ModNdR =  9; // ND-R:  9[0-15],10[0-15],11[0,1]
const UShort_t  ModNdL = 11; // ND-L: 11[2-15],12[0-15],13[2,3]
const UShort_t  ModTofR=  6; // ND-R:  6[0-15], 7[0]
const UShort_t  ModTofL=  7; // ND-L:  7[1-15], 8[0]
Int_t           ChND;
////////////////////
// Counting Rates //
////////////////////
const UShort_t ModFC = 4;
const UShort_t ChFC  =12;
Double_t     **EachEvent;
Long64_t       EachTsFirst, EachTsLast;
Double_t       CountSum, DeltaTs;
vector<Double_t>           VecCrEa;          // Ea [FNum-1]
vector<Double_t>           VecCrEb;          // Eb [FNum-1]
vector< vector<Double_t> > VecCrDE;          // DE [FNum-1][NoDE]
vector<Double_t>           VecCrDETmp(16,0); // DE [NoDE]
vector<Double_t>           VecCrFC;          // FC [FNum-1]
vector< vector<Double_t> > VecCrQ;           // Q  [FNum-1][Ch]
vector<Double_t>           VecCrQTmp(4,0);   // Q  [Ch]
vector< vector<Double_t> > VecCrX;           // X  [FNum-1][Ch]
vector<Double_t>           VecCrXTmp(4,0);   // X  [Ch]
vector< vector<Double_t> > VecCrY;           // Y  [FNum-1][Ch]
vector<Double_t>           VecCrYTmp(4,0);   // Y  [Ch]
vector< vector<Double_t> > VecCrDt;          // Dt  [FNum-1][Ch] (MWPC-MWPC)
vector<Double_t>           VecCrDtTmp(4,0);  // Dt  [Ch]
vector< vector<Double_t> > VecCrMCP;         // MCP [FNum-1][Ch] (MWPC-MCP)
vector<Double_t>           VecCrMCPTmp(4,0); // MCP [Ch]
vector<Double_t>           VecCrPhR;         // PH  [FNum-1] (ND-R)
vector<Double_t>           VecCrPhL;         // PH  [FNum-1] (ND-L)
vector<Double_t>           VecCrPsdR;        // PSD [FNum-1] (ND-R)
vector<Double_t>           VecCrPsdL;        // PSD [FNum-1] (ND-L)
vector<Double_t>           VecCrTOFR;        // TOF [FNum-1] (ND-R)
vector<Double_t>           VecCrTOFL;        // TOF [FNum-1] (ND-L)
TH1D       *HistCrEa;                        // Ea
TH1D       *HistCrEb;                        // Eb
TH1D       *HistCrDE[16];                    // DE [NoDE]
TH1D       *HistCrFC;                        // FC
TH1D       *HistCrDt[4];                     // Dt (MWPC-MWPC)
TH1D       *HistCrMCP[4];                    // Dt (MWPC-MCP)
TDirectory *CdCr;
TString     NameCdCr = "CountingRates";
TLegend    *LegEaEb;
TLegend    *LegDt;
TLegend    *LegMCP;
Double_t    MaxEaEb, MinEaEb;
Double_t    MaxDE,   MinDE;
Double_t    MaxDt,   MinDt;
Double_t    MaxMCP,  MinMCP;

///////////////////////////////
// Ts Window for coincidence //
///////////////////////////////
const Long64_t TsMinEa   = -1; // -1 < TsEa <  4
const Long64_t TsMaxEa   =  4;
const Long64_t TsMinEb   = -1; // -1 < TsEb <  4
const Long64_t TsMaxEb   =  4;
const Long64_t TsMinDE   = -3; // -3 < TsDE <  3
const Long64_t TsMaxDE   =  3;
const Long64_t TsMinQ    = -5; // -5 < TsQ  < -1
const Long64_t TsMaxQ    = -1;
const Long64_t TsMinPos  = -5; // -5 < TsPos<  0
const Long64_t TsMaxPos  =  0;
const Long64_t TsMinDt   =-12; //-12 < TsDt < -1
const Long64_t TsMaxDt   = -1;
const Long64_t TsMinMCP  = -5; // -5 < TsMCP< -1
const Long64_t TsMaxMCP  = -1;
const Long64_t TsMinPH   =  0; //  0 < TsPH < 5
const Long64_t TsMaxPH   =  5;
const Long64_t TsMinPSD  =  0; //  0 < TsPSD< 4
const Long64_t TsMaxPSD  =  4;
const Long64_t TsMinTOF  = -6; // -6 < TsTOF< 1
const Long64_t TsMaxTOF  =  1;
//////////////////////////
// Parameters for TTree //
//////////////////////////
TFile *FTree;
TTree *OutTree;
const TString NameOutTree = "SortSi02"; // use same selector for SortSi01.C
// Branch
//Long64_t DataTT;     // Ts for Trigger event
UShort_t NumDE;      // Number of Si DE events          in Trigger Ts Window
UShort_t NumEa;      // Number of strips in Si Total-Ea in Trigger Ts Window
UShort_t NumEb;      // Number of strips in Si Total-Eb in Trigger Ts Window
UShort_t NumDt;      // Number of Dt events             in Trigger Ts Window
UShort_t AdcDE[16];  // ADC values for DE
UShort_t AdcEa[16];  // ADC values for Ea
UShort_t AdcEb[16];  // ADC values for Eb
UShort_t QdcQ[  4];  // ADC values for Q (MWPC)
UShort_t TdcX[  4];  // TDC values for X (MWPC)
UShort_t TdcY[  4];  // TDC values for Y (MWPC)
UShort_t TdcDt[ 4];  // TDC values for Dt(MWPC)
UShort_t TdcMCP[4];  // TDC values for Dt(MWPC-MCP)
UShort_t PhNdR[17];  // Pulse Height of ND (R)
UShort_t PhNdL[17];  // Pulse Height of ND (L)
UShort_t PsdNdR[17]; // PSD of ND (R)
UShort_t PsdNdL[17]; // PSD of ND (L)
UShort_t TofNdR[17]; // TOF of ND (R)
UShort_t TofNdL[17]; // TOF of ND (L)

#endif
