#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TApplication.h"
#include <iostream>
#include <vector>
using namespace std;

#include "jshort.h"

void run_draw_pid(int runNo=528)
{
    int selectdE = 7;

    DetectorSetting* ds = DetectorSetting::GetDetectorSetting();

    bool pid_is_set = false;
    TCutG* pid_cut[17][8];
    auto file_pid = new TFile(Form("out/pid_de%d.root",selectdE),"read");
    if (file_pid->IsZombie()==false) {
        pid_is_set = true;
        pid_cut[selectdE][1] = (TCutG*) file_pid -> Get("pid1");
        pid_cut[selectdE][2] = (TCutG*) file_pid -> Get("pid2");
        pid_cut[selectdE][3] = (TCutG*) file_pid -> Get("pid3");
        pid_cut[selectdE][4] = (TCutG*) file_pid -> Get("pid4");
        pid_cut[selectdE][5] = (TCutG*) file_pid -> Get("pid5");
        pid_cut[selectdE][1] -> SetLineColor(kOrange);
        pid_cut[selectdE][2] -> SetLineColor(kOrange);
        pid_cut[selectdE][3] -> SetLineColor(kOrange);
        pid_cut[selectdE][4] -> SetLineColor(kOrange);
        pid_cut[selectdE][5] -> SetLineColor(kOrange);
    }

    TString inputFileName = Form("out/RUN%d.summary.root",runNo);
    TFile* file_summary = new TFile(inputFileName);
    TTree* tree = (TTree*) file_summary -> Get("event");

    TClonesArray *channelArray = nullptr;
    double de, ee;
    tree -> SetBranchAddress("de",&de);
    tree -> SetBranchAddress("ee",&ee);
    tree -> SetBranchAddress("ch",&channelArray);

    TH1D* histdECount = new TH1D("histdECount",inputFileName + ";dE channel number;count",17,0,17);
    TH1D* histS1Count = new TH1D("histS1Count",inputFileName + ";S1J channel number;count",33,0,33);
    histdECount -> SetFillColor(29);
    histS1Count -> SetFillColor(29);
    histdECount -> GetXaxis() -> SetNdivisions(117);
    histS1Count -> GetXaxis() -> SetNdivisions(217);

    vector<int> dEChannelArray = {7,8,9,11};
    TH2D* histPID[17];
    for (auto iCh=0; iCh<dEChannelArray.size(); ++iCh) {
        auto dEChannel = dEChannelArray[iCh];
        TString title = Form("dE(%d);de+s1;de",dEChannel);
        histPID[dEChannel] = new TH2D(Form("histPID_dE%d",dEChannel) ,title,200,0,30,200,0,15);
    }

    TH2D* histEVSAngle[7];
    for (auto ipid=0; ipid<7; ++ipid)
        histEVSAngle[ipid] = new TH2D(Form("histEVSAngle_pid%d",ipid),Form("PID cut %d;angle (deg);de+s1",ipid),25,20,45,200,0,30);

    Long64_t numEvents = tree -> GetEntries();
    for (Long64_t iEvent=0; iEvent<numEvents; ++iEvent)
    {
        tree -> GetEntry(iEvent);
        int numChannels = channelArray -> GetEntries();
        int founddE = -1;
        int foundS1 = -1;
        if (de>0 && ee>0)
        {
            for (int iChannel=0; iChannel<numChannels; ++iChannel)
            {
                ChannelData* channel = (ChannelData*) channelArray -> At(iChannel);

                if (channel->det==kdE) {
                    histdECount -> Fill(channel->dch);
                    for (auto iCh=0; iCh<dEChannelArray.size(); ++iCh) {
                        auto dEChannel = dEChannelArray[iCh];
                        if (dEChannel==channel->dch) {
                            if (founddE==-1) founddE = dEChannel;
                            else if (founddE>=0) founddE = -2;
                        }
                    }
                }
                if (channel->det==kS1J) {
                    histS1Count -> Fill(channel->dch);
                    if (foundS1==-1) foundS1 = channel->dch;
                    else if (foundS1>=0) foundS1 = -2;
                }
            }

            if (foundS1>0 && founddE>0)
            {
                histPID[founddE] -> Fill(ee,de);
                double angle1 = ds->S1Ch_Angle1(foundS1);
                double angle2 = ds->S1Ch_Angle2(foundS1);
                histEVSAngle[0] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                if (pid_is_set) {
                    if (pid_cut[selectdE][1]->IsInside(ee,de)) histEVSAngle[1] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                    if (pid_cut[selectdE][2]->IsInside(ee,de)) histEVSAngle[2] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                    if (pid_cut[selectdE][3]->IsInside(ee,de)) histEVSAngle[3] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                    if (pid_cut[selectdE][4]->IsInside(ee,de)) histEVSAngle[4] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                    if (pid_cut[selectdE][5]->IsInside(ee,de)) histEVSAngle[5] -> Fill(gRandom->Uniform(angle1,angle2),ee);
                }
            }
        }
    }

    auto DrawHistBar = [](TH1D* hist)
    {
        TAxis *ax = hist -> GetXaxis();
        auto nx = ax -> GetNbins();
        for (auto binx=1; binx<nx; ++binx)
        {
            auto val = hist -> GetBinContent(binx);
            if (val>0)
            {
                auto x2 = ax -> GetBinUpEdge(binx);
                auto line = new TLine(x2,0,x2,val);
                line -> SetLineColor(hist->GetLineColor());
                line -> Draw("samel");
            }
        }
    };

    auto cvs = new TCanvas("cvs","",1200,500);
    cvs -> Divide(2,1);
    cvs -> cd(1); histdECount -> Draw(""); DrawHistBar(histdECount);
    cvs -> cd(2); histS1Count -> Draw(""); DrawHistBar(histS1Count);
    cvs -> SaveAs("figures/figure_channel_count.png");

    auto cvs_de = new TCanvas("cvs_de","",1200,1000);
    cvs_de -> cd(1); histPID[selectdE] -> Draw("colz");
    if (pid_is_set) {
        pid_cut[selectdE][1] -> Draw("samel");
        pid_cut[selectdE][2] -> Draw("samel");
        pid_cut[selectdE][3] -> Draw("samel");
        pid_cut[selectdE][4] -> Draw("samel");
        pid_cut[selectdE][5] -> Draw("samel");
    }
    cvs_de -> SaveAs(Form("figures/figure_dE%d_pid.png",selectdE));

    auto cvs_ea = new TCanvas("cvs_ea","",1800,1000);
    cvs_ea -> Divide(3,2);
    cvs_ea -> cd(1); histEVSAngle[0] -> Draw("colz");
    cvs_ea -> cd(2); histEVSAngle[1] -> Draw("colz");
    cvs_ea -> cd(3); histEVSAngle[2] -> Draw("colz");
    cvs_ea -> cd(4); histEVSAngle[3] -> Draw("colz");
    cvs_ea -> cd(5); histEVSAngle[4] -> Draw("colz");
    cvs_ea -> cd(6); histEVSAngle[5] -> Draw("colz");
    cvs_ea -> SaveAs(Form("figures/figure_dE%d_E_vs_Angle.png",selectdE));

    TString outputFileName = Form("out/RUN%d.analysis.root",runNo);
    auto file_out = new TFile(outputFileName,"recreate");
    cvs -> Write();
    cvs_de -> Write();
    cvs_ea -> Write();
    histdECount -> Write();
    histS1Count -> Write();
    pid_cut[selectdE][1] -> Write();
    pid_cut[selectdE][2] -> Write();
    pid_cut[selectdE][3] -> Write();
    pid_cut[selectdE][4] -> Write();
    pid_cut[selectdE][5] -> Write();
    histEVSAngle[0] -> Write();
    histEVSAngle[1] -> Write();
    histEVSAngle[2] -> Write();
    histEVSAngle[3] -> Write();
    histEVSAngle[4] -> Write();
    histEVSAngle[5] -> Write();
    cout << outputFileName << endl;
}
