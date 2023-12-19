#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "THn.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
#include <vector>

#endif

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

struct StJetTreeStruct
{
    float d0z;
    float d0mass;
    float d0pt;
    float d0r;
    float pt;
    float ptcorr;
    float phi;
    float eta;
    float centrality;
    float lambda[4];
};

const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

const Int_t nCentralityBins = 3;

map<Int_t, Int_t> centralityMap = {
    {8, 0}, // 8->   0-10%
    {7, 1}, // 7->  10-20%
    {6, 1}, // 6->  20-30%
    {5, 1}, // 5->  30-40%
    {4, 2}, // 4->  40-50%
    {3, 2}, // 3->  50-60%
    {2, 2}, // 2->  60-70%
    {1, 2}  // 1->  70-80%
};

vector<Double_t> getAxisVector(TAxis *axis)
{
    vector<Double_t> axisVector;
    for (Int_t i = 1; i <= axis->GetNbins(); i++)
    {
        axisVector.push_back(axis->GetBinLowEdge(i));
    }
    axisVector.push_back(axis->GetBinUpEdge(axis->GetNbins()));
    return axisVector;
}

TH2D *invertAxis(TH2D *hist)
{

    vector<Double_t> xAxisVector = getAxisVector(hist->GetXaxis());
    vector<Double_t> yAxisVector = getAxisVector(hist->GetYaxis());

    Int_t nBinsX = xAxisVector.size() - 1;
    Int_t nBinsY = yAxisVector.size() - 1;

    TH2D *histInv = new TH2D(hist->GetName(), hist->GetTitle(), nBinsY, &yAxisVector[0], nBinsX, &xAxisVector[0]);
    for (Int_t iBinX = 1; iBinX <= nBinsX; iBinX++)
    {
        for (Int_t iBinY = 1; iBinY <= nBinsY; iBinY++)
        {
            histInv->SetBinContent(iBinY, iBinX, hist->GetBinContent(iBinX, iBinY));
            histInv->SetBinError(iBinY, iBinX, hist->GetBinError(iBinX, iBinY));
        }
    }
    return histInv;
}

void NormalizeByBinWidth(TH1D *hist, const Int_t color)
{
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        if (hist->GetBinWidth(i) != 0)
        {
            hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
            hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
        }
    }
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
}

double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order
Int_t getCentralityBin(const Float_t &centrality)
{
    for (Int_t i = 0; i < nCentralityBins; i++)
    {
        if (centrality >= centBins[i] && centrality < centBins[i + 1])
            return i;
    }
    cout << "Error: centrality not in range" << endl;
    return -1;
}

// https://www.star.bnl.gov/protected/lfsupc/tdrk/Centrality/Run19AuAu200/top20_tables/table_Ncoll_vs_centrality_systematicerror.txt
const Double_t Ncoll[nCentralityBins] = {952., 397., 58.};

TString centralityTitles[nCentralityBins] = {"0-10%", "10-40%", "40-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_40", "40_80"};

void plotComparison(TCanvas *can, TH2D *hUnfolded, TH2D *hRealData, TH2D *hMc, TH2D *hMcMeasured, const Int_t &iCent, TString var)
{

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    can->Clear();

    can->Divide(2, 2);
    TH1D *hUnfoldedProjX = (TH1D *)hUnfolded->ProjectionX("hUnfoldedProjX");
    TH1D *hRealDataProjX = (TH1D *)hRealData->ProjectionX("hRealDataProjX");
    TH1D *hMcProjX = (TH1D *)hMc->ProjectionX("hMcProjX");
    TH1D *hMcMeasuredProjX = (TH1D *)hMcMeasured->ProjectionX("hMcMeasuredProjX");
    NormalizeByBinWidth(hUnfoldedProjX, 2001);
    NormalizeByBinWidth(hRealDataProjX, 2002);
    NormalizeByBinWidth(hMcProjX, 2003);
    NormalizeByBinWidth(hMcMeasuredProjX, 2004);

    TLegend *leg1 = new TLegend(0.2, 0.8, 0.3, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    TLegend *leg2 = new TLegend(0.8, 0.8, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg1->AddEntry(hUnfoldedProjX, "Unfolded", "l");
    leg2->AddEntry(hRealDataProjX, "Real Data", "l");
    leg1->AddEntry(hMcProjX, "Mc", "l");
    leg2->AddEntry(hMcMeasuredProjX, "Mc Reco", "l");
    // first draw and compare x projections of hUnfolded and hMc
    can->cd(1);
    gPad->SetLogy();
    hUnfoldedProjX->GetYaxis()->SetTitle("dN/dp_{t}");
    hUnfoldedProjX->Draw();
    hMcProjX->Scale(hUnfoldedProjX->GetMaximum() / hMcProjX->GetMaximum());
    hMcProjX->Draw("same");

    // first draw and compare x projections of hRealData and hMcMeasured
    can->cd(2);
    gPad->SetLogy();

    hRealDataProjX->GetYaxis()->SetTitle("dN/dp_{t}");
    hRealDataProjX->GetXaxis()->SetTitle("p_{t}, GeV/c");
    hRealDataProjX->Draw();
    hMcMeasuredProjX->Scale(hRealDataProjX->GetMaximum() / hMcMeasuredProjX->GetMaximum());
    hMcMeasuredProjX->Draw("same");

    TH1D *hUnfoldedProjY = (TH1D *)hUnfolded->ProjectionY("hUnfoldedProjY");
    TH1D *hRealDataProjY = (TH1D *)hRealData->ProjectionY("hRealDataProjY");
    TH1D *hMcProjY = (TH1D *)hMc->ProjectionY("hMcProjY");
    TH1D *hMcMeasuredProjY = (TH1D *)hMcMeasured->ProjectionY("hMcMeasuredProjY");
    NormalizeByBinWidth(hUnfoldedProjY, 2001);
    NormalizeByBinWidth(hRealDataProjY, 2002);
    NormalizeByBinWidth(hMcProjY, 2003);
    NormalizeByBinWidth(hMcMeasuredProjY, 2004);
    // first draw and compare y projections of hUnfolded and hMc
    can->cd(3);
    gPad->SetLogy();

    hUnfoldedProjY->GetYaxis()->SetTitle("dN/d" + var);
    hUnfoldedProjY->Draw();
    hMcProjY->Scale(hUnfoldedProjY->GetMaximum() / hMcProjY->GetMaximum());
    hMcProjY->Draw("same");

    // first draw and compare y projections of hRealData and hMcMeasured
    can->cd(4);
    gPad->SetLogy();
    hRealDataProjY->GetYaxis()->SetTitle("dN/d" + var);
    hRealDataProjY->GetXaxis()->SetTitle(var);
    hRealDataProjY->Draw();
    hMcMeasuredProjY->Scale(hRealDataProjY->GetMaximum() / hMcMeasuredProjY->GetMaximum());
    hMcMeasuredProjY->Draw("same");

    can->cd();
    leg1->Draw("same");
    leg2->Draw("same");
    tex->DrawLatex(0.2, 0.65, "p_{t}");

    tex->DrawLatex(0.2, 0.25, centralityTitles[iCent] + "  " + var);

    can->SaveAs("rcp.pdf");
}

void assignTree(TTree *jetTree, StJetTreeStruct &jet);

void calculateRcp()
{
    Bool_t readTree = kFALSE;
    gSystem->Load("libRooUnfold");

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    // Read response matrix from file
    TFile *responseFile = new TFile("../unfold/responseMc5Reco10.root", "READ");
    responseFile->cd();

    RooUnfoldResponse *response[nCentralityBins];
    RooUnfoldResponse *responseAngularity[nCentralityBins][nAngularities];

    TH2D *hMeasured[nCentralityBins];
    TH2D *hUnfolded[nCentralityBins];
    TH2D *hMc[nCentralityBins];
    TH2D *hMcAngularity[nCentralityBins][nAngularities];

    TH2D *hMcMeasured[nCentralityBins];
    TH2D *hMcMeasuredAngularity[nCentralityBins][nAngularities];

    TH2D *hMeasuredAngularity[nCentralityBins][nAngularities];
    TH2D *hUnfoldedAngularity[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        response[iCent] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent]);
        hMcMeasured[iCent] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/MeasTest" + centralityNames[iCent]);
        hMc[iCent] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/TrueTest" + centralityNames[iCent]);
        // hMeasured[iCent] = new TH2D("Meas" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsDetector.size() - 1, &ptBinsDetector[0], zBinsDetector.size() - 1, &zBinsDetector[0]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            responseAngularity[iCent][iLambda] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent] + angularityNames[iLambda]);
            hMcMeasuredAngularity[iCent][iLambda] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/MeasTest" + centralityNames[iCent] + angularityNames[iLambda]);
            hMcAngularity[iCent][iLambda] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/TrueTest" + centralityNames[iCent] + angularityNames[iLambda]);

            //  hMeasuredAngularity[iCent][iLambda] = new TH2D("Meas" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsDetector.size() - 1, &ptBinsDetector[0], angularityBinsDetector[iLambda].size() - 1, &angularityBinsDetector[iLambda][0]);
        }
    }

    // TFile *treeFile; // Open the file containing the tree.
    //  treeFile = new TFile("../D0_jets_2014_231030.root", "READ");
    // //treeFile = new TFile("../output_jets.root", "READ");
    // if (!treeFile || treeFile->IsZombie())
    // {
    //     return;
    // }
    // TTree *jetTree = (TTree *)treeFile->Get("Jets");
    // StJetTreeStruct jet;
    // assignTree(jetTree, jet);
    // Long_t nEntries = jetTree->GetEntries() / 1.;

    // cout << "nEntries = " << (Float_t)nEntries / 1000. << "k" << endl
    //      << endl;

    // for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    // {
    //     Float_t progress = 0.;
    //     progress = (Float_t)iEntry / (Float_t)nEntries;
    //     if (iEntry % 1000 == 0)
    //     {
    //         cout << "\r (" << (progress * 100.0) << "%)" << std::flush;
    //     }
    //     jetTree->GetEntry(iEntry);
    //     // if (jet.d0mass < 1.82054 || jet.d0mass > 1.90946)
    //     //     continue;
    //     // if (jet.centrality <= 1)
    //     //     continue;   // skip 80-100% centrality
    //     // Int_t centBin =  centralityMap[(Int_t)jet.centrality];
    //     Int_t centBin = getCentralityBin(jet.centrality);
    //     if (jet.d0mass == 0)
    //         continue;

    //     hMeasured[centBin]->Fill(jet.ptcorr, jet.d0z);
    //     for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    //     {
    //         hMeasuredAngularity[centBin][iLambda]->Fill(jet.ptcorr, jet.lambda[iLambda]);
    //     }

    // } // end of loop over train entries

    // Read measured histograms from file

    TFile *realDataFile = new TFile("../splot_2D_histogram_samePT.root", "READ");

    TString centralityNamesTemp[nCentralityBins] = {"010", "1040", "4080"};
    TString temp = "_hist_60";
    TString varNamesTemp[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

    TCanvas *can = new TCanvas("can", "", 1200, 1200);
    can->SaveAs("rcp.pdf[");
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    can->cd();
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        hMeasured[iCent] = (TH2D *)realDataFile->Get("z" + temp + centralityNamesTemp[iCent] + "_2D");
        hMeasured[iCent] = invertAxis(hMeasured[iCent]);
        hMeasured[iCent]->Draw("colz");
        tex->DrawLatex(0.2, 0.8, centralityTitles[iCent]);
        can->SaveAs("rcp.pdf");
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularity[iCent][iLambda] = (TH2D *)realDataFile->Get(varNamesTemp[iLambda] + temp + centralityNamesTemp[iCent] + "_2D");
            hMeasuredAngularity[iCent][iLambda] = invertAxis(hMeasuredAngularity[iCent][iLambda]);
            hMeasuredAngularity[iCent][iLambda]->Draw("colz");
            tex->DrawLatex(0.2, 0.8, centralityTitles[iCent] + " " + angularityTitle[iLambda]);
            can->SaveAs("rcp.pdf");
        }
    }
    //  can->SaveAs("rcp.pdf]");
    cout << "reading finished" << endl;

    // Create RooUnfoldBayes object and run the unfolding
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        RooUnfoldBayes unfolding(response[iCent], hMeasured[iCent], 5);
        hUnfolded[iCent] = (TH2D *)unfolding.Hunfold();

        plotComparison(can, hUnfolded[iCent], hMeasured[iCent], hMc[iCent], hMcMeasured[iCent], iCent, "z");

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            RooUnfoldBayes unfoldingAngularity(responseAngularity[iCent][iLambda], hMeasuredAngularity[iCent][iLambda], 5);
            hUnfoldedAngularity[iCent][iLambda] = (TH2D *)unfoldingAngularity.Hunfold();

            plotComparison(can, hUnfoldedAngularity[iCent][iLambda], hMeasuredAngularity[iCent][iLambda], hMcAngularity[iCent][iLambda], hMcMeasuredAngularity[iCent][iLambda], iCent, angularityTitle[iLambda]);
        }
    }

    // Draw closure test check

    TString RcpTitles[3] = {"0-10/40-80", "0-10/10-40", "10-40/40-80"};
    // Create all ratios betweein centrality histograms
    // 0 - 10 / 60 - 80
    TH2D *hRcp[nCentralityBins];
    TH2D *hRcpAngularity[nCentralityBins][nAngularities];

    hRcp[0] = (TH2D *)hUnfolded[0]->Clone("hRcp" + RcpTitles[0]);
    hRcp[0]->Divide(hUnfolded[nCentralityBins - 1]);
    hRcp[0]->Scale(Ncoll[nCentralityBins - 1] / Ncoll[0]);

    can->cd();

    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hRcpAngularity[0][iLambda] = (TH2D *)hUnfoldedAngularity[0][iLambda]->Clone("hRcp" + angularityNames[iLambda] + RcpTitles[0]);
        hRcpAngularity[0][iLambda]->Divide(hUnfoldedAngularity[nCentralityBins - 1][iLambda]);
        hRcpAngularity[0][iLambda]->Scale(Ncoll[nCentralityBins - 1] / Ncoll[0]);
    }

    // for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    // {
    //     hRcp[iCent] = (TH2D *)hUnfolded[iCent]->Clone("hRcp" + centralityNames[iCent]);
    //     hRcp[iCent]->Divide(hUnfolded[iCent - 1]);
    //     hRcp[iCent]->Scale(Ncoll[iCent - 1] / Ncoll[iCent]);

    //     for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    //     {
    //         hRcpAngularity[iCent][iLambda] = (TH2D *)hUnfoldedAngularity[iCent][iLambda]->Clone("hRcp" + centralityNames[iCent] + angularityNames[iLambda]);
    //         hRcpAngularity[iCent][iLambda]->Divide(hUnfoldedAngularity[iCent - 1][iLambda]);
    //         hRcpAngularity[iCent][iLambda]->Scale(Ncoll[iCent - 1] / Ncoll[iCent]);
    //     }
    // }

    // Draw Rcp
    can->cd();
    can->Clear();
    can->Divide(2, 1);
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    const Int_t colors[5] = {2001, 2002, 2003, 2004, kMagenta};

    for (Int_t iCent = 0; iCent < 1; iCent++)
    {
        can->cd(2);
        TH1D *hz = (TH1D *)hRcp[iCent]->ProjectionY("hz" + centralityNames[iCent]);
        hz->GetYaxis()->SetTitle("R_{cp}");
        hz->SetLineColor(colors[iCent]);
        hz->SetMarkerColor(colors[iCent]);
        hz->SetMarkerStyle(20);
        // hz->GetYaxis()->SetRangeUser(0.0, 2.0);
        hz->DrawClone(iCent == 0 ? "" : "same");
        leg->AddEntry(hz, RcpTitles[iCent], "l");

        can->cd(1);
        TH1D *hptz = (TH1D *)hRcp[iCent]->ProjectionX("hptz" + centralityNames[iCent]);
        hptz->GetYaxis()->SetTitle("R_{cp}");
        hptz->SetLineColor(colors[iCent]);
        hptz->SetMarkerColor(colors[iCent]);
        hptz->SetMarkerStyle(20);
        // hz->GetYaxis()->SetRangeUser(0.0, 2.0);
        hptz->DrawClone(iCent == 0 ? "" : "same");
    }
    can->cd();
    leg->Draw("same");
    tex->DrawLatex(0.2, 0.8, "z R_{cp}");
    can->SaveAs("rcp.pdf");

    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {

        can->Clear();
        can->Divide(2, 1);
        can->cd(2);
        for (Int_t iCent = 0; iCent < 1; iCent++)
        {

            TH1D *hLambda = (TH1D *)hRcpAngularity[iCent][iLambda]->ProjectionY("hLambda" + centralityNames[iCent] + angularityNames[iLambda]);
            hLambda->SetLineColor(colors[iCent]);
            hLambda->SetMarkerColor(colors[iCent]);
            hLambda->SetMarkerStyle(20);
            hLambda->GetYaxis()->SetTitle("R_{cp}");
            // hLambda->GetYaxis()->SetRangeUser(0.0, 2.0);
            hLambda->DrawClone(iCent == 0 ? "" : "same");
        }

        can->cd(1);

        for (Int_t iCent = 0; iCent < 1; iCent++)
        {

            TH1D *hptLambda = (TH1D *)hRcpAngularity[iCent][iLambda]->ProjectionX("hptLambda" + centralityNames[iCent] + angularityNames[iLambda]);
            hptLambda->SetLineColor(colors[iCent]);
            hptLambda->SetMarkerColor(colors[iCent]);
            hptLambda->SetMarkerStyle(20);
            hptLambda->GetYaxis()->SetTitle("R_{cp}");
            // hptLambda->GetYaxis()->SetRangeUser(0.0, 2.0);
            hptLambda->DrawClone(iCent == 0 ? "" : "same");
        }

        can->cd();

        leg->Draw("same");
        tex->DrawLatex(0.2, 0.8, angularityTitle[iLambda] + " R_{cp}");
        can->SaveAs("rcp.pdf");
    }
    can->SaveAs("rcp.pdf]");
}

//==============================================================================

// void assignTree(TTree *jetTree, StJetTreeStruct &jet)
// {
//     cout << "Reading tree" << endl;

//     jetTree->SetBranchAddress("centrality", &jet.centrality);
//     jetTree->SetBranchAddress("jet_pt", &jet.pt);
//     jetTree->SetBranchAddress("jet_pt_corr", &jet.ptcorr);
//     jetTree->SetBranchAddress("jet_phi", &jet.phi);
//     jetTree->SetBranchAddress("jet_eta", &jet.eta);
//     jetTree->SetBranchAddress("D0mass", &jet.d0mass);
//     jetTree->SetBranchAddress("D0_r", &jet.d0r);
//     jetTree->SetBranchAddress("D0_pT", &jet.d0pt);
//     jetTree->SetBranchAddress("lambda_1_1", &jet.lambda[0]);
//     jetTree->SetBranchAddress("lambda_1_1half", &jet.lambda[1]);
//     jetTree->SetBranchAddress("lambda_1_2", &jet.lambda[2]);
//     jetTree->SetBranchAddress("lambda_1_3", &jet.lambda[3]);
//     jetTree->SetBranchAddress("z", &jet.d0z);
// }

void assignTree(TTree *jetTree, StJetTreeStruct &jet)
{
    cout << "Reading tree" << endl;

    jetTree->SetBranchAddress("Centrality", &jet.centrality);

    // jetTree->SetBranchAddress("McD0Pt", &jet.d0pt);
    // jetTree->SetBranchAddress("McD0Eta", &jet.d0eta);
    // jetTree->SetBranchAddress("McD0Phi", &jet.d0phi);
    // jetTree->SetBranchAddress("McPionPt", &jet.pionpt);
    // jetTree->SetBranchAddress("McPionEta", &jet.pioneta);
    // jetTree->SetBranchAddress("McPionPhi", &jet.pionphi);
    // jetTree->SetBranchAddress("McKaonPt", &jet.kaonpt);
    // jetTree->SetBranchAddress("McKaonEta", &jet.kaoneta);
    // jetTree->SetBranchAddress("McKaonPhi", &jet.kaonphi);

    // jetTree->SetBranchAddress("McJetPt", &jet.jetpt);
    // jetTree->SetBranchAddress("McJetEta", &jet.jeteta);
    // jetTree->SetBranchAddress("McJetPhi", &jet.jetphi);
    // jetTree->SetBranchAddress("McJetArea", &jet.jetarea);
    // jetTree->SetBranchAddress("McJetE", &jet.jetenergy);
    // jetTree->SetBranchAddress("McJetNConst", &jet.numberofconstituents);
    // jetTree->SetBranchAddress("McJetLambda_1_1", &jet.lambda[0]);
    // jetTree->SetBranchAddress("McJetLambda_1_1half", &jet.lambda[1]);
    // jetTree->SetBranchAddress("McJetLambda_1_2", &jet.lambda[2]);
    // jetTree->SetBranchAddress("McJetLambda_1_3", &jet.lambda[3]);
    // jetTree->SetBranchAddress("McJetD0Z", &jet.d0z);

    jetTree->SetBranchAddress("RecoJetPt", &jet.pt);
    jetTree->SetBranchAddress("RecoJetPt", &jet.ptcorr);
    jetTree->SetBranchAddress("RecoJetPhi", &jet.phi);
    jetTree->SetBranchAddress("RecoJetEta", &jet.eta);
    jetTree->SetBranchAddress("RecoJetNConst", &jet.d0mass);
    // jetTree->SetBranchAddress("D0_r", &jet.d0r);
    jetTree->SetBranchAddress("RecoD0Pt", &jet.d0pt);
    jetTree->SetBranchAddress("RecoJetLambda_1_1", &jet.lambda[0]);
    jetTree->SetBranchAddress("RecoJetLambda_1_1half", &jet.lambda[1]);
    jetTree->SetBranchAddress("RecoJetLambda_1_2", &jet.lambda[2]);
    jetTree->SetBranchAddress("RecoJetLambda_1_3", &jet.lambda[3]);
    jetTree->SetBranchAddress("RecoJetD0Z", &jet.d0z);
    // &jet.d0mass = 1.865;

    // jetTree->SetBranchAddress("McJetPt", &jet.pt);
    // jetTree->SetBranchAddress("McJetPt", &jet.ptcorr);
    // jetTree->SetBranchAddress("McJetPhi", &jet.phi);
    // jetTree->SetBranchAddress("McJetEta", &jet.eta);

    // jetTree->SetBranchAddress("McD0Pt", &jet.d0pt);
    // jetTree->SetBranchAddress("McJetLambda_1_1", &jet.lambda[0]);
    // jetTree->SetBranchAddress("McJetLambda_1_1half", &jet.lambda[1]);
    // jetTree->SetBranchAddress("McJetLambda_1_2", &jet.lambda[2]);
    // jetTree->SetBranchAddress("McJetLambda_1_3", &jet.lambda[3]);
    // jetTree->SetBranchAddress("McJetD0Z", &jet.d0z);
}