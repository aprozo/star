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

const vector<Float_t> ptBinsDetector = {-20, -10, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30};

const vector<Float_t> ptBinsParticle = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                        12, 14, 16, 18, 20, 30};

// const vector<Float_t> ptBinsParticle = {0, 15, 30};

const vector<Float_t> zBinParticle = {0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001};

const vector<Float_t> zBinsDetector = {-3., -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 3.};

const vector<vector<Float_t>> angularityBinsParticle = {
    {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7},
    {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4},
    {0, 0.025, 0.05, 0.10, 0.15, 0.2},
    {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06}};

const vector<vector<Float_t>> angularityBinsDetector{
    {-3, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3},
    {-1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2},
    {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6},
    {-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5}};

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

const Int_t nCentralityBins = 5;

map<Int_t, Int_t> centralityMap = {
    {8, 0}, // 0-10%
    {7, 1}, // 10-20%
    {6, 2}, // 20-40%
    {5, 2}, // 20-40%
    {4, 3}, // 40-60%
    {3, 3}, // 40-60%
    {2, 4}, // 60-80%
    {1, 4}  // 60-80%
};

double centBins[nCentralityBins + 1] = {0, 10, 20, 40, 60, 80}; // in icreasing order
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

const Double_t Ncoll[nCentralityBins] = {952., 599., 297., 93., 21.};

TString centralityTitles[nCentralityBins] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_20", "20_40", "40_60", "60_80"};

void assignTree(TTree *jetTree, StJetTreeStruct &jet);

void calculateRcp()
{
    Bool_t readTree = kFALSE;
    gSystem->Load("libRooUnfold");

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Read response matrix from file
    TFile *responseFile = new TFile("../unfold/response.root", "READ");
    responseFile->cd();

    RooUnfoldResponse *response[nCentralityBins];
    RooUnfoldResponse *responseAngularity[nCentralityBins][nAngularities];

    TH2D *hMeasured[nCentralityBins];
    TH2D *hUnfolded[nCentralityBins];

    TH2D *hMeasuredAngularity[nCentralityBins][nAngularities];
    TH2D *hUnfoldedAngularity[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        response[iCent] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent]);
        hMeasured[iCent] = new TH2D("Meas" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsDetector.size() - 1, &ptBinsDetector[0], zBinsDetector.size() - 1, &zBinsDetector[0]);
        hUnfolded[iCent] = new TH2D("Unfolded" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsParticle.size() - 1, &ptBinsParticle[0], zBinParticle.size() - 1, &zBinParticle[0]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            responseAngularity[iCent][iLambda] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent] + angularityNames[iLambda]);

            hMeasuredAngularity[iCent][iLambda] = new TH2D("Meas" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsDetector.size() - 1, &ptBinsDetector[0], angularityBinsDetector[iLambda].size() - 1, &angularityBinsDetector[iLambda][0]);
            hUnfoldedAngularity[iCent][iLambda] = new TH2D("Unfolded" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsParticle.size() - 1, &ptBinsParticle[0], angularityBinsParticle[iLambda].size() - 1, &angularityBinsParticle[iLambda][0]);
        }
    }

    TFile *treeFile; // Open the file containing the tree.
    // treeFile = new TFile("../D0_jets_2014_231030.root", "READ");
    treeFile = new TFile("../output_jets.root", "READ");
    if (!treeFile || treeFile->IsZombie())
    {
        return;
    }
    TTree *jetTree = (TTree *)treeFile->Get("Jets");
    StJetTreeStruct jet;
    assignTree(jetTree, jet);
    Long_t nEntries = jetTree->GetEntries() / 1.;

    cout << "nEntries = " << (Float_t)nEntries / 1000. << "k" << endl
         << endl;

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        Float_t progress = 0.;
        progress = (Float_t)iEntry / (Float_t)nEntries;
        if (iEntry % 1000 == 0)
        {
            cout << "\r (" << (progress * 100.0) << "%)" << std::flush;
        }
        jetTree->GetEntry(iEntry);
        // if (jet.d0mass < 1.82054 || jet.d0mass > 1.90946)
        //     continue;
        // if (jet.centrality <= 1)
        //     continue;   // skip 80-100% centrality
        // Int_t centBin =  centralityMap[(Int_t)jet.centrality];
        Int_t centBin = getCentralityBin(jet.centrality);
        if (jet.d0mass == 0)
            continue;

        hMeasured[centBin]->Fill(jet.ptcorr, jet.d0z);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularity[centBin][iLambda]->Fill(jet.ptcorr, jet.lambda[iLambda]);
        }

    } // end of loop over train entries

    cout << "reading finished" << endl;

    // Create RooUnfoldBayes object and run the unfolding
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        RooUnfoldBayes unfolding(response[iCent], hMeasured[iCent], 5);
        hUnfolded[iCent] = (TH2D *)unfolding.Hunfold();

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            RooUnfoldBayes unfoldingAngularity(responseAngularity[iCent][iLambda], hMeasuredAngularity[iCent][iLambda], 5);
            hUnfoldedAngularity[iCent][iLambda] = (TH2D *)unfoldingAngularity.Hunfold();
        }
    }

    // Draw closure test check
    TCanvas *can = new TCanvas("can", "", 1200, 1200);
    can->SaveAs("rcp.pdf[");
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    TString RcpTitles[5] = {"0-10/60-80", "0-10/10-20", "10-20/20-40", "20-40/40-60", "40-60/60-80"};
    // Create all ratios betweein centrality histograms
    // 0 - 10 / 60 - 80
    TH2D *hRcp[nCentralityBins];
    TH2D *hRcpAngularity[nCentralityBins][nAngularities];

    hRcp[0] = (TH2D *)hUnfolded[0]->Clone("hRcp" + RcpTitles[0]);
    hRcp[0]->Divide(hUnfolded[nCentralityBins - 1]);
    hRcp[0]->Scale(Ncoll[nCentralityBins - 1] / Ncoll[0]);

    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hRcpAngularity[0][iLambda] = (TH2D *)hUnfoldedAngularity[0][iLambda]->Clone("hRcp" + angularityNames[iLambda] + RcpTitles[0]);
        hRcpAngularity[0][iLambda]->Divide(hUnfoldedAngularity[nCentralityBins - 1][iLambda]);
        hRcpAngularity[0][iLambda]->Scale(Ncoll[nCentralityBins - 1] / Ncoll[0]);
    }

    for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    {
        hRcp[iCent] = (TH2D *)hUnfolded[iCent]->Clone("hRcp" + centralityNames[iCent]);
        hRcp[iCent]->Divide(hUnfolded[iCent - 1]);
        hRcp[iCent]->Scale(Ncoll[iCent - 1] / Ncoll[iCent]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hRcpAngularity[iCent][iLambda] = (TH2D *)hUnfoldedAngularity[iCent][iLambda]->Clone("hRcp" + centralityNames[iCent] + angularityNames[iLambda]);
            hRcpAngularity[iCent][iLambda]->Divide(hUnfoldedAngularity[iCent - 1][iLambda]);
            hRcpAngularity[iCent][iLambda]->Scale(Ncoll[iCent - 1] / Ncoll[iCent]);
        }
    }

    // Draw Rcp
    can->cd();
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    const Int_t colors[5] = {kBlack, kRed, kBlue, kGreen, kMagenta};

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        TH1D *hz = (TH1D *)hRcp[iCent]->ProjectionY("hz" + centralityNames[iCent]);
        hz->GetYaxis()->SetTitle("R_{cp}");
        hz->SetLineColor(colors[iCent]);
        hz->SetMarkerColor(colors[iCent]);
        hz->SetMarkerStyle(20);
        hz->GetYaxis()->SetRangeUser(0.0, 2.0);
        hz->DrawClone(iCent == 0 ? "" : "same");
        leg->AddEntry(hz, RcpTitles[iCent], "l");
    }
    leg->Draw("same");
    tex->DrawLatex(0.2, 0.8, "z R_{cp}");
    can->SaveAs("rcp.pdf");

    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            TH1D *hLambda = (TH1D *)hRcpAngularity[iCent][iLambda]->ProjectionY("hLambda" + centralityNames[iCent] + angularityNames[iLambda]);
            hLambda->SetLineColor(colors[iCent]);
            hLambda->SetMarkerColor(colors[iCent]);
            hLambda->SetMarkerStyle(20);
            hLambda->GetYaxis()->SetTitle("R_{cp}");
            hLambda->GetYaxis()->SetRangeUser(0.0, 2.0);
            hLambda->DrawClone(iCent == 0 ? "" : "same");
        }
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