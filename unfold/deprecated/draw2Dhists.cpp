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
#include "TPaletteAxis.h"
#include <vector>

#endif

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

const vector<Float_t> ptBinsReco = {-20, -10, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30};

const vector<Float_t> ptBinsMc = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                  12, 14, 16, 18, 20, 30};

// const vector<Float_t> ptBinsMc = {0, 15, 30};

const vector<Float_t> zBinsMc = {0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001};

const vector<Float_t> zBinsReco = {-3., -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 3.};

const vector<vector<Float_t>> angularityBinsMc = {
    {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7},
    {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4},
    {0, 0.025, 0.05, 0.10, 0.15, 0.2},
    {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06}};

const vector<vector<Float_t>> angularityBinsReco{
    {-3, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3},
    {-1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2},
    {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6},
    {-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5}};

const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

const Int_t nCentralityBins = 5;
double centBins[nCentralityBins + 1] = {0, 10, 20, 40, 60, 80}; // in icreasing order
TString centralityTitles[nCentralityBins] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
Int_t getCentralityBin(const Float_t &centrality)
{
    for (Int_t i = 0; i < nCentralityBins; i++)
    {
        if (centrality >= centBins[i] && centrality < centBins[i + 1])
            return (i);
    }
    cout << "Error: centrality not in range" << endl;
    return -1;
}

struct StJetTreeStruct
{
    float refmult;
    float grefmult;
    float centrality;
    float refcorr2;
    float mcrefmult;
    float weight;
    float jetpt;
    float jetcorrectedpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetradius;
    float jetenergy;
    float fRhoValforjet;
    int numberofconstituents;
    float d0z;
    float d0mass;
    float d0pt;
    float d0phi;
    float d0eta;
    float pionpt;
    float pioneta;
    float pionphi;
    float kaonpt;
    float kaoneta;
    float kaonphi;
    float lambda[4];
};

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet, StJetTreeStruct &mcRecoJet);

void draw2Dhists()
{

    gStyle->SetHistFillStyle(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Create histograms

    TH2D *hMeasuredTest[nCentralityBins];
    TH2D *hTruthTest[nCentralityBins];

    TH2D *hMeasuredAngularityTest[nCentralityBins][nAngularities];
    TH2D *hTruthAngularityTest[nCentralityBins][nAngularities];

    const Int_t nBins = 300;

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {

        hMeasuredTest[iCent] = new TH2D("MeasTest" + centralityNames[iCent], "Reco;p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", nBins, -20, 40, nBins, -10, 10);
        hTruthTest[iCent] = new TH2D("TrueTest" + centralityNames[iCent], "Mc;p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", nBins, 0, 30, nBins, 0, 1.001);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularityTest[iCent][iLambda] = new TH2D("MeasTest" + centralityNames[iCent] + angularityNames[iLambda], "Reco;p_{t}, GeV/c;" + angularityTitle[iLambda], nBins, -20, 40, nBins, -15, 15);
            hTruthAngularityTest[iCent][iLambda] = new TH2D("TrueTest" + centralityNames[iCent] + angularityNames[iLambda], "Mc;p_{t}, GeV/c;" + angularityTitle[iLambda], nBins, 0, 30, nBins, 0, 0.7);
        }
    }

    TFile *treeFile; // Open the file containing the tree.
    treeFile = new TFile("../output_jets.root", "READ");
    if (!treeFile || treeFile->IsZombie())
    {
        return;
    }
    TTree *jetTree = (TTree *)treeFile->Get("Jets");
    StJetTreeStruct mcJet, recoJet, mcRecoJet;
    assignTree(jetTree, mcJet, recoJet, mcRecoJet);
    Double_t nEntries = jetTree->GetEntries();
    nEntries /= 1.;
    cout << "nEntries = " << (Float_t)nEntries / 1000. << "k" << endl
         << endl;

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        jetTree->GetEntry(iEntry);
        Int_t centBin = getCentralityBin(recoJet.centrality);
        Bool_t isRecoJetFound = (recoJet.numberofconstituents != 0);
        hTruthTest[centBin]->Fill(mcJet.jetpt, mcJet.d0z);
        if (isRecoJetFound)
            hMeasuredTest[centBin]->Fill(recoJet.jetpt, recoJet.d0z);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hTruthAngularityTest[centBin][iLambda]->Fill(mcJet.jetpt, mcJet.lambda[iLambda]);
            if (isRecoJetFound)
                hMeasuredAngularityTest[centBin][iLambda]->Fill(recoJet.jetpt, recoJet.lambda[iLambda]);
        }

    } // end of loop over test entries

    // Draw closure test check
    TCanvas *can = new TCanvas("can", "Closure Test Check", 1200, 1000);
    TString outPdf = "draw2dhists.pdf";
    can->SaveAs(outPdf + "[");

    can->Divide(2, 2);
    can->SetLogz();
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        can->cd(1);
        gPad->SetLogz();
        hTruthTest[iCent]->Draw("colz");
        can->cd(3);
        gPad->SetLogz();
        hMeasuredTest[iCent]->Draw("colz");
        can->cd(2);
        gPad->SetLogz();
        hTruthAngularityTest[iCent][0]->Draw("colz");
        can->cd(4);
        gPad->SetLogz();
        hMeasuredAngularityTest[iCent][0]->Draw("colz");

        tex->DrawLatex(0.45, 0.45, centralityTitles[iCent]);
        // tex->DrawLatex(0.2, 0.6, "z");
        // tex->DrawLatex(0.6, 0.6, "#lambda_{1}^{1}");

        can->SaveAs(outPdf);
    }
    can->SaveAs(outPdf + "]");
}

//==============================================================================

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet, StJetTreeStruct &mcRecoJet)
{
    cout << "Reading tree" << endl;

    jetTree->SetBranchAddress("Centrality", &recoJet.centrality);
    // jetTree->SetBranchAddress("Weight", &recoJet.weight);
    // jetTree->SetBranchAddress("RefMult", &recoJet.refmult);
    // jetTree->SetBranchAddress("gRefMult", &recoJet.grefmult);
    // jetTree->SetBranchAddress("RefCorr2", &recoJet.refcorr2);
    // jetTree->SetBranchAddress("McRefMult", &mcJet.refmult);
    // jetTree->SetBranchAddress("RecoRefMult", &mcJet.grefmult);

    // jetTree->SetBranchAddress("McD0Pt", &mcJet.d0pt);
    // jetTree->SetBranchAddress("McD0Eta", &mcJet.d0eta);
    // jetTree->SetBranchAddress("McD0Phi", &mcJet.d0phi);
    // jetTree->SetBranchAddress("McPionPt", &mcJet.pionpt);
    // jetTree->SetBranchAddress("McPionEta", &mcJet.pioneta);
    // jetTree->SetBranchAddress("McPionPhi", &mcJet.pionphi);
    // jetTree->SetBranchAddress("McKaonPt", &mcJet.kaonpt);
    // jetTree->SetBranchAddress("McKaonEta", &mcJet.kaoneta);
    // jetTree->SetBranchAddress("McKaonPhi", &mcJet.kaonphi);

    // jetTree->SetBranchAddress("RecoD0Pt", &recoJet.d0pt);
    // jetTree->SetBranchAddress("RecoD0Eta", &recoJet.d0eta);
    // jetTree->SetBranchAddress("RecoD0Phi", &recoJet.d0phi);
    // jetTree->SetBranchAddress("RecoPionPt", &recoJet.pionpt);
    // jetTree->SetBranchAddress("RecoPionEta", &recoJet.pioneta);
    // jetTree->SetBranchAddress("RecoPionPhi", &recoJet.pionphi);
    // jetTree->SetBranchAddress("RecoKaonPt", &recoJet.kaonpt);
    // jetTree->SetBranchAddress("RecoKaonEta", &recoJet.kaoneta);
    // jetTree->SetBranchAddress("RecoKaonPhi", &recoJet.kaonphi);

    jetTree->SetBranchAddress("McJetPt", &mcJet.jetpt);
    // jetTree->SetBranchAddress("McJetEta", &mcJet.jeteta);
    // jetTree->SetBranchAddress("McJetPhi", &mcJet.jetphi);
    // jetTree->SetBranchAddress("McJetArea", &mcJet.jetarea);
    // jetTree->SetBranchAddress("McJetE", &mcJet.jetenergy);
    // jetTree->SetBranchAddress("McJetNConst", &mcJet.numberofconstituents);
    jetTree->SetBranchAddress("McJetLambda_1_1", &mcJet.lambda[0]);
    jetTree->SetBranchAddress("McJetLambda_1_1half", &mcJet.lambda[1]);
    jetTree->SetBranchAddress("McJetLambda_1_2", &mcJet.lambda[2]);
    jetTree->SetBranchAddress("McJetLambda_1_3", &mcJet.lambda[3]);
    jetTree->SetBranchAddress("McJetD0Z", &mcJet.d0z);

    // jetTree->SetBranchAddress("McRecoJetPt", &mcRecoJet.jetpt);
    // jetTree->SetBranchAddress("McRecoJetEta", &mcRecoJet.jeteta);
    // jetTree->SetBranchAddress("McRecoJetPhi", &mcRecoJet.jetphi);
    // jetTree->SetBranchAddress("McRecoJetE", &mcRecoJet.jetenergy);
    // jetTree->SetBranchAddress("McRecoJetArea", &mcRecoJet.jetarea);
    // jetTree->SetBranchAddress("McRecoJetNConst", &mcRecoJet.numberofconstituents);
    // jetTree->SetBranchAddress("McRecoJetLambda_1_1", &mcRecoJet.lambda[0]);
    // jetTree->SetBranchAddress("McRecoJetLambda_1_1half", &mcRecoJet.lambda[1]);
    // jetTree->SetBranchAddress("McRecoJetLambda_1_2", &mcRecoJet.lambda[2]);
    // jetTree->SetBranchAddress("McRecoJetLambda_1_3", &mcRecoJet.lambda[3]);
    // jetTree->SetBranchAddress("McRecoJetD0Z", &mcRecoJet.d0z);

    jetTree->SetBranchAddress("RecoJetPt", &recoJet.jetpt);
    // jetTree->SetBranchAddress("RecoJetEta", &recoJet.jeteta);
    // jetTree->SetBranchAddress("RecoJetPhi", &recoJet.jetphi);
    // jetTree->SetBranchAddress("RecoJetArea", &recoJet.jetarea);
    // jetTree->SetBranchAddress("RecoJetE", &recoJet.jetenergy);
    // jetTree->SetBranchAddress("RecoJetRhoVal", &recoJet.fRhoValforjet);
    jetTree->SetBranchAddress("RecoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress("RecoJetLambda_1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress("RecoJetLambda_1_1half", &recoJet.lambda[1]);
    jetTree->SetBranchAddress("RecoJetLambda_1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress("RecoJetLambda_1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress("RecoJetD0Z", &recoJet.d0z);
}