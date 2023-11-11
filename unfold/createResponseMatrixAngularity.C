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

const Float_t TrainToTestRatio = 0.95;

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

void drawClosure(TH2 *hUnfolded, TH2 *hTruth, TH2 *hMeasured, TCanvas *can)
{
    can->Clear();
    can->Divide(2, 2);
    can->cd(1);
    gPad->SetLogy();

    TH1D *xTruth = (TH1D *)hTruth->ProjectionX();
    TH1D *yTruth = (TH1D *)hTruth->ProjectionY();

    TH1D *xUnfolded = (TH1D *)hUnfolded->ProjectionX();
    TH1D *yUnfolded = (TH1D *)hUnfolded->ProjectionY();

    TH1D *xMeasured = (TH1D *)hMeasured->ProjectionX();
    TH1D *yMeasured = (TH1D *)hMeasured->ProjectionY();

    TH1D *xRatio = (TH1D *)xUnfolded->Clone("xRatio");
    xRatio->Divide(xTruth);
    xRatio->GetYaxis()->SetTitle("Unfolded/True");

    TH1D *yRatio = (TH1D *)yUnfolded->Clone("yRatio");
    yRatio->Divide(yTruth);
    yRatio->GetYaxis()->SetTitle("Unfolded/True");

    xUnfolded->SetLineColor(kAzure + 1);
    xUnfolded->SetMarkerColor(kAzure + 1);
    xUnfolded->SetMarkerStyle(20);
    xMeasured->GetYaxis()->SetRangeUser(1, xUnfolded->GetMaximum() * 1.2);
    xMeasured->Draw("hist");
    xTruth->SetMarkerStyle(20);
    xTruth->SetLineColor(kRed);
    xTruth->Draw("hist same");
    xMeasured->SetMarkerStyle(21);
    xMeasured->SetLineColor(kTeal - 1);
    xUnfolded->Draw("PE same");

    can->cd(3);
    gPad->SetLogy(0);
    xRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    xRatio->SetLineColor(kAzure + 1);
    xRatio->SetMarkerStyle(20);
    xRatio->SetMarkerColor(kAzure + 1);
    xRatio->Draw("PE");

    can->cd(2);
    gPad->SetLogy();
    yUnfolded->SetLineColor(kAzure + 1);
    yUnfolded->SetMarkerColor(kAzure + 1);
    yUnfolded->SetMarkerStyle(20);
    yMeasured->GetYaxis()->SetRangeUser(1, yUnfolded->GetMaximum() * 1.2);
    yMeasured->Draw("hist");
    yTruth->SetMarkerStyle(20);
    yTruth->SetLineColor(kRed);
    yTruth->Draw("hist same");
    yMeasured->SetMarkerStyle(21);
    yMeasured->SetLineColor(kTeal - 1);
    yUnfolded->Draw("PE same");

    can->cd(4);
    gPad->SetLogy(0);
    yRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    yRatio->SetLineColor(kAzure + 1);
    yRatio->SetMarkerStyle(20);
    yRatio->SetMarkerColor(kAzure + 1);
    yRatio->Draw("PE");

    can->cd();
    // Add legend
    TLegend *legend = new TLegend(0.2, 0.3, 0.3, 0.4);
    legend->AddEntry(xTruth, "True", "l");
    legend->AddEntry(xMeasured, "Measured", "l");
    legend->AddEntry(xUnfolded, "Unfolded", "pl");
    legend->Draw();
}

const Int_t colors[6] = {kMagenta, kGreen + 1, kOrange + 9, kAzure, kOrange - 9};

void plotIterations(TCanvas *can, RooUnfoldResponse *response, TH2D *hTruth, TH2D *hMeasured, const Int_t &nIter)
{

    can->Clear();
    can->Divide(2, 2);

    TH1D *xTruth = (TH1D *)hTruth->ProjectionX();
    TH1D *yTruth = (TH1D *)hTruth->ProjectionY();

    TH1D *xMeasured = (TH1D *)hMeasured->ProjectionX();
    TH1D *yMeasured = (TH1D *)hMeasured->ProjectionY();

    TLegend *legend = new TLegend(0.15, 0.32, 0.25, 0.44);
    legend->AddEntry(xTruth, "True", "l");
    legend->AddEntry(xMeasured, "Measured", "l");

    can->cd(1);

    gPad->SetLogy();
    xMeasured->GetYaxis()->SetRangeUser(1, hTruth->GetMaximum() * 1.2);
    xMeasured->Draw("hist");
    xTruth->SetMarkerStyle(20);
    xTruth->SetLineColor(kRed);
    xTruth->Draw("hist same");
    xMeasured->SetMarkerStyle(21);
    xMeasured->SetLineColor(kTeal - 1);

    can->cd(3);
    gPad->SetLogy(0);

    can->cd(2);
    gPad->SetLogy();
    yMeasured->GetYaxis()->SetRangeUser(1, hTruth->GetMaximum() * 1.2);
    yMeasured->Draw("hist");
    yTruth->SetMarkerStyle(20);
    yTruth->SetLineColor(kRed);
    yTruth->Draw("hist same");
    yMeasured->SetMarkerStyle(21);
    yMeasured->SetLineColor(kTeal - 1);

    can->cd(4);
    gPad->SetLogy(0);

    for (Int_t iter = 1; iter < nIter; iter++)
    {
        RooUnfoldBayes unfolding(response, hMeasured, iter);
        TH2D *hUnfolded = (TH2D *)unfolding.Hunfold();
        TH1D *xUnfolded = (TH1D *)hUnfolded->ProjectionX(Form("xUnfolded%i", iter));
        TH1D *yUnfolded = (TH1D *)hUnfolded->ProjectionY(Form("yUnfolded%i", iter));

        xUnfolded->SetLineColor(colors[iter]);
        xUnfolded->SetMarkerStyle(20);
        xUnfolded->SetMarkerColor(colors[iter]);

        yUnfolded->SetLineColor(colors[iter]);
        yUnfolded->SetMarkerStyle(20);
        yUnfolded->SetMarkerColor(colors[iter]);

        legend->AddEntry(xUnfolded, Form("Iter%i", iter), "pel");

        TH1D *xRatio = (TH1D *)xUnfolded->Clone("xRatio");
        xRatio->Divide(xTruth);
        xRatio->GetYaxis()->SetTitle("Unfolded/True");

        TH1D *yRatio = (TH1D *)yUnfolded->Clone("yRatio");
        yRatio->Divide(yTruth);
        yRatio->GetYaxis()->SetTitle("Unfolded/True");

        can->cd(1);
        xUnfolded->Draw("PE same");

        can->cd(3);
        xRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
        xRatio->Draw(iter == 0 ? "PE" : "PE same");

        can->cd(2);
        yUnfolded->Draw("PE same");

        can->cd(4);
        yRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
        yRatio->Draw(iter == 0 ? "PE" : "PE same");
    }
    can->cd();
    legend->Draw();
}

/// get pearson coeffs from covariance matrix
TH2D *getPearsonCoeffs(const TMatrixD &covMatrix)
{

    Int_t nrows = covMatrix.GetNrows();
    Int_t ncols = covMatrix.GetNcols();

    TH2D *PearsonCoeffs = new TH2D("PearsonCoeffs", "Pearson Coefficients;", nrows, 0, nrows, ncols, 0, ncols);
    for (Int_t row = 0; row < nrows; row++)
    {
        for (Int_t col = 0; col < ncols; col++)
        {
            Double_t pearson = 0.;
            if (covMatrix(row, row) != 0. && covMatrix(col, col) != 0.)
                pearson = covMatrix(row, col) / TMath::Sqrt(covMatrix(row, row) * covMatrix(col, col));
            PearsonCoeffs->SetBinContent(row + 1, col + 1, pearson);
        }
    }

    PearsonCoeffs->GetZaxis()->SetRangeUser(-1, 1);
    return PearsonCoeffs;
}

const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

const Int_t nCentralityBins = 6;
double centBins[nCentralityBins] = {0, 10, 20, 40, 60, 80}; // in icreasing order
TString centralityTitles[nCentralityBins] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
TString centralityNames[nCentralityBins] = {"0_80", "0_10", "10_20", "20_40", "40_60", "60_80"};
Int_t getCentralityBin(const Float_t &centrality)
{
    for (Int_t i = 0; i < nCentralityBins - 1; i++)
    {
        if (centrality >= centBins[i] && centrality < centBins[i + 1])
            return (i + 1);
    }
    cout << "Error: centrality not in range" << endl;
    return -1;
}

// const Int_t nPtBinsReco = 10;
// const Float_t ptBinsReco[nPtBinsReco + 1] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};
// const Int_t nPtBinsMc = 19;
// const Float_t ptBinsMc[nPtBinsMc + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
// const Int_t nZBinsMc = 10;
// const Float_t zBinsMc[nZBinsMc + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
// const Int_t nZBinsReco = 9;
// const Float_t zBinsReco[nZBinsReco + 1] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2., 2.5, 3.};

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet, StJetTreeStruct &mcRecoJet);

void createResponseMatrixAngularity()
{
    Bool_t readTree = kFALSE;
    gSystem->Load("libRooUnfold");

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Create histograms
    TFile *responseFile = new TFile("response.root", "READ");
    if (!responseFile || responseFile->IsZombie())
    {
        readTree = kTRUE;
        responseFile = new TFile("response.root", "RECREATE");
    }
    responseFile->cd();

    TH2D *hMeasured[nCentralityBins];
    TH2D *hTruth[nCentralityBins];
    TH2D *hMeasuredTest[nCentralityBins];
    TH2D *hTruthTest[nCentralityBins];

    // Create RooUnfoldResponse object
    RooUnfoldResponse *response[nCentralityBins];

    TH2D *hMeasuredAngularity[nCentralityBins][nAngularities];
    TH2D *hTruthAngularity[nCentralityBins][nAngularities];
    TH2D *hMeasuredAngularityTest[nCentralityBins][nAngularities];
    TH2D *hTruthAngularityTest[nCentralityBins][nAngularities];

    RooUnfoldResponse *responseAngularity[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {

        hMeasured[iCent] = new TH2D("Meas" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsReco.size() - 1, &ptBinsReco[0], zBinsReco.size() - 1, &zBinsReco[0]);
        hTruth[iCent] = new TH2D("True" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsMc.size() - 1, &ptBinsMc[0], zBinsMc.size() - 1, &zBinsMc[0]);

        hMeasuredTest[iCent] = new TH2D("MeasTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsReco.size() - 1, &ptBinsReco[0], zBinsReco.size() - 1, &zBinsReco[0]);
        hTruthTest[iCent] = new TH2D("TrueTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsMc.size() - 1, &ptBinsMc[0], zBinsMc.size() - 1, &zBinsMc[0]);

        response[iCent] = new RooUnfoldResponse("response" + centralityNames[iCent], centralityNames[iCent] + "response");
        response[iCent]->Setup(hMeasured[iCent], hTruth[iCent]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularity[iCent][iLambda] = new TH2D("Meas" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsReco.size() - 1, &ptBinsReco[0], angularityBinsReco[iLambda].size() - 1, &angularityBinsReco[iLambda][0]);
            hTruthAngularity[iCent][iLambda] = new TH2D("True" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsMc.size() - 1, &ptBinsMc[0], angularityBinsMc[iLambda].size() - 1, &angularityBinsMc[iLambda][0]);

            hMeasuredAngularityTest[iCent][iLambda] = new TH2D("MeasTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsReco.size() - 1, &ptBinsReco[0], angularityBinsReco[iLambda].size() - 1, &angularityBinsReco[iLambda][0]);
            hTruthAngularityTest[iCent][iLambda] = new TH2D("TrueTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsMc.size() - 1, &ptBinsMc[0], angularityBinsMc[iLambda].size() - 1, &angularityBinsMc[iLambda][0]);

            responseAngularity[iCent][iLambda] = new RooUnfoldResponse("response" + centralityNames[iCent] + angularityNames[iLambda], "response" + centralityNames[iCent] + angularityNames[iLambda]);
            responseAngularity[iCent][iLambda]->Setup(hMeasuredAngularity[iCent][iLambda], hTruthAngularity[iCent][iLambda]);
        }
    }

    if (readTree)
    {
        TFile *treeFile; // Open the file containing the tree.
        treeFile = new TFile("../output_jets.root", "READ");
        if (!treeFile || treeFile->IsZombie())
        {
            return;
        }
        TTree *jetTree = (TTree *)treeFile->Get("Jets");
        StJetTreeStruct mcJet, recoJet, mcRecoJet;
        assignTree(jetTree, mcJet, recoJet, mcRecoJet);
        Long_t nEntries = jetTree->GetEntries();
        nEntries /= 1.;
        cout << "nEntries = " << (Float_t)nEntries / 1000. << "k" << endl
             << endl;

        cout << "Train Sample = " << TrainToTestRatio * (Float_t)nEntries / 1000. << "k" << endl;
        cout << "Test Sample = " << (1 - TrainToTestRatio) * (Float_t)nEntries / 1000. << "k" << endl;

        for (Int_t iEntry = 0; iEntry < nEntries * TrainToTestRatio; iEntry++)
        {
            Float_t progress = 0.;
            progress = (Float_t)iEntry / (TrainToTestRatio * nEntries);
            if (iEntry % 10000 == 0)
            {
                cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;
            }
            jetTree->GetEntry(iEntry);
            Double_t weight = recoJet.weight;
            Int_t centBin = getCentralityBin(recoJet.centrality);
            Bool_t isRecoJetFound = (recoJet.numberofconstituents != 0);

            if (isRecoJetFound)
            {
                response[centBin]->Fill(recoJet.jetpt, recoJet.d0z, mcJet.jetpt, mcJet.d0z);
                response[0]->Fill(recoJet.jetpt, recoJet.d0z, mcJet.jetpt, mcJet.d0z);

                for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
                {
                    responseAngularity[centBin][iLambda]->Fill(recoJet.jetpt, recoJet.lambda[iLambda], mcJet.jetpt, mcJet.lambda[iLambda]);
                    responseAngularity[0][iLambda]->Fill(recoJet.jetpt, recoJet.lambda[iLambda], mcJet.jetpt, mcJet.lambda[iLambda]);
                }
            }

            else
            {
                response[centBin]->Miss(mcJet.jetpt, mcJet.d0z);
                response[0]->Miss(mcJet.jetpt, mcJet.d0z);

                for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
                {
                    responseAngularity[centBin][iLambda]->Miss(mcJet.jetpt, mcJet.lambda[iLambda]);
                    responseAngularity[0][iLambda]->Miss(mcJet.jetpt, mcJet.lambda[iLambda]);
                }
            }

        } // end of loop over train entries
        for (Int_t iEntry = nEntries * TrainToTestRatio; iEntry < nEntries * (1. - TrainToTestRatio); iEntry++)
        {
            jetTree->GetEntry(iEntry);
            Double_t weight = recoJet.weight;
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

        for (Int_t iCent = 1; iCent < nCentralityBins; iCent++) // make a unified centrality bin 0-80%
        {
            hMeasuredTest[0]->Add(hMeasuredTest[iCent]);
            hTruthTest[0]->Add(hTruthTest[iCent]);

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                hMeasuredAngularityTest[0][iLambda]->Add(hMeasuredAngularityTest[iCent][iLambda]);
                hTruthAngularityTest[0][iLambda]->Add(hTruthAngularityTest[iCent][iLambda]);
            }
        }
    } // if readTree

    else // if not readTree
    {
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            response[iCent] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent]);
            hMeasuredTest[iCent] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/MeasTest" + centralityNames[iCent]);
            hTruthTest[iCent] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/TrueTest" + centralityNames[iCent]);

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                responseAngularity[iCent][iLambda] = (RooUnfoldResponse *)responseFile->Get(centralityNames[iCent] + "/response" + centralityNames[iCent] + angularityNames[iLambda]);

                hMeasuredAngularityTest[iCent][iLambda] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/MeasTest" + centralityNames[iCent] + angularityNames[iLambda]);
                hTruthAngularityTest[iCent][iLambda] = (TH2D *)responseFile->Get(centralityNames[iCent] + "/TrueTest" + centralityNames[iCent] + angularityNames[iLambda]);
            }
        }
    }

    // Draw closure test check
    TCanvas *can = new TCanvas("can", "Closure Test Check", 1200, 1200);

    can->SaveAs("closure_check.pdf[");
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    // Create RooUnfoldBayes object and run the unfolding
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << "centrality  " << centralityTitles[iCent] << endl;
        RooUnfoldBayes unfolding(response[iCent], hMeasuredTest[iCent], 5);

        TH2D *hUnfolded = (TH2D *)unfolding.Hunfold();
        hUnfolded->SetName("Unfolded" + centralityNames[iCent]);

        drawClosure(hUnfolded, hTruthTest[iCent], hMeasuredTest[iCent], can);
        can->cd();
        tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + " z");
        can->SaveAs("closure_check.pdf");

        plotIterations(can, response[iCent], hTruthTest[iCent], hMeasuredTest[iCent], 5);
        can->cd();
        tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + " z");
        can->SaveAs("closure_check.pdf");
        can->SetLogz(1);

        // TH2 *hResponse = response[iCent]->Hresponse();
        // hResponse->SetTitle("Hresponse");
        // hResponse->Draw("colz");
        // tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  z");
        // tex->DrawLatex(0.5, 0.35, "Hresponse");
        // can->SaveAs("closure_check.pdf");

        // TMatrixD MresponseMatrix = response[iCent]->Mresponse(true);

        // TH2D *hMresponse = new TH2D(MresponseMatrix);
        // hMresponse->SetTitle("Mresponse");
        // hMresponse->Draw("colz");
        // tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  z");
        // tex->DrawLatex(0.5, 0.35, "Mresponse");
        // can->SaveAs("closure_check.pdf");
        // delete hMresponse;

        TMatrixD covMatrix = unfolding.UnfoldingMatrix();

        TH2D hUnfoldingMatrix = (TH2D)unfolding.UnfoldingMatrix();
        hUnfoldingMatrix.SetTitle((TString) "Unfolding Matrix;bin = z_{det}*" + ptBinsReco.size() + "+pt_{det};bin =z_{reco}*" + ptBinsMc.size() + "+pt_{reco}");
        hUnfoldingMatrix.Draw("colz");
        tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  z");
        tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
        can->SaveAs("closure_check.pdf");

        TH2D *fPearsonCoeffs = getPearsonCoeffs(unfolding.Eunfold(RooUnfold::kCovariance));
        fPearsonCoeffs->SetTitle((TString) "Pearson Coefficients;bin = z*" + ptBinsMc.size() + "+pt;bin =z*" + ptBinsMc.size() + "+pt");
        can->cd();
        can->SetLogz(0);
        fPearsonCoeffs->Draw("colz");
        tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  z");
        tex->DrawLatex(0.5, 0.35, "Pearson Coefficients");
        can->SaveAs("closure_check.pdf");
        delete fPearsonCoeffs;
        can->SetLogz();

        TH2D hCovariance = (TH2D)unfolding.Eunfold(RooUnfold::kCovariance);
        hCovariance.SetTitle((TString) "Covariance Matrix;bin = z*" + ptBinsMc.size() + "+pt;bin =z*" + ptBinsMc.size() + "+pt");
        hCovariance.Draw("colz");
        tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  z");
        tex->DrawLatex(0.5, 0.35, "Covariance Matrix");
        can->SaveAs("closure_check.pdf");

        TH2D *hUnfoldedAngularity[nAngularities];

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            RooUnfoldBayes unfoldingAngularity(responseAngularity[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], 5);
            hUnfoldedAngularity[iLambda] = (TH2D *)unfoldingAngularity.Hunfold();
            hUnfoldedAngularity[iLambda]->SetName("Unfolded" + centralityNames[iCent] + angularityNames[iLambda]);

            drawClosure(hUnfoldedAngularity[iLambda], hTruthAngularityTest[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], can);
            can->cd();
            tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            can->SaveAs("closure_check.pdf");

            plotIterations(can, responseAngularity[iCent][iLambda], hTruthAngularityTest[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], 5);
            can->cd();
            tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            can->SaveAs("closure_check.pdf");

            TH2D *fPearsonCoeffs = getPearsonCoeffs(unfoldingAngularity.Eunfold(RooUnfold::kCovariance));
            fPearsonCoeffs->SetTitle((TString) "Pearson Coefficients;bin = #lambda*" + ptBinsMc.size() + "+pt;bin =#lambda*" + ptBinsMc.size() + "+pt");
            can->cd();
            can->SetLogz(0);
            fPearsonCoeffs->Draw("colz");
            tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            tex->DrawLatex(0.5, 0.35, "Pearson Coefficients");
            can->SaveAs("closure_check.pdf");
            delete fPearsonCoeffs;
            can->SetLogz();

            TH2D hCovarianceAngularity = (TH2D)unfoldingAngularity.Eunfold(RooUnfold::kCovariance);
            hCovarianceAngularity.SetTitle((TString) "Covariance Matrix;bin = #lambda*" + ptBinsMc.size() + "+pt;bin =#lambda*" + ptBinsMc.size() + "+pt");
            hCovarianceAngularity.Draw("colz");
            tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            tex->DrawLatex(0.5, 0.35, "Covariance Matrix");
            can->SaveAs("closure_check.pdf");

            TH2D hUnfoldingMatrixAngularity = (TH2D)unfoldingAngularity.UnfoldingMatrix();
            hUnfoldingMatrixAngularity.SetTitle((TString) "Unfolding Matrix;bin = #lambda_{det}*" + ptBinsReco.size() + "+pt_{det};bin =#lambda_{reco}*" + ptBinsMc.size() + "+pt_{reco}");
            hUnfoldingMatrixAngularity.Draw("colz");
            tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
            can->SaveAs("closure_check.pdf");
        }

        if (readTree)
        {
            responseFile->cd();
            TDirectory *dir = responseFile->mkdir(centralityNames[iCent]);
            dir->cd();

            hMeasuredTest[iCent]->Write();
            hTruthTest[iCent]->Write();

            hUnfolded->Write();
            response[iCent]->Write();

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {

                hMeasuredAngularityTest[iCent][iLambda]->Write();
                hTruthAngularityTest[iCent][iLambda]->Write();

                hUnfoldedAngularity[iLambda]->Write();
                responseAngularity[iCent][iLambda]->Write();
            }
        }
    } // end of loop over centrality bins
    can->SaveAs("closure_check.pdf]");

    responseFile->Save();
}

//==============================================================================

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet, StJetTreeStruct &mcRecoJet)
{
    cout << "Reading tree" << endl;

    jetTree->SetBranchAddress("Centrality", &recoJet.centrality);
    jetTree->SetBranchAddress("Weight", &recoJet.weight);
    jetTree->SetBranchAddress("RefMult", &recoJet.refmult);
    jetTree->SetBranchAddress("gRefMult", &recoJet.grefmult);
    jetTree->SetBranchAddress("RefCorr2", &recoJet.refcorr2);
    jetTree->SetBranchAddress("McRefMult", &mcJet.refmult);
    jetTree->SetBranchAddress("RecoRefMult", &mcJet.grefmult);

    jetTree->SetBranchAddress("McD0Pt", &mcJet.d0pt);
    jetTree->SetBranchAddress("McD0Eta", &mcJet.d0eta);
    jetTree->SetBranchAddress("McD0Phi", &mcJet.d0phi);
    jetTree->SetBranchAddress("McPionPt", &mcJet.pionpt);
    jetTree->SetBranchAddress("McPionEta", &mcJet.pioneta);
    jetTree->SetBranchAddress("McPionPhi", &mcJet.pionphi);
    jetTree->SetBranchAddress("McKaonPt", &mcJet.kaonpt);
    jetTree->SetBranchAddress("McKaonEta", &mcJet.kaoneta);
    jetTree->SetBranchAddress("McKaonPhi", &mcJet.kaonphi);

    jetTree->SetBranchAddress("RecoD0Pt", &recoJet.d0pt);
    jetTree->SetBranchAddress("RecoD0Eta", &recoJet.d0eta);
    jetTree->SetBranchAddress("RecoD0Phi", &recoJet.d0phi);
    jetTree->SetBranchAddress("RecoPionPt", &recoJet.pionpt);
    jetTree->SetBranchAddress("RecoPionEta", &recoJet.pioneta);
    jetTree->SetBranchAddress("RecoPionPhi", &recoJet.pionphi);
    jetTree->SetBranchAddress("RecoKaonPt", &recoJet.kaonpt);
    jetTree->SetBranchAddress("RecoKaonEta", &recoJet.kaoneta);
    jetTree->SetBranchAddress("RecoKaonPhi", &recoJet.kaonphi);

    jetTree->SetBranchAddress("McJetPt", &mcJet.jetpt);
    jetTree->SetBranchAddress("McJetEta", &mcJet.jeteta);
    jetTree->SetBranchAddress("McJetPhi", &mcJet.jetphi);
    jetTree->SetBranchAddress("McJetArea", &mcJet.jetarea);
    jetTree->SetBranchAddress("McJetE", &mcJet.jetenergy);
    jetTree->SetBranchAddress("McJetNConst", &mcJet.numberofconstituents);
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
    jetTree->SetBranchAddress("RecoJetEta", &recoJet.jeteta);
    jetTree->SetBranchAddress("RecoJetPhi", &recoJet.jetphi);
    jetTree->SetBranchAddress("RecoJetArea", &recoJet.jetarea);
    jetTree->SetBranchAddress("RecoJetE", &recoJet.jetenergy);
    jetTree->SetBranchAddress("RecoJetRhoVal", &recoJet.fRhoValforjet);
    jetTree->SetBranchAddress("RecoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress("RecoJetLambda_1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress("RecoJetLambda_1_1half", &recoJet.lambda[1]);
    jetTree->SetBranchAddress("RecoJetLambda_1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress("RecoJetLambda_1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress("RecoJetD0Z", &recoJet.d0z);
}