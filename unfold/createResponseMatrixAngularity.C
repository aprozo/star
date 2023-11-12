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

vector<Float_t> ptMcBinsVec[nCentralityBins];
vector<Float_t> zMcBinsVec[nCentralityBins];
vector<Float_t> angularityMcBinsVec[nCentralityBins][nAngularities];

vector<Float_t> ptRecoBinsVec[nCentralityBins];
vector<Float_t> zRecoBinsVec[nCentralityBins];
vector<Float_t> angularityRecoBinsVec[nCentralityBins][nAngularities];

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

void NormalizeByBinWidth(TH1D *hist)
{
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
        hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
    }
}

// void drawClosure(TH2 *hUnfolded, TH2 *hTruth, TH2 *hMeasured, TCanvas *can)
// {
//     can->Clear();
//     can->Divide(2, 2);
//     can->cd(1);
//     gPad->SetLogy();

//     TH1D *xTruth = (TH1D *)hTruth->ProjectionX();
//     TH1D *yTruth = (TH1D *)hTruth->ProjectionY();
//     NormalizeByBinWidth(xTruth);
//     NormalizeByBinWidth(yTruth);

//     TH1D *xUnfolded = (TH1D *)hUnfolded->ProjectionX();
//     TH1D *yUnfolded = (TH1D *)hUnfolded->ProjectionY();
//     NormalizeByBinWidth(xUnfolded);
//     NormalizeByBinWidth(yUnfolded);

//     TH1D *xMeasured = (TH1D *)hMeasured->ProjectionX();
//     TH1D *yMeasured = (TH1D *)hMeasured->ProjectionY();
//     NormalizeByBinWidth(xMeasured);
//     NormalizeByBinWidth(yMeasured);

//     TH1D *xRatio = (TH1D *)xUnfolded->Clone("xRatio");
//     xRatio->Divide(xTruth);
//     xRatio->GetYaxis()->SetTitle("Unfolded/True");

//     TH1D *yRatio = (TH1D *)yUnfolded->Clone("yRatio");
//     yRatio->Divide(yTruth);
//     yRatio->GetYaxis()->SetTitle("Unfolded/True");

//     xUnfolded->SetLineColor(kAzure + 1);
//     xUnfolded->SetMarkerColor(kAzure + 1);
//     xUnfolded->SetMarkerStyle(20);
//     xMeasured->GetYaxis()->SetRangeUser(1, xUnfolded->GetMaximum() * 1.2);
//     xMeasured->Draw("hist");
//     xTruth->SetMarkerStyle(20);
//     xTruth->SetLineColor(kRed);
//     xTruth->Draw("hist same");
//     xMeasured->SetMarkerStyle(21);
//     xMeasured->SetLineColor(kTeal - 1);
//     xUnfolded->Draw("PE same");

//     can->cd(3);
//     gPad->SetLogy(0);
//     xRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
//     xRatio->SetLineColor(kAzure + 1);
//     xRatio->SetMarkerStyle(20);
//     xRatio->SetMarkerColor(kAzure + 1);
//     xRatio->Draw("PE");

//     can->cd(2);
//     gPad->SetLogy();
//     yUnfolded->SetLineColor(kAzure + 1);
//     yUnfolded->SetMarkerColor(kAzure + 1);
//     yUnfolded->SetMarkerStyle(20);
//     yMeasured->GetYaxis()->SetRangeUser(yMeasured->GetMinimum()*0.8, yUnfolded->GetMaximum() * 1.2);
//     yMeasured->GetXaxis()->SetRangeUser(-2, 2);
//     yMeasured->Draw("hist");
//     yTruth->SetMarkerStyle(20);
//     yTruth->SetLineColor(kRed);
//     yTruth->Draw("hist same");
//     yMeasured->SetMarkerStyle(21);
//     yMeasured->SetLineColor(kTeal - 1);
//     yUnfolded->Draw("PE same");

//     can->cd(4);
//     gPad->SetLogy(0);
//     yRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
//     yRatio->SetLineColor(kAzure + 1);
//     yRatio->SetMarkerStyle(20);
//     yRatio->SetMarkerColor(kAzure + 1);
//     yRatio->Draw("PE");

//     can->cd();
//     // Add legend
//     TLegend *legend = new TLegend(0.2, 0.3, 0.3, 0.4);
//     legend->AddEntry(xTruth, "True", "l");
//     legend->AddEntry(xMeasured, "Measured", "l");
//     legend->AddEntry(xUnfolded, "Unfolded", "pl");
//     legend->Draw();
// }

/// get pearson coeffs from covariance matrix
TH2D *getPearsonCoeffs(const TMatrixD &covMatrix)
{

    Int_t nrows = covMatrix.GetNrows();
    Int_t ncols = covMatrix.GetNcols();

    TH2D *PearsonCoeffs = new TH2D((TString) "PearsonCoeffs" + covMatrix.GetName(), "Pearson Coefficients;", nrows, 0, nrows, ncols, 0, ncols);
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

void plotIterations(TCanvas *can, TString outPdf, RooUnfoldResponse *response, TH2D *hTruth, TH2D *hMeasured, const Int_t &iCent, TString var)
{
    const vector<Int_t> plotIterations = {1, 2, 3, 4, 10, 15, 20};
    const Int_t nIter = plotIterations.size();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    can->Clear();
    can->Divide(2, 2);

    TH1D *xTruth = (TH1D *)hTruth->ProjectionX();
    TH1D *yTruth = (TH1D *)hTruth->ProjectionY();
    NormalizeByBinWidth(xTruth);
    NormalizeByBinWidth(yTruth);

    TH1D *xMeasured = (TH1D *)hMeasured->ProjectionX();
    TH1D *yMeasured = (TH1D *)hMeasured->ProjectionY();
    NormalizeByBinWidth(xMeasured);
    NormalizeByBinWidth(yMeasured);

    TLegend *legend = new TLegend(0.15, 0.32, 0.25, 0.44);
    legend->AddEntry(xTruth, "True", "l");
    legend->AddEntry(xMeasured, "Measured", "l");

    can->cd(1);

    gPad->SetLogy();
    xMeasured->GetYaxis()->SetRangeUser(xMeasured->GetMinimum() * 0.8, xTruth->GetMaximum() * 1.5);
    xMeasured->GetYaxis()->SetTitle("dN/dp_{t}");
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
    yMeasured->GetYaxis()->SetRangeUser(yMeasured->GetMinimum() * 0.8, yTruth->GetMaximum() * 1.5);
  
    yMeasured->GetYaxis()->SetTitle("dN/d" + var);
    yMeasured->GetXaxis()->SetRangeUser(-1.5,1.5);
    yMeasured->Draw("hist");
    yTruth->SetMarkerStyle(20);
    yTruth->SetLineColor(kRed);
    yTruth->Draw("hist same");
    yMeasured->SetMarkerStyle(21);
    yMeasured->SetLineColor(kTeal - 1);

    can->cd(4);
    gPad->SetLogy(0);

    TH2D *fPearsonCoeffs[nIter];
    for (int iter = 0; iter < nIter; iter++)
    {
        RooUnfoldBayes unfolding(response, hMeasured, plotIterations[iter]);

        TH2D *hUnfolded = (TH2D *)unfolding.Hunfold();
        TH1D *xUnfolded = (TH1D *)hUnfolded->ProjectionX(Form("xUnfolded%i", iter));
        TH1D *yUnfolded = (TH1D *)hUnfolded->ProjectionY(Form("yUnfolded%i", iter));
        fPearsonCoeffs[iter] = getPearsonCoeffs(unfolding.Eunfold(RooUnfold::kCovariance));

        NormalizeByBinWidth(xUnfolded);
        NormalizeByBinWidth(yUnfolded);

        xUnfolded->SetLineColor(2000 + iter);
        xUnfolded->SetMarkerStyle(20);
        xUnfolded->SetMarkerColor(2000 + iter);

        yUnfolded->SetLineColor(2000 + iter);
        yUnfolded->SetMarkerStyle(20);
        yUnfolded->SetMarkerColor(2000 + iter);

        legend->AddEntry(xUnfolded, Form("Iter%i", plotIterations[iter]), "pel");

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
    tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + "  " + var);
    can->SaveAs(outPdf);
    can->Clear();
    can->SetLogz(1);
    RooUnfoldBayes lastUnfolding(response, hMeasured, 4);
    TH2D *hUnfolded = (TH2D *)lastUnfolding.Hunfold();

    TH2D *hUnfoldingMatrix = new TH2D(lastUnfolding.UnfoldingMatrix());
    hUnfoldingMatrix->SetTitle("Unfolding Matrix");
    hUnfoldingMatrix->SetName("UnfoldingMatrix");
    can->cd();
    hUnfoldingMatrix->Draw("colz");
    tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + var);
    tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
    can->SaveAs(outPdf);
    delete hUnfoldingMatrix;

    can->Clear();
    can->cd();
    can->SetLogz(0);
    for (int iter = 0; iter < nIter; iter++)
    {
        fPearsonCoeffs[iter]->SetTitle((TString) "Pearson Coefficients;bin =" + var + "*10+pt;bin =" + var + "*10+pt");
        fPearsonCoeffs[iter]->Draw("colz");
        tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + var);
        tex->DrawLatex(0.2, 0.35, Form("Pearson Coefficients Iter %i", plotIterations[iter]));

        can->SaveAs(outPdf);
        delete fPearsonCoeffs[iter];
    }

    can->Clear();

    // TH2D hCovariance = (TH2D)lastUnfolding->Eunfold(RooUnfold::kCovariance);
    // hCovariance.SetTitle((TString) "Covariance Matrix;bin =" + var + "*" + ptBinsMc.size() + "+pt;bin =" + var + "*" + ptBinsMc.size() + "+pt");
    // hCovariance.Draw("colz");
    // tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + var);
    // tex->DrawLatex(0.5, 0.35, "Covariance Matrix");
    // can->SaveAs(outPdf);
    // delete lastUnfolding;
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
    gStyle->SetHistFillStyle(0);
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

    // for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    // {

    //     hMeasured[iCent] = new TH2D("Meas" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsReco.size() - 1, &ptBinsReco[0], zBinsReco.size() - 1, &zBinsReco[0]);
    //     hTruth[iCent] = new TH2D("True" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsMc.size() - 1, &ptBinsMc[0], zBinsMc.size() - 1, &zBinsMc[0]);

    //     hMeasuredTest[iCent] = new TH2D("MeasTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsReco.size() - 1, &ptBinsReco[0], zBinsReco.size() - 1, &zBinsReco[0]);
    //     hTruthTest[iCent] = new TH2D("TrueTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptBinsMc.size() - 1, &ptBinsMc[0], zBinsMc.size() - 1, &zBinsMc[0]);

    //     response[iCent] = new RooUnfoldResponse("response" + centralityNames[iCent], centralityNames[iCent] + "response");
    //     response[iCent]->Setup(hMeasured[iCent], hTruth[iCent]);

    //     for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    //     {
    //         hMeasuredAngularity[iCent][iLambda] = new TH2D("Meas" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsReco.size() - 1, &ptBinsReco[0], angularityBinsReco[iLambda].size() - 1, &angularityBinsReco[iLambda][0]);
    //         hTruthAngularity[iCent][iLambda] = new TH2D("True" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsMc.size() - 1, &ptBinsMc[0], angularityBinsMc[iLambda].size() - 1, &angularityBinsMc[iLambda][0]);

    //         hMeasuredAngularityTest[iCent][iLambda] = new TH2D("MeasTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsReco.size() - 1, &ptBinsReco[0], angularityBinsReco[iLambda].size() - 1, &angularityBinsReco[iLambda][0]);
    //         hTruthAngularityTest[iCent][iLambda] = new TH2D("TrueTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptBinsMc.size() - 1, &ptBinsMc[0], angularityBinsMc[iLambda].size() - 1, &angularityBinsMc[iLambda][0]);

    //         responseAngularity[iCent][iLambda] = new RooUnfoldResponse("response" + centralityNames[iCent] + angularityNames[iLambda], "response" + centralityNames[iCent] + angularityNames[iLambda]);
    //         responseAngularity[iCent][iLambda]->Setup(hMeasuredAngularity[iCent][iLambda], hTruthAngularity[iCent][iLambda]);
    //     }
    // }

    // vector <Float_t> ptMcBinsVec[nCentralityBins];
    // vector <Float_t> zMcBinsVec[nCentralityBins];
    // vector <vector<Float_t>> angularityMcBinsVec[nCentralityBins][nAngularities];

    // vector <Float_t> ptRecoBinsVec[nCentralityBins];
    // vector <Float_t> zRecoBinsVec[nCentralityBins];
    // vector <vector<Float_t>> angularityRecoBinsVec[nCentralityBins][nAngularities];

    // centrality0-10%
    // centrality0-10%
    ptMcBinsVec[0] = {1, 1.667, 2.102, 2.479, 2.856, 3.262, 3.697, 4.219, 4.915, 6.133, 30};
    zMcBinsVec[0] = {0, 0.557, 0.658, 0.733, 0.798, 0.854, 0.908, 1.0001};
    angularityMcBinsVec[0][0] = {0, 0.0007, 0.0616, 0.1071, 0.154, 0.2051, 0.2674, 0.3549, 0.7};
    angularityMcBinsVec[0][1] = {0, 0.0004, 0.0128, 0.03, 0.0516, 0.076, 0.106, 0.1512, 0.4};
    angularityMcBinsVec[0][2] = {0, 0.0002, 0.0032, 0.01, 0.0204, 0.0324, 0.0472, 0.0702, 0.2};
    angularityMcBinsVec[0][3] = {0, 6e-05, 0.0006, 0.00252, 0.00558, 0.0093, 0.0141, 0.03354, 0.06};
    ptRecoBinsVec[0] = {-15, -3.395, -1.745, -0.645, 0.29, 1.06, 1.83, 2.6, 3.37, 4.14, 4.91, 5.68, 6.505, 7.385, 8.32, 9.365, 10.575, 12.115, 14.26, 18.66, 40};
    zRecoBinsVec[0] = {-10, -1.94, -0.9, -0.42, 0.12, 0.18, 0.22, 0.26, 0.3, 0.36, 0.42, 0.48, 0.56, 0.66, 0.8, 1, 1.34, 2.06, 4.98, 10};
    angularityRecoBinsVec[0][0] = {-14.985, -4.55544, -2.15784, -0.8991, 0.68931, 0.86913, 0.98901, 1.10889, 1.22877, 1.34865, 1.4985, 1.64835, 1.82817, 2.06793, 2.36763, 2.78721, 3.41658, 4.52547, 7.10289, 14.985};
    angularityRecoBinsVec[0][1] = {-11.025, -2.6085, -1.19025, -0.49275, 0.34425, 0.43725, 0.507, 0.57675, 0.6465, 0.71625, 0.786, 0.879, 0.99525, 1.13475, 1.32075, 1.59975, 2.06475, 3.018, 6.34275, 12.225};
    angularityRecoBinsVec[0][2] = {-8.145, -1.48968, -0.65568, -0.27204, 0.16164, 0.22836, 0.26172, 0.29508, 0.32844, 0.3618, 0.41184, 0.46188, 0.5286, 0.612, 0.71208, 0.8622, 1.1124, 1.64616, 3.66444, 8.535};
    angularityRecoBinsVec[0][3] = {-4.515, -0.4956, -0.19992, -0.0798, 0.04956, 0.06804, 0.08652, 0.105, 0.12348, 0.14196, 0.16044, 0.18816, 0.22512, 0.28056, 0.3822, 0.63168, 2.72916, 4.725};
    // centrality10-20%
    ptMcBinsVec[1] = {1, 1.667, 2.102, 2.479, 2.856, 3.262, 3.697, 4.219, 4.915, 6.133, 30};
    zMcBinsVec[1] = {0, 0.557, 0.658, 0.733, 0.798, 0.855, 0.909, 1.0001};
    angularityMcBinsVec[1][0] = {0, 0.0007, 0.0609, 0.1064, 0.1533, 0.2044, 0.2667, 0.3535, 0.7};
    angularityMcBinsVec[1][1] = {0, 0.0004, 0.0128, 0.03, 0.0512, 0.0756, 0.1056, 0.1504, 0.4};
    angularityMcBinsVec[1][2] = {0, 0.0002, 0.0032, 0.01, 0.0204, 0.0324, 0.0472, 0.0702, 0.2};
    angularityMcBinsVec[1][3] = {0, 6e-05, 0.0006, 0.00252, 0.00564, 0.00936, 0.01422, 0.04566, 0.06};
    ptRecoBinsVec[1] = {-15, -2.295, -0.975, -0.04, 0.73, 1.39, 2.05, 2.71, 3.37, 4.03, 4.69, 5.35, 6.065, 6.835, 7.66, 8.595, 9.695, 11.07, 12.995, 17.065, 40};
    zRecoBinsVec[1] = {-10, -1.96, -0.84, 0.1, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.62, 0.72, 0.86, 1.06, 1.38, 2, 4, 10};
    angularityRecoBinsVec[1][0] = {-14.985, -3.56643, -1.40859, 0.44955, 0.71928, 0.83916, 0.92907, 1.01898, 1.10889, 1.22877, 1.34865, 1.4985, 1.67832, 1.88811, 2.18781, 2.60739, 3.2967, 4.70529, 9.50049, 14.985};
    angularityRecoBinsVec[1][1] = {-14.985, -2.12787, -0.80919, 0, 0.35964, 0.41958, 0.47952, 0.53946, 0.5994, 0.65934, 0.74925, 0.83916, 0.95904, 1.10889, 1.34865, 1.73826, 2.57742, 6.14385, 14.985};
    angularityRecoBinsVec[1][2] = {-14.145, -1.17753, -0.42327, -0.0171303, 0.18594, 0.21495, 0.24396, 0.27297, 0.30198, 0.33099, 0.38901, 0.44703, 0.50505, 0.59208, 0.73713, 0.96921, 1.5204, 5.03061, 14.865};
    angularityRecoBinsVec[1][3] = {-7.875, -0.34902, -0.12096, -0.00693002, 0.05823, 0.07452, 0.09081, 0.1071, 0.12339, 0.13968, 0.17226, 0.22113, 0.30258, 0.49806, 2.72979, 8.415};
    // centrality20-40%
    ptMcBinsVec[2] = {1, 1.667, 2.102, 2.479, 2.856, 3.262, 3.697, 4.219, 4.915, 6.133, 30};
    zMcBinsVec[2] = {0, 0.557, 0.658, 0.733, 0.798, 0.854, 0.908, 1.0001};
    angularityMcBinsVec[2][0] = {0, 0.0007, 0.0609, 0.1064, 0.1533, 0.2044, 0.266, 0.3521, 0.7};
    angularityMcBinsVec[2][1] = {0, 0.0004, 0.0128, 0.03, 0.0516, 0.076, 0.106, 0.1512, 0.4};
    angularityMcBinsVec[2][2] = {0, 0.0002, 0.0032, 0.01, 0.0204, 0.0324, 0.0472, 0.0702, 0.2};
    angularityMcBinsVec[2][3] = {0, 6e-05, 0.0006, 0.00252, 0.00564, 0.00936, 0.01416, 0.03882, 0.06};
    ptRecoBinsVec[2] = {-15, -0.865, 0.07, 0.73, 1.28, 1.83, 2.325, 2.82, 3.315, 3.81, 4.305, 4.855, 5.405, 6.01, 6.67, 7.44, 8.32, 9.475, 11.18, 15.47, 40};
    zRecoBinsVec[2] = {-10, -1.56, 0.14, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.6, 0.68, 0.76, 0.86, 1, 1.2, 1.52, 2.14, 4.14, 10};
    angularityRecoBinsVec[2][0] = {-14.985, -1.73826, 0.32967, 0.53946, 0.62937, 0.68931, 0.74925, 0.80919, 0.86913, 0.95904, 1.04895, 1.13886, 1.25874, 1.40859, 1.61838, 1.91808, 2.45754, 3.74625, 13.3666, 14.985};
    angularityRecoBinsVec[2][1] = {-14.985, -0.92907, 0.08991, 0.23976, 0.2997, 0.35964, 0.38961, 0.41958, 0.47952, 0.53946, 0.5994, 0.68931, 0.80919, 0.98901, 1.31868, 2.24775, 14.985};
    angularityRecoBinsVec[2][2] = {-13.275, -0.46938, 0.0396304, 0.12, 0.14679, 0.17358, 0.20037, 0.22716, 0.25395, 0.28074, 0.30753, 0.36111, 0.41469, 0.49506, 0.6558, 1.08444, 13.515};
    angularityRecoBinsVec[2][3] = {-6.525, -0.12531, 0.00999001, 0.03705, 0.05058, 0.06411, 0.07764, 0.09117, 0.1047, 0.13176, 0.17235, 0.26706, 1.14651, 7.005};
    // centrality40-60%
    ptMcBinsVec[3] = {1, 1.667, 2.102, 2.479, 2.856, 3.262, 3.697, 4.219, 4.915, 6.133, 30};
    zMcBinsVec[3] = {0, 0.557, 0.658, 0.733, 0.798, 0.854, 0.908, 1.0001};
    angularityMcBinsVec[3][0] = {0, 0.0007, 0.0616, 0.1071, 0.154, 0.2051, 0.2674, 0.3542, 0.7};
    angularityMcBinsVec[3][1] = {0, 0.0004, 0.0128, 0.03, 0.0516, 0.076, 0.106, 0.1512, 0.4};
    angularityMcBinsVec[3][2] = {0, 0.0002, 0.0032, 0.01, 0.0204, 0.0324, 0.047, 0.0698, 0.2};
    angularityMcBinsVec[3][3] = {0, 6e-05, 0.0006, 0.00252, 0.00564, 0.00936, 0.01416, 0.03864, 0.06};
    ptRecoBinsVec[3] = {-15, 0.455, 1.06, 1.5, 1.885, 2.27, 2.6, 2.93, 3.26, 3.59, 3.92, 4.305, 4.69, 5.13, 5.625, 6.175, 6.835, 7.715, 9.145, 40};
    zRecoBinsVec[3] = {-10, 0.26, 0.34, 0.4, 0.46, 0.52, 0.58, 0.64, 0.7, 0.76, 0.82, 0.9, 0.98, 1.08, 1.22, 1.44, 1.86, 3.28, 10};
    angularityRecoBinsVec[3][0] = {-14.985, 0.17982, 0.26973, 0.32967, 0.38961, 0.44955, 0.50949, 0.56943, 0.62937, 0.68931, 0.74925, 0.83916, 0.95904, 1.16883, 1.76823, 14.985};
    angularityRecoBinsVec[3][1] = {-11.925, 0.0624598, 0.11304, 0.13833, 0.16362, 0.18891, 0.2142, 0.23949, 0.26478, 0.29007, 0.31536, 0.36594, 0.41652, 0.49239, 0.64413, 1.78218, 13.365};
    angularityRecoBinsVec[3][2] = {-8.385, 0.0182397, 0.0533997, 0.0709797, 0.0885597, 0.10614, 0.12372, 0.1413, 0.15888, 0.17646, 0.21162, 0.26436, 0.38742, 9.195};
    angularityRecoBinsVec[3][3] = {-4.515, 0.00510018, 0.0142502, 0.0234002, 0.0325502, 0.0417002, 0.0508502, 0.0600002, 0.0783002, 0.1332, 4.635};
    // centrality60-80%
    ptMcBinsVec[4] = {1, 1.667, 2.102, 2.479, 2.856, 3.262, 3.697, 4.219, 4.915, 6.133, 30};
    zMcBinsVec[4] = {0, 0.557, 0.658, 0.733, 0.797, 0.853, 0.907, 1.0001};
    angularityMcBinsVec[4][0] = {0, 0.0007, 0.0609, 0.1064, 0.1526, 0.2037, 0.2653, 0.3514, 0.7};
    angularityMcBinsVec[4][1] = {0, 0.0004, 0.0128, 0.03, 0.0512, 0.0752, 0.1048, 0.1492, 0.4};
    angularityMcBinsVec[4][2] = {0, 0.0002, 0.0032, 0.01, 0.0204, 0.0322, 0.0468, 0.0694, 0.2};
    angularityMcBinsVec[4][3] = {0, 6e-05, 0.0006, 0.00246, 0.00552, 0.00924, 0.01398, 0.03042, 0.06};
    ptRecoBinsVec[4] = {-15, 1.225, 1.665, 1.995, 2.27, 2.545, 2.82, 3.095, 3.37, 3.645, 3.92, 4.195, 4.525, 4.855, 5.24, 5.735, 6.34, 7.22, 9.09, 40};
    zRecoBinsVec[4] = {-10, 0.38, 0.46, 0.52, 0.58, 0.62, 0.66, 0.7, 0.74, 0.78, 0.82, 0.86, 0.92, 0.98, 1.04, 1.14, 1.34, 2.46, 10};
    angularityRecoBinsVec[4][0] = {-7.245, 0.0347402, 0.0841502, 0.11709, 0.15003, 0.18297, 0.21591, 0.24885, 0.28179, 0.31473, 0.34767, 0.38061, 0.41355, 0.46296, 0.52884, 0.64413, 9.225};
    angularityRecoBinsVec[4][1] = {-4.755, 0.00926987, 0.0298499, 0.0504299, 0.0710099, 0.0915899, 0.11217, 0.13275, 0.15333, 0.17391, 0.19449, 0.22536, 0.26652, 0.47232, 5.535};
    angularityRecoBinsVec[4][2] = {-3.105, 0.00701997, 0.0138, 0.02058, 0.02736, 0.03414, 0.04092, 0.0477, 0.05448, 0.06126, 0.06804, 0.0816, 0.09516, 0.1155, 0.1494, 3.675};
    angularityRecoBinsVec[4][3] = {-1.575, 0.00135, 0.00474, 0.00813, 0.01152, 0.01491, 0.0183, 0.02169, 0.02508, 0.03186, 0.04881, 1.815};
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {

        hMeasured[iCent] = new TH2D("Meas" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0], zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0]);
        hTruth[iCent] = new TH2D("True" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptMcBinsVec[iCent].size() - 1, &ptMcBinsVec[iCent][0], zMcBinsVec[iCent].size() - 1, &zMcBinsVec[iCent][0]);

        hMeasuredTest[iCent] = new TH2D("MeasTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0], zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0]);
        hTruthTest[iCent] = new TH2D("TrueTest" + centralityNames[iCent], ";p_{t}, GeV/c; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}", ptMcBinsVec[iCent].size() - 1, &ptMcBinsVec[iCent][0], zMcBinsVec[iCent].size() - 1, &zMcBinsVec[iCent][0]);

        response[iCent] = new RooUnfoldResponse("response" + centralityNames[iCent], centralityNames[iCent] + "response");
        response[iCent]->Setup(hMeasured[iCent], hTruth[iCent]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularity[iCent][iLambda] = new TH2D("Meas" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0], angularityRecoBinsVec[iCent][iLambda].size() - 1, &angularityRecoBinsVec[iCent][iLambda][0]);
            hTruthAngularity[iCent][iLambda] = new TH2D("True" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptMcBinsVec[iCent].size() - 1, &ptMcBinsVec[iCent][0], angularityMcBinsVec[iCent][iLambda].size() - 1, &angularityMcBinsVec[iCent][iLambda][0]);

            hMeasuredAngularityTest[iCent][iLambda] = new TH2D("MeasTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0], angularityRecoBinsVec[iCent][iLambda].size() - 1, &angularityRecoBinsVec[iCent][iLambda][0]);
            hTruthAngularityTest[iCent][iLambda] = new TH2D("TrueTest" + centralityNames[iCent] + angularityNames[iLambda], ";p_{t}, GeV/c;" + angularityTitle[iLambda], ptMcBinsVec[iCent].size() - 1, &ptMcBinsVec[iCent][0], angularityMcBinsVec[iCent][iLambda].size() - 1, &angularityMcBinsVec[iCent][iLambda][0]);

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
        Double_t nEntries = jetTree->GetEntries();
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

                for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
                {
                    responseAngularity[centBin][iLambda]->Fill(recoJet.jetpt, recoJet.lambda[iLambda], mcJet.jetpt, mcJet.lambda[iLambda]);
                }
            }

            else
            {
                response[centBin]->Miss(mcJet.jetpt, mcJet.d0z);

                for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
                {
                    responseAngularity[centBin][iLambda]->Miss(mcJet.jetpt, mcJet.lambda[iLambda]);
                }
            }

        } // end of loop over train entries

        for (Int_t iEntry = nEntries * TrainToTestRatio; iEntry < nEntries; iEntry++)
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
    TCanvas *can = new TCanvas("can", "Closure Test Check", 1200, 1000);
    TString outPdf = "closure_check.pdf";
    can->SaveAs(outPdf + "[");

    // Create RooUnfoldBayes object and run the unfolding
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << "centrality  " << centralityTitles[iCent] << endl;
        // RooUnfoldBayes unfolding(response[iCent], hMeasuredTest[iCent], 20);

        // TH2D *hUnfolded = (TH2D *)unfolding.Hunfold();
        // hUnfolded->SetName("Unfolded" + centralityNames[iCent]);

        // drawClosure(hUnfolded, hTruthTest[iCent], hMeasuredTest[iCent], can);
        // can->cd();
        // tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + " z");
        // can->SaveAs(outPdf);

        plotIterations(can, outPdf, response[iCent], hTruthTest[iCent], hMeasuredTest[iCent], iCent, "z");

        // // TH2D *hUnfoldedAngularity[nAngularities];

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            // RooUnfoldResponse *resp = responseAngularity[iCent][iLambda];

            // RooUnfoldBayes unfoldingAngularity(responseAngularity[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], 5);
            // hUnfoldedAngularity[iLambda] = (TH2D *)unfoldingAngularity.Hunfold();
            // hUnfoldedAngularity[iLambda]->SetName("Unfolded" + centralityNames[iCent] + angularityNames[iLambda]);

            // drawClosure(hUnfoldedAngularity[iLambda], hTruthAngularityTest[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], can);
            // can->cd();
            // tex->DrawLatex(0.6, 0.35, centralityTitles[iCent] + "  " + angularityTitle[iLambda]);
            // can->SaveAs(outPdf);

            plotIterations(can, outPdf, responseAngularity[iCent][iLambda], hTruthAngularityTest[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], iCent, angularityTitle[iLambda]);
        }

        if (!readTree)
            continue;
        responseFile->cd();
        TDirectory *dir = responseFile->mkdir(centralityNames[iCent]);
        dir->cd();
        hMeasuredTest[iCent]->Write();
        hTruthTest[iCent]->Write();
        // hUnfolded->Write();
        response[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularityTest[iCent][iLambda]->Write();
            hTruthAngularityTest[iCent][iLambda]->Write();
            //    hUnfoldedAngularity[iLambda]->Write();
            responseAngularity[iCent][iLambda]->Write();
        }

    } // end of loop over centrality bins
    can->SaveAs(outPdf + "]");

    responseFile->Save();
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