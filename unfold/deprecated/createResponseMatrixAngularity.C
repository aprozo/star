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

vector<Double_t> ptMcBinsVec[nCentralityBins];
vector<Double_t> zMcBinsVec[nCentralityBins];
vector<Double_t> angularityMcBinsVec[nCentralityBins][nAngularities];

vector<Double_t> ptRecoBinsVec[nCentralityBins];
vector<Double_t> zRecoBinsVec[nCentralityBins];
vector<Double_t> angularityRecoBinsVec[nCentralityBins][nAngularities];

const Double_t TrainToTestRatio = 0.95;

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
        if (hist->GetBinWidth(i) != 0)
        {
            hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
            hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
        }
    }
}

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

    TLegend *legend = new TLegend(0.45, 0.62, 0.55, 0.94);
    legend->AddEntry(xTruth, "True", "l");
    legend->AddEntry(xMeasured, "Measured", "l");

    can->cd(1);

    gPad->SetLogy();

    xTruth->SetMarkerStyle(20);
    xTruth->SetLineColor(kRed);
    xTruth->GetYaxis()->SetTitle("dN/dp_{t}");
    xTruth->GetYaxis()->SetRangeUser(xMeasured->GetMinimum() * 0.8, xTruth->GetMaximum() * 1.2);

    xTruth->Draw("hist ");

    xMeasured->Draw("hist same");

    xMeasured->SetMarkerStyle(21);
    xMeasured->SetLineColor(kTeal - 1);

    can->cd(3);
    gPad->SetLogy(0);

    can->cd(2);
    gPad->SetLogy();

    yTruth->SetMarkerStyle(20);
    yTruth->SetLineColor(kRed);
    yTruth->GetYaxis()->SetTitle("dN/d" + var);
    yTruth->GetYaxis()->SetRangeUser(yMeasured->GetMinimum() * 0.8, yTruth->GetMaximum() * 1.2);
    yTruth->Draw("hist ");

    yMeasured->Draw("hist same");
    // yMeasured->GetYaxis()->SetRangeUser(yMeasured->GetMinimum() * 0.8, yTruth->GetMaximum() * 1.5);
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
    tex->DrawLatex(0.7, 0.35, centralityTitles[iCent] + "  " + var);
    can->SaveAs(outPdf);
    can->Clear();
    can->SetLogz(1);
    RooUnfoldBayes lastUnfolding(response, hMeasured, 4);
    TH2D *hUnfolded = (TH2D *)lastUnfolding.Hunfold();

    TH2D *hUnfoldingMatrix = new TH2D(lastUnfolding.UnfoldingMatrix());
    hUnfoldingMatrix->SetTitle("Bin Migration probability");
    hUnfoldingMatrix->SetName("UnfoldingMatrix");
    hUnfoldingMatrix->GetXaxis()->SetTitle("bin=N_{reco}*" + var + "+pt");
    hUnfoldingMatrix->GetYaxis()->SetTitle("bin=N_{unfold}*" + var + "+pt");
    can->cd();
    hUnfoldingMatrix->GetZaxis()->SetTitleOffset(0.7);

    hUnfoldingMatrix->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis *)hUnfoldingMatrix->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.90);
    palette->SetX2NDC(0.92);

    tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + var);
    tex->DrawLatex(0.5, 0.35, "Unfolding Matrix");
    can->SaveAs(outPdf);
    delete hUnfoldingMatrix;

    can->Clear();
    can->cd();
    can->SetLogz(0);
    for (int iter = 0; iter < nIter; iter++)
    {
        fPearsonCoeffs[iter]->SetTitle((TString) "Pearson Coefficients;bin =" + var + "*N_{bins}+pt;bin =" + var + "*N_{bins}+pt");
        fPearsonCoeffs[iter]->Draw("colz");
        gPad->Update();
        fPearsonCoeffs[iter]->GetZaxis()->SetTitleOffset(0.7);
        tex->DrawLatex(0.2, 0.65, centralityTitles[iCent] + "  " + var);
        tex->DrawLatex(0.2, 0.35, Form("Pearson Coefficients Iter %i", plotIterations[iter]));

        TPaletteAxis *palette = (TPaletteAxis *)fPearsonCoeffs[iter]->GetListOfFunctions()->FindObject("palette");
        palette->SetX1NDC(0.91);
        palette->SetX2NDC(0.93);
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

vector<Double_t> getAxis(TH1D *hist)
{
    vector<Double_t> axis;
    for (Int_t i = 1; i <= hist->GetNbinsX(); i++)
    {
        axis.push_back(hist->GetBinLowEdge(i));
    }
    axis.push_back(hist->GetBinLowEdge(hist->GetNbinsX() + 1));

    // comare last 2 bins and remove the last bins in case they are equal

    if (axis[axis.size() - 1] == axis[axis.size() - 2])
    {
        axis.pop_back();
    }

    if (axis[axis.size() - 1] == 1.)
        axis.push_back(1.00001);

    return axis;
}

void PrintComponent(const vector<Double_t> &vec)
{
    for (Int_t i = 0; i < vec.size(); i++)
    {
        cout << vec[i];
        if (i != vec.size() - 1)
        {
            cout << ", ";
        }
    }
    cout << "};" << endl;
}

void createResponseMatrixAngularity()
{

    Int_t nMcBins = 5;
    Int_t nRecoBins = 10;

    TFile *binSizes = new TFile(Form("binningMc%iReco%i.root", nMcBins, nRecoBins), "read");
    if (!binSizes || binSizes->IsZombie())
    {
        cout << "Error: file with bin sizes not found" << endl;
        return;
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        TString dir = Form("cent_%s/", centralityNames[iCent].Data());

        ///=======================================
        TH1D *hPtMcRebinned = (TH1D *)binSizes->Get(dir + "hPtMcRebinned" + centralityNames[iCent]);
        TH1D *hZMcRebinned = (TH1D *)binSizes->Get(dir + "hZMcRebinned" + centralityNames[iCent]);

        ptMcBinsVec[iCent] = getAxis(hPtMcRebinned);
        zMcBinsVec[iCent] = getAxis(hZMcRebinned);

        TH1D *hAngularityMcRebinned[nAngularities];

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hAngularityMcRebinned[iLambda] = (TH1D *)binSizes->Get(dir + Form("hAngularity_%iMcRebinned", iLambda) + centralityNames[iCent]);
            angularityMcBinsVec[iCent][iLambda] = getAxis(hAngularityMcRebinned[iLambda]);
        }

        cout << "ptMcBinsVecpt[" << iCent << "] = {";
        PrintComponent(ptMcBinsVec[iCent]);
        cout << "zMcBinsVec[" << iCent << "] = {";
        PrintComponent(zMcBinsVec[iCent]);
        for (Int_t iLambda = 0; iLambda < 4; iLambda++)
        {
            cout << "angularityMcBinsVec[" << iCent << "][" << iLambda << "] = {";
            PrintComponent(angularityMcBinsVec[iCent][iLambda]);
        }

        ///=======================================
        TH1D *hPtRecoRebinned = (TH1D *)binSizes->Get(dir + "hPtRecoRebinned" + centralityNames[iCent]);
        TH1D *hZRecoRebinned = (TH1D *)binSizes->Get(dir + "hZRecoRebinned" + centralityNames[iCent]);

        ptRecoBinsVec[iCent] = getAxis(hPtRecoRebinned);
        zRecoBinsVec[iCent] = getAxis(hZRecoRebinned);

        TH1D *hAngularityRecoRebinned[nAngularities];
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hAngularityRecoRebinned[iLambda] = (TH1D *)binSizes->Get(dir + Form("hAngularity_%iRecoRebinned", iLambda) + centralityNames[iCent]);
            angularityRecoBinsVec[iCent][iLambda] = getAxis(hAngularityRecoRebinned[iLambda]);
        }

        cout << "ptRecoBinsVec[" << iCent << "] = {";
        PrintComponent(ptRecoBinsVec[iCent]);
        cout << "zRecoBinsVec[" << iCent << "] = {";
        PrintComponent(zRecoBinsVec[iCent]);
        for (Int_t iLambda = 0; iLambda < 4; iLambda++)
        {
            cout << "angularityRecoBinsVec[" << iCent << "][" << iLambda << "] = {";
            PrintComponent(angularityRecoBinsVec[iCent][iLambda]);
        }
    }

    Bool_t readTree = kFALSE;
    gSystem->Load("libRooUnfold");
    gStyle->SetHistFillStyle(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Create histograms
    TFile *responseFile = new TFile(Form("responseMc%iReco%i.root", nMcBins, nRecoBins), "READ");
    if (!responseFile || responseFile->IsZombie())
    {
        readTree = kTRUE;
        responseFile = new TFile(Form("responseMc%iReco%i.root", nMcBins, nRecoBins), "RECREATE");
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

    // // centrality0 - 10 %
    // ptMcBinsVec[0] = {1, 2.102, 2.856, 3.668, 4.828, 30};
    // zMcBinsVec[0] = {0, 0.660007, 0.801008, 0.912009, 1.00001};
    // angularityMcBinsVec[0][0] = {0, 0.0007, 0.1064, 0.2037, 0.3507, 0.7};
    // angularityMcBinsVec[0][1] = {0, 0.0004, 0.0296, 0.0748, 0.148, 0.4};
    // angularityMcBinsVec[0][2] = {0, 0.0002, 0.01, 0.0322, 0.0694, 0.2};
    // angularityMcBinsVec[0][3] = {0, 6e-05, 0.0024, 0.00912, 0.02772, 0.06};
    // ptRecoBinsVec[0] = {-15, -1.745, 0.235, 1.775, 3.26, 4.745, 6.285, 8.045, 10.19, 13.435, 40};
    // zRecoBinsVec[0] = {-10, -0.9, 0.12, 0.22, 0.3, 0.4, 0.52, 0.7, 1.06, 2.32, 10};
    // angularityRecoBinsVec[0][0] = {-14.985, -2.15784, 0.68931, 0.95904, 1.16883, 1.40859, 1.70829, 2.15784, 2.96703, 5.12487, 14.985};
    // angularityRecoBinsVec[0][1] = {-11.025, -1.19025, 0.321, 0.48375, 0.6, 0.7395, 0.90225, 1.158, 1.64625, 3.1575, 12.225};
    // angularityRecoBinsVec[0][2] = {-8.145, -0.65568, 0.16164, 0.26172, 0.32844, 0.41184, 0.51192, 0.67872, 1.029, 2.68032, 8.535};
    // angularityRecoBinsVec[0][3] = {-4.515, -0.19992, 0.04956, 0.07728, 0.105, 0.13272, 0.16968, 0.23436, 0.40068, 4.60488, 4.725};
    // // centrality10-20%
    // ptMcBinsVec[1] = {1, 2.102, 2.856, 3.668, 4.857, 30};
    // zMcBinsVec[1] = {0, 0.660007, 0.801008, 0.912009, 1.00001};
    // angularityMcBinsVec[1][0] = {0, 0.0007, 0.1064, 0.2037, 0.3514, 0.7};
    // angularityMcBinsVec[1][1] = {0, 0.0004, 0.0296, 0.0748, 0.1484, 0.4};
    // angularityMcBinsVec[1][2] = {0, 0.0002, 0.0098, 0.032, 0.069, 0.2};
    // angularityMcBinsVec[1][3] = {0, 6e-05, 0.0024, 0.00912, 0.02748, 0.06};
    // ptRecoBinsVec[1] = {-15, -0.975, 0.675, 1.995, 3.26, 4.525, 5.9, 7.44, 9.365, 12.335, 40};
    // zRecoBinsVec[1] = {-10, -0.84, 0.18, 0.28, 0.38, 0.48, 0.62, 0.84, 1.3, 3.28, 10};
    // angularityRecoBinsVec[1][0] = {-14.985, -1.43856, 0.71928, 0.92907, 1.10889, 1.31868, 1.58841, 2.00799, 2.84715, 5.87412, 14.985};
    // angularityRecoBinsVec[1][1] = {-14.985, -0.80919, 0.35964, 0.47952, 0.56943, 0.68931, 0.83916, 1.07892, 1.61838, 4.43556, 14.985};
    // angularityRecoBinsVec[1][2] = {-14.145, -0.42327, 0.18594, 0.24396, 0.30198, 0.36, 0.44703, 0.59208, 0.9402, 3.98625, 14.865};
    // angularityRecoBinsVec[1][3] = {-7.875, -0.12096, 0.05823, 0.07452, 0.09081, 0.12339, 0.17226, 0.27, 1.01934, 8.415};
    // // centrality20-40%
    // ptMcBinsVec[2] = {1, 2.102, 2.856, 3.668, 4.857, 30};
    // zMcBinsVec[2] = {0, 0.661007, 0.802008, 0.912009, 1.00001};
    // angularityMcBinsVec[2][0] = {0, 0.0007, 0.1064, 0.2037, 0.3507, 0.7};
    // angularityMcBinsVec[2][1] = {0, 0.0004, 0.0296, 0.0748, 0.1484, 0.4};
    // angularityMcBinsVec[2][2] = {0, 0.0002, 0.01, 0.0322, 0.0694, 0.2};
    // angularityMcBinsVec[2][3] = {0, 6e-05, 0.0024, 0.00912, 0.0276, 0.06};
    // ptRecoBinsVec[2] = {-15, 0.015, 1.225, 2.27, 3.26, 4.25, 5.295, 6.505, 8.045, 10.52, 40};
    // zRecoBinsVec[2] = {-10, 0.14, 0.28, 0.38, 0.48, 0.6, 0.74, 0.96, 1.38, 2.98, 10};
    // angularityRecoBinsVec[2][0] = {-14.985, 0.32967, 0.5994, 0.71928, 0.83916, 0.98901, 1.16883, 1.43856, 1.97802, 4.07592, 14.985};
    // angularityRecoBinsVec[2][1] = {-14.985, 0.08991, 0.2997, 0.35964, 0.41958, 0.50949, 0.62937, 0.80919, 1.25874, 14.8951, 14.985};
    // angularityRecoBinsVec[2][2] = {-13.275, 0.0128404, 0.14679, 0.20037, 0.25395, 0.30753, 0.3879, 0.57543, 2.07567, 13.515};
    // angularityRecoBinsVec[2][3] = {-6.525, 0.00999001, 0.03705, 0.05058, 0.06411, 0.07764, 0.1047, 0.14529, 0.32118, 7.005};
    // // centrality40-60%
    // ptMcBinsVec[3] = {1, 2.102, 2.856, 3.668, 4.857, 30};
    // zMcBinsVec[3] = {0, 0.661007, 0.802008, 0.912009, 1.00001};
    // angularityMcBinsVec[3][0] = {0, 0.0007, 0.1064, 0.2037, 0.3507, 0.7};
    // angularityMcBinsVec[3][1] = {0, 0.0004, 0.0296, 0.0748, 0.148, 0.4};
    // angularityMcBinsVec[3][2] = {0, 0.0002, 0.01, 0.0322, 0.0692, 0.2};
    // angularityMcBinsVec[3][3] = {0, 6e-05, 0.0024, 0.00912, 0.02748, 0.06};
    // ptRecoBinsVec[3] = {-15, 1.005, 1.83, 2.545, 3.205, 3.865, 4.58, 5.405, 6.505, 8.32, 40};
    // zRecoBinsVec[3] = {-10, 0.34, 0.46, 0.56, 0.66, 0.76, 0.88, 1.06, 1.4, 2.88, 10};
    // angularityRecoBinsVec[3][0] = {-14.985, 0.26973, 0.38961, 0.47952, 0.56943, 0.65934, 0.77922, 0.98901, 1.97802, 14.985};
    // angularityRecoBinsVec[3][1] = {-11.925, 0.11304, 0.16362, 0.2142, 0.26478, 0.31536, 0.39123, 0.54297, 13.365};
    // angularityRecoBinsVec[3][2] = {-8.385, 0.0533997, 0.0885597, 0.10614, 0.1413, 0.17646, 0.2292, 0.45774, 9.195};
    // angularityRecoBinsVec[3][3] = {-4.515, 0.0142502, 0.0234002, 0.0325502, 0.0417002, 0.0600002, 0.10575, 4.635};
    // // centrality60-80%
    // ptMcBinsVec[4] = {1, 2.102, 2.856, 3.668, 4.857, 30};
    // zMcBinsVec[4] = {0, 0.661007, 0.801008, 0.911009, 1.00001};
    // angularityMcBinsVec[4][0] = {0, 0.0007, 0.1064, 0.2037, 0.3507, 0.7};
    // angularityMcBinsVec[4][1] = {0, 0.0004, 0.0296, 0.0748, 0.148, 0.4};
    // angularityMcBinsVec[4][2] = {0, 0.0002, 0.01, 0.0322, 0.0694, 0.2};
    // angularityMcBinsVec[4][3] = {0, 6e-05, 0.0024, 0.00912, 0.02766, 0.06};
    // ptRecoBinsVec[4] = {-15, 1.61, 2.215, 2.71, 3.205, 3.7, 4.25, 4.91, 5.79, 7.275, 40};
    // zRecoBinsVec[4] = {-10, 0.46, 0.56, 0.64, 0.72, 0.8, 0.88, 0.98, 1.14, 2.18, 10};
    // angularityRecoBinsVec[4][0] = {-7.245, 0.0841502, 0.15003, 0.19944, 0.24885, 0.29826, 0.34767, 0.41355, 0.4959, 0.79236, 9.225};
    // angularityRecoBinsVec[4][1] = {-4.755, 0.0195599, 0.0504299, 0.0812999, 0.11217, 0.14304, 0.17391, 0.21507, 0.31797, 5.535};
    // angularityRecoBinsVec[4][2] = {-3.105, 0.00701997, 0.02058, 0.03414, 0.0477, 0.06126, 0.0816, 0.10872, 3.675};
    // angularityRecoBinsVec[4][3] = {-1.575, 0.00135, 0.00474, 0.00813, 0.01152, 0.01491, 0.02169, 0.03186, 1.815};

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
            // if (mcJet.d0z == 1.)
            //     mcJet.d0z -= 0.000001;

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
            // if (mcJet.d0z == 1.)
            //     mcJet.d0z -= 0.000001;
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
    TString outPdf = Form("closure_check_Mc%iReco%i.pdf", nMcBins, nRecoBins);
    can->SaveAs(outPdf + "[");

    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.2);

    // Create RooUnfoldBayes object and run the unfolding
    for (Int_t iCent = 0; iCent < 1; iCent++)
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