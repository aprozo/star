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

const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

const Int_t nCentralityBins = 3;
double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order
TString centralityTitles[nCentralityBins] = {"0-10%", "10-40%", "40-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_40", "40_80"};
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
    const vector<Int_t> plotIterations = {1, 2, 3, 4, 5, 10, 15, 20};
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
}

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

void createResponseMatrixAngularity3CentralityBins()
{

    Int_t nMcBins = 5;
    Int_t nRecoBins = 10;

    Bool_t readTree = kFALSE;
    gSystem->Load("libRooUnfold");
    gStyle->SetHistFillStyle(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // reco bins

    TFile *realDataFile = new TFile("../splot_2D_histogram_samePT.root", "READ");
    if (!realDataFile || realDataFile->IsZombie())
    {
        cout << "Error: file with real data not found" << endl;
        return;
    }

    TString centralityNamesTemp[nCentralityBins] = {"010", "1040", "4080"};
    TString temp = "_hist_60";
    TString varNamesTemp[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

    TH2D *hDataZ[nCentralityBins];
    TH2D *hDataLambda[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        hDataZ[iCent] = (TH2D *)realDataFile->Get("z" + temp + centralityNamesTemp[iCent] + "_2D");
        zRecoBinsVec[iCent] = getAxisVector(hDataZ[iCent]->GetXaxis());
        ptRecoBinsVec[iCent] = getAxisVector(hDataZ[iCent]->GetYaxis());
        PrintComponent(zRecoBinsVec[iCent]);
        PrintComponent(ptRecoBinsVec[iCent]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hDataLambda[iCent][iLambda] = (TH2D *)realDataFile->Get(varNamesTemp[iLambda] + temp + centralityNamesTemp[iCent] + "_2D");
            angularityRecoBinsVec[iCent][iLambda] = getAxisVector(hDataLambda[iCent][iLambda]->GetXaxis());
            PrintComponent(angularityRecoBinsVec[iCent][iLambda]);
        }
    }

    // zRecoBinsVec[0] = {-5, -0.38, 0.115, 0.17, 0.247, 0.335, 0.445, 0.61, 0.962, 2.326, 6};
    // zRecoBinsVec[1] = {-5, -0.732, 0.159, 0.247, 0.335, 0.434, 0.577, 0.786, 1.149, 2.249, 6};
    // zRecoBinsVec[2] = {-5, 0.335, 0.478, 0.599, 0.731, 0.863, 1.039, 1.292, 1.666, 2.645, 6};

    // ptRecoBinsVec[0] = {-15, -0.425, 1.115, 2.38, 3.535, 4.745, 6.395, 8.485, 11.345, 16.68, 40};
    // ptRecoBinsVec[1] = {-15, -0.425, 0.62, 1.555, 2.435, 3.425, 4.47, 5.79, 7.605, 10.85, 40};
    // ptRecoBinsVec[2] = {-15, 0.565, 1.06, 1.5, 1.94, 2.435, 2.985, 3.645, 4.525, 6.065, 40};

    // angularityRecoBinsVec[0][0] = {-5, 0.475, 0.805, 1.1575, 1.4575, 1.7875, 2.1325, 2.605, 3.3325, 4.8475, 10};
    // angularityRecoBinsVec[1][0] = {-5, 0.3925, 0.7975, 1.045, 1.2625, 1.4875, 1.765, 2.1775, 2.77, 4.1875, 10};
    // angularityRecoBinsVec[2][0] = {-5, 0.295, 0.475, 0.5875, 0.7, 0.82, 0.94, 1.15, 1.4125, 2.395, 10};

    // angularityRecoBinsVec[0][1] = {-5, 0.355, 0.6325, 0.9325, 1.1875, 1.4575, 1.735, 2.155, 2.7775, 4.39, 10};
    // angularityRecoBinsVec[1][1] = {-5, -0.14, 0.6025, 0.8125, 0.9925, 1.1875, 1.4275, 1.7575, 2.2675, 3.4675, 10};
    // angularityRecoBinsVec[2][1] = {-5, 0.16, 0.295, 0.4075, 0.4975, 0.6025, 0.7, 0.8575, 1.075, 1.7575, 10};

    // angularityRecoBinsVec[0][2] = {-5, 0.235, 0.52, 0.76, 0.9925, 1.2325, 1.48, 1.8475, 2.4175, 3.8875, 10};
    // angularityRecoBinsVec[1][2] = {-5, -0.26, 0.4825, 0.6625, 0.82, 0.9925, 1.195, 1.495, 1.9525, 3.07, 10};
    // angularityRecoBinsVec[2][2] = {-5, 0.1, 0.2125, 0.31, 0.3775, 0.46, 0.55, 0.67, 0.8425, 1.345, 10};

    // angularityRecoBinsVec[0][3] = {-5, 0.0175, 0.385, 0.565, 0.745, 0.94, 1.1425, 1.48, 1.9525, 3.34, 10};
    // angularityRecoBinsVec[1][3] = {-5, -0.2825, 0.325, 0.4675, 0.595, 0.7375, 0.895, 1.135, 1.51, 2.425, 10};
    // angularityRecoBinsVec[2][3] = {-5, 0.04, 0.1225, 0.1825, 0.2425, 0.3025, 0.3775, 0.4675, 0.6175, 1.0675, 10};

    TFile *binSizes = new TFile(Form("binningMc%iReco%i.root", nMcBins, nRecoBins), "read");
    if (!binSizes || binSizes->IsZombie())
    {
        cout << "Error: file with bin sizes not found" << endl;
        return;
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        TString dir = Form("cent_%s/", centralityNames[0].Data());

        ///=======================================
        TH1D *hPtMcRebinned = (TH1D *)binSizes->Get(dir + "hPtMcRebinned" + centralityNames[0]);
        TH1D *hZMcRebinned = (TH1D *)binSizes->Get(dir + "hZMcRebinned" + centralityNames[0]);

        ptMcBinsVec[iCent] = getAxis(hPtMcRebinned);
        zMcBinsVec[iCent] = getAxis(hZMcRebinned);

        TH1D *hAngularityMcRebinned[nAngularities];

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hAngularityMcRebinned[iLambda] = (TH1D *)binSizes->Get(dir + Form("hAngularity_%iMcRebinned", iLambda) + centralityNames[0]);
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
    }

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
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << "centrality  " << centralityTitles[iCent] << endl;

        plotIterations(can, outPdf, response[iCent], hTruthTest[iCent], hMeasuredTest[iCent], iCent, "z");

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            plotIterations(can, outPdf, responseAngularity[iCent][iLambda], hTruthAngularityTest[iCent][iLambda], hMeasuredAngularityTest[iCent][iLambda], iCent, angularityTitle[iLambda]);
        }

        if (!readTree)
            continue;
        responseFile->cd();
        TDirectory *dir = responseFile->mkdir(centralityNames[iCent]);
        dir->cd();
        hMeasuredTest[iCent]->Write();
        hTruthTest[iCent]->Write();
        response[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMeasuredAngularityTest[iCent][iLambda]->Write();
            hTruthAngularityTest[iCent][iLambda]->Write();
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