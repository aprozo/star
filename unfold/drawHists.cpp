
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TTree.h"
#include "TLatex.h"
#include "TColor.h"

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

Int_t getCentralityBin(const Float_t &centrality)
{
    if (centrality < 10)
        return 1;
    else if (centrality < 20)
        return 2;
    else if (centrality < 40)
        return 3;
    else if (centrality < 60)
        return 4;
    else if (centrality < 80)
        return 5;
    else
        return 6;
}

const Int_t nCentralityBins = 6;
TString centralityTitles[nCentralityBins] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
TString centralityNames[nCentralityBins] = {"0_80", "0_10", "10_20", "20_40", "40_60", "60_80"};

const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

void drawDescription(const Int_t &iCent, const TString &observable)
{
    TString equation;
    if (observable == "z")
        equation = "z = #frac{#vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} }{|#vec{p}_{T, jet}|^{2}}";
    else if (observable == "jet_pt")
        equation = "p_{T, jet Mc} || p_{T, jet}^{sub} = p_{T, jet}^{raw}  - #rho A_{jet}  ";
    else if (observable == "jet_pt_corr")
        equation = "p_{T, jet}^{sub} = p_{T, jet}^{raw}  - #rho A_{jet}";
    else if (observable == "lambda_1_1")
        equation = "#lambda_{1}^{1} = #sum_{i_{track} #in jet}#left(#frac{p_{T,i}}{p_{T,jet full}}#right)#left(#frac{#Delta R_{jet,i}}{R}#right)";
    else if (observable == "lambda_1_1half")
        equation = "#lambda_{1.5}^{1} = #sum_{i_{track} #in jet}#left(#frac{p_{T,i}}{p_{T,jet full}}#right)#left(#frac{#Delta R_{jet,i}}{R}#right)^{1.5}";
    else if (observable == "lambda_1_2")
        equation = "#lambda_{2}^{1} = #sum_{i_{track} #in jet}#left(#frac{p_{T,i}}{p_{T,jet full}}#right)#left(#frac{#Delta R_{jet,i}}{R}#right)^{2}";
    else if (observable == "lambda_1_3")
        equation = "#lambda_{3}^{1} = #sum_{i_{track} #in jet}#left(#frac{p_{T,i}}{p_{T,jet full}}#right)#left(#frac{#Delta R_{jet,i}}{R}#right)^{3}";
    else
        equation = "";

    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.SetTextAlign(12);
    Float_t yShift = 0.85;

    latex.DrawLatexNDC(0.58, yShift, "run14, Au-Au, #sqrt{s_{NN}} = 200 GeV");
    latex.DrawLatexNDC(0.58, yShift - 0.05, centralityTitles[iCent] + ", R = 0.4, |#eta| < 1 - R");
    latex.DrawLatexNDC(0.6, yShift - 0.10, "p_{T,D^{0}} > 1 GeV/c");
    if (!(observable == "lambdas"))
        latex.DrawLatexNDC(0.59, yShift - 0.20, equation);
}

void getBorders(const Double_t cut, TH1D *histPassed)
{
    TH1D *hist = (TH1D *)histPassed->Clone();
    hist->Rebin(2);
    hist->Scale(1.0 / hist->Integral());
    Int_t bin1 = hist->FindFirstBinAbove(cut);
    Int_t bin2 = hist->FindLastBinAbove(cut);

    Double_t x1 = hist->GetBinCenter(bin1);
    Double_t x2 = hist->GetBinCenter(bin2);

    cout << "{" << x1 << ", " << x2 << "}," << endl;
}

void drawHists(TString inputFileName = "../output_jets.root", TString outputFileName = "jetTreePlots")
{
    Bool_t readHistograms = kTRUE;

    gStyle->SetOptStat(0);
    gStyle->SetOptDate(0);

    // Double_t red[9] = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
    // Double_t green[9] = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
    // Double_t blue[9] = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
    // Double_t stops[9] = {0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    // Int_t nb = 255;
    // TColor::CreateGradientColorTable(9, stops, red, green, blue, nb);

    Int_t colors[3] = {kBlack, kRed, kBlue};
    Int_t colors2[6] = {kMagenta, kTeal - 1, kOrange + 9, kAzure, kOrange - 9, kGreen + 1};

    Float_t angularitiesMin[nAngularities] = {-0.001, -0.001, -0.001, -0.001};
    Float_t angularitiesMax[nAngularities] = {0.8, 0.4, 0.2, 0.06};

    Float_t hiAngularityMin = -15;
    Float_t hiAngularityMax = 15;

    Double_t cutBorders = 2E-5;

    TH1D *hHIRecoJetPt[nCentralityBins];
    TH1D *hMcJetPt[nCentralityBins];
    TH1D *hMcRecoJetPt[nCentralityBins];

    TH2D *hHIRecoJetPtRecoD0Pt[nCentralityBins];
    TH2D *hHIRecoJetPtMcJetPt[nCentralityBins];
    TH2D *hHIRecoJetZMcJetZ[nCentralityBins];
    TH2D *hMcRecoJetZMcJetZ[nCentralityBins];
    TH2D *hHIRecoJetZMcRecoJetZ[nCentralityBins];
    TH1D *hCentralityNTracksRatio[nCentralityBins];

    TH1D *hHIRecoJetZ[nCentralityBins];
    TH1D *hMcRecoJetZ[nCentralityBins];
    TH1D *hMcJetZ[nCentralityBins];

    TH2D *hMcJetPtMcJetZ[nCentralityBins];
    TH2D *hHIRecoJetPtHIRecoJetZ[nCentralityBins];

    const Int_t nAngularities = 4;
    TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
    TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

    TH1D *hHIRecoJetLambda[nCentralityBins][nAngularities];
    TH1D *hMcJetLambda[nCentralityBins][nAngularities];
    TH1D *hMcRecoJetLambda[nCentralityBins][nAngularities];
    TH2D *hHIRecoJetLambdaMcJetLambda[nCentralityBins][nAngularities];
    TH2D *hMcRecoJetLambdaMcJetLambda[nCentralityBins][nAngularities];
    TH2D *hHIRecoJetLambdaMcRecoJetLambda[nCentralityBins][nAngularities];

    TFile *outputFile;
    outputFile = new TFile(outputFileName + ".root", "READ");

    if (outputFile->IsZombie() || !readHistograms)
    { // if file does not exist, create it
        cout << "File " << outputFileName << ".root does not exist. Creating it." << endl;
        readHistograms = kFALSE;
        outputFile = new TFile(outputFileName + ".root", "RECREATE");
    }

    if (readHistograms)
    {

        outputFile->cd();
        cout << "Reading histograms from file " << outputFileName << ".root" << endl;
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            // hCentralityNTracksRatio[iCent] = (TH1D *)outputFile->Get("hCentralityNTracksRatio" + centralityNames[iCent]);
            hHIRecoJetPt[iCent] = (TH1D *)outputFile->Get("hHIRecoJetPt" + centralityNames[iCent]);
            if (!hHIRecoJetPt[iCent])
                cout << "hHIRecoJetPt" + centralityNames[iCent] << endl;
            hMcJetPt[iCent] = (TH1D *)outputFile->Get("hMcJetPt" + centralityNames[iCent]);
            if (!hMcJetPt[iCent])
                cout << "hMcJetPt" + centralityNames[iCent] << endl;
            hMcRecoJetPt[iCent] = (TH1D *)outputFile->Get("hMcRecoJetPt" + centralityNames[iCent]);
            if (!hMcRecoJetPt[iCent])
                cout << "hMcRecoJetPt" + centralityNames[iCent] << endl;
            hHIRecoJetPtRecoD0Pt[iCent] = (TH2D *)outputFile->Get("hHIRecoJetPtRecoD0Pt" + centralityNames[iCent]);
            if (!hHIRecoJetPtRecoD0Pt[iCent])
                cout << "hHIRecoJetPtRecoD0Pt" + centralityNames[iCent] << endl;
            hHIRecoJetPtMcJetPt[iCent] = (TH2D *)outputFile->Get("hHIRecoJetPtMcJetPt" + centralityNames[iCent]);
            if (!hHIRecoJetPtMcJetPt[iCent])
                cout << "hHIRecoJetPtMcJetPt" + centralityNames[iCent] << endl;
            hHIRecoJetZMcJetZ[iCent] = (TH2D *)outputFile->Get("hHIRecoJetZMcJetZ" + centralityNames[iCent]);
            if (!hHIRecoJetZMcJetZ[iCent])
                cout << "hHIRecoJetZMcJetZ" + centralityNames[iCent] << endl;
            hMcRecoJetZMcJetZ[iCent] = (TH2D *)outputFile->Get("hMcRecoJetZMcJetZ" + centralityNames[iCent]);
            if (!hMcRecoJetZMcJetZ[iCent])
                cout << "hMcRecoJetZMcJetZ" + centralityNames[iCent] << endl;
            hHIRecoJetZMcRecoJetZ[iCent] = (TH2D *)outputFile->Get("hHIRecoJetZMcRecoJetZ" + centralityNames[iCent]);
            if (!hHIRecoJetZMcRecoJetZ[iCent])
                cout << "hHIRecoJetZMcRecoJetZ" + centralityNames[iCent] << endl;
            hMcJetPtMcJetZ[iCent] = (TH2D *)outputFile->Get("hMcJetPtMcJetZ" + centralityNames[iCent]);
            if (!hMcJetPtMcJetZ[iCent])
                cout << "hMcJetPtMcJetZ" + centralityNames[iCent] << endl;
            hHIRecoJetPtHIRecoJetZ[iCent] = (TH2D *)outputFile->Get("hHIRecoJetPtHIRecoJetZ" + centralityNames[iCent]);
            if (!hHIRecoJetPtHIRecoJetZ[iCent])
                cout << "hHIRecoJetPtHIRecoJetZ" + centralityNames[iCent] << endl;
            hHIRecoJetZ[iCent] = (TH1D *)outputFile->Get("hHIRecoJetZ" + centralityNames[iCent]);
            if (!hHIRecoJetZ[iCent])
                cout << "hHIRecoJetZ" + centralityNames[iCent] << endl;
            hMcRecoJetZ[iCent] = (TH1D *)outputFile->Get("hMcRecoJetZ" + centralityNames[iCent]);
            if (!hMcRecoJetZ[iCent])
                cout << "hMcRecoJetZ" + centralityNames[iCent] << endl;
            hMcJetZ[iCent] = (TH1D *)outputFile->Get("hMcJetZ" + centralityNames[iCent]);
            if (!hMcJetZ[iCent])
                cout << "hMcJetZ" + centralityNames[iCent] << endl;

            for (Int_t i = 0; i < 4; i++)
            {
                TString nameAdd = centralityNames[iCent] + "_" + angularityNames[i];
                hHIRecoJetLambda[iCent][i] = (TH1D *)outputFile->Get("hHIRecoJetLambda" + nameAdd);
                if (!hHIRecoJetLambda[iCent][i])
                    cout << "hHIRecoJetLambda" + nameAdd << endl;
                hMcJetLambda[iCent][i] = (TH1D *)outputFile->Get("hMcJetLambda" + nameAdd);
                if (!hMcJetLambda[iCent][i])
                    cout << "hMcJetLambda" + nameAdd << endl;
                hMcRecoJetLambda[iCent][i] = (TH1D *)outputFile->Get("hMcRecoJetLambda" + nameAdd);
                if (!hMcRecoJetLambda[iCent][i])
                    cout << "hMcRecoJetLambda" + nameAdd << endl;
                hHIRecoJetLambdaMcJetLambda[iCent][i] = (TH2D *)outputFile->Get("hHIRecoJetLambdaMcJetLambda" + nameAdd);
                if (!hHIRecoJetLambdaMcJetLambda[iCent][i])
                    cout << "hHIRecoJetLambdaMcJetLambda" + nameAdd << endl;
                hMcRecoJetLambdaMcJetLambda[iCent][i] = (TH2D *)outputFile->Get("hMcRecoJetLambdaMcJetLambda" + nameAdd);
                if (!hMcRecoJetLambdaMcJetLambda[iCent][i])
                    cout << "hMcRecoJetLambdaMcJetLambda" + nameAdd << endl;
                hHIRecoJetLambdaMcRecoJetLambda[iCent][i] = (TH2D *)outputFile->Get("hHIRecoJetLambdaMcRecoJetLambda" + nameAdd);
                if (!hHIRecoJetLambdaMcRecoJetLambda[iCent][i])
                    cout << "hHIRecoJetLambdaMcRecoJetLambda" + nameAdd << endl;
            }
        }
    }

    else
    {
        cout << "Reading input tree" << endl;
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            hCentralityNTracksRatio[iCent] = new TH1D("hCentralityNTracksRatio" + centralityNames[iCent], "; (N_{HI} + N_{Mc})/ N_{HI}", 300, 1, 2.2);
            hCentralityNTracksRatio[iCent]->SetLineColor(colors2[iCent]);

            hHIRecoJetPt[iCent] = new TH1D("hHIRecoJetPt" + centralityNames[iCent], "; p_{t}, GeV/c; Counts", 600, -30, 40);
            hMcJetPt[iCent] = new TH1D("hMcJetPt" + centralityNames[iCent], "; p_{t}, GeV/c; Counts", 200, 0, 30);
            hMcRecoJetPt[iCent] = new TH1D("hMcRecoJetPt" + centralityNames[iCent], "; p_{t}, GeV/c; Counts", 200, 0, 30);

            hHIRecoJetZ[iCent] = new TH1D("hHIRecoJetZ" + centralityNames[iCent], "; z;  Counts", 800, -10, 10);
            hMcRecoJetZ[iCent] = new TH1D("hMcRecoJetZ" + centralityNames[iCent], ";  z;  Counts", 400, -0.1, 1.1);
            hMcJetZ[iCent] = new TH1D("hMcJetZ" + centralityNames[iCent], ";  z;  Counts", 400, -0.1, 1.1);

            hHIRecoJetPtRecoD0Pt[iCent] = new TH2D("hHIRecoJetPtRecoD0Pt" + centralityNames[iCent], "; Jet p_{t}, GeV/c;  D_{0} p_{t}, GeV/c", 400, 0, 30, 400, 0, 30);
            hHIRecoJetPtMcJetPt[iCent] = new TH2D("hHIRecoJetPtMcJetPt" + centralityNames[iCent], "; HI Reco Jet p_{t}, GeV/c;  Mc Jet p_{t}, GeV/c", 400, -1, 40, 400, 0, 30);
            hHIRecoJetZMcJetZ[iCent] = new TH2D("hHIRecoJetZMcJetZ" + centralityNames[iCent], "; HI Reco Jet z;  Mc Jet z", 400, -2, 3, 400, -0.1, 1.1);
            hMcRecoJetZMcJetZ[iCent] = new TH2D("hMcRecoJetZMcJetZ" + centralityNames[iCent], "; Mc Reco Jet z;  Mc Jet z", 400, -0.1, 1.1, 400, -0.1, 1.1);
            hHIRecoJetZMcRecoJetZ[iCent] = new TH2D("hHIRecoJetZMcRecoJetZ" + centralityNames[iCent], "; HI Reco Jet z;  Mc Reco Jet z", 400, -2, 3, 400, -0.1, 1.1);
            hMcJetPtMcJetZ[iCent] = new TH2D("hMcJetPtMcJetZ" + centralityNames[iCent], ";   Jet p_{t}, GeV/c; z; Counts", 400, 0, 30, 400, -0.1, 1.1);
            hHIRecoJetPtHIRecoJetZ[iCent] = new TH2D("hHIRecoJetPtHIRecoJetZ" + centralityNames[iCent], ";   Jet p_{t}, GeV/c; z; Counts", 400, -1, 40, 400, -2, 3);

            for (Int_t i = 0; i < 4; i++)
            {
                TString nameAdd = centralityNames[iCent] + "_" + angularityNames[i];
                hHIRecoJetLambda[iCent][i] = new TH1D("hHIRecoJetLambda" + nameAdd, Form("; %s; Counts", angularityTitle[i].Data()), 2000, hiAngularityMin, hiAngularityMax);
                hHIRecoJetLambda[iCent][i]->SetLineColor(colors2[i]);
                hMcJetLambda[iCent][i] = new TH1D("hMcJetLambda" + nameAdd, Form("; %s; Counts", angularityTitle[i].Data()), 500, angularitiesMin[i], angularitiesMax[i]);
                hMcJetLambda[iCent][i]->SetLineColor(colors2[i]);
                hMcRecoJetLambda[iCent][i] = new TH1D("hMcRecoJetLambda" + nameAdd, Form("; %s; Counts", angularityTitle[i].Data()), 500, angularitiesMin[i], angularitiesMax[i]);
                hMcRecoJetLambda[iCent][i]->SetLineColor(colors2[i]);
                hHIRecoJetLambdaMcJetLambda[iCent][i] = new TH2D("hHIRecoJetLambdaMcJetLambda" + nameAdd, Form(";HI Reco Jet %s; Mc Jet %s", angularityTitle[i].Data(), angularityTitle[i].Data()), 500, hiAngularityMin, hiAngularityMax, 500, angularitiesMin[i], angularitiesMax[i]);
                hMcRecoJetLambdaMcJetLambda[iCent][i] = new TH2D("hMcRecoJetLambdaMcJetLambda" + nameAdd, Form(";Mc Reco Jet %s; Mc Jet %s", angularityTitle[i].Data(), angularityTitle[i].Data()), 500, angularitiesMin[i], angularitiesMax[i], 500, angularitiesMin[i], angularitiesMax[i]);
                hHIRecoJetLambdaMcRecoJetLambda[iCent][i] = new TH2D("hHIRecoJetLambdaMcRecoJetLambda" + nameAdd, Form("; HI Reco Jet %s; Mc Reco Jet %s", angularityTitle[i].Data(), angularityTitle[i].Data()), 500, hiAngularityMin, hiAngularityMax, 500, angularitiesMin[i], angularitiesMax[i]);
            }
        }

        StJetTreeStruct mcJet, recoJet, mcRecoJet;
        TFile *myFile; // Open the file containing the tree.
        myFile = new TFile(inputFileName, "READ");
        if (!myFile || myFile->IsZombie())
        {
            return;
        }

        TTree *jetTree = (TTree *)myFile->Get("Jets");
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

        jetTree->SetBranchAddress("McRecoJetPt", &mcRecoJet.jetpt);
        jetTree->SetBranchAddress("McRecoJetEta", &mcRecoJet.jeteta);
        jetTree->SetBranchAddress("McRecoJetPhi", &mcRecoJet.jetphi);
        jetTree->SetBranchAddress("McRecoJetE", &mcRecoJet.jetenergy);
        jetTree->SetBranchAddress("McRecoJetArea", &mcRecoJet.jetarea);
        jetTree->SetBranchAddress("McRecoJetNConst", &mcRecoJet.numberofconstituents);
        jetTree->SetBranchAddress("McRecoJetLambda_1_1", &mcRecoJet.lambda[0]);
        jetTree->SetBranchAddress("McRecoJetLambda_1_1half", &mcRecoJet.lambda[1]);
        jetTree->SetBranchAddress("McRecoJetLambda_1_2", &mcRecoJet.lambda[2]);
        jetTree->SetBranchAddress("McRecoJetLambda_1_3", &mcRecoJet.lambda[3]);
        jetTree->SetBranchAddress("McRecoJetD0Z", &mcRecoJet.d0z);
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

        for (Int_t iEntry = 0; iEntry < (Int_t)jetTree->GetEntries(); iEntry++)
        {
            jetTree->GetEntry(iEntry);
            Int_t centBin = getCentralityBin(recoJet.centrality);

            hMcJetPt[centBin]->Fill(mcJet.jetpt);
            hMcJetZ[centBin]->Fill(mcJet.d0z);
            hMcJetPtMcJetZ[centBin]->Fill(mcJet.jetpt, mcJet.d0z);

            for (Int_t iLambda = 0; iLambda < 4; iLambda++)
            {
                hMcJetLambda[centBin][iLambda]->Fill(mcJet.lambda[iLambda]);
            }

            if (mcRecoJet.numberofconstituents != 0)
            {
                hMcRecoJetPt[centBin]->Fill(mcRecoJet.jetpt);
                hMcRecoJetZMcJetZ[centBin]->Fill(mcRecoJet.d0z, mcJet.d0z);
                hMcRecoJetZ[centBin]->Fill(mcRecoJet.d0z);
            }

            if (recoJet.numberofconstituents == 0)
                continue;

            hHIRecoJetPtHIRecoJetZ[centBin]->Fill(recoJet.jetpt, recoJet.d0z);
            hCentralityNTracksRatio[centBin]->Fill((recoJet.grefmult + mcJet.grefmult) / recoJet.grefmult);
            hHIRecoJetPt[centBin]->Fill(recoJet.jetpt);

            hHIRecoJetPtRecoD0Pt[centBin]->Fill(recoJet.jetpt, recoJet.d0pt);
            hHIRecoJetPtMcJetPt[centBin]->Fill(recoJet.jetpt, mcJet.jetpt);
            hHIRecoJetZMcJetZ[centBin]->Fill(recoJet.d0z, mcJet.d0z);
            if (mcRecoJet.numberofconstituents != 0)
                hHIRecoJetZMcRecoJetZ[centBin]->Fill(recoJet.d0z, mcRecoJet.d0z);

            hHIRecoJetZ[centBin]->Fill(recoJet.d0z);

            for (Int_t iLambda = 0; iLambda < 4; iLambda++)
            {
                hHIRecoJetLambda[centBin][iLambda]->Fill(recoJet.lambda[iLambda]);
                hMcRecoJetLambda[centBin][iLambda]->Fill(mcRecoJet.lambda[iLambda]);
                hHIRecoJetLambdaMcJetLambda[centBin][iLambda]->Fill(recoJet.lambda[iLambda], mcJet.lambda[iLambda]);
                hMcRecoJetLambdaMcJetLambda[centBin][iLambda]->Fill(mcRecoJet.lambda[iLambda], mcJet.lambda[iLambda]);
                hHIRecoJetLambdaMcRecoJetLambda[centBin][iLambda]->Fill(recoJet.lambda[iLambda], mcRecoJet.lambda[iLambda]);
            }
        }

        for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
        {
            hHIRecoJetPt[0]->Add(hHIRecoJetPt[iCent]);
            hMcJetPt[0]->Add(hMcJetPt[iCent]);
            hMcRecoJetPt[0]->Add(hMcRecoJetPt[iCent]);
            hHIRecoJetPtRecoD0Pt[0]->Add(hHIRecoJetPtRecoD0Pt[iCent]);
            hHIRecoJetPtMcJetPt[0]->Add(hHIRecoJetPtMcJetPt[iCent]);
            hHIRecoJetZMcJetZ[0]->Add(hHIRecoJetZMcJetZ[iCent]);
            hMcRecoJetZMcJetZ[0]->Add(hMcRecoJetZMcJetZ[iCent]);
            hHIRecoJetZMcRecoJetZ[0]->Add(hHIRecoJetZMcRecoJetZ[iCent]);
            hCentralityNTracksRatio[0]->Add(hCentralityNTracksRatio[iCent]);

            hHIRecoJetZ[0]->Add(hHIRecoJetZ[iCent]);
            hMcRecoJetZ[0]->Add(hMcRecoJetZ[iCent]);
            hMcJetZ[0]->Add(hMcJetZ[iCent]);

            hMcJetPtMcJetZ[0]->Add(hMcJetPtMcJetZ[iCent]);
            hHIRecoJetPtHIRecoJetZ[0]->Add(hHIRecoJetPtHIRecoJetZ[iCent]);
            for (Int_t iLambda = 0; iLambda < 4; iLambda++)
            {
                hHIRecoJetLambda[0][iLambda]->Add(hHIRecoJetLambda[iCent][iLambda]);
                hMcJetLambda[0][iLambda]->Add(hMcJetLambda[iCent][iLambda]);
                hMcRecoJetLambda[0][iLambda]->Add(hMcRecoJetLambda[iCent][iLambda]);
                hHIRecoJetLambdaMcJetLambda[0][iLambda]->Add(hHIRecoJetLambdaMcJetLambda[iCent][iLambda]);
                hMcRecoJetLambdaMcJetLambda[0][iLambda]->Add(hMcRecoJetLambdaMcJetLambda[iCent][iLambda]);
                hHIRecoJetLambdaMcRecoJetLambda[0][iLambda]->Add(hHIRecoJetLambdaMcRecoJetLambda[iCent][iLambda]);
            }
        }
    }

    TCanvas *can = new TCanvas("can", "", 800, 600);
    can->SetLogz();
    can->SaveAs(outputFileName + ".pdf[");
    TLegend *leg = new TLegend(0.12, 0.7, 0.3, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.055);
    latex->SetTextFont(42);
    latex->SetTextAlign(12);

    can->cd();

    //     hCentralityNTracksRatio[1]->Draw("hist");
    // hCentralityNTracksRatio[0]->Draw("same hist");
    // hCentralityNTracksRatio[2]->Draw("same hist");

    // for (Int_t iCent = 2; iCent >= 0; iCent--)
    // {
    //     leg->AddEntry(hCentralityNTracksRatio[iCent], centralityTitles[iCent], "l");
    // }
    // leg->Draw();
    // latex->DrawLatexNDC(0.32, 0.95, "#frac{N_{HI tracks} + N_{Mc reco tracks}}{ N_{HI tracks}}");
    // can->SaveAs(outputFileName + ".pdf");
    // leg->Clear();
    can->SetLogy();
    for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    {
        hHIRecoJetPt[iCent]->GetYaxis()->SetRangeUser(1, hHIRecoJetPt[nCentralityBins - 1]->GetMaximum() * 1.2);
        hHIRecoJetPt[iCent]->SetLineColor(colors2[iCent - 1]);
        hHIRecoJetPt[iCent]->DrawNormalized(iCent == 1 ? "hist" : "same hist");
        leg->AddEntry(hHIRecoJetPt[iCent], centralityTitles[iCent], "l");
    }

    leg->Draw();
    latex->DrawLatexNDC(0.32, 0.95, "HI Reco Jet p_{t}");
    can->SaveAs(outputFileName + ".pdf");
    leg->Clear();

    for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    {
        hMcJetPt[iCent]->SetLineColor(colors2[iCent - 1]);
        hMcJetPt[iCent]->DrawNormalized(iCent == 1 ? "hist" : "same hist");
        leg->AddEntry(hMcJetPt[iCent], centralityTitles[iCent], "l");
    }
    leg->Draw();
    latex->DrawLatexNDC(0.32, 0.95, "Mc Jet p_{t}");
    can->SaveAs(outputFileName + ".pdf");
    leg->Clear();

    for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    {
        hHIRecoJetZ[iCent]->SetLineColor(colors2[iCent - 1]);
        hHIRecoJetZ[iCent]->DrawNormalized(iCent == 1 ? "hist" : "same hist");
        leg->AddEntry(hHIRecoJetZ[iCent], centralityTitles[iCent], "l");
    }

    leg->Draw();
    latex->DrawLatexNDC(0.32, 0.95, "HI Reco Jet z");

    can->SaveAs(outputFileName + ".pdf");
    leg->Clear();
    for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
    {
        hMcJetZ[iCent]->SetLineColor(colors2[iCent - 1]);
        hMcJetZ[iCent]->DrawNormalized(iCent == 1 ? "hist" : "same hist");
        leg->AddEntry(hMcJetZ[iCent], centralityTitles[iCent], "l");
    }
    leg->Draw();
    latex->DrawLatexNDC(0.42, 0.95, "Mc Jet z");
    can->SaveAs(outputFileName + ".pdf");
    leg->Clear();

    TLine *line = new TLine();
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(10);

    for (Int_t iLambda = 0; iLambda <= 3; iLambda++)
    {
        for (Int_t iCent = 1; iCent < nCentralityBins; iCent++)
        {
            hHIRecoJetLambda[iCent][iLambda]->GetYaxis()->SetRangeUser(1, hHIRecoJetLambda[nCentralityBins - 1][iLambda]->GetMaximum() * 1.2);
            hHIRecoJetLambda[iCent][iLambda]->SetLineColor(colors2[iCent - 1]);
            hHIRecoJetLambda[iCent][iLambda]->DrawNormalized(iCent == 1 ? "hist" : "same hist");
            leg->AddEntry(hHIRecoJetLambda[iCent][iLambda], centralityTitles[iCent], "l");
        }

        line->DrawLine(hiAngularityMin, cutBorders, hiAngularityMax, cutBorders);
        latex->DrawLatexNDC(0.42, 0.95, "HI " + angularityTitle[iLambda]);
        leg->Draw();
        can->SaveAs(outputFileName + ".pdf");
        leg->Clear();
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        can->SetLogy();
        hHIRecoJetPt[iCent]->SetLineColor(colors[0]);
        hMcJetPt[iCent]->SetLineColor(colors[1]);
        hMcRecoJetPt[iCent]->SetLineColor(colors[2]);
        hHIRecoJetZ[iCent]->SetLineColor(colors[0]);
        hMcJetZ[iCent]->SetLineColor(colors[1]);
        hMcRecoJetZ[iCent]->SetLineColor(colors[2]);

        hHIRecoJetPt[iCent]->GetYaxis()->SetRangeUser(1, hMcJetPt[iCent]->GetMaximum() * 1.05);
        hHIRecoJetPt[iCent]->Draw("hist");
        hMcJetPt[iCent]->Draw("same hist");
        hMcRecoJetPt[iCent]->Draw("same hist");

        leg->AddEntry(hMcJetPt[iCent], "Mc Jet", "l");
        leg->AddEntry(hMcRecoJetPt[iCent], "Mc Reco Jet", "l");
        leg->AddEntry(hHIRecoJetPt[iCent], "HI Reco Jet", "l");
        leg->Draw();
        drawDescription(iCent, "jet_pt");
        can->SaveAs(outputFileName + ".pdf");
        leg->Clear();

        can->SetLogy(0);
        hHIRecoJetPtMcJetPt[iCent]->Draw("colz");
        drawDescription(iCent, "jet_pt");
        can->SaveAs(outputFileName + ".pdf");

        can->SetLogy(1);
        hHIRecoJetZ[iCent]->GetYaxis()->SetRangeUser(1, hMcJetZ[iCent]->GetMaximum() * 1.05);
        hHIRecoJetZ[iCent]->Draw(" hist");
        hMcJetZ[iCent]->Draw("same hist");
        hMcRecoJetZ[iCent]->Draw("same hist");

        leg->AddEntry(hMcJetZ[iCent], "Mc Jet", "l");
        leg->AddEntry(hMcRecoJetZ[iCent], "Mc Reco Jet", "l");
        leg->AddEntry(hHIRecoJetZ[iCent], "HI Reco Jet", "l");
        leg->Draw();
        drawDescription(iCent, "z");
        can->SaveAs(outputFileName + ".pdf");
        leg->Clear();

        can->SetLogy(0);
        hHIRecoJetZMcJetZ[iCent]->Draw("colz");
        drawDescription(iCent, "z");
        can->SaveAs(outputFileName + ".pdf");

        hMcRecoJetZMcJetZ[iCent]->Draw("colz");
        drawDescription(iCent, "z");
        can->SaveAs(outputFileName + ".pdf");

        hHIRecoJetZMcRecoJetZ[iCent]->Draw("colz");
        drawDescription(iCent, "z");
        can->SaveAs(outputFileName + ".pdf");

        //======================================== Angularities==============================================//
        //===================================================================================================//
        can->SetLogy(1);
        for (Int_t iLambda = 0; iLambda <= 3; iLambda++)
        {
            hHIRecoJetLambda[iCent][iLambda]->GetYaxis()->SetRangeUser(1, hHIRecoJetLambda[iCent][3]->GetMaximum() * 1.05);
            hHIRecoJetLambda[iCent][iLambda]->SetLineColor(colors2[iLambda]);
            hHIRecoJetLambda[iCent][iLambda]->Draw(iLambda == 0 ? "hist" : "same hist");
            leg->AddEntry(hHIRecoJetLambda[iCent][iLambda], angularityTitle[iLambda], "l");
        }
        leg->Draw();
        latex->DrawLatexNDC(0.32, 0.95, "HI Reco Jet Lambda");
        drawDescription(iCent, "lambdas");
        can->SaveAs(outputFileName + ".pdf");

        leg->Clear();
        for (Int_t iLambda = 0; iLambda <= 3; iLambda++)
        {
            hMcJetLambda[iCent][iLambda]->GetYaxis()->SetRangeUser(1, hMcJetLambda[iCent][3]->GetMaximum() * 1.05);
            hMcJetLambda[iCent][iLambda]->SetLineColor(colors2[iLambda]);
            hMcJetLambda[iCent][iLambda]->Draw(iLambda == 0 ? "hist" : "same hist");
            leg->AddEntry(hMcJetLambda[iCent][iLambda], angularityTitle[iLambda], "l");
        }
        leg->Draw();
        latex->DrawLatexNDC(0.32, 0.95, "Mc Jet Lambda");
        drawDescription(iCent, "lambdas");
        can->SaveAs(outputFileName + ".pdf");

        leg->Clear();
        for (Int_t iLambda = 0; iLambda <= 3; iLambda++)
        {
            hMcRecoJetLambda[iCent][iLambda]->GetYaxis()->SetRangeUser(1, hMcRecoJetLambda[iCent][3]->GetMaximum() * 1.05);
            hMcRecoJetLambda[iCent][iLambda]->SetLineColor(colors2[iLambda]);
            hMcRecoJetLambda[iCent][iLambda]->Draw(iLambda == 0 ? "hist" : "same hist");
            leg->AddEntry(hMcRecoJetLambda[iCent][iLambda], angularityTitle[iLambda], "l");
        }
        leg->Draw();
        latex->DrawLatexNDC(0.32, 0.95, "Mc Reco Jet Lambda");
        drawDescription(iCent, "lambdas");
        can->SaveAs(outputFileName + ".pdf");

        leg->Clear();
        for (Int_t iLambda = 0; iLambda <= 3; iLambda++)
        {
            hHIRecoJetLambda[iCent][iLambda]->GetYaxis()->SetRangeUser(1, hMcJetLambda[iCent][3]->GetMaximum() * 1.05);

            hHIRecoJetLambda[iCent][iLambda]->SetLineColor(colors[0]);
            hHIRecoJetLambda[iCent][iLambda]->Draw("hist ");

            hMcJetLambda[iCent][iLambda]->SetLineColor(colors[1]);
            hMcJetLambda[iCent][iLambda]->Draw("hist same");
            hMcRecoJetLambda[iCent][iLambda]->SetLineColor(colors[2]);
            hMcRecoJetLambda[iCent][iLambda]->Draw("hist same");

            leg->AddEntry(hMcJetLambda[iCent][iLambda], "Mc", "l");
            leg->AddEntry(hMcRecoJetLambda[iCent][iLambda], "Mc Reco", "l");
            leg->AddEntry(hHIRecoJetLambda[iCent][iLambda], "HI Reco", "l");

            leg->Draw();
            drawDescription(iCent, angularityNames[iLambda]);
            can->SaveAs(outputFileName + ".pdf");
            leg->Clear();
        }
        can->SetLogy(0);

        for (Int_t iLambda = 0; iLambda < 4; iLambda++)
        {
            hHIRecoJetLambdaMcJetLambda[iCent][iLambda]->Draw("colz");
            drawDescription(iCent, angularityNames[iLambda]);
            can->SaveAs(outputFileName + ".pdf");
            hMcRecoJetLambdaMcJetLambda[iCent][iLambda]->Draw("colz");
            drawDescription(iCent, angularityNames[iLambda]);
            can->SaveAs(outputFileName + ".pdf");
            hHIRecoJetLambdaMcRecoJetLambda[iCent][iLambda]->Draw("colz");
            drawDescription(iCent, angularityNames[iLambda]);
            can->SaveAs(outputFileName + ".pdf");
        }
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << "{ //" << centralityTitles[iCent] << endl;
        for (Int_t iLambda = 0; iLambda < 4; iLambda++)
        {

            getBorders(cutBorders, hHIRecoJetLambda[iCent][iLambda]);
        }
        cout << "}," << endl;
    }

    cout << "  pt borders: " << endl;
    cout << "{" << endl;
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << " //" << centralityTitles[iCent] << endl;
        getBorders(10E-6, hHIRecoJetPt[iCent]);
    }
    cout << "}," << endl;

    can->SaveAs(outputFileName + ".pdf]");
    if (!readHistograms)
    {
        outputFile->cd();
        outputFile->GetList()->Write();
        outputFile->Close();
    }
}