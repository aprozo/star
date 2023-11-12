#include <iostream>
#include <iomanip>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"

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

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet);

const vector<vector<pair<Float_t, Float_t>>> angularityBinsRecoBorders{
    {
        // 0-80%
        {-14.985, 14.985},
        {-11.025, 12.225},
        {-8.145, 8.535},
        {-4.515, 4.725},
    },
    {
        // 0-10%
        {-14.985, 14.985},
        {-14.985, 14.985},
        {-14.145, 14.865},
        {-7.875, 8.415},
    },
    {
        // 10-20%
        {-14.985, 14.985},
        {-14.985, 14.985},
        {-13.275, 13.515},
        {-6.525, 7.005},
    },
    {
        // 20-40%
        {-14.985, 14.985},
        {-11.925, 13.365},
        {-8.385, 9.195},
        {-4.515, 4.635},
    },
    {
        // 40-60%
        {-7.245, 9.225},
        {-4.755, 5.535},
        {-3.105, 3.675},
        {-1.575, 1.815},
    },
    {
        // 60-80%
        {-0.015, 2.085},
        {-0.165, 1.125},
        {-0.195, 0.675},
        {-0.075, 0.225},
    }};

const vector<pair<Float_t, Float_t>> angularityBinsMcBorders{
    {

        {0, 0.7},
        {0, 0.4},
        {0, 0.2},
        {0, 0.06},
    }};

// vector<pair<Float_t, Float_t>> ptBinsReco =
//     {
//         // 0-80%
//         {-11.6833, 28.6833},
//         // 0-10%
//         {-13.7833, 31.7167},
//         // 10-20%
//         {-10.5167, 29.3833},
//         // 20-40%
//         {-7.01667, 25.65},
//         // 40-60%
//         {-3.05, 22.6167},
//         // 60-80%
//         {-0.716667, 20.9833}};

vector<pair<Float_t, Float_t>> ptBinsReco =
    {
        // 0-80%
        {-15, 40},
        // 0-10%
        {-15, 40},
        // 10-20%
        {-15, 40},
        // 20-40%
        {-15, 40},
        // 40-60%
        {-15, 40},
        // 60-80%
        {-15, 40}};

const Int_t nCentralityBins = 5;
double centBins[nCentralityBins + 1] = {0, 10, 20, 40, 60, 80}; // in icreasing order
TString centralityTitles[nCentralityBins] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

vector<Double_t> getEqualBining(TH1D *hist, const Int_t &nBins)
{
    Int_t nBinsHist = hist->GetNbinsX();
    Double_t integral = hist->Integral();
    Double_t sum = 0;
    vector<Double_t> binEdges;
    binEdges.push_back(hist->GetBinLowEdge(1));
    for (Int_t i = 1; i <= nBinsHist; i++)
    {
        sum += hist->GetBinContent(i);
        if (sum > integral / nBins)
        {
            binEdges.push_back(hist->GetBinLowEdge(i + 1));
            sum = 0;
        }
    }
    binEdges.push_back(hist->GetBinLowEdge(nBinsHist + 1));
    return binEdges;
}

struct Sizes;

struct Hists
{
    Hists(){}; // default constructor
    Hists(TString _name, const Int_t &nBins, const Int_t &iCent);
    Hists(TString _name, const Sizes &sizes, const Int_t &iCent);
    void Fill(StJetTreeStruct &jet);
    void Write();
    TH1D *hPt;
    TH1D *hZ;
    TH1D *hAngularity[nAngularities];
    TString name;
};
Hists::Hists(TString _name, const Int_t &nBins, const Int_t &iCent)
{
    name = _name;
    Double_t xMinPt, xMaxPt, xMinZ, xMaxZ, xMinAngularity, xMaxAngularity;
    if (name.Contains("Mc"))
    {
        xMinPt = 1;
        xMaxPt = 30;
        xMinZ = 0;
        xMaxZ = 1;
    }
    else
    {
        xMinPt = ptBinsReco[iCent].first;
        xMaxPt = ptBinsReco[iCent].second;
        xMinZ = -10;
        xMaxZ = 10;
    }

    hPt = new TH1D("hPt" + name + centralityNames[iCent], ";pt", nBins, xMinPt, xMaxPt);
    hZ = new TH1D("hZ" + name + centralityNames[iCent], ";z", nBins, xMinZ, xMaxZ);
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        if (name.Contains("Mc"))
        {

            xMinAngularity = 0;
            xMaxAngularity = angularityBinsMcBorders[iLambda].second;
        }
        else
        {
            xMinAngularity = angularityBinsRecoBorders[iCent][iLambda].first;
            xMaxAngularity = angularityBinsRecoBorders[iCent][iLambda].second;
        }
        hAngularity[iLambda] = new TH1D(Form("hAngularity_%d", iLambda) + name + centralityNames[iCent], ";angularity", nBins, xMinAngularity, xMaxAngularity);
    }
};

void Hists::Fill(StJetTreeStruct &jet)
{

    if (name.Contains("Reco") && jet.numberofconstituents == 0)
    {
        return;
    }

    hPt->Fill(jet.jetpt);
    if (jet.d0z == 1.)
        jet.d0z -= 0.0001;
    hZ->Fill(jet.d0z);

    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hAngularity[iLambda]->Fill(jet.lambda[iLambda]);
    }
};
void Hists::Write()
{
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hAngularity[iLambda]->Write();
    }
    hPt->Write();
    hZ->Write();
};

struct Sizes
{
    vector<Double_t> pt;
    vector<Double_t> z;
    vector<Double_t> lambda[4];

    void calculateBinning(const Hists &hists, const Int_t &nBins)
    {
        pt = getEqualBining(hists.hPt, nBins);
        z = getEqualBining(hists.hZ, nBins);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            lambda[iLambda] = getEqualBining(hists.hAngularity[iLambda], nBins);
        }
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
    void Print(const TString name, const Int_t &iCent)
    {
        cout << "pt" + name + "[" << iCent << "] = {";
        PrintComponent(pt);
        cout << "z" + name + "[" << iCent << "] = {";
        PrintComponent(z);
        for (Int_t iLambda = 0; iLambda < 4; iLambda++)
        {
            cout << "angularity" + name + "[" << iCent << "][" << iLambda << "] = {";
            PrintComponent(lambda[iLambda]);
        }
    }
};

Hists::Hists(const TString _name, const Sizes &sizes, const Int_t &iCent)
{
    name = _name;

    hPt = new TH1D("hPt" + name + centralityNames[iCent], ";pt", sizes.pt.size() - 1, &sizes.pt[0]);
    hZ = new TH1D("hZ" + name + centralityNames[iCent], ";z", sizes.z.size() - 1, &sizes.z[0]);
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hAngularity[iLambda] = new TH1D(Form("hAngularity_%d", iLambda) + name + centralityNames[iCent], ";angularity", sizes.lambda[iLambda].size() - 1, &sizes.lambda[iLambda][0]);
    }
};

Int_t getCentralityBin(const Float_t &centrality)
{
    for (Int_t i = 0; i < nCentralityBins; i++)
    {
        if (centrality >= centBins[i] && centrality < centBins[i + 1])
            return (i);
    }
    cout << "Error: centrality not in range" << centrality << endl;
    return -1;
}

void FillFromTree(TTree *jetTree, Hists *histsMc, Hists *histsReco)
{
    StJetTreeStruct mcJet, recoJet;
    assignTree(jetTree, mcJet, recoJet);
    Double_t nEntries = jetTree->GetEntries();
    nEntries /= 1.;
    cout << "nEntries = " << (Float_t)nEntries / 1000. << "k" << endl
         << endl;

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        Float_t progress = 0.;
        progress = (Float_t)iEntry / nEntries;
        if (iEntry % 10000 == 0)
        {
            cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;
        }
        jetTree->GetEntry(iEntry);
        Int_t centBin = getCentralityBin(recoJet.centrality);
        if (centBin < 0)
        {
            continue;
        }

        histsMc[centBin].Fill(mcJet);
        histsReco[centBin].Fill(recoJet);
    } // end of loop over  entries
}

void fillTestHistsImproved()
{
    const Int_t nRecoBins = 20;
    const Int_t nMcBins = 10;

    Hists Mc[nCentralityBins];
    Hists Reco[nCentralityBins];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        Mc[iCent] = Hists("Mc", 1000, iCent);
        Reco[iCent] = Hists("Reco", 1000, iCent);
    }

    TFile *treeFile = new TFile("../output_jets.root", "READ");
    if (treeFile->IsZombie())
    {
        return;
    }
    TTree *jetTree = (TTree *)treeFile->Get("Jets");

    FillFromTree(jetTree, Mc, Reco);

    Sizes sizesMc[nCentralityBins];
    Sizes sizesReco[nCentralityBins];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        sizesMc[iCent].calculateBinning(Mc[iCent], nMcBins);
        sizesReco[iCent].calculateBinning(Reco[iCent], nRecoBins);
        cout << "// centrality" << centralityTitles[iCent] << endl;
        sizesMc[iCent].Print("McBinsVec", iCent);
        sizesReco[iCent].Print("RecoBinsVec", iCent);
    }

    Hists McRebinned[nCentralityBins];
    Hists RecoRebinned[nCentralityBins];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        McRebinned[iCent] = Hists("McRebinned", sizesMc[iCent], iCent);
        RecoRebinned[iCent] = Hists("RecoRebinned", sizesReco[iCent], iCent);
    }

    // fill them with tree data

    FillFromTree(jetTree, McRebinned, RecoRebinned);
    TFile *outFile = new TFile("testHists.root", "RECREATE");
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        outFile->cd();
        TDirectory *centDir = outFile->mkdir(Form("cent_%s", centralityNames[iCent].Data()));
        centDir->cd();
        Mc[iCent].Write();
        Reco[iCent].Write();
        McRebinned[iCent].Write();
        RecoRebinned[iCent].Write();
    }
    outFile->Save();
    outFile->Close();
}

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet)
{
    cout << "Reading tree" << endl;

    jetTree->SetBranchAddress("Centrality", &recoJet.centrality);

    jetTree->SetBranchAddress("McJetPt", &mcJet.jetpt);
    jetTree->SetBranchAddress("McJetLambda_1_1", &mcJet.lambda[0]);
    jetTree->SetBranchAddress("McJetLambda_1_1half", &mcJet.lambda[1]);
    jetTree->SetBranchAddress("McJetLambda_1_2", &mcJet.lambda[2]);
    jetTree->SetBranchAddress("McJetLambda_1_3", &mcJet.lambda[3]);
    jetTree->SetBranchAddress("McJetD0Z", &mcJet.d0z);

    jetTree->SetBranchAddress("RecoJetPt", &recoJet.jetpt);
    jetTree->SetBranchAddress("RecoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress("RecoJetLambda_1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress("RecoJetLambda_1_1half", &recoJet.lambda[1]);
    jetTree->SetBranchAddress("RecoJetLambda_1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress("RecoJetLambda_1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress("RecoJetD0Z", &recoJet.d0z);
}