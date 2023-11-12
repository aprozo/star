// Description: This file fills the histograms for the test data

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

void fillTestHists()
{

    const Int_t nRecoBins = 20;
    const Int_t nMcBins = 10;

    const Int_t nAngularities = 4;
    TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
    TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

    TH1D *hMcPt[nCentralityBins];
    TH1D *hMcZ[nCentralityBins];
    TH1D *hMcAngularity[nCentralityBins][nAngularities];

    TH1D *hRecoPt[nCentralityBins];
    TH1D *hRecoZ[nCentralityBins];
    TH1D *hRecoAngularity[nCentralityBins][nAngularities];
    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        hMcPt[iCent] = new TH1D("hMcPt" + centralityNames[iCent], ";pt", 10000, 1, 30);
        hMcZ[iCent] = new TH1D("hMcZ" + centralityNames[iCent], ";z", 10000, 0, 1.001);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularity[iCent][iLambda] = new TH1D(Form("hMcAngularity_%d", iLambda) + centralityNames[iCent], ";angularity", 10000, 0, angularityBinsMcBorders[iLambda].second);
        }
        hRecoPt[iCent] = new TH1D("hRecoPt" + centralityNames[iCent], ";pt", 10000, ptBinsReco[iCent].first, ptBinsReco[iCent].second);
        hRecoZ[iCent] = new TH1D("hRecoZ" + centralityNames[iCent], ";z", 10000, -10, 10);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hRecoAngularity[iCent][iLambda] = new TH1D(Form("hRecoAngularity_%d", iLambda) + centralityNames[iCent], ";angularity", 10000, angularityBinsRecoBorders[iCent][iLambda].first, angularityBinsRecoBorders[iCent][iLambda].second);
        }
    }

    TFile *treeFile = new TFile("../output_jets.root", "READ");
    if (treeFile->IsZombie())
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
        Float_t progress = 0.;
        progress = (Float_t)iEntry / nEntries;
        if (iEntry % 10000 == 0)
        {
            cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;
        }
        jetTree->GetEntry(iEntry);

        Int_t centBin = getCentralityBin(recoJet.centrality);
        Bool_t isRecoJetFound = (recoJet.numberofconstituents != 0);

        hMcPt[centBin]->Fill(mcJet.jetpt);
        hMcZ[centBin]->Fill(mcJet.d0z);

        if (isRecoJetFound)
        {
            hRecoPt[centBin]->Fill(recoJet.jetpt);
            hRecoZ[centBin]->Fill(recoJet.d0z);
        }

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularity[centBin][iLambda]->Fill(mcJet.lambda[iLambda]);
            if (isRecoJetFound)
                hRecoAngularity[centBin][iLambda]->Fill(recoJet.lambda[iLambda]);
        }
    } // end of loop over  entries

    vector<Double_t> ptMcBinsVec[nCentralityBins];
    vector<Double_t> zMcBinsVec[nCentralityBins];
    vector<Double_t> angularityMcBinsVec[nCentralityBins][nAngularities];

    vector<Double_t> ptRecoBinsVec[nCentralityBins];
    vector<Double_t> zRecoBinsVec[nCentralityBins];
    vector<Double_t> angularityRecoBinsVec[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        vector<Double_t> ptMcBins = getEqualBining(hMcPt[iCent], nMcBins);
        vector<Double_t> zMcBins = getEqualBining(hMcZ[iCent], nMcBins);
        vector<Double_t> angularityMcBins[nAngularities];
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            angularityMcBins[iLambda] = getEqualBining(hMcAngularity[iCent][iLambda], nMcBins);
        }

        vector<Double_t> ptRecoBins = getEqualBining(hRecoPt[iCent], nRecoBins);
        vector<Double_t> zRecoBins = getEqualBining(hRecoZ[iCent], nRecoBins);
        vector<Double_t> angularityRecoBins[nAngularities];
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            angularityRecoBins[iLambda] = getEqualBining(hRecoAngularity[iCent][iLambda], nRecoBins);
        }

        ptMcBinsVec[iCent] = ptMcBins;
        zMcBinsVec[iCent] = zMcBins;
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            angularityMcBinsVec[iCent][iLambda] = angularityMcBins[iLambda];
        }

        ptRecoBinsVec[iCent] = ptRecoBins;
        zRecoBinsVec[iCent] = zRecoBins;
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            angularityRecoBinsVec[iCent][iLambda] = angularityRecoBins[iLambda];
        }
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        cout << "// centrality" << centralityTitles[iCent] << endl;

        cout << "ptMcBinsVec[" << iCent << "] = {";
        for (Int_t i = 0; i < ptMcBinsVec[iCent].size(); i++)
        {
            cout << ptMcBinsVec[iCent][i];
            if (i != ptMcBinsVec[iCent].size() - 1)
            {
                cout << ", ";
            }
        }
        cout << "};" << endl;

        cout << "zMcBinsVec[" << iCent << "] = {";
        for (Int_t i = 0; i < zMcBinsVec[iCent].size(); i++)
        {
            cout << zMcBinsVec[iCent][i];
            if (i != zMcBinsVec[iCent].size() - 1)
            {
                cout << ", ";
            }
        }
        cout << "};" << endl;

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            cout << "angularityMcBinsVec[" << iCent << "][" << iLambda << "] = {";
            for (Int_t i = 0; i < angularityMcBinsVec[iCent][iLambda].size(); i++)
            {
                cout << angularityMcBinsVec[iCent][iLambda][i];
                if (i != angularityMcBinsVec[iCent][iLambda].size() - 1)
                {
                    cout << ", ";
                }
            }
            cout << "};" << endl;
        }
    }

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {

        cout << "// centrality" << centralityTitles[iCent] << endl;

        cout << "ptRecoBinsVec[" << iCent << "] = {";
        for (Int_t i = 0; i < ptRecoBinsVec[iCent].size(); i++)
        {
            cout << ptRecoBinsVec[iCent][i];
            if (i != ptRecoBinsVec[iCent].size() - 1)
            {
                cout << ", ";
            }
        }
        cout << "};" << endl;

        cout << "zRecoBinsVec[" << iCent << "] = {";
        for (Int_t i = 0; i < zRecoBinsVec[iCent].size(); i++)
        {
            cout << zRecoBinsVec[iCent][i];
            if (i != zRecoBinsVec[iCent].size() - 1)
            {
                cout << ", ";
            }
        }
        cout << "};" << endl;

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            cout << "angularityRecoBinsVec[" << iCent << "][" << iLambda << "] = {";
            for (Int_t i = 0; i < angularityRecoBinsVec[iCent][iLambda].size(); i++)
            {
                cout << angularityRecoBinsVec[iCent][iLambda][i];
                if (i != angularityRecoBinsVec[iCent][iLambda].size() - 1)
                {
                    cout << ", ";
                }
            }
            cout << "};" << endl;
        }
    }

    // Now fill the histograms with the new binning

    TH1D *hMcPtRebinned[nCentralityBins];
    TH1D *hMcZRebinned[nCentralityBins];
    TH1D *hMcAngularityRebinned[nCentralityBins][nAngularities];

    TH1D *hRecoPtRebinned[nCentralityBins];
    TH1D *hRecoZRebinned[nCentralityBins];
    TH1D *hRecoAngularityRebinned[nCentralityBins][nAngularities];

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        hMcPtRebinned[iCent] = new TH1D("hMcPtRebinned" + centralityNames[iCent], ";pt", ptMcBinsVec[iCent].size() - 1, &ptMcBinsVec[iCent][0]);
        hMcZRebinned[iCent] = new TH1D("hMcZRebinned" + centralityNames[iCent], ";z", zMcBinsVec[iCent].size() - 1, &zMcBinsVec[iCent][0]);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularityRebinned[iCent][iLambda] = new TH1D(Form("hMcAngularityRebinned_%d", iLambda) + centralityNames[iCent], ";angularity", angularityMcBinsVec[iCent][iLambda].size() - 1, &angularityMcBinsVec[iCent][iLambda][0]);
        }
        hRecoPtRebinned[iCent] = new TH1D("hRecoPtRebinned" + centralityNames[iCent], ";pt", ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0]);
        hRecoZRebinned[iCent] = new TH1D("hRecoZRebinned" + centralityNames[iCent], ";z", zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0]);

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hRecoAngularityRebinned[iCent][iLambda] = new TH1D(Form("hRecoAngularityRebinned_%d", iLambda) + centralityNames[iCent], ";angularity", angularityRecoBinsVec[iCent][iLambda].size() - 1, &angularityRecoBinsVec[iCent][iLambda][0]);
        }
    }

    // fill them with tree data

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++)
    {
        Float_t progress = 0.;
        progress = (Float_t)iEntry / nEntries;
        if (iEntry % 10000 == 0)
        {
            cout << "New binning : \r (" << (progress * 100.0) << "%)" << std::flush;
        }
        jetTree->GetEntry(iEntry);

        Int_t centBin = getCentralityBin(recoJet.centrality);
        Bool_t isRecoJetFound = (recoJet.numberofconstituents != 0);

        hMcPtRebinned[centBin]->Fill(mcJet.jetpt);
        hMcZRebinned[centBin]->Fill(mcJet.d0z);

        if (isRecoJetFound)
        {
            hRecoPtRebinned[centBin]->Fill(recoJet.jetpt);
            hRecoZRebinned[centBin]->Fill(recoJet.d0z);
        }

        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularityRebinned[centBin][iLambda]->Fill(mcJet.lambda[iLambda]);
            if (isRecoJetFound)
                hRecoAngularityRebinned[centBin][iLambda]->Fill(recoJet.lambda[iLambda]);
        }
    } // end of loop over  entries

    TFile *outFile = new TFile("testHists.root", "RECREATE");

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        outFile->cd();
        TDirectory *centDir = outFile->mkdir(Form("cent_%s", centralityNames[iCent].Data()));
        centDir->cd();

        hMcPt[iCent]->Write();
        hMcZ[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularity[iCent][iLambda]->Write();
        }
        hRecoPt[iCent]->Write();
        hRecoZ[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hRecoAngularity[iCent][iLambda]->Write();
        }

        hMcPtRebinned[iCent]->Write();
        hMcZRebinned[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hMcAngularityRebinned[iCent][iLambda]->Write();
        }
        hRecoPtRebinned[iCent]->Write();
        hRecoZRebinned[iCent]->Write();
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hRecoAngularityRebinned[iCent][iLambda]->Write();
        }
    }
    outFile->Save();
    outFile->Close();
}