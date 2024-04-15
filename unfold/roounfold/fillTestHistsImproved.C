#include <iostream>
#include <iomanip>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"

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

const Int_t nCentralityBins = 3;
double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order
TString centralityTitles[nCentralityBins] = {"0-10%", "10-40%", "40-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_40", "40_80"};
const Int_t nAngularities = 4;
TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}^{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}"};
TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

void check(TH1D *h, TString name)
{
    if (h == nullptr)
    {
        cout << "Error: " << name << " not found" << endl;
    }
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
void NormalizeByBinWidth(TH1D *hist)
{
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        if (hist->GetBinWidth(i) == 0)
            continue;
        hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
        hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
    }
}

struct Sizes;

struct Hists
{
    Hists(){}; // default constructor
    Hists(TString _name, const Int_t &nBins, const Int_t &iCent);
    Hists(TString _name, const Sizes &sizes, const Int_t &iCent);
    Sizes getSizes(void);
    void SetColor(const Int_t &color)
    {
        hPt->SetLineColor(color);
        hZ->SetLineColor(color);
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
        {
            hAngularity[iLambda]->SetLineColor(color);
        }
    };
    void Fill(StJetTreeStruct &jet);
    void Write();
    Hists *getNormalized(void);
    TH1D *getPt(void) { return hPt; }
    TH1D *getZ(void) { return hZ; }
    TH1D *getAngularity(Int_t iLambda) { return hAngularity[iLambda]; }

    TH1D *hPt;
    TH1D *hZ;
    TH1D *hAngularity[nAngularities];
    TString name;
};
Hists::Hists(TString _name, const Int_t &nBins, const Int_t &iCent)
{
    name = _name;
    Double_t xMinPt, xMaxPt, xMinZ, xMaxZ, xMinAngularity, xMaxAngularity;
    Int_t color;
    if (name.Contains("Mc"))
    {
        xMinPt = 1;
        xMaxPt = 30;
        xMinZ = 0;
        xMaxZ = 1.00001;
    }
    else
    {
        xMinPt = ptBinsReco[iCent].first;
        xMaxPt = ptBinsReco[iCent].second;
        xMinZ = -10;
        xMaxZ = 10;
    }

    hPt = new TH1D("hPt" + name + centralityNames[iCent], ";p_{t}, GeV/c; dN/dp_{t}", nBins, xMinPt, xMaxPt);
    hZ = new TH1D("hZ" + name + centralityNames[iCent], ";z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2};dN/dz", nBins, xMinZ, xMaxZ);
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
        hAngularity[iLambda] = new TH1D(Form("hAngularity_%d", iLambda) + name + centralityNames[iCent], ";#lambda_{#alpha}^{1}; dN/d#lambda", nBins, xMinAngularity, xMaxAngularity);
    }
};

void Hists::Fill(StJetTreeStruct &jet)
{

    if (name.Contains("Reco") && jet.numberofconstituents == 0)
    {
        return;
    }

    hPt->Fill(jet.jetpt);
    // if (jet.d0z == 1.)
    //     jet.d0z -= 0.000001;
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

Hists *Hists::getNormalized(void)
{
    Hists *ret = new Hists();
    ret->name = name;
    ret->hPt = (TH1D *)hPt->Clone((TString)hPt->GetName() + "_normalized");
    ret->hZ = (TH1D *)hZ->Clone((TString)hZ->GetName() + "_normalized");
    ret->hPt->GetYaxis()->SetTitleOffset(1.2);
    ret->hZ->GetYaxis()->SetTitleOffset(1.2);
    NormalizeByBinWidth(ret->hPt);
    NormalizeByBinWidth(ret->hZ);
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        ret->hAngularity[iLambda] = (TH1D *)hAngularity[iLambda]->Clone((TString)hAngularity[iLambda]->GetName() + "_normalized");

        NormalizeByBinWidth(ret->hAngularity[iLambda]);
        ret->hAngularity[iLambda]->GetYaxis()->SetTitleOffset(1.2);
    }
    return ret;
}

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

vector<Double_t> getAxis(TH1D *hist)
{
    vector<Double_t> axis;
    for (Int_t i = 1; i <= hist->GetNbinsX(); i++)
    {
        axis.push_back(hist->GetBinLowEdge(i));
    }
    axis.push_back(hist->GetBinLowEdge(hist->GetNbinsX() + 1));
    return axis;
}

Sizes Hists::getSizes(void)
{
    Sizes sizes;
    sizes.pt = getAxis(hPt);
    sizes.z = getAxis(hZ);
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        sizes.lambda[iLambda] = getAxis(hAngularity[iLambda]);
    }
    return sizes;
}

Hists::Hists(const TString _name, const Sizes &sizes, const Int_t &iCent)
{
    name = _name;

    hPt = new TH1D("hPt" + name + centralityNames[iCent], ";p_{t}, GeV/c; dN/dp_{t}", sizes.pt.size() - 1, &sizes.pt[0]);
    hZ = new TH1D("hZ" + name + centralityNames[iCent], ";z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2};dN/dz", sizes.z.size() - 1, &sizes.z[0]);
    for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
    {
        hAngularity[iLambda] = new TH1D(Form("hAngularity_%d", iLambda) + name + centralityNames[iCent], ";#lambda_{#alpha}^{1}; dN/d#lambda", sizes.lambda[iLambda].size() - 1, &sizes.lambda[iLambda][0]);
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

    const Int_t nRecoBins = 10;
    const Int_t nMcBins = 5;

    Bool_t createFile = kFALSE;
    TFile *outFile = new TFile(Form("binningMc%iReco%i.root", nMcBins, nRecoBins), "read");
    if (outFile->IsZombie())
    {
        cout << "Creating file" << endl;
        createFile = kTRUE;
    }

    Hists Mc[nCentralityBins];
    Hists Reco[nCentralityBins];
    Hists McRebinned[nCentralityBins];
    Hists RecoRebinned[nCentralityBins];

    if (createFile)
    {
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            Mc[iCent] = Hists("Mc", 10000, iCent);
            Reco[iCent] = Hists("Reco", 10000, iCent);
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

        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            McRebinned[iCent] = Hists("McRebinned", sizesMc[iCent], iCent);
            RecoRebinned[iCent] = Hists("RecoRebinned", sizesReco[iCent], iCent);
        }

        // fill them with tree data

        FillFromTree(jetTree, McRebinned, RecoRebinned);
        outFile = new TFile(Form("binningMc%iReco%i.root", nMcBins, nRecoBins), "RECREATE");
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            outFile->cd();
            TDirectory *centDir = outFile->mkdir(Form("cent_%s", centralityNames[iCent].Data()));
            centDir->cd();
            Mc[iCent].SetColor(2000);
            Mc[iCent].Write();
            Reco[iCent].SetColor(2002);
            Reco[iCent].Write();
            McRebinned[iCent].SetColor(2001);
            McRebinned[iCent].Write();
            RecoRebinned[iCent].SetColor(2003);
            RecoRebinned[iCent].Write();
        }
    }

    if (!createFile)
    {
        cout << "Reading file" << endl;
        outFile->cd();
        for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
        {
            TString dir = Form("cent_%s/", centralityNames[iCent].Data());

            ///=======================================
            Mc[iCent].hPt = (TH1D *)outFile->Get(dir + "hPtMc" + centralityNames[iCent]);
            check(Mc[iCent].hPt, dir + "hPtMc" + centralityNames[iCent]);
            Mc[iCent].hZ = (TH1D *)outFile->Get(dir + "hZMc" + centralityNames[iCent]);
            check(Mc[iCent].hZ, dir + "hZMc" + centralityNames[iCent]);
            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                Mc[iCent].hAngularity[iLambda] = (TH1D *)outFile->Get(dir + Form("hAngularity_%iMc", iLambda) + centralityNames[iCent]);
                check(Mc[iCent].hAngularity[iLambda], dir + Form("hAngularity_%iiMc", iLambda) + centralityNames[iCent]);
            }
            ///=======================================
            Reco[iCent].hPt = (TH1D *)outFile->Get(dir + "hPtReco" + centralityNames[iCent]);
            check(Reco[iCent].hPt, dir + "hPtReco" + centralityNames[iCent]);
            Reco[iCent].hZ = (TH1D *)outFile->Get(dir + "hZReco" + centralityNames[iCent]);
            check(Reco[iCent].hZ, dir + "hZReco" + centralityNames[iCent]);

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                Reco[iCent].hAngularity[iLambda] = (TH1D *)outFile->Get(dir + Form("hAngularity_%iReco", iLambda) + centralityNames[iCent]);
                check(Reco[iCent].hAngularity[iLambda], dir + Form("hAngularity_%iReco", iLambda) + centralityNames[iCent]);
            }
            ///=======================================
            McRebinned[iCent].hPt = (TH1D *)outFile->Get(dir + "hPtMcRebinned" + centralityNames[iCent]);
            check(McRebinned[iCent].hPt, dir + "hPtMcRebinned" + centralityNames[iCent]);
            McRebinned[iCent].hZ = (TH1D *)outFile->Get(dir + "hZMcRebinned" + centralityNames[iCent]);
            check(McRebinned[iCent].hZ, dir + "hZMcRebinned" + centralityNames[iCent]);

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                McRebinned[iCent].hAngularity[iLambda] = (TH1D *)outFile->Get(dir + Form("hAngularity_%iMcRebinned", iLambda) + centralityNames[iCent]);
                check(McRebinned[iCent].hAngularity[iLambda], dir + Form("hAngularity_%iMcRebinned", iLambda) + centralityNames[iCent]);
            }
            ///=======================================
            RecoRebinned[iCent].hPt = (TH1D *)outFile->Get(dir + "hPtRecoRebinned" + centralityNames[iCent]);
            check(RecoRebinned[iCent].hPt, dir + "hPtRecoRebinned" + centralityNames[iCent]);
            RecoRebinned[iCent].hZ = (TH1D *)outFile->Get(dir + "hZRecoRebinned" + centralityNames[iCent]);
            check(RecoRebinned[iCent].hZ, dir + "hZRecoRebinned" + centralityNames[iCent]);

            for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++)
            {
                RecoRebinned[iCent].hAngularity[iLambda] = (TH1D *)outFile->Get(dir + Form("hAngularity_%iRecoRebinned", iLambda) + centralityNames[iCent]);
                check(RecoRebinned[iCent].hAngularity[iLambda], dir + Form("hAngularity_%iRecoRebinned", iLambda) + centralityNames[iCent]);
            }

            Sizes Mc = McRebinned[iCent].getSizes();
            Sizes Reco = RecoRebinned[iCent].getSizes();
            Mc.Print("McBinsVec", iCent);
            Reco.Print("RecoBinsVec", iCent);
        }
    }
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("can", "can", 900, 600);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    TString outPdf = Form("binningMc%iReco%i.pdf", nMcBins, nRecoBins);
    can->SaveAs(outPdf + "[");

    gStyle->SetPadRightMargin(0.01);
    gStyle->SetPadLeftMargin(0.15);

    for (Int_t iCent = 0; iCent < nCentralityBins; iCent++)
    {
        can->Clear();
        can->Divide(3, 2);

        cout << "Drawing " << centralityTitles[iCent] << endl;

        Hists *normMc = Mc[iCent].getNormalized();
        Hists *normReco = Reco[iCent].getNormalized();
        Hists *normMcRebinned = McRebinned[iCent].getNormalized();
        McRebinned[iCent].SetColor(kBlack);
        Hists *normRecoRebinned = RecoRebinned[iCent].getNormalized();
        RecoRebinned[iCent].SetColor(kBlack);

        can->cd(1);
        gPad->SetLogy();

        cout << "hPt" << endl;
        normMc->hPt->Draw("hist");
        normMcRebinned->hPt->Draw("hist same");
        McRebinned[iCent].hPt->Draw("hist same");
        can->cd(2);
        gPad->SetLogy();

        cout << "hZ" << endl;
        normMc->hZ->Draw("hist");
        normMcRebinned->hZ->Draw("hist same");
        McRebinned[iCent].hZ->Draw("hist same");

        can->cd(3);
        gPad->SetLogy();

        cout << "hAngularity" << endl;
        normMc->hAngularity[0]->Draw("hist");
        normMcRebinned->hAngularity[0]->Draw("hist same");
        McRebinned[iCent].hAngularity[0]->Draw("hist same");
        can->cd(4);
        gPad->SetLogy();

        cout << "hPt Rebinned" << endl;
        normReco->hPt->Draw("hist");
        normRecoRebinned->hPt->Draw("hist same");
        RecoRebinned[iCent].hPt->Draw("hist same");
        can->cd(5);
        gPad->SetLogy();

        cout << "hZ Rebinned" << endl;
        normReco->hZ->Draw("hist");
        normRecoRebinned->hZ->Draw("hist same");
        RecoRebinned[iCent].hZ->Draw("hist same");
        can->cd(6);
        gPad->SetLogy();

        cout << "hAngularity Rebinned" << endl;
        normReco->hAngularity[0]->Draw("hist");
        normRecoRebinned->hAngularity[0]->Draw("hist same");
        RecoRebinned[iCent].hAngularity[0]->Draw("hist same");

        can->cd(0);
        tex->DrawLatex(0.42, 0.85, centralityTitles[iCent]);

        tex->DrawLatex(0.44, 0.96, Form("Mc N_{bins}=%i", nMcBins));
        tex->DrawLatex(0.44, 0.47, Form("Reco N_{bins}=%i", nRecoBins));

        TLegend *leg = new TLegend(0.82, 0.85, 0.97, 0.99);
        leg->AddEntry(normMc->hPt, "Mc", "l");
        leg->AddEntry(normMcRebinned->hPt, "McRebinned", "l");
        leg->AddEntry(McRebinned[iCent].hPt, "Not normalized", "l");
        leg->AddEntry(normReco->hPt, "Reco", "l");
        leg->AddEntry(normRecoRebinned->hPt, "RecoRebinned", "l");
        leg->Draw("same");

        can->SaveAs(outPdf + "");
    }
    can->SaveAs(outPdf + "]");

    outFile->Save();
    // outFile->Close();
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