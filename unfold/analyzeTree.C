void analyzeTree()
{

    gStyle->SetOptStat(0);
    TFile *treeFile = new TFile("~/dev/star/output_jets.root");
    TTree *Jets = (TTree *)treeFile->Get("Jets");
    Double_t centBins[4] = {0, 10, 40, 80};
    Int_t nCentrality = 3;
    Double_t xbins[14] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 30};
    TH2D *h2 = new TH2D("h2", "", 13, xbins, nCentrality, centBins);
    TH2D *h1 = new TH2D("h1", "", 13, xbins, nCentrality, centBins);
    TH2D *h3 = new TH2D("JetFinderEfficiency", "", 13, xbins, nCentrality, centBins);

    Jets->Draw("centrality:McRecoJetPt>>h2", "RecoJetNConst==0");

    Jets->Draw("centrality:McRecoJetPt>>h1", "");

    h2->Divide(h1);

    for (int xBin = 1; xBin <= h3->GetNbinsX(); xBin++)
    {
        for (int yBin = 1; yBin <= h3->GetNbinsY(); yBin++)
        {
            h3->SetBinContent(xBin, yBin, 1.);
        }
    }
    h3->Add(h2, -1.);

    h3->SetTitle("Efficiency of Jet Finder; MC jet p_{t},GeV/c; Centrality, % ; Efficiency");
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    h3->SetMarkerStyle(20);
    h3->SetMarkerColor(2002);
    // h3->GetZaxis()->SetRangeUser(0.84, 1.01);
    h3->Draw("colz");
    c1->SaveAs("JetFinderEfficiency.pdf");

    TFile *f = new TFile("JetFinderEfficiency.root", "RECREATE");
    h3->Write();
    f->Close();
}