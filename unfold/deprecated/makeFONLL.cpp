#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include "TH1D.h"
#include "TFile.h"

int makeFONLL()
{
    TH1D *hist = new TH1D("FONLL", ";p_{t}, GeV/c; d#sigma/dp_{t}", 145, 1, 30);

    std::ifstream file("FONLL.txt");
    std::vector<double> y_values;

    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            if (line[0] == '#')
            {
                continue;
            }
            double pt, sigma;
            std::istringstream iss(line);
            iss >> pt >> sigma;
            hist->SetBinContent(hist->FindBin(pt + 0.0001), sigma);
        }
        file.close();
    }

    // Write histogram to output file
    TFile *outfile = new TFile("FONLL.root", "RECREATE");
    hist->Write();
    outfile->Close();

    return 0;
}
