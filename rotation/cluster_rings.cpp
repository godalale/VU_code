#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TIterator.h>
#include <TList.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <TStyle.h>

void cluster_rings() {
    // Open a canvas
    TCanvas *can = new TCanvas("can", "Histograms", 800, 600);
    TString outfilename = "pdfs/cluster_rings.pdf";

    // Create a PDF file to store the histograms
    can->Print(outfilename + "[");

    std::vector<int> r_position = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270};

    Int_t bin_width = 10;

    // Iterate over the specified files
    for (auto i : r_position) {
        TString filename = TString::Format("processed_eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", i);

        // Open the .root file
        TFile *infile = TFile::Open("hists/" + filename);
        if (!infile || infile->IsZombie()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }

        // Get the histogram
        TH1D *hist = nullptr;
        infile->GetObject(Form("hRadius_%dcm", i), hist);
        hist->SetTitle(TString::Format("R = %i cm", i));
        if (!hist) {
            std::cerr << "Histogram not found in file: " << filename << std::endl;
            infile->Close();
            continue;
        }

        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
            //std::cout << "current bin " << bin << std::endl;
            Double_t r_in = r_position[bin-1]*1.0;
            Double_t r_out = r_in + bin_width*1.0;
            //std::cout << "r_in " << r_in << " r_out " << r_out << std::endl;
            Double_t area = TMath::Pi() * (r_out * r_out - r_in * r_in);
            //std::cout << "area " << area << std::endl;
            Double_t binContent = hist->GetBinContent(bin);
            //std::cout << "binContent before " << binContent << std::endl;
            binContent = binContent/area;
            //std::cout << "binContent after " << binContent << std::endl;
            //std::cout << std::endl;
            hist->SetBinContent(bin, binContent);
        }

        // Draw the histogram on the canvas
        gStyle->SetOptStat(0); // Hide the statistics box
        hist->GetYaxis()->SetTitle("N_{clusters}/A_{ring}");
        hist->Draw();

        // Save the canvas to the PDF file
        can->Print(outfilename);

        // Close the file
        infile->Close();
    }

    // Close the PDF file
    can->Print(outfilename + "]");
}
