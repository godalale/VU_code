#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TIterator.h>
#include <TList.h>
#include <iostream>
#include <vector>
#include <TStyle.h>

void cluster_R() {
    // Open a canvas
    TCanvas *can = new TCanvas("can", "Histograms", 800, 600);
    TString outfilename = "pdfs/clusters_R.pdf";

    // Create a PDF file to store the histograms
    can->Print(outfilename + "[");

    std::vector<int> r_position = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270};
    // Iterate over the specified files
    for (auto i : r_position) {
        TString filename = TString::Format("processed_eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", i);

        // Open the .root file
        TFile *file = TFile::Open("hists/" + filename);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }

        // Get the histogram
        TH1D *hist = nullptr;
        file->GetObject(Form("hnClusters_%dcm", i), hist);
        hist->SetTitle(TString::Format("R = %i cm", i));
        if (!hist) {
            std::cerr << "Histogram not found in file: " << filename << std::endl;
            file->Close();
            continue;
        }

        // Draw the histogram on the canvas
        gStyle->SetOptStat(0); // Hide the statistics box
        hist->GetYaxis()->SetTitle("");
        hist->Draw();

        // Save the canvas to the PDF file
        can->Print(outfilename);

        // Close the file
        file->Close();
    }

    // Close the PDF file
    can->Print(outfilename + "]");
}
