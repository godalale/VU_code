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
#include <TStyle.h>

void cluster_x_y() {
    // Open a canvas
    TCanvas *c = new TCanvas("c", "Histograms", 800, 600);
    TString outfilename = "pdfs/cluster_x_y.pdf";

    // Create a PDF file to store the histograms
    c->Print(outfilename + "[");

    std::vector<int> r_position = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270};
    for (auto r : r_position){
        TString filename = TString::Format("hists/processed_eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", r);

        // Open the .root file
        TFile *file = TFile::Open(filename);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }

        // Get the histogram
        TH2D *hist = nullptr;
        file->GetObject("HCal_Reco/hPosition_HCal_Reco", hist);
        if (!hist) {
            std::cerr << "Histogram not found in file: " << filename << std::endl;
            file->Close();
            continue;
        }

        // Draw the histogram and add it to the PDF
        hist->SetTitle(Form("Cluster x-y position for r = %d cm", r));

        hist->GetYaxis()->SetTitleOffset(1.2);

        hist->Draw();
        gStyle->SetOptStat(0);
        c->Print(outfilename);

        // Close the file
        file->Close();
    }

    // Close the PDF file
    c->Print(outfilename + "]");
}
