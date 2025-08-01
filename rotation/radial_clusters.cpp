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

void averageClusters() {
    // Open a canvas
    TCanvas *can = new TCanvas("can", "Histograms", 800, 600);
    TString outfilename = "pdfs/radial_clusters.pdf";

    // Create a PDF file to store the histograms
    can->Print(outfilename + "[");

    // Create a histogram to store the average number of clusters
    TH1D *hAverage = new TH1D("hAverage", "", 300, 0, 300);
    hAverage->GetXaxis()->SetTitle("R [cm]");
    hAverage->GetYaxis()->SetTitle("#LT N_{clusters} #GT");
    hAverage->SetTitle("");

    // Set the top margin smaller
    can->SetTopMargin(0.05);
    can->SetRightMargin(0.05);

    std::vector<int> r_position = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270};
    for (auto i : r_position) {
        TString filename = TString::Format("processed_eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", i);

        // Open the .root file
        TFile *file = TFile::Open("hists/"+filename);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }

        // Get the histogram
        TH1D *hist = nullptr;
        file->GetObject(Form("hRadius_%dcm", i), hist);
        if (!hist) {
            std::cerr << "Histogram not found in file: " << filename << std::endl;
            file->Close();
            continue;
        }

        // Get mean from each histogram
        Double_t mean = hist->GetMean();
        Double_t sigma = hist->GetStdDev();
        hAverage->SetBinContent(i, mean);
        hAverage->SetBinError(i, sigma);

        // Draw only the top point of the histogram
        gStyle->SetOptStat(0); // Hide the statistics box

        hAverage->SetMarkerStyle(5); // Set marker style
        hAverage->SetMarkerSize(2.0); // Set marker size
        hAverage->SetMarkerColor(kBlack); // Set marker color
        hAverage->Draw("SAME PE"); // Draw only the points on the same canvas

        // Close the file
        file->Close();
    }

    // Save the canvas to the PDF file
    can->Print(outfilename);

    // Close the PDF file
    can->Print(outfilename + "]");
}