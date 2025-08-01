#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <string>
#include <vector>

Double_t zGlobal = 395;

Int_t getMaximumEnergyIndex(TTreeReaderArray<Float_t> &energyVector)
{
    Int_t maxIndex = -1;
    Float_t maxEnergy = 0.;

    for (unsigned int i = 0; i < energyVector.GetSize(); ++i)
    {
        if (energyVector[i] > maxEnergy)
        {
            maxEnergy = energyVector[i];
            maxIndex = i;
        }
    }
    return maxIndex;
}

struct Cluster
{
    Cluster(TString _type) : energy(0), x(0), y(0), r(0) { type = _type; };
    Cluster() : energy(0), x(0), y(0), r(0), type("None") {};

    Double_t energy;
    Double_t x;
    Double_t y;
    Double_t r;

    Double_t xPosition;
    Double_t yPosition;

    TString type;
};

struct ClusterHists
{
    const Double_t xrange = 10.;

    ClusterHists(TString _name)
    {
        hEnergy = new TH1D("hEnergy_" + _name, _name + " Cluster energy; E [GeV]; Entries", 1000, 0, 25);
        hPos = new TH2D("hPosition_" + _name, _name + " Cluster position x,y; x [mm]; y [mm]; Entries", 5000, -3000, 3000, 500, -3000, 3000);
        name = _name;
    }

    void Fill(const Cluster &cluster)
    {
        hEnergy->Fill(cluster.energy);
        hPos->Fill(cluster.xPosition, cluster.yPosition);
    }

    void Write(TDirectory *output)
    {
        std::cout << "Writing histograms " << name << std::endl;
        output->cd();
        TDirectory *dir = output->mkdir(name);
        dir->cd();
        hEnergy->Write();
        hPos->Write();
    }

    TString name;
    TH1D *hEnergy;
    TH2D *hPos;
};

void analyzeTrees(TString inFileName, TString outFileName, int radial, TH1D *hAllRadii)
{
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleOffset(.85, "X");
    gStyle->SetTitleOffset(.85, "Y");
    gStyle->SetTitleSize(.04, "X");
    gStyle->SetTitleSize(.04, "Y");
    gStyle->SetLabelSize(.04, "X");
    gStyle->SetLabelSize(.04, "Y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);

    TFile *file = new TFile(inFileName);
    TTreeReader myReader("events", file);

    TTreeReaderArray<Int_t> status(myReader, "MCParticles.generatorStatus");
    TTreeReaderArray<Float_t> hcal_reco_E(myReader, "HcalEndcapNClusters.energy");
    TTreeReaderArray<Float_t> hcal_reco_x(myReader, "HcalEndcapNClusters.position.x");
    TTreeReaderArray<Float_t> hcal_reco_y(myReader, "HcalEndcapNClusters.position.y");

    TCanvas *can = new TCanvas("can", "can", 1200, 1000);
    can->SetMargin(0.09, 0.1, 0.1, 0.06);

    TFile *output = new TFile(outFileName, "RECREATE");
    TH1D *hnClusters = new TH1D(Form("hnClusters_%icm", radial), " Number of Clusters; Number of Clusters; Entries", 10, 0, 10);
    TH1D *hRadius = new TH1D(Form("hRadius_%icm", radial), "Cluster Radius; r [cm]; Entries", 27, 20, 280); 

    ClusterHists hcalRecoHist("HCal_Reco");

    Int_t nEvents = myReader.GetEntries();
    std::cout << "Corresponding to radius :" << radial << std::endl;
    std::cout << "Total Events: " << nEvents << std::endl;

    for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        myReader.SetEntry(iEvent);
        if (iEvent % 100000 == 0)
            std::cout << "Event No: " << iEvent << std::endl;

        Int_t indexHClusterReco = getMaximumEnergyIndex(hcal_reco_E);

        Int_t nClusters = hcal_reco_E.GetSize();
        hnClusters->Fill(nClusters);

        Cluster hcalReco("hcalReco");

        if (indexHClusterReco >= 0)
        {
            hcalReco.energy = hcal_reco_E[indexHClusterReco];
            hcalReco.xPosition = hcal_reco_x[indexHClusterReco];
            hcalReco.yPosition = hcal_reco_y[indexHClusterReco];

            // Compute radius and fill global histogram
            hcalReco.r = TMath::Sqrt(hcalReco.xPosition * hcalReco.xPosition + hcalReco.yPosition * hcalReco.yPosition)/10.; //in cm !
            hRadius->Fill(hcalReco.r);
            hAllRadii->Fill(hcalReco.r);

            hcalRecoHist.Fill(hcalReco);
        }
    }

    output->cd();
    hcalRecoHist.Write(output);
    output->cd("../");
    hnClusters->Write();
    hRadius->Write();
    output->Save();
    output->Close();

    delete can;
    delete output;
    delete file;
}

void analyzeMultipleFiles()
{
    std::vector<int> r_position = {20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270};

    // Global histogram for combined radii
    TH1D *hAllRadii = new TH1D("hAllRadii", "Cluster Radius Distribution; r [cm]; Entries", 27, 20, 280);  

    for (auto radial : r_position)
    {
        TString inFileName = TString::Format("output/eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", radial);
        TString outFileName = Form("hists/processed_eicrecon_neutron_500000events_r%icm_p2gev_phi45.edm4eic.root", radial);

        std::cout << std::endl;
        std::cout << "Analyzing file: " << inFileName << std::endl;
        analyzeTrees(inFileName, outFileName, radial, hAllRadii);
    }

    // Save combined histogram after all files processed
    TFile *finalOutput = new TFile("hists/combined_radius_distribution.root", "RECREATE");
    hAllRadii->Write();
    finalOutput->Close();

    delete hAllRadii;
}