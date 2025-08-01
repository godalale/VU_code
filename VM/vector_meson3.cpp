#include <TChain.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <vector>
#include <cmath>
#include <edm4eic/ClusterData.h>
#include <edm4eic/ReconstructedParticleData.h>
#include <podio/ObjectID.h>
#include <fstream>
#include <iostream>

TChain* getChainedTree(const TString& listFilePath, const TString& treeName, Int_t maxFiles, Int_t startLine) {
  std::ifstream infile(listFilePath.Data());
  std::string line;
  TChain* chain = new TChain(treeName);
  Int_t count = 0;
  Int_t lineNumber = 0;

  while (std::getline(infile, line) && count < maxFiles) {
    ++lineNumber;
    if (lineNumber < startLine) continue;
    chain->Add(line.c_str());
    ++count;
  }

  std::cout << "Added " << count << " files to the TChain, starting from line " << startLine << "." << std::endl;
  return chain;
}

void vector_meson3(const Int_t numberOfFiles = 2, const Int_t firstLine = 1) {
  Int_t kaons = 0;
  Int_t kaonsFromPhi = 0;
  Double_t eta_min = -4.14;
  Double_t eta_max = -1.18;

  TString listFile = "old_files.txt";
  TChain* chain = getChainedTree(listFile, "events", numberOfFiles, firstLine);

  TTreeReader reader(chain);

  TTreeReaderValue<std::vector<edm4eic::ClusterData>> hcalClusters(reader, "HcalEndcapNClusters");
  TTreeReaderValue<std::vector<edm4eic::ReconstructedParticleData>> reconstructedParticles(reader, "ReconstructedParticles");
  TTreeReaderValue<std::vector<podio::ObjectID>> clusterToRecoIndex(reader, "_HcalEndcapNClusterAssociations_rec");
  TTreeReaderValue<std::vector<podio::ObjectID>> recoToClusterIndex(reader, "_ReconstructedParticles_clusters");
  TTreeReaderValue<std::vector<podio::ObjectID>> recoToMCIndex(reader, "_ReconstructedParticleAssociations_rec");
  TTreeReaderValue<std::vector<podio::ObjectID>> mcFromRecoIndex(reader, "_ReconstructedParticleAssociations_sim");

  TTreeReaderArray<Int_t> mcPDG(reader, "MCParticles.PDG");
  TTreeReaderArray<UInt_t> mcParentStart(reader, "MCParticles.parents_begin");
  TTreeReaderArray<UInt_t> mcParentEnd(reader, "MCParticles.parents_end");
  TTreeReaderArray<Int_t> mcParentIndex(reader, "_MCParticles_parents.index");

  const Float_t r_Hcal = 2670.0;
  const Int_t nbins = 100;
  TH2F* hClusterXY = new TH2F("hClusterXY", "Cluster position of #phi->K^{+}K^{-} daughters;x [mm];y [mm]", nbins, -r_Hcal, r_Hcal, nbins, -r_Hcal, r_Hcal);

  TH1D* hNumerator = new TH1D("hNumerator", "Kaons from #phi with cluster in acceptance; p_{T} [GeV]; Count", 100, 0., 50.);
  TH1D* hDenominator = new TH1D("hDenominator", "All kaons from #phi in acceptance; p_{T} [GeV]; Count", 100, 0., 50.);
  TH1D* hEfficiency = new TH1D("hEfficiency", "Cluster efficiency for kaons from #phi; p_{T} [GeV]; Efficiency", 100, 0., 50.);

  Long64_t nEntries = chain->GetEntries();
  std::cout << "\n Processing " << nEntries << " events... \n" << std::endl;

  Long64_t i = 0;
  while (reader.Next()) {
    if (i % 1000 == 0) {
      std::cout << "Currently processing event: " << i << "          " << (Int_t)((1. * i / nEntries) * 100) << "% done" << std::endl;
    }

    for (size_t recoIndex = 0; recoIndex < reconstructedParticles->size(); ++recoIndex) {
      const auto& particle = reconstructedParticles->at(recoIndex);
      if (std::abs(particle.PDG) != 321) continue; // only kaons+- from now on

      int mcIndex = -1;
      for (size_t j = 0; j < recoToMCIndex->size(); ++j) {
        if ((Int_t)recoToMCIndex->at(j).index == recoIndex) {
          mcIndex = mcFromRecoIndex->at(j).index;
          break;
        }
      }
      if (mcIndex < 0 || mcIndex >= mcPDG.GetSize()) continue;
      if (std::abs(mcPDG[mcIndex]) != 321) continue;


      // determination, if is a phi->KK kaon
      UInt_t parentStart = mcParentStart[mcIndex];
      UInt_t parentEnd = mcParentEnd[mcIndex];

      Bool_t fromPhi = false;
      for (UInt_t ip = parentStart; ip < parentEnd; ++ip) {
        Int_t parentIdx = mcParentIndex[ip];
        if (parentIdx >= 0 && parentIdx < mcPDG.GetSize() && mcPDG[parentIdx] == 333) {
          fromPhi = true;
          break;
        }
      }
      if (!fromPhi) continue; // only phi->KK kaons from now on

      Double_t px = particle.momentum.x;
      Double_t py = particle.momentum.y;
      Double_t pz = particle.momentum.z;
      Double_t p  = TMath::Sqrt(px*px + py*py + pz*pz);
      Double_t eta = 0.5 * TMath::Log((p + pz) / (p - pz));

      if (eta < eta_min || eta > eta_max) continue;
      hDenominator->Fill(std::hypot(particle.momentum.x, particle.momentum.y));

      int clusterIndex = -1;
      if (recoIndex < (int)recoToClusterIndex->size()) {
        clusterIndex = recoToClusterIndex->at(recoIndex).index;
      }
      if (clusterIndex >= 0 && clusterIndex < (int)hcalClusters->size()) {
        const auto& cluster = hcalClusters->at(clusterIndex);
        hNumerator->Fill(std::hypot(particle.momentum.x, particle.momentum.y));
        hClusterXY->Fill(cluster.position.x, cluster.position.y);
      }
    }
    ++i;
  }

  hNumerator->Sumw2();
  hDenominator->Sumw2();
  hEfficiency->Sumw2();
  hEfficiency->Divide(hNumerator, hDenominator);

  TString pdfName = "pdfs/vector_meson_plots_3.pdf";

  TCanvas c1("c1", "Cluster position", 800, 800);
  c1.Print(pdfName + "[");
  gStyle->SetPalette(55);
  hClusterXY->SetMinimum(1);
  hClusterXY->Draw("COLZ");
  gPad->SetLogz();
  c1.Print(pdfName);

  TCanvas c2("c2", "Efficiency", 800, 800);
  hEfficiency->Draw("E");
  c2.Print(pdfName);
  c2.Print(pdfName + "]");
}

// num = th1(100, 0 , 50 )
// denom = same

// if track is a daughter of phi from mc information  // via a loop
// if (track in eta aceptance and has association with nhcal)
// num->Fill(pt)
// if track in eta acceptance
// denom->Fill(pt)

// eff = num/denom
// x = pT
// y = numerator denominator (all reconstrected Kaons in nHCal acceptance)
