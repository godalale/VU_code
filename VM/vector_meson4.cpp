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

void vector_meson4(const Int_t numberOfFiles = 2, const Int_t firstLine = 1) {
  Double_t eta_min = -4.14;
  Double_t eta_max = -1.18;

  TString listFile = "old_files.txt";
  TChain* chain = getChainedTree(listFile, "events", numberOfFiles, firstLine);

  TTreeReader reader(chain);

  TTreeReaderValue<std::vector<edm4eic::ClusterData>> hcalClusters(reader, "HcalEndcapNClusters");
  TTreeReaderValue<std::vector<edm4eic::ReconstructedParticleData>> reconstructedParticles(reader, "ReconstructedParticles");
  TTreeReaderValue<std::vector<podio::ObjectID>> clusterRecIDs(reader, "_HcalEndcapNClusterAssociations_rec");
  TTreeReaderValue<std::vector<podio::ObjectID>> recoToMCIndex(reader, "_ReconstructedParticleAssociations_rec");
  TTreeReaderValue<std::vector<podio::ObjectID>> mcFromRecoIndex(reader, "_ReconstructedParticleAssociations_sim");
  TTreeReaderValue<std::vector<podio::ObjectID>> particleClusterIDs(reader, "_ReconstructedParticles_clusters");

  TTreeReaderArray<Int_t> mcPDG(reader, "MCParticles.PDG");
  TTreeReaderArray<UInt_t> mcParentStart(reader, "MCParticles.parents_begin");
  TTreeReaderArray<UInt_t> mcParentEnd(reader, "MCParticles.parents_end");
  TTreeReaderArray<Int_t> mcParentIndex(reader, "_MCParticles_parents.index");

  const Float_t r_Hcal = 3000.0;
  const Int_t nbins = 1000;
  TH2F* hClusterXY = new TH2F("hClusterXY", "Cluster position of #phi->K^{+}K^{-} daughters;x [mm];y [mm]", nbins, -r_Hcal, r_Hcal, nbins, -r_Hcal, r_Hcal);

  TH1D* hNumerator = new TH1D("hNumerator", "Kaons from #phi with cluster in acceptance; p_{T} [MeV]; Count", 200, 0., 2000.);
  TH1D* hDenominator = new TH1D("hDenominator", "All kaons from #phi in acceptance; p_{T} [MeV]; Count", 200, 0., 2000.);
  hNumerator->Sumw2();
  hDenominator->Sumw2();

  Long64_t nEntries = chain->GetEntries();
  std::cout << "\n Processing " << nEntries << " events... \n" << std::endl;

  Long64_t i = 0;
  while (reader.Next()) {
    if (i % 1000 == 0) {
      std::cout << "Currently processing event: " << i << "          " << (Int_t)((1. * i / nEntries) * 100) << "% done" << std::endl;
    }

    for (size_t clusterIndex = 0; clusterIndex < hcalClusters->size(); ++clusterIndex) {
      const auto& cluster = hcalClusters->at(clusterIndex);
      const auto& clusterAssoc = clusterRecIDs->at(clusterIndex);
      Int_t recoIndex = clusterAssoc.index;
      if (recoIndex < 0 || recoIndex >= (Int_t)reconstructedParticles->size()) continue;

      const auto& particle = reconstructedParticles->at(recoIndex);
      if (std::abs(particle.PDG) != 321) continue;
      if (recoIndex >= (Int_t)particleClusterIDs->size()) continue;
      if (particleClusterIDs->at(recoIndex).index != (Int_t)clusterIndex) continue;

      Int_t mcIndex = -1;
      for (size_t j = 0; j < recoToMCIndex->size(); ++j) {
        if ((Int_t)recoToMCIndex->at(j).index == recoIndex) {
          mcIndex = mcFromRecoIndex->at(j).index;
          break;
        }
      }
      if (mcIndex < 0 || mcIndex >= mcPDG.GetSize()) continue;
      if (std::abs(mcPDG[mcIndex]) != 321) continue;

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
      if (!fromPhi) continue;

      Double_t px = particle.momentum.x;
      Double_t py = particle.momentum.y;
      Double_t pz = particle.momentum.z;
      Double_t p  = TMath::Sqrt(px*px + py*py + pz*pz);
      Double_t eta = 0.5 * TMath::Log((p + pz) / (p - pz));
      if (eta < eta_min || eta > eta_max) continue;

      Double_t pt = TMath::Hypot(px, py) * 1000.;
      hNumerator->Fill(pt);
      hClusterXY->Fill(cluster.position.x, cluster.position.y);
    }



    // Separate loop to fill denominator with all phi->K kaons in acceptance
    for (size_t recoIndex = 0; recoIndex < reconstructedParticles->size(); ++recoIndex) {
      const auto& particle = reconstructedParticles->at(recoIndex);
      if (std::abs(particle.PDG) != 321) continue;

      Int_t mcIndex = -1;
      for (size_t j = 0; j < recoToMCIndex->size(); ++j) {
        if ((Int_t)recoToMCIndex->at(j).index == recoIndex) {
          mcIndex = mcFromRecoIndex->at(j).index;
          break;
        }
      }
      if (mcIndex < 0 || mcIndex >= mcPDG.GetSize()) continue;
      if (std::abs(mcPDG[mcIndex]) != 321) continue;

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
      if (!fromPhi) continue;

      Double_t px = particle.momentum.x;
      Double_t py = particle.momentum.y;
      Double_t pz = particle.momentum.z;
      Double_t p  = TMath::Sqrt(px*px + py*py + pz*pz);
      Double_t eta = 0.5 * TMath::Log((p + pz) / (p - pz));
      if (eta < eta_min || eta > eta_max) continue;

      Double_t pt = TMath::Hypot(px, py) * 1000.;
      hDenominator->Fill(pt);
    }
    ++i;
  }

  TH1D *hEfficiency = (TH1D*)hNumerator->Clone("hEfficiency");
  hEfficiency->Divide(hDenominator,1,1,"b");

  TString pdfName = "pdfs/vector_meson_plots_4_2.pdf";

  TCanvas c1("c1", "Cluster position", 800, 800);
  c1.Print(pdfName + "[");
  gStyle->SetPalette(55);
  hClusterXY->SetMinimum(1);
  hClusterXY->Draw("surf1");
  gPad->SetLogz();

  TEllipse *circle = new TEllipse(0., 0., 2760.);
  circle->SetLineColor(kRed);
  circle->SetLineWidth(2);
  circle->SetFillStyle(0);
  circle->Draw("same");
  c1.Print(pdfName);

  TCanvas c2("c2", "Efficiency", 800, 800);
  hEfficiency->Draw("E");
  hEfficiency->SetTitle("Detection Efficiency");
  hEfficiency->GetXaxis()->SetTitle("#p_T [MeV]");
  c2.Print(pdfName);
  c2.Print(pdfName + "]");
}
