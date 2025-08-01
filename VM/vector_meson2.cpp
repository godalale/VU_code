#include <TChain.h>
#include <TString.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
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
  int count = 0;
  int lineNumber = 0;

  while (std::getline(infile, line) && count < maxFiles) {
    ++lineNumber;
    if (lineNumber < startLine) continue;
    chain->Add(line.c_str());
    ++count;
  }

  std::cout << "Added " << count << " files to the TChain, starting from line " << startLine << "." << std::endl;
  return chain;
}

void vector_meson2(const Int_t numberOfFiles = 2, const Int_t firstLine = 1) {
  TString listFile = "old_files.txt";
  TChain* chain = getChainedTree(listFile, "events", numberOfFiles, firstLine);

  std::vector<edm4eic::ClusterData>* clusters = nullptr;
  std::vector<edm4eic::ReconstructedParticleData>* particles = nullptr;
  std::vector<podio::ObjectID>* clusterRecIDs = nullptr;
  std::vector<podio::ObjectID>* particleClusterIDs = nullptr;
  std::vector<podio::ObjectID>* recoAssoc = nullptr;
  std::vector<podio::ObjectID>* simuAssoc = nullptr;

  // Set legacy-style branch addresses
  chain->SetBranchAddress("HcalEndcapNClusters", &clusters);
  chain->SetBranchAddress("ReconstructedParticles", &particles);
  chain->SetBranchAddress("_HcalEndcapNClusterAssociations_rec", &clusterRecIDs);
  chain->SetBranchAddress("_ReconstructedParticles_clusters", &particleClusterIDs);
  chain->SetBranchAddress("_ReconstructedParticleAssociations_rec", &recoAssoc);
  chain->SetBranchAddress("_ReconstructedParticleAssociations_sim", &simuAssoc);

  // Set up TTreeReader for MC truth info
  TTreeReader reader(chain);
  TTreeReaderArray<int> mcPdg(reader, "MCParticles.PDG");
  TTreeReaderArray<unsigned int> parentsBegin(reader, "MCParticles.parents_begin");
  TTreeReaderArray<unsigned int> parentsEnd(reader, "MCParticles.parents_end");
  TTreeReaderArray<int> parentsIndex(reader, "_MCParticles_parents.index");

  // Histograms
  const float r_Hcal = 2670.0;
  const int nbins = 100;
  TH2F* hxy_assoc = new TH2F("hxy_assoc", "Cluster position of #phi->K^{+}K^{-} daughters;x [mm];y [mm]", nbins, -r_Hcal, r_Hcal, nbins, -r_Hcal, r_Hcal);
  TH2F* h_energy_assoc = new TH2F("h_energy_assoc", "Cluster Energy vs Distance (phi->K^{+}K^{-}); Cluster Energy [GeV]; Distance [mm]", 100, 0, 20, nbins, 0, r_Hcal);

  Long64_t nEntries = chain->GetEntries();
  std::cout << "Processing " << nEntries << " events..." << std::endl;

  for (Long64_t i = 0; i < nEntries; ++i) {
    reader.SetEntry(i);
    chain->GetEntry(i);

    if (i % 1000 == 0) {
      std::cout << "Currently processing event: " << i << "          " << (Int_t)((1. * i / nEntries) * 100) << "% done" << std::endl;
    }

    if (!clusters || !particles || !clusterRecIDs || !particleClusterIDs || !recoAssoc || !simuAssoc) {
      continue;
    }

    for (size_t icl = 0; icl < clusters->size(); ++icl) {
      const auto& cluster = clusters->at(icl);
      const auto& crec = clusterRecIDs->at(icl);
      int recIndex = crec.index;
      if (recIndex < 0 || recIndex >= (int)particles->size()) continue;

      const auto& part = particles->at(recIndex);
      if (std::abs(part.PDG) != 321) continue; // reconstructed kaon
      if (recIndex >= (int)particleClusterIDs->size()) continue;
      if (particleClusterIDs->at(recIndex).index != (int)icl) continue;

      // find MC matching index
      int mcIndex = -1;
      for (size_t j = 0; j < recoAssoc->size(); ++j) {
        if ((int)recoAssoc->at(j).index == recIndex) {
          mcIndex = simuAssoc->at(j).index;
          break;
        }
      }
      if (mcIndex < 0 || mcIndex >= mcPdg.GetSize()) continue;
      if (std::abs(mcPdg[mcIndex]) != 321) continue;

      // check for phi(1020) parent
      unsigned int pBeg = parentsBegin[mcIndex];
      unsigned int pEnd = parentsEnd[mcIndex];
      bool fromPhi = false;
      for (unsigned int ip = pBeg; ip < pEnd; ++ip) {
        int pidx = parentsIndex[ip];
        if (pidx >= 0 && pidx < mcPdg.GetSize() && mcPdg[pidx] == 333) {
          fromPhi = true;
          break;
        }
      }
      if (!fromPhi) continue;

      // fill histos
      hxy_assoc->Fill(cluster.position.x, cluster.position.y);
      h_energy_assoc->Fill(cluster.energy, std::hypot(cluster.position.x, cluster.position.y));
    }
  }

  // Save PDF
  TString pdfName = "pdfs/vector_meson_phi_plots.pdf";

  TCanvas c1("c1", "Cluster position", 800, 800);
  c1.Print(pdfName + "[");
  gStyle->SetPalette(55);
  hxy_assoc->SetMinimum(1);
  hxy_assoc->Draw("COLZ");
  gPad->SetLogz();
  c1.Print(pdfName);

  TCanvas c2("c2", "Cluster Energy vs Distance", 800, 800);
  h_energy_assoc->SetMinimum(1);
  h_energy_assoc->Draw("COLZ");
  gPad->SetLogz();
  c2.Print(pdfName);
  c2.Print(pdfName + "]");
}
