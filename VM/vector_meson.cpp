#include <TChain.h>
#include <TString.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <cmath>
#include <edm4eic/ClusterData.h>
#include <edm4eic/ReconstructedParticleData.h>
#include <podio/ObjectID.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <string>

TChain* getChainedTree(const TString& listFilePath, const TString& treeName, Int_t maxFiles, Int_t startLine) {
  std::ifstream infile(listFilePath.Data());
  std::string line;
  TChain* chain = new TChain(treeName);
  int count = 0;
  int lineNumber = 0;

  while (std::getline(infile, line) && count < maxFiles) {
    lineNumber++;
    if (lineNumber < startLine) {
      continue;
    }
    chain->Add(line.c_str());
    ++count;
  }

  std::cout << "Added " << count << " files to the TChain, starting from line " << startLine << "." << std::endl;
  return chain;
}

void vector_meson(const Int_t numberOfFiles = 2, const Int_t firstLine = 1) {
  TString listFile = "old_files.txt";
  TChain *chain = getChainedTree(listFile, "events", numberOfFiles, firstLine);

  std::vector<edm4eic::ClusterData>* clusters = nullptr;
  std::vector<edm4eic::ReconstructedParticleData>* particles = nullptr;

  std::vector<podio::ObjectID>* clusterRecIDs = nullptr;
  std::vector<podio::ObjectID>* particleClusterIDs = nullptr;

  chain->SetBranchAddress("HcalEndcapNClusters", &clusters);
  chain->SetBranchAddress("ReconstructedParticles", &particles);
  chain->SetBranchAddress("_HcalEndcapNClusterAssociations_rec", &clusterRecIDs);
  chain->SetBranchAddress("_ReconstructedParticles_clusters", &particleClusterIDs);

  TTreeReader tree_reader(chain);
  TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");

  TTreeReaderArray<int> parents_index(tree_reader, "_MCParticles_parents.index");
  TTreeReaderArray<unsigned int> parents_begin(tree_reader, "MCParticles.parents_begin");
  TTreeReaderArray<unsigned int> parents_end(tree_reader, "MCParticles.parents_end");

  TTreeReaderArray<int> daughters_index(tree_reader, "_MCParticles_daughters.index");
  TTreeReaderArray<unsigned int> daughters_begin(tree_reader, "MCParticles.daughters_begin");
  TTreeReaderArray<unsigned int> daughters_end(tree_reader, "MCParticles.daughters_end");

  float r_Hcal = 2670.0;
  const int nbins = 100;

  TH2F* hxy_assoc = new TH2F("hxy_assoc", "Cluster position;x [mm];y [mm]", nbins, -r_Hcal, r_Hcal, nbins, -r_Hcal, r_Hcal);
  TH2F* h_energy_assoc = new TH2F("h_energy_assoc", "Distance from origin vs. Cluster Energy; Cluster Energy [GeV]; Distance from origin [mm]", 100, 0, 20, nbins, 0, r_Hcal);

  Long64_t nEntries = chain->GetEntries();
  std::cout << "Number of events: " << nEntries << std::endl;

  for (Long64_t i = 0; i < nEntries; ++i) {
    chain->GetEntry(i);
    
    if (i % 1000 == 0) {
    	std::cout << "Currently processing event: " << i << "          " << (Int_t)((1.*i/nEntries)*100) << "% done" << std::endl;
    }

    if (!particles || !clusters || !clusterRecIDs || !particleClusterIDs) continue;
    for (size_t icluster = 0; icluster < clusters->size(); ++icluster) {
      const auto& cluster = clusters->at(icluster);
      const auto& clusterAssoc = clusterRecIDs->at(icluster); // Associate Cluster ID to Reconctructed Particle ID

      if (clusterAssoc.index < 0 || clusterAssoc.index >= (int)particles->size()) continue; // Out of bounds check
        

      const auto& particle = particles->at(clusterAssoc.index);
      if (TMath::Abs(particle.PDG) != 321) continue;  // kaons +-

      if (clusterAssoc.index >= (int)particleClusterIDs->size()) continue; // Safety check
    
      const auto& reverseAssoc = particleClusterIDs->at(clusterAssoc.index); // Verify Matching Particle ID to Cluster ID
      
      if (reverseAssoc.index != (int)icluster) continue;

      hxy_assoc->Fill(cluster.position.x, cluster.position.y);
      h_energy_assoc->Fill(cluster.energy, std::hypot(cluster.position.x, cluster.position.y));
    }
  }

  TString pdfName = "pdfs/vector_meson_plots.pdf";

  TCanvas* c1 = new TCanvas("c1", "Cluster position", 800, 800);
  c1->Print(pdfName + "[");

 
  gStyle->SetPalette(55);
  hxy_assoc->Draw("COLZ");
  hxy_assoc->SetMinimum(1);
  gPad->SetLogz();
  c1->Print(pdfName);

  TCanvas* c2 = new TCanvas("c2", "Cluster Energy vs Distance from origin", 800, 800);
  h_energy_assoc->Draw("COLZ");
  h_energy_assoc->SetMinimum(1);
  gPad->SetLogz();
  c2->Print(pdfName);

  c2->Print(pdfName + "]");
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
