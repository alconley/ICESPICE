// ROOT macro file for plotting example B4 histograms
//
// Can be run from ROOT session:
// root[0] .x plotHisto.C

{

gROOT->Reset();
gROOT->SetStyle("Plain");

// Open file filled by Geant4 simulation
TFile f("build/MiniOrange.root");

// Create a canvas and divide it into 2x2 pads
TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);

TH1D* hist1 = (TH1D*)f.Get("Esil");
hist1->Draw("HIST");

// set y axis to log scale
c1->SetLogy();

hist1->GetXaxis()->SetTitle("Energy [MeV]");
hist1->GetYaxis()->SetTitle("Counts");

// get counts in the histogram
Int_t total_counts = hist1->GetEntries();

Int_t total_counts_4pi = total_counts*2; //only sending particle in 2pi in the simulation

// get the counts between the first and the last bin
Int_t non_zero_events = hist1->Integral(2, hist1->GetNbinsX());

// divde the counts by the total counts
Double_t max_transmission_prob = (Double_t)non_zero_events / (Double_t)total_counts_4pi * 100.0; //in percent 

// print the counts to terminal
std::cout << "Total Counts (2π): " << total_counts << std::endl;
std::cout << "Total Counts (4π): " << total_counts << std::endl;
std::cout << "Non-zero Counts: " << non_zero_events << std::endl;
std::cout << "Max transmission (All counts above the first bin): " << max_transmission_prob << "%" << std::endl;

}
