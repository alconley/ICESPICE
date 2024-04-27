#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TKey.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <map>
#include <TDirectory.h>

// Forward declaration of processFile
void processFile(TString inputFileName, TFile* outputFile);

// Forward declaration
void plotData(const std::map<double, double>& data, const char* dirName, const char* plotTitleSuffix);


void processDirectory(const char* inputDir) {
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    TString outputFileName = "./AllModifiedOutputs.root";
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    // Open directory
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(inputDir)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            TString fileName = ent->d_name;
            if (fileName.EndsWith(".root")) {
                TString inputFileName = TString(inputDir) + "/" + fileName;
                processFile(inputFileName, outputFile);
            }
        }
        closedir(dir);
    } else {
        perror("Could not open directory");
        return;
    }

    outputFile->Write();
    outputFile->Close();
}

void processFile(TString inputFileName, TFile* outputFile) {
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Cannot open file: " << inputFileName << std::endl;
        return;
    }

    std::stringstream ss(inputFileName.Data());
    std::string item, lastItem;
    std::vector<std::string> tokens;
    while (getline(ss, item, '/')) { lastItem = item; }
    std::stringstream ss2(lastItem);
    while (getline(ss2, item, '_')) { tokens.push_back(item); }

    std::string histName = tokens.back();
    size_t pos = histName.find(".root");
    if (pos != std::string::npos) { histName = histName.substr(0, pos); }
    tokens.pop_back(); // Remove last item which includes ".root"

    std::string currentPath = "";
    outputFile->cd(); // Ensure starting from the root directory of outputFile
    for (auto& token : tokens) {
        currentPath += (currentPath.empty() ? "" : "/") + token;
        if (!outputFile->GetDirectory(currentPath.c_str())) {
            outputFile->mkdir(currentPath.c_str());
        }
        outputFile->cd(currentPath.c_str());
    }

    TH1D* hist = (TH1D*)inputFile->Get("Esil");
    if (!hist) {
        std::cerr << "Histogram 'Esil' not found in file: " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    hist->SetName(histName.c_str());
    hist->SetTitle(histName.c_str());
    hist->GetXaxis()->SetTitle("Energy [MeV]");
    hist->GetYaxis()->SetTitle("Counts");

    hist->Write();
    inputFile->Close();
}

// Define a function to process histograms and calculate transmission probabilities
void processHistograms(const char* fileName, const char* dirName) {
    // Open the ROOT file
    TFile* file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    // Navigate to the specified directory
    TDirectory* dir = file->GetDirectory(dirName);
    if (!dir) {
        std::cerr << "Directory not found: " << dirName << std::endl;
        file->Close();
        return;
    }

    // Prepare to loop over all objects in this directory
    TIter nextkey(dir->GetListOfKeys());
    TKey* key;
    std::map<double, double> transmissionData;
    std::map<double, double> fullEnergyTransmissionData;

    while ((key = (TKey*)nextkey())) {
        // Ensure the object is a histogram
        if (strcmp(key->GetClassName(), "TH1F") == 0 || strcmp(key->GetClassName(), "TH1D") == 0) {
            TH1* hist = (TH1*)key->ReadObj();

            // Calculate transmission probability
            double totalCounts = hist->GetEntries();
            double totalCounts_4pi = totalCounts * 2.0;

            double nonZeroCounts = hist->Integral(2, hist->GetNbinsX());
            double transmissionProb = (totalCounts > 0) ? (nonZeroCounts / totalCounts_4pi * 100.0) : 0.0;

            // Assuming histogram name is the energy
            double energy = std::stod(hist->GetName());

            // Store the result
            transmissionData[energy] = transmissionProb;

            //print energy
            std::cout << "Energy: " << energy << std::endl;
             // Finding the bin corresponding to the full energy deposition
            int energyBin = energy; //one to one binning
            //print out energy bin
            std::cout << "Energy Bin: " << energyBin << std::endl;
            double fullEnergyCount = hist->GetBinContent(energyBin);
            //print out count
            std::cout << "Full Energy Count: " << fullEnergyCount << std::endl;
            double fullEnergyProb = (totalCounts > 0) ? (fullEnergyCount / totalCounts_4pi * 100.0) : 0.0;
            fullEnergyTransmissionData[energy] = fullEnergyProb;
        }
    }

    std::cout << "Transmission Probability in " << std::string(dirName) << std::endl;
    for (const auto& data : transmissionData) {
        std::cout << "Energy: " << data.first << " keV, Transmission Probability: " << data.second << "%" << std::endl;
    }

    // Clean up
    file->Close();

     // Plot the general transmission data
    plotData(transmissionData, dirName, "General Transmission Probability");
    // Plot the full energy specific transmission data
    plotData(fullEnergyTransmissionData, dirName, "Full Energy Deposited Transmission Probability");

}

void plotData(const std::map<double, double>& data, const char* dirName, const char* plotTitleSuffix) {
    if (data.empty()) {
        std::cerr << "No data to plot for " << plotTitleSuffix << "." << std::endl;
        return;
    }

    TGraph *graph = new TGraph();
    int i = 0;
    for (const auto &entry : data) {
        graph->SetPoint(i++, entry.first, entry.second);
    }

    std::string title = std::string(dirName) + " - " + plotTitleSuffix;
    graph->SetTitle((title + ";Energy (keV);Transmission Probability (%)").c_str());
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 800, 600);
    graph->Draw("APL");
    // c->SaveAs((title + ".png").c_str());
    c->Draw();
}











    // get counts in the histogram
    // Int_t total_counts = hist->GetEntries();

    // Int_t total_counts_4pi = total_counts*2; //only sending particle in 2pi in the simulation

    // // get the counts between the first and the last bin
    // Int_t non_zero_events = hist->Integral(2, hist->GetNbinsX());

    // // divde the counts by the total counts
    // Double_t max_transmission_prob = (Double_t)non_zero_events / (Double_t)total_counts_4pi * 100.0; //in percent 

    // // Extract energy from the file name
    // std::string histName = inputFileName(inputFileName.Last('/') + 1, inputFileName.Length());
    // Double_t energy = std::stod(histName.substr(0, histName.find("_")));

    // // Store energy and transmission probability in the map
    // transmissionData[energy] = max_transmission_prob;

    // std::cout "Energy: " << energy << "keV -> " << "Transmission Prob: " << max_transmission_prob << "%" << std::endl;
