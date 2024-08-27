#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include <iostream>

void Plots(double beamEnergy);

void Plots(double beamEnergy) {
    // Open the ROOT file with update option
    TFile* file = TFile::Open("build/ICESPICE.root", "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file build/ICESPICE.root" << std::endl;
        return;
    }

    // Get the histogram named "Esil"
    TH1F* srcHist = static_cast<TH1F*>(file->Get("Esil"));
    if (!srcHist) {
        std::cerr << "Error: Could not find histogram Esil" << std::endl;
        return;
    }

    // Function for FWHM in MeV to σ conversion
    auto fwhmToSigma = [](double fwhm) {
        return (fwhm / 1000.0) / 2.3548;
    };

    // Define Gaussian smearing parameters
    double fwhm = 2.0; // FWHM in keV
    TRandom3 rand;

    // Create a new histogram for the convoluted data
    TH1F* outHist = static_cast<TH1F*>(srcHist->Clone("ESil_Folded"));
    outHist->SetTitle(Form("Folded with Gaussian (FWHM=%.2f keV)", fwhm));
    outHist->GetXaxis()->SetTitle("Energy [MeV]");
    outHist->GetYaxis()->SetTitle("Counts/keV");
    outHist->Reset(); // Clear the histogram
    outHist->SetOption("HIST"); // Set options to draw as bar charts

    // Set the title of the original histogram
    srcHist->SetTitle("Energy Deposition in PIPS");
    srcHist->GetXaxis()->SetTitle("Energy [MeV]");
    srcHist->GetYaxis()->SetTitle("Counts/keV");

    // Loop over all bins in the original histogram
    int maxFilledBin = srcHist->FindLastBinAbove();
    for (int i = 1; i <= maxFilledBin; i++) {
        double counts = srcHist->GetBinContent(i);
        double binCenter = srcHist->GetBinCenter(i);

        // Calculate the standard deviation σ at this bin energy
        double sigma = fwhmToSigma(fwhm) * binCenter;

        // Smear the content of each bin
        for (int j = 0; j < counts; j++) {
            double smearedValue = rand.Gaus(binCenter, sigma);
            outHist->Fill(smearedValue);
        }
    }

    // Calculate the transmission probability: sum counts above the first bin
    double totalCounts = 0.0;
    double totalInteractionCounts = 0.0;
    for (int i = 1; i <= srcHist->GetNbinsX(); i++) {
        double binContent = srcHist->GetBinContent(i);
        totalCounts += binContent;
        if (i > 1) { // Exclude the first bin
            totalInteractionCounts += binContent;
        }
    }
    double fullEnergyCounts = srcHist->GetBinContent(beamEnergy);
    double transmissionProbability = totalInteractionCounts / totalCounts * 100.0 / 2.0;
    double transmissionProbabilityFullEnergy = fullEnergyCounts / totalCounts * 100.0 / 2.0;
    double detectorEfficiency = fullEnergyCounts / totalInteractionCounts * 100.0;

    // Print out the statistics
    std::cout << "\tTotal Counts (2pi): " << totalCounts << std::endl;
    std::cout << "\tCounts in Detector: " << totalInteractionCounts << std::endl;
    std::cout << "\tCounts at Full Energy: " << fullEnergyCounts << std::endl;
    std::cout << "\tTransmission Probability in 4pi (Counts in Detector / (Total Counts*2) * 100) [%]: " << transmissionProbability << "%" << std::endl;
    std::cout << "\tFull Energy Transmission Probability in 4pi (Full Energy Deposited Counts/ (Total Counts*2) *100) [%]: " << transmissionProbabilityFullEnergy << "%" << std::endl;
    std::cout << "\tDetector Efficiency (Counts in Detector / Counts at Full Energy * 100) [%]: " << detectorEfficiency << "%" << std::endl;

    // Create a canvas with six pads (3x2 layout)
    TCanvas* c1 = new TCanvas("c1", "Plots", 1200, 900);
    c1->Divide(2, 3);  // Divide the canvas into 2 columns and 3 rows

    // Draw the original histogram in the first pad
    c1->cd(1);  // Select the first pad
    gPad->SetLogy();  // Set the y-axis to a logarithmic scale
    srcHist->SetLineColor(kBlue);
    srcHist->Draw();

    // Draw the convoluted histogram in the second pad
    c1->cd(2);  // Select the second pad
    gPad->SetLogy();  // Set the y-axis to a logarithmic scale
    outHist->SetLineColor(kRed);
    outHist->Draw("HIST");

    // Create and draw the transmission probability on the third pad
    c1->cd(3);  // Select the third pad
    TGraph* transProbGraph = new TGraph();
    transProbGraph->SetPoint(0, beamEnergy / 1000.0, transmissionProbability); // Convert keV to MeV for x-axis
    transProbGraph->SetMarkerStyle(21);
    transProbGraph->SetMarkerColor(kBlue);
    transProbGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    transProbGraph->GetYaxis()->SetRangeUser(0, 10); // Set y-axis range to include 0
    transProbGraph->GetXaxis()->SetTitle("Energy [MeV]");
    transProbGraph->GetYaxis()->SetTitle("Probability [%]");
    transProbGraph->SetTitle("4pi- Transmission Probability");
    transProbGraph->Draw("AP");

    // Create and draw the full energy transmission probability on the fourth pad
    c1->cd(4);  // Select the fourth pad
    TGraph* fullEnergyProbGraph = new TGraph();
    fullEnergyProbGraph->SetPoint(0, beamEnergy / 1000.0, transmissionProbabilityFullEnergy); // Convert keV to MeV for x-axis
    fullEnergyProbGraph->SetMarkerStyle(22);
    fullEnergyProbGraph->SetMarkerColor(kRed);
    fullEnergyProbGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    fullEnergyProbGraph->GetYaxis()->SetRangeUser(0, 10); // Set y-axis range to include 0
    fullEnergyProbGraph->GetXaxis()->SetTitle("Energy [MeV]");
    fullEnergyProbGraph->GetYaxis()->SetTitle("Probability [%]");
    fullEnergyProbGraph->SetTitle("4pi- Full Energy Transmission Probability");
    fullEnergyProbGraph->Draw("AP");

    // Create and draw the detector efficiency in the fifth pad
    c1->cd(5);  // Select the fifth pad
    TGraph* efficiencyGraph = new TGraph();
    efficiencyGraph->SetPoint(0, beamEnergy / 1000.0, detectorEfficiency); // Convert keV to MeV for x-axis
    efficiencyGraph->SetMarkerStyle(21);
    efficiencyGraph->SetTitle("Detector Efficiency");
    efficiencyGraph->GetXaxis()->SetLimits(0, 2); // Set x-axis range from 0 to 2 MeV
    efficiencyGraph->GetYaxis()->SetRangeUser(0, 100); // Set y-axis range to include 0
    efficiencyGraph->GetXaxis()->SetTitle("Energy [MeV]");
    efficiencyGraph->GetYaxis()->SetTitle("Detector Efficiency [%]");
    efficiencyGraph->Draw("AP");

    // Update and save the canvas
    c1->Update();
    c1->Modified();
    
    // Save the canvas to a file
    // c1->SaveAs("Plots.png");

    // Write the convoluted histogram and graphs into the same ROOT file
    file->cd();           // Ensure we are in the right directory in the file
    outHist->Write();     // Write the convoluted histogram to the file
    transProbGraph->Write("TransmissionProbability"); // Write the transmission probability graph
    fullEnergyProbGraph->Write("FullEnergyTransmissionProbability"); // Write the full energy transmission probability graph
    efficiencyGraph->Write("DetectorEfficiency"); // Write the detector efficiency graph


    // Clean up
    // file->Close();
    // delete file;
}
